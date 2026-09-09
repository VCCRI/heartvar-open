from __future__ import annotations

import json
import logging
import os
import re
import sqlite3
from pathlib import Path

import httpx  # noqa: F401  (re-exported as a monkeypatch surface; tests patch uniprot.httpx.AsyncClient)

from ..localio import run_local
from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict
from ._paths import PROJECT_ROOT

log = logging.getLogger("heartvar.uniprot")

UNIPROT_API = "https://rest.uniprot.org"

_PROJECT_ROOT = PROJECT_ROOT
_DEFAULT_DB_PATH = _PROJECT_ROOT / "data" / "uniprot.db"


def _db_path() -> Path:
    """Resolve the local UniProt DB path, honouring the UNIPROT_DB_PATH
    env override at call time (so tests can repoint it per-test)."""
    override = os.environ.get("UNIPROT_DB_PATH")
    return Path(override) if override else _DEFAULT_DB_PATH


DB_PATH = _DEFAULT_DB_PATH


def _loc(feat: dict) -> tuple[int | None, int | None]:
    loc = feat.get("location") or {}
    s = (loc.get("start") or {}).get("value")
    e = (loc.get("end") or {}).get("value")
    return s, e


_TRACK_FEATURE_TYPES = (
    "Domain",
    "Region",
    "Motif",
    "Repeat",
    "Zinc finger",
    "DNA binding",
    "Active site",
    "Binding site",
    "Transmembrane",
    "Signal",
    "Coiled coil",
)


_DISEASE_PHRASE_RE = re.compile(r"^\s*in\s+(?!dbsnp\b)([A-Za-z][A-Za-z0-9-]*)", re.IGNORECASE)
_EXPLICIT_BENIGN_RE = re.compile(r"\blikely[\s-]benign\b|\bbenign\b", re.IGNORECASE)
_EXPLICIT_PATH_RE = re.compile(r"\blikely[\s-]pathogenic\b|\bpathogenic\b", re.IGNORECASE)
_EXPLICIT_VUS_RE = re.compile(r"\buncertain significance\b|\bvus\b", re.IGNORECASE)


def _classify_natural_variant(description: str | None) -> str:
    """Bucket a UniProt natural-variant description into one of:
    "pathogenic" / "benign" / "vus" / "tolerated" / "unknown".

    Explicit ACMG-style labels in the description win over the
    disease-code heuristic so curated benign / VUS calls are not
    mis-bucketed as pathogenic on the basis of an associated disease
    code (e.g. "in CMH1; benign; dbSNP:..." → benign).
    """
    if not description:
        return "unknown"
    desc = description.strip()
    if not desc:
        return "unknown"
    if _EXPLICIT_BENIGN_RE.search(desc):
        return "benign"
    if _EXPLICIT_VUS_RE.search(desc):
        return "vus"
    if _EXPLICIT_PATH_RE.search(desc):
        return "pathogenic"
    if _DISEASE_PHRASE_RE.match(desc):
        return "pathogenic"
    return "tolerated"


def _xref_strings(feat: dict) -> dict:
    """Pull dbSNP / ClinVar IDs out of UniProt's featureCrossReferences
    array. Returns ``{"dbsnp": "rs...", "clinvar": "VCV..."}`` with
    ``None`` for any missing source."""
    xrefs = feat.get("featureCrossReferences") or []
    out: dict[str, str | None] = {"dbsnp": None, "clinvar": None}
    for ref in xrefs:
        db = (ref.get("database") or "").lower()
        rid = ref.get("id")
        if not rid:
            continue
        if db == "dbsnp" and not out["dbsnp"]:
            out["dbsnp"] = rid
        elif db == "clinvar" and not out["clinvar"]:
            out["clinvar"] = rid
    return out


def _parse_natural_variant(feat: dict) -> dict | None:
    """Convert a raw UniProt Natural-variant feature into the structured
    form consumed by ``process_natural_variants``. Returns ``None`` when
    the feature lacks a usable position."""
    pos, _end = _loc(feat)
    if not isinstance(pos, int):
        return None
    alt = feat.get("alternativeSequence") or {}
    original_aa = alt.get("originalSequence")
    alt_list = alt.get("alternativeSequences") or []
    variant_aa = alt_list[0] if alt_list else None
    description = feat.get("description") or ""
    xrefs = _xref_strings(feat)
    return {
        "position": pos,
        "original_aa": original_aa,
        "variant_aa": variant_aa,
        "description": description[:200],
        "clinical_significance": _classify_natural_variant(description),
        "dbsnp": xrefs["dbsnp"],
        "clinvar": xrefs["clinvar"],
        "uniprot_var_id": feat.get("featureId"),
    }


_AA3_TO_AA1: dict[str, str] = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}


def _aa3_to_aa1(aa: str | None) -> str | None:
    """Normalise a proband alt-AA token to one-letter. Accepts the
    3-letter HGVS form ("Leu") or an already-1-letter token ("L");
    returns None for unrecognised / missing input."""
    if not aa:
        return None
    aa = aa.strip()
    if len(aa) == 1 and aa.upper() in set(_AA3_TO_AA1.values()) | {"*"}:
        return aa.upper()
    return _AA3_TO_AA1.get(aa.capitalize())


def process_natural_variants(
    variants: list[dict],
    protein_position: int | None,
    domains: list[dict] | None = None,
    window: int = 15,
    proband_variant_aa: str | None = None,
) -> dict:
    """Slice the full natural-variant list into the position-aware views
    each ACMG criterion needs (PS1 / PM5 / PM1 / BP1).

    Inputs:
      - ``variants``: full structured list from
        ``_parse_natural_variant`` (one dict per UniProt entry).
      - ``protein_position``: the residue index of the curated variant
        (``None`` for non-missense / VEP failures — in that case the
        position-relative slices come back empty but the aggregate
        stats are still computed).
      - ``domains``: the protein's domain list as returned by
        ``fetch_uniprot`` (``[{name, start, end}, ...]``).
      - ``window``: residues either side of ``protein_position`` to
        count as "nearby" — defaults to 15, matching the working
        definition of a hotspot used by ClinGen SVI in the PM1 calibration.

    Output dict (all fields always present):
      - ``nearby_variants`` — variants within ±window of
        ``protein_position`` (sorted by position).
      - ``position_exact_match`` — variants at exactly
        ``protein_position`` (subset of nearby).
      - ``same_aa_match`` — exact-residue variants whose alt AA equals
        the proband's (``proband_variant_aa``): the ACMG-2015 PS1 surface
        ("same amino acid change as a previously established variant").
        Empty when ``proband_variant_aa`` is not supplied.
      - ``different_aa_match`` — exact-residue variants whose alt AA
        DIFFERS from the proband's: the PM5 surface ("different missense
        change at the same residue"). Falls back to the full exact list
        when ``proband_variant_aa`` is unknown (can't tell same from
        different, so don't claim any are same-AA PS1 support).
      - ``domain_variants`` — list of ``{domain, total, pathogenic,
        density_per_residue}`` for every domain that contains
        ``protein_position``.
      - ``benign_missense_stats`` — aggregate counts for BP1:
        ``{total, benign, pathogenic, tolerated, vus, unknown,
        benign_fraction}``.
      - ``full_variant_count`` — total natural-variant entries in
        the protein (so Claude knows the full denominator).
    """
    total = len(variants or [])

    by_class: dict[str, int] = {
        "pathogenic": 0,
        "benign": 0,
        "vus": 0,
        "tolerated": 0,
        "unknown": 0,
    }
    for v in variants or []:
        cls = v.get("clinical_significance") or "unknown"
        if cls not in by_class:
            by_class[cls] = 0
        by_class[cls] += 1
    benign_fraction = (by_class["benign"] / total) if total else 0.0

    nearby: list[dict] = []
    exact: list[dict] = []
    if isinstance(protein_position, int) and protein_position > 0:
        lo = protein_position - window
        hi = protein_position + window
        for v in variants or []:
            pos = v.get("position")
            if not isinstance(pos, int):
                continue
            if lo <= pos <= hi:
                nearby.append(v)
                if pos == protein_position:
                    exact.append(v)
        nearby.sort(key=lambda v: v.get("position", 0))
        exact.sort(key=lambda v: (v.get("variant_aa") or ""))

    proband_aa1 = _aa3_to_aa1(proband_variant_aa)
    same_aa_match: list[dict] = []
    different_aa_match: list[dict] = []
    if proband_aa1 is not None:
        for v in exact:
            v_aa1 = _aa3_to_aa1(v.get("variant_aa"))
            if v_aa1 is not None and v_aa1 == proband_aa1:
                same_aa_match.append(v)
            else:
                different_aa_match.append(v)
    else:
        different_aa_match = list(exact)

    domain_variants: list[dict] = []
    if isinstance(protein_position, int) and protein_position > 0:
        for d in domains or []:
            ds, de = d.get("start"), d.get("end")
            if not (isinstance(ds, int) and isinstance(de, int)):
                continue
            if not (ds <= protein_position <= de):
                continue
            in_domain = [
                v for v in (variants or [])
                if isinstance(v.get("position"), int) and ds <= v["position"] <= de
            ]
            domain_total = len(in_domain)
            domain_path = sum(
                1 for v in in_domain
                if v.get("clinical_significance") == "pathogenic"
            )
            span = max(1, de - ds + 1)
            domain_variants.append({
                "domain": d.get("name") or "(unnamed)",
                "start": ds,
                "end": de,
                "total": domain_total,
                "pathogenic": domain_path,
                "density_per_residue": round(domain_path / span, 4),
            })

    return {
        "nearby_variants": nearby,
        "position_exact_match": exact,
        "same_aa_match": same_aa_match,
        "different_aa_match": different_aa_match,
        "proband_variant_aa": proband_variant_aa,
        "domain_variants": domain_variants,
        "benign_missense_stats": {
            "total": total,
            "benign": by_class["benign"],
            "pathogenic": by_class["pathogenic"],
            "tolerated": by_class["tolerated"],
            "vus": by_class["vus"],
            "unknown": by_class["unknown"],
            "benign_fraction": round(benign_fraction, 4),
        },
        "full_variant_count": total,
        "protein_position": protein_position,
        "window": window,
    }


UNIPROT_FIELDS = (
    "accession", "id", "protein_name", "gene_names", "length",
    "sequence",
    "ft_domain", "ft_region", "ft_motif", "ft_repeat",
    "ft_zn_fing", "ft_dna_bind",
    "ft_act_site", "ft_binding",
    "ft_transmem", "ft_signal", "ft_coiled",
    "ft_variant",
)


def parse_entry(entry: dict) -> dict:
    """Parse a single UniProtKB JSON entry (as returned by either the
    ``/uniprotkb/search`` endpoint OR the ``/uniprotkb/stream`` bulk
    download — the per-entry JSON shape is identical) into the structured
    return dict ``fetch_uniprot`` exposes.

    The returned dict deliberately OMITS the ``"gene"`` key: the live
    search echoes back the *requested* symbol (which may be a synonym),
    so the caller injects ``"gene"`` after lookup. This keeps DB-backed
    and live results byte-identical, including for synonym queries.

    Build script (``scripts/build_uniprot_db.py``) calls this on every
    streamed entry so the stored payload is exactly what the live client
    would have produced.
    """
    accession = entry.get("primaryAccession")
    protein_name = (
        ((entry.get("proteinDescription") or {}).get("recommendedName") or {})
        .get("fullName", {})
        .get("value")
    )
    seq_len = (entry.get("sequence") or {}).get("length")
    seq_value = (entry.get("sequence") or {}).get("value")

    feats = entry.get("features", []) or []
    domains, active_sites, binding_sites = [], [], []
    natural_variants: list[dict] = []
    track_features: list[dict] = []
    for f in feats:
        t = f.get("type")
        s, e = _loc(f)
        desc = f.get("description") or ""
        if t == "Domain":
            domains.append({"name": desc or "(unnamed)", "start": s, "end": e})
        elif t == "Active site":
            active_sites.append({"position": s, "description": desc})
        elif t == "Binding site":
            binding_sites.append({"start": s, "end": e, "description": desc})
        elif t == "Natural variant":
            parsed = _parse_natural_variant(f)
            if parsed is not None:
                natural_variants.append(parsed)
        if t in _TRACK_FEATURE_TYPES and isinstance(s, int):
            track_features.append({
                "type": t,
                "start": s,
                "end": e if isinstance(e, int) else s,
                "description": desc,
            })

    return {
        "ok": True,
        "found": True,
        "accession": accession,
        "uniprot_id": entry.get("uniProtkbId"),
        "protein_name": protein_name,
        "length": seq_len,
        "sequence": seq_value,
        "url": f"https://www.uniprot.org/uniprotkb/{accession}" if accession else None,
        "domains": domains,
        "active_sites": active_sites,
        "binding_sites": binding_sites,
        "features": track_features,
        "natural_variant_count": len(natural_variants),
        "natural_variants": natural_variants,
    }


def gene_symbols(entry: dict) -> tuple[str | None, list[str]]:
    """Extract the (primary, [synonyms...]) gene symbols from a UniProtKB
    entry's ``genes`` array. Used by the build script to index an entry
    under its primary symbol AND every listed synonym, mirroring the live
    client's ``gene_exact`` matching (which UniProt resolves against both
    primary and synonym symbols)."""
    primary: str | None = None
    synonyms: list[str] = []
    for g in entry.get("genes") or []:
        name = (g.get("geneName") or {}).get("value")
        if name and primary is None:
            primary = name
        elif name:
            synonyms.append(name)
        for syn in g.get("synonyms") or []:
            v = syn.get("value")
            if v:
                synonyms.append(v)
    return primary, synonyms


def _query_local_sync(gene: str) -> dict | None:
    """Look the gene up in the local UniProt DB by UPPER(symbol).

    Returns the stored payload (with ``"gene"`` injected) on a hit, or
    ``None`` when the DB file is absent OR the symbol isn't present — both
    of which signal the caller to fall back to the live HTTP path.
    """
    db_path = _db_path()
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        log.warning("UniProt local DB open failed (%s); falling back to live", exc)
        return None
    try:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            "SELECT payload FROM uniprot_entry WHERE gene_symbol = ? LIMIT 1",
            (gene.upper(),),
        )
        row = cur.fetchone()
    except sqlite3.Error as exc:
        log.warning("UniProt local DB query failed for %s (%s); falling back to live", gene, exc)
        return None
    finally:
        conn.close()
    if row is None:
        return None
    try:
        payload = json.loads(row["payload"])
    except (ValueError, TypeError) as exc:
        log.warning("UniProt local DB payload decode failed for %s (%s); falling back to live", gene, exc)
        return None
    payload["gene"] = gene
    return payload


async def _fetch_uniprot_live(gene: str) -> dict:
    """Live ``/uniprotkb/search`` lookup — the graceful fallback used when
    the local DB is absent or doesn't contain the gene. This is the
    original ``fetch_uniprot`` implementation, preserved verbatim.

    Returns domains, active sites, binding sites, and the **full**
    natural-variant annotation list so callers can downstream-slice
    by curated-variant position via ``process_natural_variants``.
    Also returns a flat ``features`` list that powers the frontend
    Protein-tab track.
    """
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "fields": ",".join(UNIPROT_FIELDS),
        "size": "1",
    }
    async with make_async_client() as c:
        r = await request_with_retry(c, "GET", f"{UNIPROT_API}/uniprotkb/search", params=params, timeout=20.0, name="UniProt")
        if r is None:
            return {"ok": False, "error": "UniProt transient failure after retries"}
        if r.status_code != 200:
            return {"ok": False, "error": f"{r.status_code}: {r.text[:200]}"}
        results = r.json().get("results", []) or []
        if not results:
            return {"ok": True, "found": False, "gene": gene}

        payload = parse_entry(results[0])
        payload["gene"] = gene
        return payload


async def fetch_uniprot(gene: str) -> dict:
    """Look up the canonical human SwissProt entry for a gene symbol.

    Prefers the local SQLite cache built by
    ``scripts/build_uniprot_db.py`` and falls back to the live
    ``/uniprotkb/search`` call when the DB file is absent OR the gene is
    not found locally. The return shape is identical on both paths:
    domains, active sites, binding sites, the **full** natural-variant
    annotation list (for ``process_natural_variants``) and a flat
    ``features`` list for the Protein-tab track.
    """
    local = await run_local(_query_local_sync, gene)
    if local is not None:
        return local
    if offline_strict():
        return {"ok": True, "found": False, "gene": gene}
    return await _fetch_uniprot_live(gene)
