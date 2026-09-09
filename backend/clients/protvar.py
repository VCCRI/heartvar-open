"""ProtVar — functional, conservation, co-located, and population annotations
for missense variants.

LOCAL-FIRST (2026 rewrite): HeartVar fans out ~16 live lookups per curation
and is deploying publicly, so this source is served from the LOCAL caches the
project already ships rather than hammering EBI's ProtVar REST API on every
call:

  * ``am_pathogenicity`` / ``am_class`` — from the local AlphaMissense tabix
    table (``backend.clients.alphamissense.fetch_alphamissense``).
  * ``gene`` / ``protein_position`` / ``ref_aa`` / ``alt_aa`` /
    ``consequence`` — derived from the VEP result already in hand (gene
    symbol, ``hgvsp``, ``most_severe_consequence``).
  * ``uniprot`` (accession) / ``feature_types`` / ``function_summary`` /
    ``colocated_variants`` — from the local UniProt entry for the gene
    (``backend.clients.uniprot.fetch_uniprot`` +
    ``process_natural_variants``).
  * ``conservation_score`` / ``pocket_score`` — ``None`` (DROPPED per the
    owner's decision: no free local source; ProtVar's conservation/FoldX
    pockets came from EBI sub-endpoints we no longer call).

The original EBI implementation is preserved verbatim as
``_fetch_protvar_live`` and is used ONLY when BOTH local sources are
unavailable (AlphaMissense not configured AND no local UniProt entry for the
gene). When either local source resolves, EBI is never touched.

This client is only invoked when the VEP `most_severe_consequence` is
`missense_variant` — other variant classes render as
"Not applicable (non-missense variant)" in the frontend.

The output is CRITERIA-CRITICAL (feeds PP3 / PM1 SUPPORTING context only; the
prompt forbids it from satisfying PS3). On any ambiguity in the local path we
fail safe — a missing local annotation surfaces as ``None`` / ``[]`` rather
than a guessed value, and only a TOTAL local miss escalates to the live
fallback.
"""

from __future__ import annotations

import re
from typing import Any

import httpx

from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict
from .alphamissense import fetch_alphamissense
from .uniprot import fetch_uniprot, process_natural_variants

PROTVAR_API_BASE = "https://www.ebi.ac.uk/ProtVar/api"
PROTVAR_MAPPING = f"{PROTVAR_API_BASE}/mapping"

_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*", "Stop": "*",
}
_ONE_TO_THREE = {v: k for k, v in _THREE_TO_ONE.items() if len(v) == 1}


def _missense(vep: dict[str, Any] | None) -> bool:
    if not isinstance(vep, dict) or not vep.get("ok"):
        return False
    return (vep.get("most_severe_consequence") or "").lower() == "missense_variant"


def _variant_aa_three_upper(variant_aa: Any) -> str | None:
    """Normalise the isoform variantAA to the 3-letter UPPER form /function wants.

    ProtVar may report the variant amino acid as a 3-letter mixed-case code
    ("Tyr") or, defensively, a single-letter code ("Y"). Either is converted to
    "TYR". Returns None if it cannot be mapped.
    """
    if not variant_aa or not isinstance(variant_aa, str):
        return None
    aa = variant_aa.strip()
    if len(aa) == 1:
        three = _ONE_TO_THREE.get(aa.upper())
        return three.upper() if three else None
    title = aa[:1].upper() + aa[1:].lower()
    if title in _THREE_TO_ONE:
        return title.upper()
    return aa.upper()


_HGVSP_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


def _parse_hgvsp(hgvsp: Any) -> tuple[str | None, int | None, str | None]:
    """Parse ``(ref_aa, protein_position, alt_aa)`` from a VEP ``hgvsp``.

    Returns 3-letter ref/alt AA tokens (matching the ProtVar contract, which
    surfaced 3-letter forms) and the integer residue position. Any field that
    can't be parsed comes back ``None`` rather than a guess.
    """
    if not hgvsp or not isinstance(hgvsp, str):
        return None, None, None
    m = _HGVSP_RE.search(hgvsp.split(":")[-1])
    if not m:
        return None, None, None
    ref_aa = m.group(1)
    try:
        pos = int(m.group(2))
    except (TypeError, ValueError):
        pos = None
    alt_aa = m.group(3)
    return ref_aa, pos, alt_aa


def _feature_types_at(uni: dict, protein_position: int | None) -> list[str]:
    """Feature types from the local UniProt entry that overlap the residue.

    ProtVar's original ``feature_types`` listed the functional features the
    variant's position falls inside (domains, motifs, active sites, etc.). The
    local UniProt entry exposes a flat ``features`` track of
    ``{type, start, end, description}``; we keep the distinct types whose span
    contains ``protein_position``. When the position is unknown we fall back to
    the protein's overall feature-type vocabulary so the panel still has
    context. Deduplicated, sorted, capped at 6 (matching the live contract).
    """
    feats = uni.get("features") or []
    names: list[str] = []
    if isinstance(protein_position, int) and protein_position > 0:
        for ft in feats:
            s, e = ft.get("start"), ft.get("end")
            if not isinstance(s, int):
                continue
            end = e if isinstance(e, int) else s
            if s <= protein_position <= end:
                t = ft.get("type")
                if t:
                    names.append(t)
    else:
        for ft in feats:
            t = ft.get("type")
            if t:
                names.append(t)
    return sorted(set(names))[:6]


_NV_SIG_LABEL = {
    "pathogenic": "pathogenic",
    "benign": "benign",
    "vus": "uncertain significance",
    "tolerated": "tolerated",
    "unknown": None,
}


def _colocated_from_uniprot(proc: dict) -> list[dict]:
    """Map UniProt natural variants at the curated residue to the
    colocated-variant summary shape ``[{variant, clinical_significance,
    source_db}]``.

    Uses ``process_natural_variants``'s ``position_exact_match`` slice — the
    other natural-variant entries reported at the SAME residue. Best-effort:
    returns ``[]`` when the slice is empty or unavailable.
    """
    out: list[dict] = []
    exact = (proc or {}).get("position_exact_match") or []
    for v in exact[:8]:
        orig = v.get("original_aa")
        alt = v.get("variant_aa")
        if orig and alt:
            name = f"{orig}->{alt}"
        else:
            name = v.get("uniprot_var_id") or v.get("dbsnp") or v.get("clinvar")
        sig = _NV_SIG_LABEL.get(v.get("clinical_significance") or "unknown")
        out.append({
            "variant": name,
            "clinical_significance": sig,
            "source_db": "UniProt",
        })
    return out


async def _fetch_function(client, accession: str, position: Any,
                          variant_aa_3u: str | None) -> dict:
    """Fetch /function/{accession}/{position}. Non-fatal: returns {} on failure."""
    url = f"{PROTVAR_API_BASE}/function/{accession}/{position}"
    params = {"variantAA": variant_aa_3u} if variant_aa_3u else None
    try:
        r = await request_with_retry(
            client, "GET", url, params=params, timeout=30.0, name="ProtVar function",
        )
    except httpx.HTTPError:
        return {}
    if r is None or r.status_code != 200:
        return {}
    try:
        return r.json() or {}
    except ValueError:
        return {}


async def _fetch_population(client, accession: str, position: Any,
                            genomic_variant: str) -> dict:
    """Fetch /population/{accession}/{position}. Non-fatal: returns {} on failure."""
    url = f"{PROTVAR_API_BASE}/population/{accession}/{position}"
    try:
        r = await request_with_retry(
            client, "GET", url, params={"genomicVariant": genomic_variant},
            timeout=30.0, name="ProtVar population",
        )
    except httpx.HTTPError:
        return {}
    if r is None or r.status_code != 200:
        return {}
    try:
        return r.json() or {}
    except ValueError:
        return {}


def _function_summary(func: dict) -> str | None:
    """Extract the FUNCTION comment's first text value from /function `comments`."""
    for c in func.get("comments") or []:
        if (c.get("type") or "").upper() == "FUNCTION":
            for t in c.get("text") or []:
                val = t.get("value")
                if val:
                    return val[:600]
    return None


def _colocated_summary(pop: dict) -> list[dict]:
    """Map /population `variants[]` to the colocated-variant summary shape."""
    out: list[dict] = []
    for v in (pop.get("variants") or [])[:8]:
        sigs = v.get("clinicalSignificances") or []
        sig_str = ", ".join(
            s.get("type") or "" for s in sigs if isinstance(s, dict) and s.get("type")
        ) or None
        wt = v.get("wildType")
        alt = v.get("alternativeSequence")
        if wt and alt:
            name = f"{wt}->{alt}"
        else:
            name = v.get("ftId")
        if not name:
            for x in v.get("xrefs") or []:
                if x.get("id"):
                    name = x["id"]
                    break
        out.append({
            "variant": name or v.get("ftId"),
            "clinical_significance": sig_str,
            "source_db": v.get("sourceType"),
        })
    return out


def _input_str_and_url(vep: dict[str, Any]) -> tuple[str | None, str, str, str, Any, Any]:
    """Shared coordinate/strand-complement derivation used by BOTH the local
    path (for ``input`` / ``url``) and the live fallback (for the /mapping
    query). Returns ``(input_str, chrom, ref, alt, pos, error)`` where
    ``error`` is None on success or a dict to return verbatim on failure.

    Kept BYTE-IDENTICAL to the original derivation: VEP reports
    ``allele_string`` in TRANSCRIPT orientation for HGVS-c input, so on a
    minus-strand gene the bases are reverse-complemented relative to the
    forward genome before building the input string / url.
    """
    chrom = vep.get("seq_region_name")
    pos = vep.get("start")
    allele = (vep.get("allele_string") or "").split("/")
    if not chrom or not pos or len(allele) != 2:
        return None, "", "", "", None, {
            "ok": False, "error": "missing chromosomal coordinates in VEP data"
        }
    ref, alt = allele[0], allele[1]
    if vep.get("strand") == -1:
        _comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        ref = "".join(_comp.get(b, b) for b in ref.upper())
        alt = "".join(_comp.get(b, b) for b in alt.upper())
    input_str = f"{chrom} {pos} {ref} {alt}"
    return input_str, chrom, ref, alt, pos, None


async def _fetch_protvar_live(vep: dict[str, Any]) -> dict:
    """Live EBI ProtVar lookup — the graceful fallback used ONLY when both
    local sources (AlphaMissense + local UniProt entry) are unavailable.

    This is the original ``fetch_protvar`` implementation, preserved verbatim.
    The caller has already gated on missense + VEP-ok, so this assumes a
    usable VEP result and derives coordinates the same way.
    """
    input_str, chrom, ref, alt, pos, err = _input_str_and_url(vep)
    if err is not None:
        return err

    async with make_async_client() as c:
        try:
            r = await request_with_retry(
                c, "GET", PROTVAR_MAPPING,
                params={"q": input_str, "assembly": "AUTO"},
                timeout=30.0, name="ProtVar",
            )
        except httpx.HTTPError as e:
            return {"ok": False, "error": f"transport error: {e!r}"}
        if r is None:
            return {"ok": False, "error": "ProtVar transient failure after retries"}
        if r.status_code != 200:
            return {"ok": False, "error": f"{r.status_code}: {r.text[:200]}"}
        try:
            payload = r.json()
        except ValueError:
            return {"ok": False, "error": "non-JSON response from ProtVar"}

        content = payload.get("content") or {}
        inputs = content.get("inputs") or []
        if not inputs:
            return {"ok": True, "applicable": True, "found": False, "input": input_str}
        first = inputs[0]

        derived = first.get("derivedGenomicVariants") or []
        dgv = derived[0] if derived else {}
        genes = dgv.get("genes") or []
        if not genes:
            return {"ok": True, "applicable": True, "found": False, "input": input_str}

        g0 = genes[0]
        canonical = None
        for iso in g0.get("isoforms") or []:
            if iso.get("canonical"):
                canonical = iso
                break
        if canonical is None and g0.get("isoforms"):
            canonical = g0["isoforms"][0]
        canonical = canonical or {}

        accession = canonical.get("accession")
        protein_position = canonical.get("isoformPosition")
        variant_aa = canonical.get("variantAA")


        gv_chrom = dgv.get("chromosome") or chrom
        gv_pos = dgv.get("position") or pos
        gv_ref = dgv.get("refBase") or g0.get("refAllele") or ref
        gv_alt = dgv.get("altBase") or g0.get("altAllele") or alt
        genomic_variant = f"{gv_chrom}-{gv_pos}-{gv_ref}-{gv_alt}"

        func: dict = {}
        pop: dict = {}
        if accession and protein_position is not None:
            variant_aa_3u = _variant_aa_three_upper(variant_aa)
            func = await _fetch_function(c, accession, protein_position, variant_aa_3u)
            pop = await _fetch_population(c, accession, protein_position, genomic_variant)

    conserv_raw = func.get("conservScore")
    conservation = conserv_raw.get("score") if isinstance(conserv_raw, dict) else conserv_raw

    feature_types: list[str] = []
    for ft in func.get("features") or []:
        name = ft.get("type") or ft.get("category") or ft.get("description") or ""
        if name:
            feature_types.append(name)
    feature_types = sorted(set(feature_types))[:6]

    function_summary = _function_summary(func)

    pockets = func.get("pockets")
    pocket_score = None
    if isinstance(pockets, dict):
        pocket_score = pockets.get("pocketScore") or pockets.get("score")
    elif isinstance(pockets, list) and pockets and isinstance(pockets[0], dict):
        pocket_score = pockets[0].get("pocketScore") or pockets[0].get("score")

    colocated_summary = _colocated_summary(pop)
    colocated_total = len(pop.get("variants") or [])

    return {
        "ok": True,
        "applicable": True,
        "found": True,
        "input": input_str,
        "gene": g0.get("geneName"),
        "uniprot": accession,
        "protein_position": protein_position,
        "ref_aa": canonical.get("refAA"),
        "alt_aa": variant_aa,
        "consequence": canonical.get("consequences"),
        "conservation_score": conservation,
        "feature_types": feature_types,
        "function_summary": function_summary,
        "pocket_score": pocket_score,
        "colocated_variants": colocated_summary,
        "colocated_total": colocated_total,
        "am_pathogenicity": None,
        "am_class": None,
        "url": (
            f"https://www.ebi.ac.uk/ProtVar/query?search="
            f"{chrom}%20{pos}%20{ref}%20{alt}"
        ),
    }


async def fetch_protvar(vep: dict[str, Any] | None) -> dict:
    """Missense functional / population annotation for a curated variant.

    LOCAL-FIRST: AlphaMissense (local tabix) supplies ``am_pathogenicity`` /
    ``am_class``; the local UniProt entry supplies ``uniprot`` /
    ``feature_types`` / ``function_summary`` / ``colocated_variants``; the VEP
    result in hand supplies ``gene`` / ``protein_position`` / ``ref_aa`` /
    ``alt_aa`` / ``consequence``. ``conservation_score`` and ``pocket_score``
    are ``None`` (dropped — no free local source). The live EBI lookup
    (``_fetch_protvar_live``) runs ONLY when BOTH local sources are
    unavailable. Never raises.

    The return-dict shape/keys are byte-identical to the original EBI client.
    """
    if not isinstance(vep, dict) or not vep.get("ok"):
        return {"ok": False, "error": "VEP lookup unavailable — cannot derive coordinates"}
    if not _missense(vep):
        return {
            "ok": True,
            "applicable": False,
            "reason": "non-missense variant",
            "consequence": vep.get("most_severe_consequence"),
        }

    input_str, chrom, ref, alt, pos, err = _input_str_and_url(vep)
    if err is not None:
        return err

    am = await fetch_alphamissense(vep)
    am_available = bool(am.get("available"))
    am_pathogenicity = am.get("score") if am_available else None
    am_class = am.get("classification") if am_available else None

    gene = vep.get("gene_symbol") or vep.get("derived_gene_symbol")
    uni: dict = {}
    if gene:
        uni = await fetch_uniprot(gene)
    uni_found = bool(uni.get("ok") and uni.get("found"))

    if not am_available and not uni_found and not offline_strict():
        return await _fetch_protvar_live(vep)

    ref_aa, protein_position, alt_aa = _parse_hgvsp(vep.get("hgvsp"))
    consequence = vep.get("most_severe_consequence")
    accession = uni.get("accession") if uni_found else None

    feature_types = _feature_types_at(uni, protein_position) if uni_found else []

    function_summary = None

    colocated_variants: list[dict] = []
    colocated_total = 0
    if uni_found:
        proc = process_natural_variants(
            uni.get("natural_variants") or [],
            protein_position,
            domains=uni.get("domains"),
            proband_variant_aa=alt_aa,
        )
        colocated_variants = _colocated_from_uniprot(proc)
        colocated_total = len(proc.get("position_exact_match") or [])

    found = bool(accession or (am_pathogenicity is not None))

    return {
        "ok": True,
        "applicable": True,
        "found": found,
        "input": input_str,
        "gene": gene,
        "uniprot": accession,
        "protein_position": protein_position,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "consequence": consequence,
        "conservation_score": None,
        "feature_types": feature_types,
        "function_summary": function_summary,
        "pocket_score": None,
        "colocated_variants": colocated_variants,
        "colocated_total": colocated_total,
        "am_pathogenicity": am_pathogenicity,
        "am_class": am_class,
        "url": (
            f"https://www.ebi.ac.uk/ProtVar/query?search="
            f"{chrom}%20{pos}%20{ref}%20{alt}"
        ),
    }
