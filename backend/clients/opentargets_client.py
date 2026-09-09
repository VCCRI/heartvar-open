"""Open Targets Platform — gene-disease association evidence.

Public GraphQL API at https://api.platform.opentargets.org/api/v4/graphql,
no authentication required. We use it to surface a quantitative measure
of how well-supported the gene-disease association is — overall
association score (0-1) plus a breakdown across the eight Open Targets
datatypes (genetic_association, genetic_literature, somatic_mutation,
affected_pathway, animal_model, literature, rna_expression, clinical).

Keyed on an Ensembl gene id, which the caller does not reliably have:
VEP returns an **EntrezGene** id for RefSeq-transcript input (see
``_gene_ids``), so the id is resolved here — from the supplied value or
the gene symbol — before anything else. Without that, every NM_-input
curation reports "unavailable" while the API is perfectly healthy.

The disease side uses EFO IDs. Open Targets does not accept HPO IDs
directly on the association endpoint, so the client maps a small set of
CHD/cardiomyopathy HPO terms to their EFO equivalents locally. When the
proband HPO has no mapped EFO match (or no HPO was supplied), the client
falls back to the top-N associated diseases for the gene, with no
filtering — that path lets the curator see the gene's strongest
associations even when the proband's phenotype isn't in the local map.
"""

from __future__ import annotations

import json
import logging
import sqlite3
from pathlib import Path

import httpx

from ..localio import run_local
from ._gene_ids import resolve_ensembl_gene_id
from ._http_retry import make_async_client, request_with_retry
from ._offline import offline_strict
from ._paths import db_path

log = logging.getLogger("heartvar.opentargets")

OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"


def _opentargets_db_path() -> Path:
    """Resolve the local snapshot path at call time so an env override / test
    can point it elsewhere (mirrors the gtex/uniprot db_path pattern)."""
    return db_path("opentargets.db", "OPENTARGETS_DB_PATH")


def _query_local_sync(ensembl_id: str) -> dict | None:
    """Read the offline Open Targets snapshot for an Ensembl gene id, returning
    the SAME ``{"target": {...}}`` shape the live GraphQL ``data`` field yields
    (approvedSymbol + associatedDiseases.rows), or None when the snapshot DB is
    absent or the gene isn't in it. Never raises — a bad snapshot must not break
    a curation (the caller falls back to live unless offline-strict)."""
    p = _opentargets_db_path()
    try:
        if not p.exists():
            return None
        con = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
        try:
            row = con.execute(
                "SELECT payload FROM opentargets WHERE ensembl_id = ?",
                (ensembl_id,),
            ).fetchone()
        finally:
            con.close()
    except Exception:  # noqa: BLE001 — a bad/locked snapshot must not break curation
        return None
    if not row or not row[0]:
        return None
    try:
        return {"target": json.loads(row[0])}
    except (ValueError, TypeError):
        return None

_HPO_TO_DISEASE: dict[str, dict] = {
    "HP:0001629": {
        "efo": "EFO_0000544",
        "name_keywords": ["ventricular septal defect"],
        "alt_ids": ["MONDO_0013746"],
    },
    "HP:0001647": {
        "efo": "EFO_0004274",
        "name_keywords": ["bicuspid aortic valve"],
        "alt_ids": [],
    },
    "HP:0001631": {
        "efo": "EFO_0009820",
        "name_keywords": ["atrial septal defect"],
        "alt_ids": ["MONDO_0011938"],
    },
    "HP:0001674": {
        "efo": "EFO_0003900",
        "name_keywords": ["tetralogy of fallot"],
        "alt_ids": ["HP_0001636"],
    },
    "HP:0001680": {
        "efo": "EFO_0004713",
        "name_keywords": ["coarctation of aorta", "coarctation of the aorta"],
        "alt_ids": [],
    },
    "HP:0001644": {
        "efo": "EFO_0000400",
        "name_keywords": ["dilated cardiomyopathy"],
        "alt_ids": ["EFO_0000407"],
    },
    "HP:0001639": {
        "efo": "EFO_0000408",
        "name_keywords": ["hypertrophic cardiomyopathy"],
        "alt_ids": ["EFO_0000538"],
    },
}

_DATATYPES: tuple[str, ...] = (
    "genetic_association",
    "genetic_literature",
    "somatic_mutation",
    "affected_pathway",
    "animal_model",
    "literature",
    "rna_expression",
    "clinical",
)

_QUERY_TOP_DISEASES = """
query TopDiseases($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(page: {size: 50, index: 0}) {
      rows {
        disease { id name }
        score
        datatypeScores { id score }
      }
    }
  }
}
""".strip()


def _normalise_hpo(token: str) -> str | None:
    """Trim whitespace and uppercase the HP prefix so "hp:0001639" and
    " HP:0001639 " both match the mapping dict."""
    if not token:
        return None
    t = token.strip()
    if not t:
        return None
    if t[:3].lower() == "hp:":
        return "HP:" + t[3:]
    return t


def _parse_hpo_input(hpo_terms) -> list[str]:
    """Accept either a list of HPO tokens or the free-text ``req.hpo`` string
    the frontend sends ("HP:0001639,HP:0001644" or "HP:0001639 HP:0001644")
    and return a deduplicated list of canonicalised tokens."""
    if hpo_terms is None:
        return []
    if isinstance(hpo_terms, str):
        raw = [p for p in hpo_terms.replace(",", " ").split() if p]
    else:
        try:
            raw = list(hpo_terms)
        except TypeError:
            return []
    out: list[str] = []
    seen: set[str] = set()
    for r in raw:
        n = _normalise_hpo(r)
        if n and n not in seen:
            seen.add(n)
            out.append(n)
    return out


def _datatype_scores_dict(rows: list[dict]) -> dict[str, float]:
    """Flatten Open Targets' [{id, score}, ...] datatypeScores array into
    a {datatype_id: score} dict. When the same datatype appears in
    multiple association rows (fallback path with multiple diseases),
    the maximum score wins so the breakdown reflects the strongest
    evidence the gene has across the surfaced diseases."""
    out: dict[str, float] = {}
    for r in rows or []:
        for ds in r.get("datatypeScores") or []:
            dt_id = ds.get("id")
            score = ds.get("score")
            if not dt_id or score is None:
                continue
            try:
                fscore = float(score)
            except (TypeError, ValueError):
                continue
            prev = out.get(dt_id)
            if prev is None or fscore > prev:
                out[dt_id] = round(fscore, 3)
    for dt in _DATATYPES:
        out.setdefault(dt, 0.0)
    return out


def _build_summary(
    matched_disease_name: str | None,
    overall: float | None,
    datatype_scores: dict[str, float],
    top_diseases: list[dict],
) -> str:
    """One-liner the prompt can quote verbatim. Highlights the overall
    score band, the matched disease (or "top association" when in
    fallback), and the two strongest datatypes."""
    if overall is None:
        return "No Open Targets association data available for this gene."
    if overall >= 0.75:
        band = "Strong"
    elif overall >= 0.5:
        band = "Moderate"
    elif overall >= 0.25:
        band = "Emerging"
    else:
        band = "Weak"
    if matched_disease_name:
        anchor = f"{band.lower()} association with {matched_disease_name}"
    elif top_diseases:
        anchor = (
            f"top association {top_diseases[0]['name']} ({band.lower()})"
        )
    else:
        anchor = f"{band.lower()} associations"
    ranked = sorted(
        datatype_scores.items(), key=lambda kv: kv[1], reverse=True
    )
    top_dt = [
        f"{name.replace('_', ' ')} {score:.2f}"
        for name, score in ranked[:2]
        if score > 0
    ]
    dt_phrase = (", top datatypes: " + "; ".join(top_dt)) if top_dt else ""
    return f"{band} ({anchor}, overall score {overall:.2f}){dt_phrase}."


class OpenTargetsClient:
    """Thin async wrapper around the Open Targets GraphQL endpoint."""

    def __init__(self, endpoint: str = OPEN_TARGETS_GRAPHQL, timeout: float = 12.0):
        self.endpoint = endpoint
        self.timeout = timeout
        self.last_error: str | None = None

    async def _post(self, client: httpx.AsyncClient, query: str, variables: dict) -> dict | None:
        try:
            r = await request_with_retry(
                client, "POST", self.endpoint,
                json={"query": query, "variables": variables},
                timeout=self.timeout, name="OpenTargets", logger=log,
            )
        except httpx.HTTPError as e:
            log.warning("Open Targets POST failed: %r", e)
            self.last_error = f"Open Targets request failed: {e!r}"
            return None
        if r is None:
            self.last_error = "Open Targets unreachable after retries"
            return None
        if r.status_code != 200:
            log.warning("Open Targets non-200: %s — %s", r.status_code, r.text[:200])
            self.last_error = f"Open Targets returned HTTP {r.status_code}"
            return None
        try:
            payload = r.json()
        except ValueError:
            self.last_error = "Open Targets returned a non-JSON response"
            return None
        if "errors" in payload and payload.get("errors"):
            log.warning("Open Targets GraphQL errors: %r", payload.get("errors"))
            self.last_error = (
                f"Open Targets GraphQL error: {payload['errors']!r:.200}"
            )
            return None
        data = payload.get("data") or None
        if data is None:
            self.last_error = "Open Targets returned an empty data field"
        return data

    async def get_gene_disease_evidence(
        self,
        ensembl_id: str,
        hpo_terms=None,
        gene_symbol_hint: str | None = None,
    ) -> dict | None:
        """Look up the gene's Open Targets association profile.

        ``ensembl_id`` is whatever VEP surfaced as ``gene_id`` — an ``ENSG`` for
        Ensembl-transcript input but an **EntrezGene** id for RefSeq input. It is
        resolved to an Ensembl id first, falling back to ``gene_symbol_hint``.

        Strategy:
          0. Resolve ``ensembl_id`` to a real ``ENSG`` (offline, via HGNC).
          1. Fetch the gene's top-50 associated diseases (no filter).
          2. For each proband HPO term that maps to a known disease,
             scan the top-N rows for a match by disease ID (EFO,
             alt IDs) or by canonical-name substring.
          3. On match, the matched row provides the
             ``overall_association_score`` + ``datatype_scores`` for
             the phenotype-anchored interpretation.
          4. On no match (no HPO, no mapping, or HPO disease not in
             top-50), fall back to top-5 diseases for the gene and
             leave ``matched_disease`` null.

        Returns None on any resolution/network/parse failure so the caller
        can surface the reason recorded in ``last_error`` rather than
        tripping the whole evidence-gather phase.
        """
        resolved = resolve_ensembl_gene_id(ensembl_id, gene_symbol_hint)
        if not resolved:
            supplied = ensembl_id or "none"
            self.last_error = (
                f"No Ensembl gene ID resolved for gene_id={supplied!r} / "
                f"symbol={gene_symbol_hint or 'none'} — cannot query Open Targets"
            )
            log.warning("%s", self.last_error)
            return None
        ensembl_id = resolved

        data = await run_local(_query_local_sync, ensembl_id)
        if data is None:
            if offline_strict():
                self.last_error = (
                    f"{ensembl_id} not in the local Open Targets snapshot and "
                    "offline-strict mode forbids the live fallback — rebuild "
                    "data/opentargets.db (scripts/build_opentargets_db.py)"
                )
                log.warning("%s", self.last_error)
                return None
            async with make_async_client() as client:
                data = await self._post(
                    client,
                    _QUERY_TOP_DISEASES,
                    {"ensemblId": ensembl_id},
                )
        if not data or not data.get("target"):
            if not self.last_error:
                self.last_error = (
                    f"Open Targets has no target record for {ensembl_id}"
                )
            return None

        target = data["target"]
        gene_symbol = target.get("approvedSymbol")
        rows = (target.get("associatedDiseases") or {}).get("rows") or []
        url = f"https://platform.opentargets.org/target/{ensembl_id}"

        if not rows:
            return {
                "ok": True,
                "gene_symbol": gene_symbol,
                "ensembl_id": ensembl_id,
                "matched_hpo": None,
                "matched_disease": None,
                "overall_association_score": None,
                "datatype_scores": _datatype_scores_dict([]),
                "top_diseases": [],
                "evidence_summary": (
                    "No Open Targets associations indexed for this gene."
                ),
                "url": url,
            }

        hpo_list = _parse_hpo_input(hpo_terms)
        matched_row = None
        matched_hp: str | None = None
        for hp in hpo_list:
            entry = _HPO_TO_DISEASE.get(hp)
            if not entry:
                continue
            wanted_ids = {entry.get("efo")} | set(entry.get("alt_ids") or [])
            wanted_ids = {i for i in wanted_ids if i}
            name_keywords = [k.lower() for k in (entry.get("name_keywords") or []) if k]
            for r in rows:
                d = r.get("disease") or {}
                d_id = d.get("id") or ""
                d_name = (d.get("name") or "").lower()
                if d_id in wanted_ids:
                    matched_row = r
                    matched_hp = hp
                    break
                if name_keywords and any(kw in d_name for kw in name_keywords):
                    matched_row = r
                    matched_hp = hp
                    break
            if matched_row is not None:
                break

        if matched_row is not None:
            disease = matched_row.get("disease") or {}
            overall = matched_row.get("score")
            datatype_scores = _datatype_scores_dict([matched_row])
            summary = _build_summary(
                disease.get("name"), overall, datatype_scores, []
            )
            return {
                "ok": True,
                "gene_symbol": gene_symbol,
                "ensembl_id": ensembl_id,
                "matched_hpo": matched_hp,
                "matched_disease": {
                    "id": disease.get("id"),
                    "name": disease.get("name"),
                },
                "overall_association_score": (
                    round(float(overall), 3) if overall is not None else None
                ),
                "datatype_scores": datatype_scores,
                "top_diseases": [],
                "evidence_summary": summary,
                "url": url,
            }

        top_diseases = []
        for r in rows[:5]:
            d = r.get("disease") or {}
            score = r.get("score")
            top_diseases.append({
                "id": d.get("id"),
                "name": d.get("name"),
                "score": (
                    round(float(score), 3) if score is not None else None
                ),
            })
        top_overall = rows[0].get("score")
        datatype_scores = _datatype_scores_dict(rows[:5])
        summary = _build_summary(
            None,
            round(float(top_overall), 3) if top_overall is not None else None,
            datatype_scores,
            top_diseases,
        )
        return {
            "ok": True,
            "gene_symbol": gene_symbol,
            "ensembl_id": ensembl_id,
            "matched_hpo": None,
            "matched_disease": None,
            "overall_association_score": (
                round(float(top_overall), 3) if top_overall is not None else None
            ),
            "datatype_scores": datatype_scores,
            "top_diseases": top_diseases,
            "evidence_summary": summary,
            "url": url,
        }


async def fetch_opentargets(
    ensembl_id: str, hpo_terms=None, gene_symbol: str | None = None
) -> dict:
    """Module-level helper matching the ``fetch_*`` pattern used by the
    rest of ``backend/clients``. Wraps OpenTargetsClient so app.py can
    schedule it the same way as the other DB tasks.

    ``gene_symbol`` is the curated gene, used to resolve an Ensembl gene id
    when VEP supplied an EntrezGene one (RefSeq input) or none at all.
    """
    client = OpenTargetsClient()
    try:
        result = await client.get_gene_disease_evidence(
            ensembl_id, hpo_terms, gene_symbol_hint=gene_symbol
        )
    except Exception as e:
        log.exception("Open Targets fetch raised for %s", ensembl_id)
        return {"ok": False, "error": repr(e)}
    if result is None:
        return {
            "ok": False,
            "error": client.last_error or "Open Targets unavailable or no data",
        }
    return result
