"""AlphaFold structure-availability lookup.

Unlike the other gather sources, the AlphaFold model itself is not fetched
during curation — the protein tab streams the bundled PDB from
``/structures/<ACCESSION>.pdb`` at view time (offline, no network call; see
``build_alphafold_structures.py``). This client only reports, for the loading
panel, whether a model is bundled for the gene, so AlphaFold appears in the
"Querying databases" checklist alongside every other source.

It reads the same build manifest the frontend viewer uses
(``data/alphafold/manifest.json``), keyed by gene symbol. The payload is
display-only and never reaches ``format_evidence_block`` / the prompt — the
3-D structure is a visualization aid, never an ACMG input.
"""

from __future__ import annotations

import json

from ..localio import run_local
from ._paths import data_path

_MANIFEST: dict[str, dict] | None = None


def _load_manifest() -> dict[str, dict]:
    global _MANIFEST
    if _MANIFEST is None:
        path = data_path("alphafold/manifest.json")
        try:
            _MANIFEST = json.loads(path.read_text())
        except (OSError, ValueError):
            _MANIFEST = {}
    return _MANIFEST


def _lookup(gene: str) -> dict:
    manifest = _load_manifest()
    entry = manifest.get((gene or "").upper()) or manifest.get(gene or "")
    if entry and entry.get("status") == "downloaded":
        return {
            "ok": True,
            "available": True,
            "gene": gene,
            "accession": entry.get("accession") or entry.get("model"),
            "safe_max_residue": entry.get("safe_max_residue"),
        }
    return {
        "ok": True,
        "available": False,
        "not_applicable": True,
        "skipped": True,
        "gene": gene,
        "reason": (entry or {}).get("status") or "no model bundled for this gene",
    }


async def fetch_alphafold(gene: str) -> dict:
    """Report whether a bundled AlphaFold model exists for ``gene``.

    Local manifest read (microseconds); wrapped in a thread only to keep the
    gather machinery uniform with the network-bound sources.
    """
    return await run_local(_lookup, gene)
