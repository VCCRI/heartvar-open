#!/usr/bin/env python3
"""Bundle AlphaFold protein structures for the cardiovascular gene panel.

HeartVar's 3-D structure viewer (Protein tab) is served FULLY OFFLINE: the
deploy box serves AlphaFold ``.pdb`` files from its own ``/structures`` mount,
so the curator's browser makes no external call at view time. This script is
the BUILD STEP that populates that local cache — run it once on the build /
deploy box (same pattern as the offline VEP / gnomAD / ClinVar caches, which
are also gitignored build artifacts).

What it does, per gene in ``backend/data/cvd_gene_panel.json``:
  1. Resolve the gene symbol → UniProt accession from the LOCAL UniProt DB
     (``data/uniprot.db``; no live UniProt call needed — 671/677 panel genes
     resolve, the misses being mitochondrial tRNA genes with no protein).
  2. Ask the AlphaFold API for the current model URL (handles the version
     drift — v4 → v6 etc. — so the file URL is never hard-coded).
  3. Download the fragment-1 ``.pdb`` to ``data/alphafold/<ACCESSION>.pdb``.

The frontend requests ``/structures/<ACCESSION>.pdb`` (the accession is already
in ``db_evidence``), so the on-disk filename is just ``<ACCESSION>.pdb``.

Model selection & residue numbering: AlphaFold serves some large proteins only
as *isoform* models (e.g. CHD7 → Q9P2D1-4), whose residue numbering does NOT
match the canonical sequence a curator's HGVS p. notation is written against.
For each gene this script picks the canonical-accession model when one exists,
otherwise the isoform whose sequence shares the longest common prefix with the
canonical sequence, and records ``safe_max_residue`` — the highest residue index
the model numbers canonically. The viewer only highlights at or below that index
(and double-checks the residue identity), and otherwise shows an honest "outside
the reliable region — open in AlphaFold" message instead of a possibly-wrong pin.
Genes with no usable model (giants like TTN/RYR1; or no canonical alignment)
record ``no_model`` and the viewer falls back to the external AlphaFold link.

Usage:
    python scripts/build_alphafold_structures.py                # whole panel
    python scripts/build_alphafold_structures.py --genes MYH7,HRAS,TPM1
    python scripts/build_alphafold_structures.py --limit 25 --force

Idempotent: existing ``.pdb`` files are skipped unless ``--force``.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

import httpx

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
_PANEL = _PROJECT_ROOT / "backend" / "data" / "cvd_gene_panel.json"
_UNIPROT_DB = _PROJECT_ROOT / "data" / "uniprot.db"
_OUT_DIR = _PROJECT_ROOT / "data" / "alphafold"
_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
_UA = "HeartVar-structure-builder/1.0 (offline cache build)"
_TIMEOUT = 30
_MIN_ISOFORM_SAFE = 50


def _panel_genes() -> list[str]:
    data = json.loads(_PANEL.read_text())
    genes = data.get("genes") if isinstance(data, dict) else data
    if isinstance(genes, dict):
        return list(genes.keys())
    return list(genes or [])


def _accession_for(conn: sqlite3.Connection, gene: str) -> str | None:
    row = conn.execute(
        "SELECT accession FROM uniprot_entry WHERE gene_symbol = ? LIMIT 1",
        (gene,),
    ).fetchone()
    return row[0] if row and row[0] else None


_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O",
}
_UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def _lcp(a: str, b: str) -> int:
    """Length of the longest common prefix of two sequences."""
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _canonical_seq(client: httpx.Client, acc: str) -> str | None:
    """Canonical UniProt sequence for a bare accession — the sequence the
    curator's HGVS p. notation is numbered against. Used to align an
    isoform-only AlphaFold model back to canonical residue numbering."""
    try:
        r = client.get(_UNIPROT_FASTA.format(acc=acc))
        if r.status_code != 200:
            return None
        return "".join(s for s in r.text.splitlines() if not s.startswith(">"))
    except httpx.HTTPError:
        return None


def _pdb_ca_sequence(pdb_bytes: bytes) -> str:
    """One-letter sequence read from a PDB's CA atoms, ordered by residue
    number. AlphaFold models are gapless and numbered 1..N, so position i in
    this string is the residue the file numbers `i`."""
    res: dict[int, str] = {}
    for ln in pdb_bytes.decode("ascii", "replace").splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            try:
                res[int(ln[22:26])] = _AA3.get(ln[17:20].strip(), "X")
            except ValueError:
                pass
    return "".join(res[k] for k in sorted(res))


def _model_acc_from_pdb(pdb_bytes: bytes) -> str | None:
    """The UniProt accession a bundled PDB was built from, read from its TITLE
    (e.g. '... (Q9P2D1-4)' → 'Q9P2D1-4'). A bare accession means a canonical
    model; an isoform suffix means it was aligned to canonical."""
    title = ""
    for ln in pdb_bytes.decode("ascii", "replace").splitlines():
        if ln.startswith("TITLE"):
            title += ln[10:].rstrip()
        elif ln.startswith(("ATOM", "HETATM")):
            break
    m = re.search(r"\(([A-Za-z0-9]+(?:-\d+)?)\)\s*$", title)
    return m.group(1) if m else None


def _select_model(client: httpx.Client, acc: str):
    """Choose which AlphaFold model to bundle for `acc`, and the residue range
    over which its numbering is provably the canonical numbering.

    Returns (pdb_url, model_acc, safe_max, model_seq) or None when AlphaFold has
    no usable model. `safe_max` is the highest residue index safe to highlight:
    for a canonical-accession model that's the whole model; for an isoform-only
    model it's the longest common prefix with the canonical sequence (past it
    the isoform diverges, so the numbering can no longer be trusted)."""
    resp = client.get(_API.format(acc=acc))
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    entries = resp.json()
    if not isinstance(entries, list) or not entries:
        return None
    canon = [e for e in entries if e.get("uniprotAccession") == acc]
    if canon:
        e = min(canon, key=lambda e: e.get("uniprotStart") or 1)
        seq = e.get("uniprotSequence") or ""
        safe = min(e.get("uniprotEnd") or len(seq), len(seq))
        return e.get("pdbUrl"), acc, safe, seq
    cseq = _canonical_seq(client, acc)
    if not cseq:
        return None
    best = max(entries, key=lambda e: _lcp(cseq, e.get("uniprotSequence") or ""))
    bseq = best.get("uniprotSequence") or ""
    safe = _lcp(cseq, bseq)
    if safe < _MIN_ISOFORM_SAFE:
        return None
    return best.get("pdbUrl"), best.get("uniprotAccession"), safe, bseq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", help="comma-separated subset (default: whole panel)")
    ap.add_argument("--limit", type=int, default=0, help="cap number of genes (testing)")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    ap.add_argument("--delay", type=float, default=0.15, help="seconds between requests")
    args = ap.parse_args()

    if not _UNIPROT_DB.exists():
        print(f"ERROR: local UniProt DB not found at {_UNIPROT_DB}", file=sys.stderr)
        return 2

    genes = (
        [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        if args.genes else _panel_genes()
    )
    if args.limit:
        genes = genes[: args.limit]

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(f"file:{_UNIPROT_DB}?mode=ro", uri=True)
    client = httpx.Client(
        follow_redirects=True, timeout=_TIMEOUT, headers={"User-Agent": _UA},
    )

    prior: dict[str, dict] = {}
    mpath = _OUT_DIR / "manifest.json"
    if mpath.exists():
        try:
            prior = json.loads(mpath.read_text())
        except (ValueError, OSError):
            prior = {}

    manifest: dict[str, dict] = dict(prior)
    n_ok = n_skip = n_no_acc = n_no_model = n_err = 0
    try:
        for i, gene in enumerate(genes, 1):
            acc = _accession_for(conn, gene)
            if not acc:
                n_no_acc += 1
                manifest[gene] = {"accession": None, "status": "no_accession"}
                continue
            out = _OUT_DIR / f"{acc}.pdb"
            pinfo = prior.get(gene) or {}
            if (out.exists() and not args.force
                    and pinfo.get("status") == "downloaded"
                    and pinfo.get("model")
                    and isinstance(pinfo.get("safe_max_residue"), int)):
                n_skip += 1
                manifest[gene] = pinfo
                continue
            try:
                if out.exists() and not args.force:
                    cached = out.read_bytes()
                    if _model_acc_from_pdb(cached) == acc:
                        seq = _pdb_ca_sequence(cached)
                        n_skip += 1
                        manifest[gene] = {
                            "accession": acc, "status": "downloaded",
                            "model": acc, "safe_max_residue": len(seq),
                            "bytes": len(cached),
                        }
                        continue
                sel = _select_model(client, acc)
                if not sel:
                    if out.exists():
                        out.unlink()
                    n_no_model += 1
                    manifest[gene] = {"accession": acc, "status": "no_model"}
                    print(f"[{i}/{len(genes)}] {gene} ({acc}): no usable AlphaFold model")
                    continue
                url, model_acc, safe_max, model_seq = sel
                resp = client.get(url)
                resp.raise_for_status()
                data = resp.content
                pdb_seq = _pdb_ca_sequence(data)
                safe_max = min(safe_max, len(pdb_seq))
                if pdb_seq[:safe_max] != model_seq[:safe_max]:
                    n_err += 1
                    manifest[gene] = {"accession": acc, "status": "error: pdb/seq mismatch"}
                    print(f"[{i}/{len(genes)}] {gene} ({acc}): PDB/seq mismatch — skipped",
                          file=sys.stderr)
                    continue
                out.write_bytes(data)
                n_ok += 1
                manifest[gene] = {
                    "accession": acc, "status": "downloaded",
                    "model": model_acc, "safe_max_residue": safe_max,
                    "bytes": len(data), "source": url,
                }
                tag = "" if model_acc == acc else f"  isoform {model_acc} safe≤{safe_max}"
                print(f"[{i}/{len(genes)}] {gene} ({acc}): {len(data)//1024} KB{tag}")
                time.sleep(args.delay)
            except (httpx.HTTPError, OSError, ValueError) as e:
                n_err += 1
                manifest[gene] = {"accession": acc, "status": f"error: {e!r}"}
                print(f"[{i}/{len(genes)}] {gene} ({acc}): ERROR {e!r}", file=sys.stderr)
    finally:
        conn.close()
        client.close()

    _manifest_final = _OUT_DIR / "manifest.json"
    _manifest_tmp = _manifest_final.with_name(_manifest_final.name + ".tmp")
    _manifest_tmp.write_text(json.dumps(manifest, indent=1))
    os.replace(_manifest_tmp, _manifest_final)
    print(
        f"\nDone. downloaded={n_ok} cached={n_skip} no_accession={n_no_acc} "
        f"no_model={n_no_model} errors={n_err} (total genes={len(genes)})"
    )
    print(f"Structures → {_OUT_DIR}  (serve via the /structures mount)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
