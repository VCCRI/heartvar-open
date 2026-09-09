#!/usr/bin/env python3
"""build_hpo_labels.py — build a slim HPO term-id → label lookup.

Downloads the latest ``hp.obo`` release from the canonical Human
Phenotype Ontology GitHub release and writes a JSON map of
``HP:NNNNNNN`` → primary term label to
``backend/data/hpo_labels.json``. Used at runtime to turn curator-
submitted HPO IDs into plain-text phenotype names for downstream
search URLs (e.g. PubMed) that don't accept ontology codes.

Re-run quarterly. The HPO release line is stable; the term names
rarely change but new terms are added in every release.

Usage::

    python3 build_hpo_labels.py                 # download + rebuild
    python3 build_hpo_labels.py --force-download
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import ssl
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

URL = (
    "https://github.com/obophenotype/human-phenotype-ontology/"
    "releases/latest/download/hp.obo"
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_ROOT / "data"
CURATED_DIR = PROJECT_ROOT / "backend" / "data"
OBO_PATH = RAW_DIR / "hp.obo"
JSON_PATH = CURATED_DIR / "hpo_labels.json"
DESCENDANTS_JSON_PATH = CURATED_DIR / "cardiac_hpo_descendants.json"
SYNONYMS_JSON_PATH = CURATED_DIR / "cardiac_hpo_synonyms.json"


def _download_curl(url: str, dest: Path) -> None:
    """Fallback when Python's bundled CA store can't verify the host —
    macOS Python builds frequently can't resolve modern Let's Encrypt
    chains; system curl uses the OS keychain which has them."""
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError("curl not found on PATH; cannot bypass SSL failure")
    cmd = [curl, "--location", "--fail", "--retry", "3",
           "--silent", "--show-error", "--output", str(dest), url]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"curl exited with status {proc.returncode}")


def download_obo(force: bool) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CURATED_DIR.mkdir(parents=True, exist_ok=True)
    if OBO_PATH.exists() and not force:
        size_mb = OBO_PATH.stat().st_size / 1024 / 1024
        print(f"[hpo] using cached {OBO_PATH} ({size_mb:.1f} MB)")
        return OBO_PATH
    print(f"[hpo] downloading {URL}")
    try:
        with urllib.request.urlopen(URL, timeout=120) as resp:
            data = resp.read()
        OBO_PATH.write_bytes(data)
    except (ssl.SSLError, urllib.error.URLError) as e:
        msg = str(e)
        if "CERTIFICATE_VERIFY_FAILED" in msg or "SSL" in msg:
            print(f"[hpo] urllib TLS failed ({e}); retrying via curl")
            _download_curl(URL, OBO_PATH)
        else:
            print(f"[hpo] download failed: {e}", file=sys.stderr)
            sys.exit(1)
    size_mb = OBO_PATH.stat().st_size / 1024 / 1024
    print(f"[hpo] wrote {OBO_PATH} ({size_mb:.1f} MB)")
    return OBO_PATH


def parse_obo(path: Path) -> dict[str, str]:
    """Walk hp.obo and emit ``HP:NNNNNNN`` → primary ``name`` mapping.

    Skips any [Typedef] stanzas and entries marked ``is_obsolete: true``.
    Only the first ``id:`` line in each [Term] is honoured; alt-ids are
    ignored (callers should canonicalise via the obsolete-replacement
    fields if they need full alt-id resolution, which we don't here).
    """
    out: dict[str, str] = {}
    in_term = False
    cur_id: str | None = None
    cur_name: str | None = None
    cur_obsolete = False
    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("[Term]"):
                if in_term and cur_id and cur_name and not cur_obsolete:
                    out[cur_id] = cur_name
                in_term = True
                cur_id = None
                cur_name = None
                cur_obsolete = False
                continue
            if line.startswith("[") and line.endswith("]"):
                if in_term and cur_id and cur_name and not cur_obsolete:
                    out[cur_id] = cur_name
                in_term = False
                cur_id = None
                cur_name = None
                cur_obsolete = False
                continue
            if not in_term:
                continue
            if line.startswith("id: ") and cur_id is None:
                cur_id = line[4:].strip()
            elif line.startswith("name: ") and cur_name is None:
                cur_name = line[6:].strip()
            elif line.startswith("is_obsolete: ") and line.endswith("true"):
                cur_obsolete = True
        if in_term and cur_id and cur_name and not cur_obsolete:
            out[cur_id] = cur_name
    out = {k: v for k, v in out.items() if k.startswith("HP:")}
    return out


def parse_isa_edges(path: Path) -> dict[str, list[str]]:
    """Walk hp.obo and emit a parent→children adjacency map from ``is_a``.

    Each non-obsolete [Term] with ``id: HP:child`` and one or more
    ``is_a: HP:parent ! label`` lines contributes ``parent -> child`` edges.
    Only HP: ids are kept (hp.obo cross-references foreign ontologies).
    """
    children: dict[str, list[str]] = {}
    in_term = False
    cur_id: str | None = None
    cur_parents: list[str] = []
    cur_obsolete = False

    def flush() -> None:
        if in_term and cur_id and not cur_obsolete and cur_id.startswith("HP:"):
            for parent in cur_parents:
                if parent.startswith("HP:"):
                    children.setdefault(parent, []).append(cur_id)

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("[Term]"):
                flush()
                in_term, cur_id, cur_parents, cur_obsolete = True, None, [], False
                continue
            if line.startswith("[") and line.endswith("]"):
                flush()
                in_term, cur_id, cur_parents, cur_obsolete = False, None, [], False
                continue
            if not in_term:
                continue
            if line.startswith("id: ") and cur_id is None:
                cur_id = line[4:].strip()
            elif line.startswith("is_a: "):
                cur_parents.append(line[6:].split("!", 1)[0].strip())
            elif line.startswith("is_obsolete: ") and line.endswith("true"):
                cur_obsolete = True
        flush()
    return children


def _descendants(children: dict[str, list[str]], root: str) -> set[str]:
    """All transitive children of ``root`` (excluding ``root`` itself),
    matching JAX's /descendants semantics."""
    seen: set[str] = set()
    stack = list(children.get(root, ()))
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(children.get(node, ()))
    return seen


def build_cardiac_descendants(children: dict[str, list[str]]) -> dict[str, list[str]]:
    """Expand the 5 cardiac ancestor roots into the {HP:id -> [category…]}
    closure PanelApp uses, reproducing ``_ensure_descendant_cache`` offline:
    each descendant inherits the union of its ancestor roots' categories.
    """
    sys.path.insert(0, str(PROJECT_ROOT))
    from backend.clients.panelapp import _PARENT_TO_CATEGORIES

    closure: dict[str, set[str]] = {}
    for root, cats in _PARENT_TO_CATEGORIES.items():
        for d in _descendants(children, root):
            closure.setdefault(d, set()).update(cats)
    return {hpo_id: sorted(cats) for hpo_id, cats in sorted(closure.items())}


def parse_exact_synonyms(path: Path) -> dict[str, list[str]]:
    """Walk hp.obo and emit ``HP:NNNNNNN`` → [EXACT synonym…].

    EXACT scope ONLY. BROAD synonyms name a more general concept than the
    term itself, so accepting them would let a specific phenotype be
    recognised through a looser label and widen category matching; NARROW
    does the reverse. Neither is safe when the result gates PP4.
    """
    out: dict[str, list[str]] = {}
    in_term = False
    cur_id: str | None = None
    cur_syns: list[str] = []
    cur_obsolete = False

    def flush() -> None:
        if in_term and cur_id and cur_syns and not cur_obsolete:
            out.setdefault(cur_id, []).extend(cur_syns)

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("[") and line.endswith("]"):
                flush()
                in_term = line.startswith("[Term]")
                cur_id, cur_syns, cur_obsolete = None, [], False
                continue
            if not in_term:
                continue
            if line.startswith("id: ") and cur_id is None:
                cur_id = line[4:].strip()
            elif line.startswith("is_obsolete: ") and line.endswith("true"):
                cur_obsolete = True
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "((?:[^"\\]|\\.)*)" EXACT', line)
                if m:
                    cur_syns.append(m.group(1).replace('\\"', '"'))
        flush()
    return {k: v for k, v in out.items() if k.startswith("HP:")}


def build_cardiac_synonyms(
    labels: dict[str, str],
    synonyms: dict[str, list[str]],
    descendants: dict[str, list[str]],
) -> dict[str, str]:
    """Build the {normalised phrase -> HP:id} lookup for curator free text.

    Restricted to terms inside the cardiac closure — the index exists to
    recognise cardiovascular phenotypes, and indexing all ~20k HPO terms
    would let unrelated phenotypes resolve to an id that matches no
    category anyway.

    Includes each term's primary name plus its EXACT synonyms, which is
    where the acronyms live (DORV, TGA, D-TGA, ASD, TAPVR, HCM, SVT, PFO
    are all HPO's own).

    Keys of <= 3 characters are admitted only when the original synonym is
    all-caps, i.e. acronym-shaped. Without that guard an ontology label
    could hand us a key like "as" or "on" and turn an ordinary English word
    into a phenotype term. Curator-facing shorthand that HPO lacks is
    handled separately by the curated map in the PanelApp client, where a
    human has signed off on each entry.
    """
    sys.path.insert(0, str(PROJECT_ROOT))
    from backend.clients.panelapp import _normalise_phrase

    index: dict[str, str] = {}
    for hpo_id in descendants:
        candidates = []
        if labels.get(hpo_id):
            candidates.append(labels[hpo_id])
        candidates.extend(synonyms.get(hpo_id, []))
        for raw in candidates:
            key = _normalise_phrase(raw)
            if not key:
                continue
            if len(key) <= 3 and raw != raw.upper():
                continue
            index.setdefault(key, hpo_id)
    return dict(sorted(index.items()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--force-download",
        action="store_true",
        help="redownload hp.obo even if cached copy exists",
    )
    args = ap.parse_args()

    obo_path = download_obo(force=args.force_download)
    labels = parse_obo(obo_path)
    if not labels:
        print("[hpo] no HP: terms parsed — aborting", file=sys.stderr)
        return 1
    json_tmp = JSON_PATH.with_name(JSON_PATH.name + ".tmp")
    json_tmp.write_text(json.dumps(labels, ensure_ascii=False, sort_keys=True))
    os.replace(json_tmp, JSON_PATH)
    print(f"[hpo] wrote {JSON_PATH} ({len(labels):,} terms)")

    children = parse_isa_edges(obo_path)
    descendants = build_cardiac_descendants(children)
    if not descendants:
        print("[hpo] WARNING: cardiac descendant closure empty — skipping",
              file=sys.stderr)
    else:
        descendants_tmp = DESCENDANTS_JSON_PATH.with_name(
            DESCENDANTS_JSON_PATH.name + ".tmp")
        descendants_tmp.write_text(
            json.dumps(descendants, ensure_ascii=False, sort_keys=True))
        os.replace(descendants_tmp, DESCENDANTS_JSON_PATH)
        print(f"[hpo] wrote {DESCENDANTS_JSON_PATH} "
              f"({len(descendants):,} cardiac descendant IDs)")

        synonyms = parse_exact_synonyms(obo_path)
        index = build_cardiac_synonyms(labels, synonyms, descendants)
        if not index:
            print("[hpo] WARNING: cardiac synonym index empty — skipping",
                  file=sys.stderr)
        else:
            syn_tmp = SYNONYMS_JSON_PATH.with_name(
                SYNONYMS_JSON_PATH.name + ".tmp")
            syn_tmp.write_text(
                json.dumps(index, ensure_ascii=False, sort_keys=True))
            os.replace(syn_tmp, SYNONYMS_JSON_PATH)
            print(f"[hpo] wrote {SYNONYMS_JSON_PATH} "
                  f"({len(index):,} phrase → HP-id entries)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
