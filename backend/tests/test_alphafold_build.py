"""Unit tests for the pure helpers in scripts/build_alphafold_structures.py.

These lock in the residue-numbering correctness the 3-D viewer relies on: the
longest-common-prefix alignment that defines ``safe_max_residue``, reading a
model's sequence/accession back off a PDB, and the build-time verification that
a bundled file actually matches the aligned sequence. Network-dependent code
(``_select_model`` / downloads) is exercised live by the build script itself."""
from __future__ import annotations

import importlib.util
from pathlib import Path

_BUILD = (
    Path(__file__).resolve().parent.parent.parent
    / "scripts" / "build_alphafold_structures.py"
)
_spec = importlib.util.spec_from_file_location("af_build", _BUILD)
afb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(afb)


def _pdb(residues, *, title="ALPHAFOLD MONOMER V2.0 PREDICTION FOR X (P12883)"):
    """Tiny PDB with one CA atom per (resseq, resname-3letter) in `residues`."""
    lines = [f"TITLE     {title}"]
    for resseq, resn in residues:
        lines.append(
            f"ATOM  {resseq:>5}  CA  {resn} A{resseq:>4}      "
            "0.000   0.000   0.000  1.00 90.00           C"
        )
    return ("\n".join(lines) + "\n").encode()


def test_lcp_basic():
    assert afb._lcp("MARTK", "MARQE") == 3
    assert afb._lcp("MART", "MART") == 4
    assert afb._lcp("AART", "MART") == 0
    assert afb._lcp("MAR", "MARTKLE") == 3
    assert afb._lcp("", "MART") == 0


def test_pdb_ca_sequence_orders_by_residue_number():
    pdb = _pdb([(2, "ALA"), (1, "MET"), (3, "ARG")])
    assert afb._pdb_ca_sequence(pdb) == "MAR"


def test_pdb_ca_sequence_unknown_residue_is_x():
    assert afb._pdb_ca_sequence(_pdb([(1, "FOO")])) == "X"


def test_model_acc_canonical_vs_isoform():
    canon = _pdb([(1, "MET")], title="... PREDICTION FOR MYOSIN-7 (P12883)")
    iso = _pdb([(1, "MET")], title="... PREDICTION FOR CHD7 (Q9P2D1-4)")
    assert afb._model_acc_from_pdb(canon) == "P12883"
    assert afb._model_acc_from_pdb(iso) == "Q9P2D1-4"


def test_safe_region_sequence_identity_holds():
    canonical = "MARTKLEINS"
    isoform = "MARQQ"
    safe = afb._lcp(canonical, isoform)
    assert safe == 3
    assert isoform[:safe] == canonical[:safe]
    assert isoform[safe] != canonical[safe]
