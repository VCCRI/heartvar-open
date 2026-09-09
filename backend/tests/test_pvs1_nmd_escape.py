"""Unit tests for the PVS1 NMD-escape strength refinement.

Completes HeartVar's Abou-Tayoun / ClinGen SVI PVS1 decision tree: a
nonsense/frameshift variant predicted to ESCAPE NMD (PTC in the last exon /
last ~50 bp of the penultimate exon) is PVS1_Strong only if it removes >10% of
the protein (or a critical region), otherwise PVS1_Moderate. Previously such
variants were flat-downgraded to PVS1_Strong regardless of how little protein
was truncated (e.g. MYBPC3 p.Arg1271Ter removes ~0.2%).
"""
from __future__ import annotations

from backend.acmg.hard_coded import _pvs1_nmd_escape_strength


def test_tiny_cterminal_truncation_is_moderate():
    tier, why = _pvs1_nmd_escape_strength({"protein_start": 1271}, {"uniprot": {"length": 1274}})
    assert tier == "PVS1_Moderate"
    assert "≤10%" in why


def test_large_truncation_stays_strong():
    tier, _ = _pvs1_nmd_escape_strength({"protein_start": 500}, {"uniprot": {"length": 1274}})
    assert tier == "PVS1_Strong"


def test_boundary_just_over_ten_percent_is_strong():
    assert _pvs1_nmd_escape_strength({"protein_start": 899}, {"uniprot": {"length": 1000}})[0] == "PVS1_Strong"
    assert _pvs1_nmd_escape_strength({"protein_start": 901}, {"uniprot": {"length": 1000}})[0] == "PVS1_Moderate"


def test_exon_lengths_fallback_when_no_uniprot():
    tier, _ = _pvs1_nmd_escape_strength({"protein_start": 1271, "exon_lengths": [3822]}, {})
    assert tier == "PVS1_Moderate"


def test_unquantifiable_defaults_to_strong():
    for vep, ev in (({}, {}), ({"protein_start": 100}, {}), ({}, {"uniprot": {"length": 500}})):
        assert _pvs1_nmd_escape_strength(vep, ev)[0] == "PVS1_Strong"


def test_hgvsp_position_used_when_no_protein_start():
    tier, _ = _pvs1_nmd_escape_strength(
        {"hgvsp": "ENSP00000001.1:p.Arg1271Ter"}, {"uniprot": {"length": 1274}})
    assert tier == "PVS1_Moderate"


if __name__ == "__main__":
    import sys
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
