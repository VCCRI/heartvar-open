"""BP3 must require a repeat region that has NO known function.

ACMG/AMP 2015: BP3 = "In-frame deletions/insertions in a repetitive region
WITHOUT a known function." Two conditions, and the engine honoured neither:

  * "repetitive region" was read from ``vep["in_repeat_region"]``, a key NOTHING
    in production ever writes (grep: only tests set it). So BP3 was dead code —
    it could never fire on a real curation, and PM4's repeat carve-out never
    engaged either.
  * "without a known function" was not checked at all.

BP3 is now derived from UniProt annotation: the variant must sit inside an
annotated Repeat feature, and that repeat must not overlap any feature that
asserts a function (domain, motif, zinc finger, DNA-binding, active site,
binding site, transmembrane, signal peptide), nor be a named functional repeat
family (WD40, ankyrin, LRR, ...) which are structural folds rather than the
low-complexity filler BP3 targets.

Fail-closed throughout: absent UniProt data or an unresolvable protein position
means "without a known function" cannot be established, so BP3 stays not_met.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import _eval_pm4_bp3


def _up(features, *, ok=True, domains=None):
    return {"uniprot": {"ok": True, "found": True,
                        "features": features,
                        "domains": domains or []}} if ok else {"uniprot": {"ok": False}}


def _feat(t, start, end, desc=""):
    return {"type": t, "start": start, "end": end, "description": desc}


def _run(ev, vep, consequence="inframe_deletion"):
    return {c["code"]: c for c in _eval_pm4_bp3(vep, consequence, ev)}


VEP = {"protein_start": 100}


def test_bp3_fires_for_inframe_indel_in_a_functionless_repeat():
    ev = _up([_feat("Repeat", 90, 120, "Gln-rich")])
    d = _run(ev, VEP)
    assert d["BP3"]["status"] == "met"
    assert d["BP3"]["criteria_strength"] == "BP3_Supporting"
    assert d["PM4"]["status"] == "not_met", "PM4 and BP3 are mutually exclusive"


@pytest.mark.parametrize("ftype", [
    "Domain", "Motif", "Zinc finger", "DNA binding",
    "Active site", "Binding site", "Transmembrane", "Signal",
])
def test_bp3_blocked_when_the_repeat_overlaps_a_functional_feature(ftype):
    ev = _up([
        _feat("Repeat", 90, 120, "repeat 3"),
        _feat(ftype, 110, 115, "something functional"),
    ])
    d = _run(ev, VEP)
    assert d["BP3"]["status"] == "not_met", (
        f"BP3 fired in a repeat overlapping a {ftype}: {d['BP3'].get('evidence')!r}"
    )
    assert d["PM4"]["status"] == "met", (
        "a length change in a functionally annotated region is PM4, not nothing"
    )


@pytest.mark.parametrize("desc", [
    "WD 4", "ANK 2", "Leucine-rich repeat 7", "LRR 3", "TPR 1",
    "Armadillo 5", "Kelch 2", "EF-hand 3", "Spectrin 4",
    "Fibronectin type-III 2", "Cadherin 1", "Ig-like C2-type 3",
])
def test_bp3_blocked_for_named_functional_repeat_families(desc):
    """WD40/ankyrin/LRR repeats are functional folds, not low-complexity filler."""
    ev = _up([_feat("Repeat", 90, 120, desc)])
    d = _run(ev, VEP)
    assert d["BP3"]["status"] == "not_met", (
        f"BP3 fired on a {desc!r} repeat: {d['BP3'].get('evidence')!r}"
    )


def test_generic_region_annotation_does_not_block_bp3():
    """UniProt "Region" is generic (often "Disordered") — not a known function.

    Treating it as functional would block BP3 almost everywhere.
    """
    ev = _up([
        _feat("Repeat", 90, 120, "Pro-rich"),
        _feat("Region", 80, 130, "Disordered"),
    ])
    assert _run(ev, VEP)["BP3"]["status"] == "met"


def test_bp3_not_met_when_variant_is_outside_every_annotated_repeat():
    ev = _up([_feat("Repeat", 200, 240, "Gln-rich")])
    d = _run(ev, VEP)
    assert d["BP3"]["status"] == "not_met"
    assert d["PM4"]["status"] == "met"


def test_bp3_not_met_when_protein_has_no_repeat_annotation():
    d = _run(_up([_feat("Domain", 50, 200, "Myosin motor")]), VEP)
    assert d["BP3"]["status"] == "not_met"
    assert d["PM4"]["status"] == "met"


def test_bp3_fails_closed_when_uniprot_unavailable():
    d = _run(_up([], ok=False), VEP)
    assert d["BP3"]["status"] == "not_met"
    assert "uniprot" in d["BP3"]["evidence"].lower() or \
           "annotation" in d["BP3"]["evidence"].lower()


def test_bp3_fails_closed_when_protein_position_unresolvable():
    ev = _up([_feat("Repeat", 90, 120, "Gln-rich")])
    d = _run(ev, {})
    assert d["BP3"]["status"] == "not_met"


def test_bp3_fails_closed_with_no_evidence_at_all():
    d = _run({}, {})
    assert d["BP3"]["status"] == "not_met"
    assert d["PM4"]["status"] == "met", "in-frame indel still earns PM4 by default"


def test_legacy_in_repeat_region_flag_no_longer_drives_bp3():
    """Nothing in production ever set this key; it must not be a back door."""
    d = _run({}, {"protein_start": 100, "in_repeat_region": True})
    assert d["BP3"]["status"] == "not_met", (
        "BP3 still keys off the never-populated in_repeat_region flag"
    )


def test_stop_loss_still_earns_pm4_regardless_of_repeats():
    ev = _up([_feat("Repeat", 90, 120, "Gln-rich")])
    d = _run(ev, VEP, consequence="stop_lost")
    assert d["PM4"]["status"] == "met"
    assert d["BP3"]["status"] == "not_met", "BP3 has no stop-loss arm"


def test_missense_earns_neither():
    ev = _up([_feat("Repeat", 90, 120, "Gln-rich")])
    d = _run(ev, VEP, consequence="missense_variant")
    assert d["PM4"]["status"] == "not_met"
    assert d["BP3"]["status"] == "not_met"


def test_bp3_uses_hgvsp_position_when_protein_start_absent():
    ev = _up([_feat("Repeat", 90, 120, "Gln-rich")])
    d = _run(ev, {"hgvsp": "ENSP00000347155.1:p.Gln100del"})
    assert d["BP3"]["status"] == "met"
