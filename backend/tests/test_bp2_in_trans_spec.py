"""BP2 must mean what ACMG says it means: a SECOND pathogenic allele.

ACMG/AMP 2015 (Richards et al., Table 4) defines BP2 as:

    "Observed in trans with a pathogenic variant for a fully penetrant
     dominant gene/disorder; or observed in cis with a pathogenic variant
     in any inheritance pattern."

Both arms require a *second, pathogenic variant* in the same gene. Neither arm
says anything about who transmitted the variant.

The engine originally fired BP2 on (AD + het + inherited_unaffected) — i.e.
"inherited from an unaffected parent" — which is a different observation
entirely (a penetrance/segregation argument, BS2/BS4 territory) and never
consulted the ``in_trans_pathogenic`` input the app already collects and
already passes to PM3. Because BP2 is frequently the ONLY benign criterion to
fire, that spurious -1 was enough on its own to drop a variant to Likely
benign under the Tavtigian point system.

Regression case that motivated this: CHD7 8-60794987-A-G (hg38), DORV-TGA/VSD
proband, heterozygous, AD, full trio, inherited from an unaffected parent, no
second allele reported anywhere. The engine returned Likely benign (ClinVar:
VUS) citing BP2 "in trans / cis" while the curator had never entered a partner
variant.

These tests pin the spec-correct behaviour in both directions: BP2 fires when a
confirmed in-trans pathogenic allele IS present, and stays silent when it is
not.
"""
from __future__ import annotations

import pytest

from backend.acmg.hard_coded import compute_hard_coded_criteria
from backend.acmg.tiers import classification_for, compute_points_total


def _crit(cc: dict, gene: str = "ZZZZ9") -> dict:
    """Run the deterministic engine and return {code: entry}."""
    out = compute_hard_coded_criteria({}, cc, gene=gene)
    return {c["code"]: c for c in out}


def _bp2(**cc) -> dict:
    base = {"inheritance_input": "AD", "zygosity": "het"}
    base.update(cc)
    return _crit(base)["BP2"]


def test_bp2_not_met_when_inherited_from_unaffected_parent_only():
    """Transmission by an unaffected parent is NOT a BP2 observation.

    This is the exact reported case: the sole trigger was the inheritance
    field, with in_trans_pathogenic left blank.
    """
    entry = _bp2(trio_status="trio", denovo_status="inherited_unaffected")
    assert entry["status"] == "not_met", (
        "BP2 fired without any second pathogenic allele — it was keying off "
        f"the inheritance field: {entry.get('evidence')!r}"
    )
    assert entry.get("criteria_strength") in (None, "")


def test_chd7_reported_case_does_not_reach_likely_benign_on_bp2():
    """The reported CHD7 trio must not be pushed benign by a phantom BP2."""
    cc = {
        "inheritance_input": "AD",
        "zygosity": "het",
        "trio_status": "trio",
        "denovo_status": "inherited_unaffected",
        "in_trans_pathogenic": "",
        "sex": "male",
    }
    crit = _crit(cc, gene="CHD7")
    assert crit["BP2"]["status"] == "not_met"

    met = [c for c in crit.values() if c["status"] == "met"]
    assert not [c for c in met if c["direction"] == "benign"], (
        "no benign criterion is supported by these inputs, but the engine "
        f"applied: {[(c['code'], c.get('evidence')) for c in met]}"
    )
    pts = compute_points_total(list(crit.values()))
    assert classification_for(pts) != "Likely benign"


@pytest.mark.parametrize("in_trans", ["", "no", "unknown"])
def test_bp2_not_met_when_second_allele_absent_or_denied(in_trans):
    """Only an affirmative in-trans pathogenic allele can support BP2."""
    entry = _bp2(denovo_status="inherited_unaffected", in_trans_pathogenic=in_trans)
    assert entry["status"] == "not_met"


def test_bp2_met_when_in_trans_pathogenic_confirmed_in_dominant_gene():
    """BP2's actual first arm: het, dominant disorder, P/LP allele in trans."""
    entry = _bp2(in_trans_pathogenic="yes")
    assert entry["status"] == "met"
    assert entry["criteria_strength"] == "BP2_Supporting"
    assert entry["direction"] == "benign"
    ev_text = entry["evidence"].lower()
    assert "trans" in ev_text and "pathogenic" in ev_text


def test_bp2_met_regardless_of_who_transmitted_it():
    """Phase, not provenance, is what BP2 turns on."""
    for status in ("inherited_unaffected", "inherited_affected", "unconfirmed", ""):
        entry = _bp2(in_trans_pathogenic="yes", denovo_status=status)
        assert entry["status"] == "met", f"denovo_status={status!r} suppressed BP2"


def test_recessive_in_trans_routes_to_pm3_not_bp2():
    crit = _crit({
        "inheritance_input": "AR", "zygosity": "het", "in_trans_pathogenic": "yes",
    })
    assert crit["PM3"]["status"] == "met"
    assert crit["BP2"]["status"] == "not_met"


def test_dominant_in_trans_routes_to_bp2_not_pm3():
    crit = _crit({
        "inheritance_input": "AD", "zygosity": "het", "in_trans_pathogenic": "yes",
    })
    assert crit["BP2"]["status"] == "met"
    assert crit["PM3"]["status"] == "not_met"


def test_bp2_requires_heterozygous():
    """A homozygous call cannot be "in trans with" a different allele."""
    for zyg in ("hom", "hemi"):
        entry = _crit({
            "inheritance_input": "AD", "zygosity": zyg,
            "in_trans_pathogenic": "yes",
        })["BP2"]
        assert entry["status"] == "not_met", f"BP2 fired for zygosity={zyg!r}"
