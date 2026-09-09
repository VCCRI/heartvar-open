"""Guards on the switch from the hand-transcribed CSpec table to the registry
harvest (backend/data/vcep_criteria_spec.json).

WHY THIS FILE EXISTS. The engine read `vcep_criteria_applicability.json`,
transcribed by hand from the CSpec PDFs. It was not wrong about anything it
recorded — but it was missing 454 entries and could not represent the 193 criteria
that publish more than one strength, and FBN1's PM1 had arrived as a single merged
cell joining a Strong rule to a Moderate one. Switching tables is therefore an
improvement AND a behaviour risk, because two things read rule TEXT rather than
just applicability: the PM1 hotspot-range restriction and the PM1/PM5
co-application check. The shapes differ — a flat `text` key versus per-strength
rows — so a careless switch leaves those reading empty strings, and both fail
OPEN. PM1 would quietly stop being range-restricted, which moves calls UP.

So these tests pin the two properties that make the switch safe:
  1. applicability is unchanged everywhere it was recorded, and
  2. the PM1 ranges parsed from the new shape are identical to the old ones.
Plus the one deliberate change: FBN1's PP4, which the registry publishes and the
hand table recorded as null.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from backend.acmg.hard_coded import (
    _APPLICABILITY_GATED_CODES,
    _VCEP_CRIT_FILE,
    _crit_spec,
    _crit_text,
    _criterion_applicable,
    _parse_pm1_hotspots,
)

DATA = Path(__file__).resolve().parents[1] / "data"
SPEC = json.loads((DATA / "vcep_criteria_spec.json").read_text())["genes"]
HAND = json.loads((DATA / "vcep_criteria_applicability.json").read_text())["genes"]

POST_SWITCHOVER_GENES = {"BMPR2", "PPP1CB"}


def _codes(entry: dict) -> set[str]:
    return {k for k, v in entry.items()
            if isinstance(v, dict) and "applicability" in v}


def test_the_engine_reads_the_registry_harvest():
    assert _VCEP_CRIT_FILE.name == "vcep_criteria_spec.json", (
        "the hand-transcribed table is superseded; it is missing 454 entries and "
        "cannot represent per-strength rules"
    )


def test_both_tables_cover_the_same_genes():
    """A coverage regression here would silently stop suppressing criteria for a
    whole gene, and default-open means that fails toward APPLICABLE.

    Genes added to the harvest after the switchover are new coverage, not a
    coverage regression, so they are excluded — the direction this test guards
    is SPEC losing a gene the hand table had."""
    assert set(SPEC) - POST_SWITCHOVER_GENES == set(HAND)
    assert POST_SWITCHOVER_GENES <= set(SPEC), (
        "a gene listed as post-switchover is missing from the harvest — either "
        "the harvest lost it or the constant is stale"
    )


def test_applicability_is_unchanged_except_three_known_fbn1_rows():
    """The safety property of the switch. Three rows differ, all FBN1, all
    unknown -> applicable — and _criterion_applicable treats both as applicable,
    so no suppression changes."""
    diffs = []
    for gene in sorted(set(SPEC) | set(HAND)):
        old, new = HAND.get(gene, {}), SPEC.get(gene, {})
        for code in sorted(_codes(old) & _codes(new)):
            oa = old[code]["applicability"]
            na = new[code]["applicability"]
            if oa != na:
                diffs.append((gene, code, oa, na))
    assert diffs == [
        ("FBN1", "PM1", "unknown", "applicable"),
        ("FBN1", "PP2", "unknown", "applicable"),
        ("FBN1", "PVS1", "unknown", "applicable"),
    ], diffs
    for gene, code, _, _ in diffs:
        assert _criterion_applicable(gene, code) is True


def test_nothing_recorded_by_hand_is_missing_from_the_registry():
    missing = [(g, c) for g in HAND for c in _codes(HAND[g])
               if c not in _codes(SPEC.get(g, {}))]
    assert missing == [], missing


def test_every_enforced_not_applicable_row_is_a_code_we_chose_to_gate():
    """Replaces the old "touch no gated code" guard, which asserted the switchover
    was behaviour-NEUTRAL. It no longer is: on 2026-09-01 nine further codes were
    added to _APPLICABILITY_GATED_CODES, so the 97 not_applicable rows the
    registry added — previously known and ignored — now all suppress.

    The guard therefore changes shape rather than being deleted. It no longer
    asks "does anything suppress?" but "does anything suppress that we did not
    deliberately gate?", which is the question that still protects us: a future
    harvest marking, say, PS4 Not Applicable on a gene would show up here.
    """
    newly = [(g, c) for g in SPEC
             for c in sorted(_codes(SPEC[g]) - _codes(HAND.get(g, {})))
             if SPEC[g][c]["applicability"] == "not_applicable"]
    assert newly, "expected the registry to add not_applicable entries"

    DELIBERATE = {
        "PP2", "BP1", "PP4", "PM1", "PVS1",
        "PM3", "PM6", "PS3",
        "BP3", "BP5", "BP6", "BS2", "BS3", "PP5",
        "BP2",
    }
    surprises = sorted({
        (g, c) for g, c in newly
        if c in _APPLICABILITY_GATED_CODES and c not in DELIBERATE
    })
    assert surprises == [], (
        f"a not_applicable row suppresses a criterion nobody reviewed: {surprises}"
    )

    enforced = [(g, c) for g, c in newly if c in _APPLICABILITY_GATED_CODES]
    assert len(enforced) == 114, (
        f"{len(enforced)} enforced not_applicable rows, expected 114 — the "
        f"harvest or the gated-code list changed"
    )
    kcnq1 = sorted(c for g, c in enforced if g == "KCNQ1")
    assert kcnq1 == ["BP2", "BP3", "BP6", "PP5"], kcnq1
    ppp1cb = sorted(c for g, c in enforced if g == "PPP1CB")
    assert ppp1cb == ["BP3", "BP6", "BS3", "PM1", "PM3", "PP4", "PP5",
                      "PS3",
                      "PVS1"], ppp1cb


def test_ppp1cb_is_a_deliberate_scoring_change():
    """PPP1CB / GN128, added 2026-09-02. 8 of the 852 cardiac-VCEP records.

    Until now PPP1CB had no rule set loaded, so all 8 were scored against
    generic house-rule defaults. Loading GN128 brings nine Not Applicable rows
    into force, and the direction is MIXED — which is why it is pinned here
    rather than assumed benign:

      * pathogenic-side suppressions, direction DOWN (toward benign):
        PM1, PM3, PP4, PS3, PVS1
      * benign-side suppressions, direction UP (toward pathogenic):
        BP3, BP6, BS3, PP5

    ⚠ The net effect on the 852 must be MEASURED before any conformance number
    is quoted. Nothing here asserts a direction on the corpus.
    """
    spec = SPEC.get("PPP1CB")
    assert spec, "PPP1CB missing from the harvest"
    assert spec["gn"] == "GN128", spec["gn"]
    assert "PPP1CB" not in HAND, "PPP1CB should be new coverage, not a changed row"

    na = sorted(c for c in _codes(spec)
                if spec[c]["applicability"] == "not_applicable")
    assert na == ["BP3", "BP6", "BS3", "PM1", "PM3", "PP4", "PP5",
                  "PS3", "PVS1"], na
    assert all(c in _APPLICABILITY_GATED_CODES for c in na), na
    applicable = sorted(c for c in _codes(spec)
                        if spec[c]["applicability"] == "applicable")
    assert len(applicable) >= 15, applicable
    for code in ("PM2", "PS4", "PP1", "BS1"):
        assert code in applicable, (code, applicable)


def test_bmpr2_is_a_deliberate_scoring_change():
    """BMPR2 is NOT behaviour-neutral, unlike the switchover itself, and this
    test is the record of what changed and which way each effect pushes.

    All three come from the rendered GN125 rule set, read on 2026-09-01. Before
    this, BMPR2 had no entry at all, so every gated criterion defaulted OPEN.
    """
    bm = SPEC["BMPR2"]

    assert bm["BP1"]["applicability"] == "not_applicable"
    assert _criterion_applicable("BMPR2", "BP1") is False

    assert bm["PP4"]["applicability"] == "not_applicable"
    assert _criterion_applicable("BMPR2", "PP4") is False

    ranges, _ = _parse_pm1_hotspots(_crit_text("BMPR2", "PM1"))
    assert (33, 131) in ranges and (203, 504) in ranges, ranges

    assert bm["PVS1"]["applicability"] == "applicable"
    assert _criterion_applicable("BMPR2", "PVS1") is True


def test_pm1_hotspot_ranges_survive_the_shape_change():
    """The real hazard. _parse_pm1_hotspots reads rule TEXT, the two tables store
    it differently, and an empty string fails OPEN — PM1 would stop being
    range-restricted, which moves calls toward pathogenic without a word in the
    log.

    Post-switchover genes are excluded: they have no hand-table text to compare
    against, so the comparison is vacuous. BMPR2's own ranges are asserted in
    test_bmpr2_is_a_deliberate_scoring_change instead."""
    for gene in sorted((set(SPEC) | set(HAND)) - POST_SWITCHOVER_GENES):
        old_text = (HAND.get(gene, {}).get("PM1") or {}).get("text") or ""
        assert _parse_pm1_hotspots(_crit_text(gene, "PM1")) == \
            _parse_pm1_hotspots(old_text), gene


def test_crit_text_joins_the_strength_rows():
    """MYH7's codon range lives in the Moderate row; FBN1's cbEGF rule is Strong
    and its EGF/TB rules Moderate. Callers scanning for a range or a caveat want
    the union of the rows, not one of them."""
    myh7 = _crit_text("MYH7", "PM1")
    assert "167-931" in myh7
    assert "phenotype other than HCM" in myh7
    fbn1 = _crit_text("FBN1", "PM1")
    assert "cbEGF" in fbn1 and "TB domain" in fbn1, fbn1


def test_crit_text_is_empty_for_an_unknown_gene():
    """Default-open depends on this being falsy rather than raising."""
    assert _crit_text("NOTAGENE", "PM1") == ""
    assert _crit_spec(None, "PM1") == {}


def test_fbn1_pp4_is_published_by_the_registry():
    """The one deliberate behaviour change. The hand table recorded FBN1 PP4 as
    null, which is where the claim "GN022 has no PP4 key at all" came from. The
    registry publishes it, so FBN1's PP4 rests on spec rather than only on the
    evidence allow-list. Outcome unchanged — experts applied PP4 on 63 of FBN1's
    114 eRepo records — but the justification is now the right one."""
    assert (HAND["FBN1"].get("PP4")) is None
    pp4 = _crit_spec("FBN1", "PP4")
    assert pp4.get("applicability") == "applicable", pp4
    assert "Ghent" in _crit_text("FBN1", "PP4")


KCNQ1_WITHDRAWN = ("BP1", "BP2", "BP3", "BS2", "PP2")


def test_kcnq1_gn112_withdrawals_are_recorded_and_enforced():
    """Quoted from the rendered GN112 rule set:

      PP2  "Not applicable due to presence of benign variation throughout the
            KCNQ1 gene (since the missense constraint Z-score in gnomAD is
            1.83, lower than 3)."
      BP1  "Not applicable, as pathogenic KCNQ1 variants are not limited to
            truncating variants, but can be missense as well."
      BP2  "Not applicable to KCNQ1 due to biallelic cases (Jervell and
            Lange-Nielsen syndrome)"
      BP3  "Not applicable to KCNQ1"
      BS2  "Not applicable due to incomplete penetrance. Please note that
            hearing loss and other phenotypes are not completely penetrant in
            homozygotes (PMID: 23392653 ...)"

    Recorded, gated, AND the text kept — a flag with the rule text discarded
    would not be auditable against the source."""
    spec = SPEC["KCNQ1"]
    assert spec["gn"] == "GN112", spec["gn"]
    for code in KCNQ1_WITHDRAWN:
        assert spec[code]["applicability"] == "not_applicable", code
        assert _criterion_applicable("KCNQ1", code) is False, code
        assert code in _APPLICABILITY_GATED_CODES, code
        assert "not applicable" in _crit_text("KCNQ1", code).lower(), code


def test_kcnq1_bp5_is_limited_not_withdrawn():
    """The discriminator's hard case, and why it has to be positional. KCNQ1
    BP5 contains "not applicable" too, but as a caveat that NARROWS BP5 rather
    than withdrawing it: "BP5 is only applicable when the phenotypes match
    another form of LQTS". Matching the phrase anywhere in the row would have
    withdrawn a criterion the VCEP publishes."""
    assert SPEC["KCNQ1"]["BP5"]["applicability"] == "applicable"
    assert _criterion_applicable("KCNQ1", "BP5") is True
    assert "only applicable when the phenotypes match" in _crit_text(
        "KCNQ1", "BP5")


def test_kcnq1_keeps_everything_gn112_does_publish():
    """A parse that withdrew too much would look like a clean fix. PP4 and PM1
    matter most: KCNQ1 is one of only two affirmative PP4 rows in the whole
    table (see _PP4_EVIDENCE_ALLOW_LIST), and its PM1 pore-helix range comes
    from the same strength rows the withdrawal check now reads."""
    spec = SPEC["KCNQ1"]
    for code in ("PVS1", "PS1", "PS4", "PM1", "PM2", "PM5", "PP1", "PP3",
                 "PP4", "BA1", "BS1", "BP4", "BP7"):
        assert spec[code]["applicability"] == "applicable", code
    ranges, _ = _parse_pm1_hotspots(_crit_text("KCNQ1", "PM1"))
    assert (300, 320) in ranges, ranges


def test_the_withdrawal_discriminator_separates_five_rows_from_thirty_two():
    """The scope boundary, proved rather than asserted.

    37 gene+criterion rows in the harvest carry the phrase "not applicable" in
    their rule text. Only the five KCNQ1 rows are withdrawals. The other 32 are
    CAVEATS on criteria that DO apply: PM1 "Not applicable to specific amino
    acid residues (see PM5)" on 13 RASopathy genes, the PM4 and PM5 comparator
    caveats on 8 cardiomyopathy genes each, BMPR2's BP4 and BP7 conditions, and
    KCNQ1's own BP5. A phrase-anywhere rule would have withdrawn all 37."""
    phrase = [(g, c) for g in sorted(SPEC) for c in sorted(_codes(SPEC[g]))
              if "not applicable" in _crit_text(g, c).lower()]
    withdrawn = [(g, c) for g, c in phrase
                 if SPEC[g][c]["applicability"] == "not_applicable"]
    COMMENT_WITHDRAWN = [("SOS1", "PM1"), ("SOS2", "PM1")]
    assert withdrawn == sorted(
        [("KCNQ1", c) for c in KCNQ1_WITHDRAWN] + COMMENT_WITHDRAWN
    ), withdrawn
    caveats = [(g, c) for g, c in phrase if (g, c) not in withdrawn]
    assert len(caveats) == 30, caveats
    assert ("KCNQ1", "BP5") in caveats
    assert len([1 for g, c in caveats if c == "PM1"]) == 11, caveats
    assert len([1 for g, c in caveats if c == "PM4"]) == 8, caveats
    assert len([1 for g, c in caveats if c == "PM5"]) == 8, caveats


def test_harvester_reads_the_first_bullet_not_the_whole_row():
    """Guards the PARSE, not only the shipped file — a re-harvest has to keep
    recognising the withdrawal. Fixtures use the registry's real markup shape:
    a <p> restating the generic ACMG definition, then <ul><li> rules."""
    import importlib.util

    path = Path(__file__).resolve().parents[2] / "scripts" / "build_vcep_criteria.py"
    modspec = importlib.util.spec_from_file_location("_bvc_test", path)
    bvc = importlib.util.module_from_spec(modspec)
    modspec.loader.exec_module(bvc)

    def row(strength: str, body: str) -> str:
        return (
            '<div class="row strength ">'
            f'<div class="col strength-label"> {strength}</div>'
            f'<div class="strength-text-cmpts col strength-text">{body}</div>'
        )

    withdrawal = row(
        "Supporting",
        "<p>Observed in trans with a pathogenic variant for a fully penetrant "
        "dominant gene/disorder.</p><ul><li>Not applicable to <em>KCNQ1</em> "
        "due to biallelic cases (Jervell and Lange-Nielsen syndrome)</li></ul>",
    )
    assert bvc.parse_criterion(withdrawal)["applicability"] == "not_applicable"
    assert "biallelic" in " ".join(
        bvc.parse_criterion(withdrawal)["strengths"].values())

    caveat = row(
        "Supporting",
        "<p>Variant found in a case with an alternate molecular basis for "
        "disease</p><ul><li>Caveat: BP5 is not applicable when the phenotypes "
        "indicate <em>KCNQ1</em> as the cause.</li></ul>",
    )
    assert bvc.parse_criterion(caveat)["applicability"] == "applicable"

    trailing = row(
        "Moderate",
        "<p>Located in a mutational hot spot.</p><ul><li>Applicable only to "
        "the domains listed in the supplement.</li><li>Not applicable to "
        "specific amino acid residues (see PM5).</li></ul>",
    )
    assert bvc.parse_criterion(trailing)["applicability"] == "applicable"

    mixed = row("Strong", "<p>Def.</p><ul><li>Not applicable at Strong.</li></ul>") \
        + row("Moderate", "<p>Def.</p><ul><li>Met by a rare variant in the "
                          "pore helix (amino acids 300 to 320).</li></ul>")
    assert bvc.parse_criterion(mixed)["applicability"] == "applicable"

    plain = row("Supporting", "<p>Same amino acid change as a previously "
                              "established pathogenic variant.</p>")
    assert bvc.parse_criterion(plain)["applicability"] == "applicable"


def test_per_strength_rules_are_actually_represented():
    """The capability the hand table lacked, and the reason FBN1 PM1 could not be
    implemented from it: a criterion publishing more than one strength."""
    multi = [(g, c) for g in SPEC for c, e in SPEC[g].items()
             if isinstance(e, dict) and len(e.get("strengths") or {}) > 1]
    assert len(multi) > 100, len(multi)
    fbn1_pm1 = _crit_spec("FBN1", "PM1")["strengths"]
    assert set(fbn1_pm1) == {"Strong", "Moderate"}, fbn1_pm1
    assert "cbEGF-like domains" in fbn1_pm1["Strong"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_ps3_ceiling_is_read_from_the_harvest_not_a_gene_list():
    """14 of the 26 spec'd genes publish NO Strong row for PS3 — the whole
    RASopathy VCEP, which is 403 of the 852 cardiac-VCEP records. The engine
    scored every met PS3 at the 2015 base Strong (+4) because it never read
    these rows, so it could award +4 where the VCEP allows at most +2.

    Reading the ceiling from the harvest means a re-harvest updates it with no
    code change, and a gene the harvest does not cover is left alone."""
    from backend.acmg.hard_coded import _ps3_spec_ceiling
    capped = [g for g in SPEC if _ps3_spec_ceiling(g) == "Moderate"]
    assert sorted(capped) == [
        "BRAF", "HRAS", "KRAS", "LZTR1", "MAP2K1", "MAP2K2", "MRAS", "NRAS",
        "PTPN11", "RAF1", "RIT1", "RRAS2", "SOS1", "SOS2",
    ], sorted(capped)
    for gene in ("MYH7", "MYBPC3", "FBN1", "KCNQ1", "BMPR2"):
        assert _ps3_spec_ceiling(gene) == "Strong", gene
    for gene in ("SHOC2", "HEY2", "TTN", None):
        assert _ps3_spec_ceiling(gene) is None, gene


def test_ps3_cap_lowers_but_never_promotes_and_never_invents_a_ceiling():
    from backend.acmg.hard_coded import apply_cross_criterion_exclusions as _ax

    def _ps3(strength="PS3", ev="Functional assay shows loss of function (PMID: 12345678)"):
        return [{"code": "PS3", "status": "met", "criteria_strength": strength,
                 "direction": "pathogenic", "evidence": ev, "source": "ai"}]

    def _one(gene, **kw):
        crits = _ps3(**kw)
        crits[0]["facts"] = {"assay_lab_controls": True}
        out, _ = _ax(crits, gene=gene)
        return next(c for c in out if c["code"] == "PS3")

    for gene in ("BRAF", "MYH7", "HEY2"):
        assert _one(gene)["criteria_strength"] == "PS3_Supporting", gene

    def _with_controls(gene, n):
        crits = _ps3()
        crits[0]["facts"] = {"assay_variant_controls": n}
        out, _ = _ax(crits, gene=gene)
        return next(c for c in out if c["code"] == "PS3")

    assert _with_controls("BRAF", 40)["criteria_strength"] == "PS3_Moderate"
    assert _with_controls("MYH7", 40)["criteria_strength"] == "PS3_Moderate"

    assert _one("BRAF", strength="PS3_Supporting")["criteria_strength"] == "PS3_Supporting"

    sup = _one("BRAF", ev="ClinVar lists this variant as pathogenic")
    assert sup["status"] == "not_met" and sup["criteria_strength"] is None


def test_a_splicing_only_assay_is_not_ps3_evidence():
    """Walker et al. 2023, ClinGen SVI Splicing Subgroup (PMID 37352859): PS3/BS3
    apply "only for well-established assays that measure functional impact not
    directly captured by RNA-splicing assays". Splicing results belong on PVS1 or
    BP7.

    The prompt previously said the OPPOSITE — "only a published RNA / minigene /
    transcript assay of this variant satisfies PS3" — so this is a reversal, and
    it narrows PS3. The evidence is not lost: see
    test_a_splicing_assay_is_routed_to_pvs1_not_dropped for where it goes."""
    from backend.acmg.hard_coded import (
        _is_splicing_only_assay, apply_cross_criterion_exclusions as _ax)

    assert _is_splicing_only_assay("minigene splicing assay showed exon skipping")
    assert _is_splicing_only_assay("rt-pcr of patient rna showed aberrant transcript")
    assert not _is_splicing_only_assay("patch-clamp showed loss of current")
    assert not _is_splicing_only_assay("minigene assay plus western blot of protein")

    def _one(ev, facts=None):
        c = {"code": "PS3", "status": "met", "criteria_strength": "PS3",
             "direction": "pathogenic", "evidence": ev, "source": "ai",
             "facts": {"assay_lab_controls": True}}
        if facts:
            c["facts"] = {**c["facts"], **facts}
        out, _ = _ax([c], gene="MYH7")
        return next(x for x in out if x["code"] == "PS3")

    sup = _one("Minigene splicing assay shows exon skipping (PMID: 12345678)")
    assert sup["status"] == "not_met"
    assert "Walker 2023" in sup["evidence"]
    assert "applied to PVS1 instead" in sup["evidence"]

    sup_facts = _one(
        "Splicing assay, aberrant transcript (PMID: 12345678)",
        facts={"assay_pmids": [12345678], "assay_type": "minigene splicing assay"},
    )
    assert sup_facts["status"] == "not_met", sup_facts["evidence"]

    assert _one("In vitro functional assay of ATPase activity (PMID: 12345678)"
                )["status"] == "met"


def test_the_gate_enforces_not_applicable_for_the_pathogenic_codes():
    """A VCEP specification marking a criterion Not Applicable is binding for
    that gene — that is what a gene-specific specification IS. But the gate only
    enforced five codes, so every other not_applicable row was known and
    ignored: `_criterion_applicable` returned False and nothing acted on it.

    PM3 was the exposure. It is "detected in trans with a pathogenic variant", a
    recessive criterion that 24 of the 26 spec'd genes mark inapplicable, and it
    fired anyway for +2 on all of them.

    Note WHY this went unnoticed: PM3's not_applicable rows were already in the
    hand-transcribed table, so they were never "newly added" and the guard above
    — which only inspects rows the registry ADDED — could not see them."""
    from backend.acmg.hard_coded import _gate_criteria_applicability

    def _fires(gene, code):
        c = {"code": code, "status": "met", "criteria_strength": code,
             "direction": "pathogenic", "evidence": "e", "source": "ai"}
        out = _gate_criteria_applicability([c], gene)
        return next(x for x in out if x["code"] == code)["status"] == "met"

    assert not _fires("MYH7", "PM3"), "24 of 26 specs mark PM3 Not Applicable"
    assert not _fires("MYBPC3", "PM3")
    assert not _fires("BMPR2", "PM6")
    assert not _fires("SHOC2", "PS3")

    assert _fires("MYH7", "PS3")
    assert _fires("HEY2", "PM3"), "no spec covers HEY2, so nothing to enforce"

    for gene, code in (("BRAF", "BS3"), ("HRAS", "BS3"), ("SOS1", "BS3"),
                       ("MYH7", "BP3"), ("MYH7", "PP5"), ("MYH7", "BP6")):
        assert not _fires(gene, code), (
            f"{gene}/{code} is marked Not Applicable and must be suppressed"
        )

    assert _fires("MYH7", "BS3"), "MYH7's spec marks BS3 applicable"
    assert _fires("MYH7", "BS1"), "BS1 is applicable for MYH7"


def _route(*entries, gene="MYH7"):
    from backend.acmg.hard_coded import apply_cross_criterion_exclusions as _ax
    out, _ = _ax(list(entries), gene=gene)
    return {c["code"]: c for c in out}


def _ps3(ev):
    """`facts` carries laboratory controls so the Brnich strength ladder does not
    suppress PS3 before the rule under test is reached."""
    return {"code": "PS3", "status": "met", "criteria_strength": "PS3",
            "direction": "pathogenic", "evidence": ev, "source": "ai",
            "facts": {"assay_lab_controls": True}}


def _bs3(ev):
    return {"code": "BS3", "status": "met", "criteria_strength": "BS3",
            "direction": "benign", "evidence": ev, "source": "ai"}


def _stub(code, direction, **flags):
    c = {"code": code, "status": "not_met", "criteria_strength": None,
         "direction": direction, "evidence": "declined", "source": "engine"}
    c.update(flags)
    return c


def test_a_splicing_assay_is_routed_to_pvs1_not_dropped():
    """Walker 2023: "repurposing the PVS1_Strength code to capture splicing assay
    data that provide experimental evidence for variants resulting in RNA
    transcript(s) with loss of function."

    _eval_pvs1 already declined non-canonical splice variants with the reason
    "canonical +/-1,2 or RNA-confirmed only" — the RNA-confirmed arm was never
    implemented. A suppressed splicing-only PS3 IS that confirmation.

    PVS1_Strong (+4), not full PVS1 (+8): the variant is not at a canonical site
    and the transcript's frame/NMD competence cannot be read from the assay."""
    r = _route(
        _ps3("Minigene assay shows aberrant splicing (PMID: 12345678)"),
        _stub("PVS1", "pathogenic", _pvs1_rna_routable=True),
    )
    assert r["PS3"]["status"] == "not_met"
    assert r["PVS1"]["status"] == "met"
    assert r["PVS1"]["criteria_strength"] == "PVS1_Strong"
    assert "Walker 2023" in r["PVS1"]["evidence"]


def test_the_pvs1_route_needs_the_engine_to_have_flagged_the_variant():
    """The flag is only set where _eval_pvs1 declined for the RNA-confirmable
    reason. Without it — a missense, a non-LoF gene, anything else — a splicing
    assay must not conjure PVS1. Absent flag means do not route."""
    r = _route(
        _ps3("Minigene assay shows aberrant splicing (PMID: 12345678)"),
        _stub("PVS1", "pathogenic"),
    )
    assert r["PVS1"]["status"] == "not_met"


def test_a_protein_function_assay_does_not_trigger_the_pvs1_route():
    """Only SPLICING evidence routes. A protein-function assay keeps its PS3 and
    must leave PVS1 alone, or every assayed variant would collect +4."""
    r = _route(
        _ps3("Patch-clamp shows loss of current (PMID: 12345678)"),
        _stub("PVS1", "pathogenic", _pvs1_rna_routable=True),
    )
    assert r["PS3"]["status"] == "met"
    assert r["PVS1"]["status"] == "not_met"


def test_an_rna_negative_result_is_routed_to_bp7():
    """The benign half: "BP7 may be used to capture RNA results demonstrating no
    splicing impact for intronic and synonymous variants."

    Note the direction. Suppressing BS3 (-4) and applying BP7 (-1) moves the
    variant UP by 3 points, so this route is not conservative — it is what the
    recommendation says, and it is better than dropping the evidence entirely."""
    r = _route(
        _bs3("RT-PCR of patient RNA shows normal splicing (PMID: 12345678)"),
        _stub("BP7", "benign", _bp7_rna_routable=True),
    )
    assert r["BS3"]["status"] == "not_met"
    assert r["BP7"]["status"] == "met"
    assert r["BP7"]["criteria_strength"] == "BP7_Supporting"


def test_bp7_is_not_applied_to_a_variant_class_it_does_not_cover():
    """BP7 is a synonymous / non-coding code. Walker scopes the RNA route to
    intronic and synonymous variants, so a missense must not collect BP7 however
    the assay read out."""
    r = _route(
        _bs3("RT-PCR shows normal splicing (PMID: 12345678)"),
        _stub("BP7", "benign", _bp7_rna_routable=False),
    )
    assert r["BP7"]["status"] == "not_met"


def test_canonical_splice_offset_parsing_is_restricted_to_plus_minus_1_and_2():
    """Scope is what makes this rule implementable. At a canonical +/-1,2
    position, two different changes abolish the SAME splice site, so their
    predicted RNA-splicing effects are equivalent BY CONSTRUCTION and no
    predictor comparison is needed.

    Deeper intronic and exonic splice-region variants are excluded on purpose:
    there, "similar predicted effect" is a genuine prediction that needs SpliceAI
    for the COMPARATOR variant too, and we hold predictions only for the variant
    under assessment."""
    from backend.clients.clinvar import parse_canonical_splice_offset as _p
    assert _p("NM_000257.4:c.732+1G>A") == ("732+1", "G>A")
    assert _p("c.732+2T>C") == ("732+2", "T>C")
    assert _p("NM_000257.4:c.732-1G>A") == ("732-1", "G>A")
    assert _p("c.732+12A>G") is None
    assert _p("c.611G>A") is None
    assert _p(None) is None


def test_ps1_fires_on_a_different_change_at_the_same_canonical_splice_site():
    """MYH7 c.732+1 is the worked example: G>A, G>C, G>T and del are all P/LP at
    >=2*, so a proband carrying one of them has established comparators at the
    same site. The proband's own record is excluded, exactly as on the missense
    PS1 route, or a variant would confirm itself."""
    from .conftest import require_db
    require_db("clinvar")
    import asyncio
    from backend.clients.clinvar import get_pm5_evidence
    from backend.acmg.hard_coded import _clinvar_ps1_criterion

    def _ps1(hgvs):
        ev = asyncio.run(get_pm5_evidence("MYH7", None, proband_hgvs_c=hgvs))
        return _clinvar_ps1_criterion(ev, "MYH7")

    hit = _ps1("NM_000257.4:c.732+1G>A")
    assert hit["status"] == "met"
    assert hit["criteria_strength"] == "PS1"
    assert "Walker 2023" in hit["evidence"]
    assert "c.732+1" in hit["evidence"]
    assert "own ClinVar record" in hit["evidence"], "self-exclusion must be shown"

    miss = _ps1("NM_000257.4:c.5655+1G>C")
    assert miss["status"] == "not_met"
    assert "splice route" in miss["evidence"]

    oos = _ps1("NM_000257.4:c.732+12A>G")
    assert oos["status"] == "not_met"
    assert "not evaluated" in oos["evidence"]


def test_the_splice_ps1_bar_matches_the_missense_ps1_bar():
    """Established P/LP at >=2* review, and "Conflicting classifications of
    pathogenicity" must NOT count even though the string contains "pathogenic"."""
    from .conftest import require_db
    require_db("clinvar")
    from backend.clients.clinvar import _splice_ps1_sync
    r = _splice_ps1_sync("MYH7", "732+1", "G>A")
    assert r["ok"] and r["count"] >= 2
    assert all(c["stars"] >= 2 for c in r["candidates"])
    assert all("conflicting" not in c["significance"].lower() for c in r["candidates"])
    assert all("benign" not in c["significance"].lower() for c in r["candidates"])
    assert any("c.732+1G>A" in c["name"] for c in r["self_excluded"])


def test_only_kcnq1_needs_a_gene_specific_pvs1_tree():
    """The survey that bounds this problem, asserted rather than trusted.

    AutoPVS1 (Xiang et al., Hum Mutat 2020, doi:10.1002/humu.24051), the
    reference implementation of the SVI PVS1 tree at 95% concordance, does not
    parse specifications — it implements the generic tree and hardcodes the few
    genes with published gene-specific PVS1 recommendations. This test pins why
    the same approach is sufficient here.

    PVS1 is applicable for only 5 of 26 genes. Three of those defer to the
    generic tree in their own words, FBN1's rules ARE the generic tree restated,
    and KCNQ1 is the only one with numeric content of its own."""
    applicable = sorted(
        g for g in SPEC if (SPEC[g].get("PVS1") or {}).get("applicability") == "applicable"
    )
    assert applicable == ["BMPR2", "FBN1", "KCNQ1", "LZTR1", "MYBPC3"], applicable

    def _rows(gene):
        return " ".join((SPEC[gene]["PVS1"].get("strengths") or {}).values()).lower()

    assert "pvs1 decision tree guide" in _rows("BMPR2")
    assert "svi guidance" in _rows("MYBPC3") or "abou tayoun" in _rows("MYBPC3")
    assert "null variant in a gene where loss of function" in _rows("LZTR1")

    kcnq1 = _rows("KCNQ1")
    for token in ("1-581", "582 and 620", "621 and 676", "589-620"):
        assert token in kcnq1, token


def test_kcnq1_pvs1_bands_match_its_published_codon_ranges():
    """GN112: VeryStrong for codons 1-581 (NMD predicted), Moderate for 582-620
    (NMD not predicted but the SAD at 589-620 is removed, so the channel cannot
    tetramerize), Supporting for 621-676 (SAD retained, may still yield
    functional channels)."""
    from backend.acmg.hard_coded import _pvs1_kcnq1_strength

    def _tier(codon, consequence="stop_gained"):
        r = _pvs1_kcnq1_strength({"most_severe_consequence": consequence,
                                  "protein_position": str(codon),
                                  "gene_symbol": "KCNQ1"})
        return r[0] if r else None

    assert _tier(1) == "PVS1" and _tier(581) == "PVS1"
    assert _tier(582) == "PVS1_Moderate" and _tier(620) == "PVS1_Moderate"
    assert _tier(621) == "PVS1_Supporting" and _tier(676) == "PVS1_Supporting"

    assert _tier(677) is None
    assert _tier(700) is None

    assert _tier(300, consequence="missense_variant") is None
    assert _pvs1_kcnq1_strength(
        {"most_severe_consequence": "stop_gained", "protein_position": None}) is None
    assert _tier("589-591", consequence="frameshift_variant") == "PVS1_Moderate"


def test_no_criterion_base_strength_exceeds_its_published_max():
    """The guard that makes the PS3 and PM2 defects a class, not two incidents.

    Both were the engine emitting a criterion at the ACMG 2015 BASE strength on
    genes whose VCEP publishes a lower maximum: PS3 at Strong where 14 genes cap
    it at Moderate, and PM2 at Moderate where ALL 26 publish Supporting only.

    `cap_to_spec_strength` now caps every deterministic criterion against the
    harvest, so today's answer is that nothing over-calls. This test states that
    as a property: if a re-harvest lowers some criterion's published maximum
    below our base strength, and the cap were ever removed or bypassed, this
    fails rather than the engine quietly over-calling."""
    from backend.acmg.hard_coded import (
        _BARE_CODE_POINTS, _SPEC_STRENGTH_RANK, spec_strength_ceiling,
    )
    offenders = []
    for code, pts in _BARE_CODE_POINTS.items():
        base = abs(pts)
        for gene in SPEC:
            ceiling = spec_strength_ceiling(gene, code)
            if ceiling and base > _SPEC_STRENGTH_RANK[ceiling]:
                offenders.append((gene, code, base, ceiling))
    codes = sorted({c for _, c, _, _ in offenders})
    assert codes == ["BS3", "PM1", "PM2", "PS3"], (
        f"the set of criteria whose base exceeds a gene's published maximum has "
        f"changed: {sorted(set(offenders))}"
    )
    from backend.acmg.hard_coded import cap_to_spec_strength
    for gene, code, _base, ceiling in offenders:
        e = {"code": code, "status": "met", "criteria_strength": code,
             "direction": "benign" if code.startswith("B") else "pathogenic",
             "evidence": "x", "source": "hard_coded"}
        got = cap_to_spec_strength(e, gene)["criteria_strength"]
        assert got == f"{code}_{ceiling}", (gene, code, got, ceiling)


def test_the_cap_lowers_pm2_to_the_supporting_every_spec_publishes():
    from backend.acmg.hard_coded import cap_to_spec_strength, spec_strength_ceiling

    def _pm2(gene):
        e = {"code": "PM2", "status": "met", "criteria_strength": "PM2",
             "direction": "pathogenic", "evidence": "rare", "source": "hard_coded"}
        return cap_to_spec_strength(e, gene)

    for gene in ("MYH7", "BRAF", "BMPR2", "FBN1", "KCNQ1"):
        assert spec_strength_ceiling(gene, "PM2") == "Supporting", gene
        got = _pm2(gene)
        assert got["criteria_strength"] == "PM2_Supporting", gene
        assert "publishes PM2 at Supporting only" in got["evidence"]

    assert spec_strength_ceiling("HEY2", "PM2") is None
    assert _pm2("HEY2")["criteria_strength"] == "PM2"

    weak = {"code": "PM2", "status": "met", "criteria_strength": "PM2_Supporting",
            "direction": "pathogenic", "evidence": "keep me", "source": "hard_coded"}
    out = cap_to_spec_strength(weak, "MYH7")
    assert out["criteria_strength"] == "PM2_Supporting" and out["evidence"] == "keep me"

    nm = {"code": "PM2", "status": "not_met", "criteria_strength": None,
          "direction": "pathogenic", "evidence": "absent", "source": "hard_coded"}
    assert cap_to_spec_strength(nm, "MYH7")["criteria_strength"] is None


def test_the_bp7_rna_route_does_not_waive_bp7s_conservation_requirement():
    """A false-benign risk I introduced and then fixed.

    Per ACMG/AMP, BP7 needs THREE conjuncts: an eligible consequence, no splice
    impact, AND the nucleotide not highly conserved (phyloP100way <= 2.0 — the
    MM-VCEP operationalisation _eval_bp7 documents). Walker 2023's RNA route
    substitutes an experimental result for the SPLICE-IMPACT half only.

    The first version of the routing marker tested the consequence alone, so a
    splicing-only BS3 could set BP7 met on a HIGHLY CONSERVED nucleotide, or
    where phyloP was unavailable — exactly the two cases _eval_bp7 fails closed
    on because firing ungated "risks a false-benign call"."""
    from backend.acmg.hard_coded import compute_hard_coded_criteria

    def _routable(phylop, consequence="intron_variant"):
        ev = {"vep": {"ok": True, "most_severe_consequence": consequence,
                      "gene_symbol": "MYH7", "phylop100way": phylop}}
        out = compute_hard_coded_criteria(ev, {"inheritance_input": "AD"}, "MYH7")
        bp7 = next(c for c in out if c["code"] == "BP7")
        return bp7.get("_bp7_rna_routable")

    assert _routable(0.5) is True
    assert _routable(2.0) is True, "2.0 is the inclusive boundary"

    assert _routable(2.1) is False
    assert _routable(8.0) is False

    assert _routable(None) is False

    assert _routable(0.5, consequence="missense_variant") is False
