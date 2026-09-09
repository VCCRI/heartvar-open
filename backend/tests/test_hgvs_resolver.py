"""Guards on clients/hgvs_resolver.py — local HGVS → genomic resolution.

WHAT THIS MODULE EXISTS TO FIX, so the tests are read in context: `vep --offline`
refuses `--format hgvs` (VEP 113 Parser/HGVS.pm:104, unconditional, verified in
the official image), so the offline cache could only serve COORDINATE input while
curators type HGVS. Every HGVS curation went to rest.ensembl.org — the shared
quota whose exhaustion would take the site down for everyone.

The resolver's contract is narrow and these tests hold it to it:
  * it either resolves or returns None — it NEVER raises;
  * it emits TWO representations, because they have different consumers;
  * it never silently substitutes a transcript the curator did not ask for.

The `hgvs`/`cdot` stack is heavy and needs ~104 MB of transcript JSON, so the
arithmetic-level tests stub the mapper. The one test that exercises the real
library is skipped unless HEARTVAR_CDOT_FILES points at real data — the full
end-to-end parity gate lives outside the suite (it needs the VEP cache too).
"""
from __future__ import annotations

import os
import pathlib
import re

import pytest

from backend.clients import hgvs_resolver as hr


def test_the_vep_row_is_tab_delimited():
    """⚠ THE SILENT FAILURE THIS SHARES WITH fetch_region. Parser/VCF.pm splits
    on TAB. A space-separated row gives exit 0, ZERO bytes of stdout and seven
    "uninitialized value $ref" warnings — measured — and _run only logs stderr on
    a non-zero status, so it is invisible. The row is built inside Resolved so
    no caller can get it wrong."""
    r = hr.Resolved(chrom="14", vep_pos=23429278, vep_ref="C", vep_alt="T",
                    vcf_pos=23429278, vcf_ref="C", vcf_alt="T",
                    accession="NM_000257.4")
    row = r.vep_row()
    assert row.split("\t") == ["14", "23429278", ".", "C", "T", ".", ".", "."]
    assert " " not in row, f"a space survives in the VEP row: {row!r}"


def test_an_empty_allele_becomes_a_dash_for_vep_not_an_empty_field():
    """VEP's own convention for a deletion is `CCG/-`. An empty field would shift
    every column after it and be parsed as a malformed row."""
    r = hr.Resolved(chrom="14", vep_pos=23429277, vep_ref="CCG", vep_alt="",
                    vcf_pos=23429276, vcf_ref="CCCG", vcf_alt="C",
                    accession="NM_000257.4")
    assert r.vep_row().split("\t")[4] == "-"


def test_the_gnomad_key_uses_the_ANCHORED_form_not_the_vep_form():
    """THE WHOLE REASON Resolved carries two tuples. gnomAD stores anchored
    left-aligned keys, and the VEP CLI does NOT emit vcf_string (REST-only —
    verified: with and without --var_synonyms the top-level keys are only
    allele_string, assembly_name, end, id, input, most_severe_consequence,
    seq_region_name, start, strand, transcript_consequences). So if the gnomAD key
    were built from the VEP form, every indel would get a key that per
    ensembl_vep.gnomad_variant_id_with_provenance "can never match" — and PM2
    reads absence as evidence, so a meaningless miss is worse than no lookup."""
    r = hr.Resolved(chrom="14", vep_pos=23429277, vep_ref="CCG", vep_alt="",
                    vcf_pos=23429276, vcf_ref="CCCG", vcf_alt="C",
                    accession="NM_000257.4")
    assert r.gnomad_id() == "14-23429276-CCCG-C"
    assert "23429277" not in r.gnomad_id()


class _Edit:
    def __init__(self, name, ref=None, alt=None):
        self.__class__ = type(name, (_Edit,), {})
        self.ref = ref
        if alt is not None:
            self.alt = alt


def _edit(name, ref=None, alt=None):
    obj = _Edit.__new__(type(name, (object,), {}))
    obj.ref = ref
    if alt is not None:
        obj.alt = alt
    return obj


@pytest.mark.parametrize("kind,ref,alt,expected", [
    ("NARefAlt", "C", "T", ("C", "T")),
    ("NARefAlt", "CCG", "", ("CCG", "")),
    ("NARefAlt", "CCG", "TTT", ("CCG", "TTT")),
    ("Dup", "G", None, ("G", "GG")),
    ("Inv", "CCG", None, ("CCG", "CGG")),
])
def test_every_edit_class_yields_a_ref_and_alt(kind, ref, alt, expected):
    """⚠ THE PROTOTYPE HANDLED ONE CLASS. It read edit.ref/edit.alt directly,
    which does not exist on Dup or Inv: `c.1208dup` and `c.1207_1209inv` both
    raise `AttributeError: 'Dup' object has no attribute 'alt'` — measured
    against the real library. It went unnoticed because the inputs it had been
    exercised on are all NARefAlt, so neither crashing class was ever
    reached. A crash is not acceptable when the design rule is "degrade, never
    block"."""
    assert hr._ref_alt_from_edit(_edit(kind, ref, alt)) == expected


def test_an_unknown_edit_class_returns_None_rather_than_raising():
    """The module's contract is resolve-or-None. A future hgvs release adding an
    edit class must degrade to REST, not take the curation down."""
    assert hr._ref_alt_from_edit(_edit("Repeat", "CAG")) is None


def test_inversion_is_reverse_complemented_not_merely_complemented():
    """An inversion is revcomp, and getting it wrong is silent: the alt is still
    the right length and still valid DNA."""
    assert hr._revcomp("CCG") == "CGG"
    assert hr._revcomp("AAAT") == "ATTT"


@pytest.mark.parametrize("acc,expected", [
    ("NC_000014.9", "14"),
    ("NC_000001.11", "1"),
    ("NC_000023.11", "X"),
    ("NC_000024.10", "Y"),
    ("NC_012920.1", "MT"),
])
def test_nc_accessions_map_to_the_chromosome_names_the_cache_uses(acc, expected):
    """The VEP cache and parse_variant_input both use bare 1-22/X/Y/MT (the
    FASTA header is `>14`), so the mapping has to land there exactly."""
    assert hr._nc_to_chrom(acc) == expected


@pytest.mark.parametrize("acc", [
    "NW_009646201.1",
    "NT_187633.1",
    "NC_000025.1",
    "", None, "garbage",
])
def test_a_non_chromosomal_alignment_is_a_failure_not_a_guess(acc):
    """A transcript that aligns to a scaffold must resolve to None. Returning a
    plausible-looking chromosome would put ACMG on the wrong locus."""
    assert hr._nc_to_chrom(acc) is None


class _FakeFasta:
    """Reference sequence for one contig, 1-based inclusive coordinates mapped to
    pysam's 0-based half-open fetch."""

    def __init__(self, seq: str, offset: int = 1):
        self._seq, self._offset = seq, offset

    def fetch(self, _chrom, start0, end0):
        i = start0 - (self._offset - 1)
        j = end0 - (self._offset - 1)
        if i < 0 or j > len(self._seq):
            return ""
        return self._seq[i:j]


def test_a_substitution_is_already_anchored_and_is_left_alone():
    fa = _FakeFasta("ACGTACGT", offset=100)
    assert hr._left_anchor(fa, "14", 103, "T", "A") == (103, "T", "A")


def test_a_deletion_gains_the_preceding_base_and_moves_left():
    """HGVS deletes bases; VCF anchors on the base before. This is the
    representation difference behind any left-anchored deletion, and it shifts
    the POSITION — so an unanchored key misses gnomAD entirely."""
    fa = _FakeFasta("ACCGT", offset=100)
    assert hr._left_anchor(fa, "14", 101, "CCG", "") == (100, "ACCG", "A")


def test_an_insertion_gains_the_preceding_base():
    fa = _FakeFasta("ACCGT", offset=100)
    assert hr._left_anchor(fa, "14", 102, "", "TT") == (101, "C", "CTT")


def test_a_shared_suffix_is_trimmed_and_the_result_is_fully_left_aligned():
    """Untrimmed alleles produce a different key for the same variant, and
    half-trimmed ones produce a different key again.

    Reference ACCGT at 100-104. `101 CCGT>CGT` deletes one C from the CC at
    101-102. Three representations all describe that:
        101 CCGT CGT   (as supplied, untrimmed)
        101 CC   C     (suffix trimmed, not shifted)
        100 AC   A     (fully left-aligned + anchored)  <- the only one gnomAD has
    Only the last matches a gnomAD key, which is the entire point of computing it.
    """
    fa = _FakeFasta("ACCGT", offset=100)
    assert hr._left_anchor(fa, "14", 101, "CCGT", "CGT") == (100, "AC", "A")


def test_anchoring_at_the_contig_start_fails_rather_than_reading_off_the_end():
    fa = _FakeFasta("ACCGT", offset=1)
    assert hr._left_anchor(fa, "14", 1, "A", "") is None


def test_no_reference_means_no_anchored_key():
    """Without the FASTA an indel cannot get a meaningful gnomAD key, and the
    resolver must say so rather than emit the unanchored form."""
    assert hr._left_anchor(None, "14", 101, "CCG", "") is None


def test_disabled_by_default(monkeypatch):
    monkeypatch.delenv("HEARTVAR_HGVS_RESOLVER", raising=False)
    assert hr.resolver_enabled() is False
    assert hr.resolve("c.1208G>A", "NM_000257.4") is None


@pytest.mark.parametrize("val", ["1", "true", "TRUE", "yes", "on"])
def test_enabled_by_truthy_values(monkeypatch, val):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", val)
    assert hr.resolver_enabled() is True


def test_an_empty_hgvs_never_resolves(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    assert hr.resolve("", "NM_000257.4", "MYH7") is None


def test_a_bare_gene_hgvs_resolves_against_MANE_SELECT(monkeypatch):
    """A designation, not a guess. Verified against the failing production log:
    CHD7's local MANE Select is ENST00000423902.7, the same transcript REST
    itself selected before it started 500ing."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "/fake/cdot.db")
    from backend.clients import cdot_sqlite
    monkeypatch.setattr(
        cdot_sqlite, "mane_select_for_gene",
        lambda _db, _gene, *a: cdot_sqlite.ManeSelect("ENST00000423902.7",
                                                      "NM_017780.4"))

    accession, note = hr._mane_select_for_gene("CHD7")
    assert accession == "ENST00000423902.7"
    assert "MANE Select" in note and "CHD7" in note


def test_the_ensembl_accession_is_preferred_over_the_refseq_one(monkeypatch):
    """⚠ A COMPATIBILITY CONSTRAINT, not a clinical preference. REST queried
    bare-gene HGVS with VEP_PARAMS, not _REFSEQ_PARAMS, so it returned
    ENST/ENSP. vep_offline keys `refseq=` off this accession, so handing back
    NM_ here would flip the transcript set on every transcript-less curation —
    a change measured once before at hgvsp 89/101 mismatched, which feeds
    PS1/PM5 residue matching."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "/fake/cdot.db")
    from backend.clients import cdot_sqlite
    monkeypatch.setattr(
        cdot_sqlite, "mane_select_for_gene",
        lambda _db, _gene, *a: cdot_sqlite.ManeSelect("ENST00000423902.7",
                                                      "NM_017780.4"))
    assert hr._mane_select_for_gene("CHD7")[0].startswith("ENST")


def test_the_refseq_accession_is_used_when_there_is_no_ensembl_one(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "/fake/cdot.db")
    from backend.clients import cdot_sqlite
    monkeypatch.setattr(
        cdot_sqlite, "mane_select_for_gene",
        lambda _db, _gene, *a: cdot_sqlite.ManeSelect(None, "NM_017780.4"))
    assert hr._mane_select_for_gene("CHD7")[0] == "NM_017780.4"


def test_a_gene_with_no_MANE_SELECT_still_goes_to_REST(monkeypatch):
    """The fallback that made the old behaviour safe is retained for exactly
    the case it was reasoning about."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "/fake/cdot.db")
    from backend.clients import cdot_sqlite
    monkeypatch.setattr(
        cdot_sqlite, "mane_select_for_gene",
        lambda _db, _gene, *a: cdot_sqlite.ManeSelect(None, None))
    assert hr._mane_select_for_gene("OBSCURE1") == (None, None)
    assert hr.resolve("c.1208G>A", None, "OBSCURE1") is None


def test_no_cdot_database_means_REST_not_a_crash(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "")
    assert hr._mane_select_for_gene("CHD7") == (None, None)
    assert hr.resolve("c.1208G>A", None, "CHD7") is None


def test_a_lookup_that_explodes_falls_back_to_REST(monkeypatch):
    """Contract: resolve-or-None. A broken mount must not take the curation
    down — it must produce the pre-2026-08-31 behaviour."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.setattr(hr, "_db_path", lambda: "/fake/cdot.db")
    from backend.clients import cdot_sqlite

    def _boom(*_a, **_k):
        raise RuntimeError("mount dropped")

    monkeypatch.setattr(cdot_sqlite, "mane_select_for_gene", _boom)
    assert hr._mane_select_for_gene("CHD7") == (None, None)
    assert hr.resolve("c.1208G>A", None, "CHD7") is None


def test_a_bare_gene_with_no_gene_at_all_resolves_to_nothing(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    assert hr._mane_select_for_gene(None) == (None, None)
    assert hr._mane_select_for_gene("") == (None, None)
    assert hr.resolve("c.1208G>A", None, None) is None


def test_resolve_never_raises_when_the_mapper_explodes(monkeypatch):
    """Contract: resolve-or-None. A library exception is a REST fallback, not a
    500 and not a dead curation."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")

    class _Boom:
        def parse_hgvs_variant(self, _v):
            raise RuntimeError("kaboom")

    monkeypatch.setattr(hr, "_providers", lambda: (None, _Boom(), None))
    assert hr.resolve("c.1208G>A", "NM_000257.4") is None


def test_missing_cdot_files_degrade_to_REST(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.delenv("HEARTVAR_CDOT_DB", raising=False)
    monkeypatch.setenv("HEARTVAR_CDOT_FILES", "/nonexistent/a.json.gz")
    hr._providers.cache_clear()
    try:
        assert hr._providers() is None
        assert hr.resolve("c.1208G>A", "NM_000257.4") is None
    finally:
        hr._providers.cache_clear()


def test_a_version_substitution_is_reported_not_silent(monkeypatch):
    """A curator typing `NM_000257` gets resolved against a specific version, and
    that choice changes the coordinates. It has to be visible: the REST path
    deliberately strips versions ("RefSeq versions are frequently stale"), and
    cdot REQUIRES one — `NM_000257:c.1208G>A` raises HGVSDataNotAvailableError
    while the stale `NM_000257.3` resolves fine, both measured."""
    monkeypatch.setattr(hr, "_versions_by_accession",
                        lambda: {"NM_000257": ["NM_000257.3", "NM_000257.4"]})
    acc, note = hr._pin_version("NM_000257")
    assert acc == "NM_000257.4", "must take the newest version available"
    assert note and "without a version" in note and "NM_000257.4" in note


def test_an_explicit_version_is_never_second_guessed(monkeypatch):
    """A stale explicit version still resolves correctly (measured: NM_000257.3
    gives the right coordinates), so silently upgrading it would change the
    transcript the curator asked about."""
    monkeypatch.setattr(hr, "_versions_by_accession",
                        lambda: {"NM_000257": ["NM_000257.3", "NM_000257.4"]})
    assert hr._pin_version("NM_000257.3") == ("NM_000257.3", None)


def test_refseq_detection_matches_the_rest_clients_notion(monkeypatch):
    for acc in ("NM_000257.4", "NR_1.1", "XM_1.1", "XR_1.1"):
        assert hr._is_refseq(acc) is True
    for acc in ("ENST00000355349.4", "", None):
        assert hr._is_refseq(acc) is False


@pytest.mark.skipif(not (os.environ.get("HEARTVAR_CDOT_DB")
                         or os.environ.get("HEARTVAR_CDOT_FILES")),
                    reason="needs real cdot transcript data")
def test_the_canonical_variant_resolves_against_real_data(monkeypatch):
    """MYH7 R403Q, the canonical HCM founder missense and the one variant every
    cardiac tool is expected to get right. 14:23,429,278 C>T on GRCh38."""
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    hr._providers.cache_clear()
    try:
        r = hr.resolve("c.1208G>A", "NM_000257.4", "MYH7")
    finally:
        hr._providers.cache_clear()
    assert r is not None
    assert (r.chrom, r.vep_pos, r.vep_ref, r.vep_alt) == ("14", 23429278, "C", "T")
    assert r.gnomad_id() == "14-23429278-C-T"


def test_the_db_is_preferred_over_the_json_files(monkeypatch):
    """⚠ THE JSON ROUTE COSTS 4480 MB RSS AGAINST A 4 GiB CONTAINER.

    Measured 2026-08-28: cdot's JSONDataProvider parses the whole transcript set
    (906,754 transcripts) into memory, so the resolver was not deployable at all.
    The SQLite provider serves the same records a row at a time — 46 MB, the same
    0.19 ms resolution, and the SAME FULL COVERAGE.

    Panel-scoping the JSON was rejected: it got memory to 498 MB but dropped
    any off-panel variant to REST. A gene being outside our panel is not a
    reason to annotate it worse.
    """
    monkeypatch.setenv("HEARTVAR_CDOT_DB", "/some/cdot.db")
    monkeypatch.setenv("HEARTVAR_CDOT_FILES", "/a.json.gz,/b.json.gz")
    assert hr._db_path() == "/some/cdot.db"

    calls = {}

    def _fake_make(path):
        calls["db"] = path
        return object()

    import backend.clients.cdot_sqlite as cs
    monkeypatch.setattr(cs, "make_provider", _fake_make)
    monkeypatch.setattr(hr, "_versions_by_accession", lambda: {})
    hr._providers.cache_clear()
    try:
        hr._providers()
    finally:
        hr._providers.cache_clear()
    assert calls.get("db") == "/some/cdot.db", (
        "the JSON path was taken even though HEARTVAR_CDOT_DB is set; that is a "
        "4.5 GB memory regression"
    )


def test_no_cdot_source_at_all_degrades_to_REST(monkeypatch):
    monkeypatch.setenv("HEARTVAR_HGVS_RESOLVER", "1")
    monkeypatch.delenv("HEARTVAR_CDOT_DB", raising=False)
    monkeypatch.delenv("HEARTVAR_CDOT_FILES", raising=False)
    hr._providers.cache_clear()
    try:
        assert hr._providers() is None
        assert hr.resolve("c.1208G>A", "NM_000257.4") is None
    finally:
        hr._providers.cache_clear()


def test_the_sqlite_provider_refuses_region_queries_rather_than_faking_them():
    """get_tx_for_region needs an interval tree over every transcript, which is
    the resident-memory cost the provider exists to avoid. Returning an empty
    tree would report "no transcripts in this region" as a fact — nothing in
    HeartVar calls it, so raising is the honest answer."""
    import inspect

    from backend.clients import cdot_sqlite

    src = inspect.getsource(cdot_sqlite)
    assert "_get_contig_interval_tree" in src
    assert "NotImplementedError" in src, (
        "the provider fakes region queries instead of refusing them"
    )


def test_the_rss_ceiling_is_shared_with_the_builder():
    """The number has to live in ONE place. The builder asserts it in the real
    container against the real mount; this file asserts it in the dev loop. Two
    copies would drift, and the first sign would be a 4 GiB container OOMing."""
    from backend.clients import cdot_sqlite

    assert isinstance(cdot_sqlite.RSS_CEILING_MB, int)
    assert cdot_sqlite.RSS_CEILING_MB > 0
    builder = (_ROOT_PATH / "scripts" / "build_cdot_transcripts.sh").read_text(
        encoding="utf-8")
    assert "RSS_CEILING_MB" in builder and "current_rss_mb" in builder, (
        "the builder no longer asserts the resident-memory ceiling — that is the "
        "only check that runs with the real data at the real memory limit"
    )


def test_ru_maxrss_units_are_handled_per_platform():
    """ru_maxrss is BYTES on Darwin and KILOBYTES on Linux. Getting it backwards
    makes the guard either never fire or always fire — and 'never fire' is
    indistinguishable from 'no problem', which is how the 4.5 GB shipped."""
    from backend.clients import cdot_sqlite

    rss = cdot_sqlite.current_rss_mb()
    assert 1 < rss < 20000, f"current_rss_mb() returned {rss}, units look wrong"


@pytest.mark.skipif(not os.environ.get("HEARTVAR_CDOT_DB"),
                    reason="needs the real cdot database")
def test_resolving_with_real_data_stays_under_the_ceiling():
    """THE CHECK WHOSE ABSENCE LET AN UNDEPLOYABLE RESOLVER MERGE.

    #41 was correct and could not run: cdot's JSONDataProvider loads all 906,754
    transcripts into RAM — 4480 MB measured, against a 4 GiB container. Every
    test passed, because the unit tests stub the mapper and the real-data test
    ran in a process with the whole machine to itself. Nothing in the suite had
    any notion of resource cost.

    Verified to FIRE on the regression: forcing the JSON provider gives 4358 MB
    against this 512 MB ceiling.
    """
    from backend.clients import hgvs_resolver as hr
    from backend.clients.cdot_sqlite import RSS_CEILING_MB, current_rss_mb

    _prior = os.environ.get("HEARTVAR_HGVS_RESOLVER")
    os.environ["HEARTVAR_HGVS_RESOLVER"] = "1"
    hr._providers.cache_clear()
    try:
        r = hr.resolve("c.1208G>A", "NM_000257.4", "MYH7")
    finally:
        hr._providers.cache_clear()
        if _prior is None:
            os.environ.pop("HEARTVAR_HGVS_RESOLVER", None)
        else:
            os.environ["HEARTVAR_HGVS_RESOLVER"] = _prior
    assert r is not None, "the real database did not resolve MYH7 R403Q"
    rss = current_rss_mb()
    assert rss <= RSS_CEILING_MB, (
        f"resolving with real data used {rss:.0f} MB against a "
        f"{RSS_CEILING_MB} MB ceiling. ~4.5 GB means the in-memory JSON provider "
        f"was used — check HEARTVAR_CDOT_DB points at an existing .db"
    )


_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_ROOT_PATH = pathlib.Path(_ROOT)


def _requirements() -> str:
    with open(os.path.join(_ROOT, "requirements.txt"), encoding="utf-8") as fh:
        return fh.read()


def _dockerfile() -> str:
    with open(os.path.join(_ROOT, "Dockerfile"), encoding="utf-8") as fh:
        return fh.read()


def test_hgvs_and_cdot_are_installed_no_deps_and_not_via_requirements():
    """⚠ MOVING EITHER ONE INTO requirements.txt BREAKS THE IMAGE BUILD.

    `hgvs` requires `psycopg2` BY NAME, and `psycopg2-binary` does not satisfy
    that requirement — measured: `pip install hgvs` with psycopg2-binary already
    present still tries to BUILD psycopg2 and fails with
    "Failed to build 'psycopg2' when getting requirements to build wheel"
    (it wants pg_config/libpq). `cdot` depends on `hgvs`, so it drags the same
    resolution in.

    The arrangement: both --no-deps in the Dockerfile, every real dependency
    pinned in requirements.txt. This test is what stops someone "tidying" them
    into requirements.txt and discovering it at deploy time.
    """
    req, dockerfile = _requirements(), _dockerfile()
    assert "--no-deps hgvs==" in dockerfile and "cdot==" in dockerfile, (
        "the Dockerfile no longer installs hgvs/cdot with --no-deps"
    )
    for pkg in ("hgvs", "cdot"):
        for line in req.splitlines():
            line = line.split("#")[0].strip()
            if line.lower().startswith((f"{pkg}==", f"{pkg}>=", f"{pkg}~=")):
                raise AssertionError(
                    f"{pkg} is pinned in requirements.txt; pip will then resolve "
                    "its dependency on psycopg2 (by name) and try to build it "
                    "from source. Install it --no-deps in the Dockerfile instead."
                )


def test_psycopg2_binary_is_used_and_plain_psycopg2_is_not():
    """psycopg2-binary provides the importable module without needing libpq or a
    compiler. Plain psycopg2 would require both, to talk to a database this app
    never connects to (UTA is unused — the transcript data is cdot JSON)."""
    req = _requirements()
    assert "psycopg2-binary==" in req
    for line in req.splitlines():
        line = line.split("#")[0].strip()
        assert not line.startswith("psycopg2=="), (
            "plain psycopg2 pulls in a source build; use psycopg2-binary"
        )


def test_every_real_dependency_of_the_hgvs_stack_is_pinned():
    """Because they are installed --no-deps, pip will NOT fetch these for us. A
    missing one is an ImportError at first resolve, in production, on the one
    code path that was supposed to make the site more robust."""
    req = _requirements()
    for dep in ("attrs", "bioutils", "configparser", "importlib_resources",
                "parsley", "biocommons.seqrepo", "intervaltree",
                "more-itertools", "lazy", "msgspec", "requests", "ipython"):
        assert f"{dep}==" in req, (
            f"{dep} is a real runtime dependency of the hgvs/cdot stack but is "
            "not pinned in requirements.txt, and --no-deps means pip will not "
            "add it"
        )


def _builder_dockerfile() -> str:
    with open(os.path.join(_ROOT, "Dockerfile.builder"), encoding="utf-8") as fh:
        return fh.read()


def _no_deps_pins(dockerfile: str) -> dict[str, str]:
    """The hgvs/cdot versions a Dockerfile installs --no-deps, as {pkg: version}."""
    out: dict[str, str] = {}
    for m in re.finditer(r"--no-deps\s+([^\n\\]+)", dockerfile):
        for tok in m.group(1).split():
            if "==" in tok:
                pkg, _, ver = tok.partition("==")
                out[pkg.strip()] = ver.strip()
    return out


def test_the_builder_image_installs_the_hgvs_stack_too():
    """⚠ THE BUILDER RUNS THE RESOLVER, SO IT MUST BE ABLE TO IMPORT IT.

    scripts/build_cdot_transcripts.sh ends by calling
    ``hgvs_resolver.resolve()`` on MYH7 R403Q and then asserting the RSS
    ceiling. Both run in the DATA-BUILDER image, not the runtime one. Measured
    2026-08-29 in heartvar-builder: the step built the database correctly
    (906,754 transcripts, 608 MB) and then

        could not build the HGVS resolver (ModuleNotFoundError:
        No module named 'hgvs')
            RESOLVER RETURNED None for MYH7 R403Q
        SCRIPT_EXIT=1

    Two consequences, and the second is the one that matters:
      * `cdot` reports FAILED on every build even though its output is fine.
      * THE RSS GUARD CAN NEVER RUN. It is placed in that script precisely
        because that is "the only place that has the real data, in the real
        container, at the real memory limit" — but resolve() returns None
        first, so the guard added to catch the 4.5 GB regression exits before
        reaching it. A guard that cannot execute is not a guard.
    """
    builder = _builder_dockerfile()
    assert "--no-deps" in builder and "hgvs==" in builder and "cdot==" in builder, (
        "Dockerfile.builder does not install the hgvs/cdot stack, so "
        "build_cdot_transcripts.sh cannot verify the resolver or measure its "
        "RSS — see backend/tests/test_hgvs_resolver.py for why that matters"
    )


def test_builder_and_runtime_pin_the_same_hgvs_stack():
    """A skew here would have the builder VERIFY one resolver and production RUN
    another. The RSS measurement and the R403Q coordinate check are only
    evidence about the deployment if both images hold the same code."""
    runtime_pins = _no_deps_pins(_dockerfile())
    builder_pins = _no_deps_pins(_builder_dockerfile())
    for pkg in ("hgvs", "cdot"):
        assert pkg in runtime_pins, f"Dockerfile no longer pins {pkg} --no-deps"
        assert pkg in builder_pins, f"Dockerfile.builder no longer pins {pkg} --no-deps"
        assert runtime_pins[pkg] == builder_pins[pkg], (
            f"{pkg} is {runtime_pins[pkg]} in the runtime image but "
            f"{builder_pins[pkg]} in the builder image. The builder's resolver "
            "check would then be evidence about a different resolver than the "
            "one production runs."
        )


def test_the_image_verifies_the_stack_imports_at_build_time():
    """A --no-deps install that is missing something fails at first USE, which
    for this module is a live curation. The Dockerfile imports the stack in the
    same layer so the build fails instead."""
    dockerfile = _dockerfile()
    assert "import hgvs.parser" in dockerfile and "psycopg2" in dockerfile, (
        "the Dockerfile no longer smoke-imports the hgvs stack, so a broken "
        "--no-deps install would ship"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
