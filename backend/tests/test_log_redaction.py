"""The curated gene and variant must not reach a log record.

The About page tells curators their input is "not written to a database, logged
with their content, or retained after the request completes". Sixty log call
sites across eighteen modules used to interpolate the gene, the HGVS string, the
gnomAD id or the coordinates, and the deployment log retains every record for
the life of the server process. These tests hold the fix in place.
"""
import io
import logging

import pytest

from backend.logredact import (
    REDACTED,
    RedactCurationIdentifiers,
    add_values,
    begin_curation,
    current_id,
    end_curation,
    install,
    unfiltered_handlers,
)


@pytest.fixture
def cap():
    """A 'heartvar' logger with the filter installed, and a buffer to read."""
    lg = logging.getLogger("heartvar")
    saved_handlers, saved_filters = lg.handlers[:], lg.filters[:]
    saved_level = lg.level
    lg.handlers, lg.filters = [], []
    lg.setLevel(logging.DEBUG)
    buf = io.StringIO()
    lg.addHandler(logging.StreamHandler(buf))
    install("heartvar")
    try:
        yield buf
    finally:
        end_curation()
        lg.handlers, lg.filters = saved_handlers, saved_filters
        lg.setLevel(saved_level)


def test_gene_and_hgvs_are_redacted_from_a_child_logger(cap):
    """Every client logs through a child logger ("heartvar.clinvar"), which is
    the case the first implementation got wrong."""
    begin_curation("MYH7", "NM_000257.4:c.1208G>A")
    logging.getLogger("heartvar.clinvar").exception(
        "ClinVar local DB query failed for %s %s", "MYH7", "NM_000257.4:c.1208G>A")
    out = cap.getvalue()
    assert "MYH7" not in out
    assert "c.1208G>A" not in out
    assert REDACTED in out


def test_the_filter_must_be_on_handlers_not_the_logger():
    """A filter on a LOGGER never sees records propagated from a child logger.
    Pinned because it is invisible: the code looks right and leaks everything.
    """
    lg = logging.getLogger("hv_placement_probe")
    lg.setLevel(logging.DEBUG)
    lg.handlers, lg.filters = [], []
    buf = io.StringIO()
    lg.addHandler(logging.StreamHandler(buf))
    lg.addFilter(RedactCurationIdentifiers())
    begin_curation("MYH7")
    try:
        logging.getLogger("hv_placement_probe.child").warning("gene %s", "MYH7")
        assert "MYH7" in buf.getvalue(), (
            "logger-level filters now catch propagated records; if this ever "
            "starts passing, the handler-level install is no longer required")
    finally:
        end_curation()


def test_derived_identifiers_can_be_added_later(cap):
    """Coordinates, the gnomAD id and the MANE HGVS are produced after the
    request starts and are just as identifying."""
    begin_curation("MYH7", "NM_000257.4:c.1208G>A")
    add_values(47348490, "11-47348490-C-T", "11-47348490", "11:47348490")
    logging.getLogger("heartvar.gnomad").warning(
        "gnomAD freq payload for %s is not valid JSON", "11-47348490-C-T")
    out = cap.getvalue()
    assert "47348490" not in out


def test_timings_and_counts_survive(cap):
    """The redaction must not corrupt the operational content of a log line.
    A chromosome is "11", and a naive word-boundary rule replaces the "11"
    inside a timing like "0.11"."""
    begin_curation("MYH7")
    add_values(47348490, "11-47348490")
    log = logging.getLogger("heartvar.evidence")
    log.info("[timing] clinvar total=%.3fs waited=%.3fs", 0.117, 0.011)
    log.info("gathered %d of %d sources in %.2fs", 11, 25, 3.11)
    out = cap.getvalue()
    assert "0.117" in out and "0.011" in out
    assert "11 of 25" in out


def test_short_numeric_values_are_never_registered(cap):
    """A bare chromosome identifies nothing and collides with ordinary numbers,
    so it is deliberately not redacted."""
    begin_curation("11", "7", "X")
    logging.getLogger("heartvar.x").info("processed %d records in %.2fs", 11, 7.0)
    assert "11 records" in cap.getvalue()


def test_nothing_is_redacted_outside_a_curation(cap):
    """Startup and admin logging must be unaffected."""
    end_curation()
    logging.getLogger("heartvar.app").info("Deploy log: %s", "2026-09-08T00-00-00Z.log")
    assert "2026-09-08T00-00-00Z.log" in cap.getvalue()


def test_scope_closes(cap):
    begin_curation("MYH7")
    end_curation()
    logging.getLogger("heartvar.x").info("gene %s", "MYH7")
    assert "MYH7" in cap.getvalue()


def test_correlation_id_is_random_and_not_derived_from_content():
    """A hash of the gene would be brute-forced in milliseconds — there are
    about 20,000 symbols. The id must carry no information about the variant."""
    a = begin_curation("MYH7", "NM_000257.4:c.1208G>A")
    b = begin_curation("MYH7", "NM_000257.4:c.1208G>A")
    end_curation()
    assert a != b, "same input gave the same id — the id is content-derived"
    assert len(a) == 8 and all(c in "0123456789abcdef" for c in a)


def test_records_carry_the_correlation_id(cap):
    seen = {}

    class Grab(logging.Handler):
        def emit(self, record):
            seen["curation"] = getattr(record, "curation", None)

    lg = logging.getLogger("heartvar")
    h = Grab()
    h.addFilter(RedactCurationIdentifiers())
    lg.addHandler(h)
    try:
        cid = begin_curation("MYH7")
        logging.getLogger("heartvar.x").info("something happened")
        assert seen["curation"] == cid == current_id()
    finally:
        lg.removeHandler(h)
        end_curation()


def test_a_broken_format_string_does_not_stop_logging():
    """A record whose own format string is wrong must still pass the filter.
    Exercised against the filter directly: routing it through a logger makes
    pytest's logging plugin raise on the malformed record before the filter is
    reached, which tests pytest rather than this code."""
    filt = RedactCurationIdentifiers()
    record = logging.LogRecord(
        name="heartvar.x", level=logging.INFO, pathname=__file__, lineno=1,
        msg="two placeholders %s %s", args=("only-one",), exc_info=None,
    )
    begin_curation("MYH7")
    try:
        assert filt.filter(record) is True
        assert record.msg == "two placeholders %s %s"
        assert record.args == ("only-one",)
    finally:
        end_curation()


def test_install_is_idempotent():
    lg = logging.getLogger("hv_idem")
    lg.handlers = []
    lg.addHandler(logging.StreamHandler(io.StringIO()))
    first = install("hv_idem")
    second = install("hv_idem")
    assert first == 1 and second == 0
    assert unfiltered_handlers("hv_idem") == []


def test_the_app_installs_it_on_its_own_handlers():
    """Importing the app must leave no unfiltered handler on 'heartvar' or root,
    or the deployment log leaks."""
    import backend.app  # noqa: F401
    assert unfiltered_handlers("heartvar", "") == []
