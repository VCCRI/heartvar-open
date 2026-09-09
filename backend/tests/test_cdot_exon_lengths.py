"""Exon lengths served from cdot, so offline VEP can apply the FULL PVS1 NMD rule.

WHAT THIS CLOSES. clients/vep_offline set ``exon_lengths: None`` and computed
``nmd_escape`` from the exon "N/total" string alone, so the third ClinGen PVS1
rule — a PTC in the LAST 50 nt of the PENULTIMATE exon escapes NMD — could never
fire offline: offline said False where REST said True. The direction matters —
nmd_escape=False keeps PVS1 at full strength, so offline OVER-CALLED PVS1 for
exactly that class.

⚠ THE OFF-BY-ONE THAT WOULD HAVE MADE THIS SILENTLY WRONG. cdot stores exons as
    [genomic_start, genomic_end, transcript_index, cdna_start, cdna_end, gap]
with genomic starts 0-BASED HALF-OPEN and cDNA 1-based inclusive. Measured on
NM_000257.4: every exon's ``end - start + 1`` is exactly 1 MORE than its cDNA
span, and the sums differ by exactly the exon count (6067 vs 6027 over 40
exons). Using the genomic form would have shifted the cumulative cDNA offset by
40 bases on MYH7 and mis-called the 50 nt boundary.

The cDNA spans reproduce Ensembl's own numbers exactly. REST for
ENST00000355349: [41, 56, 209, 144, 157, ...] summing to 6027 — and cdot gives
41, 56, ..., 132 for the same transcript.
"""
from __future__ import annotations

import asyncio
import gzip
import json
import sqlite3

import pytest

from backend.clients import cdot_sqlite


def _db(tmp_path, records: dict[str, dict]):
    path = tmp_path / "cdot.db"
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE transcripts (accession TEXT PRIMARY KEY, "
                 "gene TEXT, data BLOB)")
    for acc, rec in records.items():
        conn.execute("INSERT INTO transcripts VALUES (?,?,?)",
                     (acc, rec.get("gene_name"),
                      gzip.compress(json.dumps(rec).encode())))
    conn.commit()
    conn.close()
    return str(path)


_MYH7_TAIL = {
    "gene_name": "MYH7",
    "genome_builds": {"GRCh38": {"strand": "-", "contig": "NC_000014.9", "exons": [
        [23412739, 23412871, 3, 5896, 6027, None],
        [23413758, 23413893, 2, 5761, 5895, None],
        [23414006, 23414102, 1, 5665, 5760, None],
        [23435619, 23435660, 0, 1, 41, None],
    ]}},
}


def test_lengths_are_cdna_spans_in_transcript_order(tmp_path):
    """Transcript order, from the exon INDEX — not the stored list order, which
    is ascending genomic and therefore reversed on the minus strand."""
    db = _db(tmp_path, {"NM_000257.4": _MYH7_TAIL})
    assert cdot_sqlite.exon_cdna_lengths(db, "NM_000257.4") == [41, 96, 135, 132]


def test_the_genomic_span_form_is_not_used(tmp_path):
    """⚠ REGRESSION GUARD ON THE OFF-BY-ONE. `end - start + 1` would give
    [42, 97, 136, 133] — one too many per exon, 4 too many in the cumulative
    sum here and 40 too many on the real 40-exon transcript."""
    db = _db(tmp_path, {"NM_000257.4": _MYH7_TAIL})
    got = cdot_sqlite.exon_cdna_lengths(db, "NM_000257.4")
    assert got != [42, 97, 136, 133], (
        "exon lengths were computed from the 0-based half-open GENOMIC span; "
        "cdot's cDNA columns are the authoritative 1-based inclusive lengths"
    )
    assert got == [41, 96, 135, 132]
    assert sum(got) == 404


def test_a_versionless_accession_resolves_to_the_newest_version(tmp_path):
    """The VEP CLI reports `ENST00000355349` without a version, while cdot keys
    on `ENST00000355349.4`. A miss here would silently reinstate the old
    exon-number-only behaviour."""
    db = _db(tmp_path, {"ENST00000355349.4": _MYH7_TAIL,
                        "ENST00000355349.3": {"genome_builds": {"GRCh38": {
                            "strand": "-", "exons": [[1, 11, 0, 1, 10, None]]}}}})
    assert cdot_sqlite.exon_cdna_lengths(db, "ENST00000355349") == [41, 96, 135, 132]
    assert cdot_sqlite.exon_cdna_lengths(db, "ENST00000355349.3") == [10]


def test_missing_or_unusable_returns_none(tmp_path):
    db = _db(tmp_path, {"NM_1.1": {"genome_builds": {}}})
    assert cdot_sqlite.exon_cdna_lengths(db, "NM_1.1") is None
    assert cdot_sqlite.exon_cdna_lengths(db, "NM_NOPE.9") is None
    assert cdot_sqlite.exon_cdna_lengths(db, None) is None
    assert cdot_sqlite.exon_cdna_lengths("/no/such.db", "NM_1.1") is None


def test_transcript_exons_are_one_based_like_ensembl(tmp_path):
    """⚠ THE SECOND OFF-BY-ONE IN THE SAME RECORD. cdot's genomic starts are
    0-BASED HALF-OPEN; Ensembl's /lookup/id Exon array is 1-based inclusive,
    and codon_genomic_positions compares genomic positions against these spans
    (``host[0] <= pos <= host[1]``). Verified against live REST for
    ENST00000355349: ``start + 1`` reproduces Ensembl's 40 exons EXACTLY
    (first 23435620-23435660, last 23412740-23412871), and without the shift
    they differ."""
    db = _db(tmp_path, {"NM_000257.4": _MYH7_TAIL})
    got = cdot_sqlite.transcript_exons(db, "NM_000257.4")
    assert got["ok"] is True
    assert got["strand"] == -1, "cdot stores '-'/'+'; callers expect -1/+1"
    assert got["transcript_id"] == "NM_000257.4"
    assert [(e["rank"], e["start"], e["end"]) for e in got["exons"]] == [
        (1, 23435620, 23435660),
        (2, 23414007, 23414102),
        (3, 23413759, 23413893),
        (4, 23412740, 23412871),
    ]
    lengths = cdot_sqlite.exon_cdna_lengths(db, "NM_000257.4")
    assert [e["end"] - e["start"] + 1 for e in got["exons"]] == lengths


def test_transcript_exons_missing_returns_none(tmp_path):
    db = _db(tmp_path, {"NM_1.1": {"genome_builds": {}}})
    assert cdot_sqlite.transcript_exons(db, "NM_1.1") is None
    assert cdot_sqlite.transcript_exons(db, "NOPE.1") is None
    assert cdot_sqlite.transcript_exons(None, "NM_1.1") is None


def test_no_cdot_import_is_needed(tmp_path):
    """Deliberately stdlib-only. make_provider drags in cdot -> hgvs ->
    psycopg2 -> ipython; exon lengths are a blob read and must not pay for it,
    and must still work if cdot is not importable at all."""
    import inspect
    for fn in (cdot_sqlite.exon_cdna_lengths, cdot_sqlite.transcript_exons):
        src = inspect.getsource(fn)
        assert "import cdot" not in src and "LocalDataProvider" not in src


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_fetch_transcript_exons_prefers_cdot_and_makes_no_network_call(
        tmp_path, monkeypatch):
    """⚠ THE ~46% OF A WARM CURATION THAT WAS A CALL THAT COULD NOT SUCCEED.

    /lookup/id resolves ENSEMBL ids only, so for the RefSeq transcript the HGVS
    resolver now picks, the round trip is guaranteed to 400 — measured 1.3-1.8 s
    — and evidence.py then spends a SECOND one on the MANE fallback. Warm-cache
    curations measured 2026-08-29: MYH7 3.25 s total with transcript_exons=1.50 s,
    MYBPC3 2.84 s with 1.30 s.

    The tripwire is the point: if this ever reaches httpx again for a transcript
    cdot holds, the test fails rather than the latency quietly returning."""
    import backend.clients.ensembl_vep as ev

    db = _db(tmp_path, {"NM_000257.4": _MYH7_TAIL})
    monkeypatch.setenv("HEARTVAR_CDOT_DB", db)
    monkeypatch.setenv("HEARTVAR_VEP_OFFLINE", "1")
    monkeypatch.setattr(ev, "_TRANSCRIPT_EXONS_CACHE", {})

    class _Boom:
        def __init__(self, *a, **k): raise AssertionError(
            "fetch_transcript_exons went to the network for a transcript cdot holds")
    monkeypatch.setattr(ev.httpx, "AsyncClient", _Boom)

    got = asyncio.run(ev.fetch_transcript_exons("NM_000257.4"))
    assert got is not None and got["ok"] is True
    assert got["transcript_id"] == "NM_000257.4"
    assert got["strand"] == -1
    assert got["exons"][0] == {"start": 23435620, "end": 23435660, "rank": 1}


def test_fetch_transcript_exons_keeps_rest_when_offline_vep_is_off(
        tmp_path, monkeypatch):
    """The ladder must come from the SAME RELEASE as the annotation. cdot is
    matched to the offline cache (113); rest.ensembl.org is on 116. With the
    offline flag off the annotation is REST's, so the ladder must be too —
    mixing them would be the release skew this pairing exists to avoid."""
    import backend.clients.ensembl_vep as ev

    db = _db(tmp_path, {"NM_000257.4": _MYH7_TAIL})
    monkeypatch.setenv("HEARTVAR_CDOT_DB", db)
    monkeypatch.delenv("HEARTVAR_VEP_OFFLINE", raising=False)
    monkeypatch.setattr(ev, "_TRANSCRIPT_EXONS_CACHE", {})

    reached = {"network": False}

    class _Marker:
        def __init__(self, *a, **k): reached["network"] = True
        async def __aenter__(self): raise RuntimeError("stop here")
        async def __aexit__(self, *a): return False
    monkeypatch.setattr(ev.httpx, "AsyncClient", _Marker)

    asyncio.run(ev.fetch_transcript_exons("NM_000257.4"))
    assert reached["network"], (
        "with offline VEP off the ladder must still come from Ensembl, so it "
        "matches the release the annotation came from"
    )


def _versioned_db(tmp_path, accessions):
    import gzip, json, sqlite3
    path = tmp_path / "cdot.db"
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE transcripts (accession TEXT PRIMARY KEY, "
                "gene TEXT, data BLOB NOT NULL)")
    for acc in accessions:
        con.execute("INSERT INTO transcripts VALUES (?,?,?)",
                    (acc, "MYH7", gzip.compress(json.dumps({"id": acc}).encode())))
    con.commit(); con.close()
    return path


def test_the_version_fallback_uses_the_index(tmp_path):
    """A plan containing SCAN means the fallback reverted to LIKE."""
    import sqlite3
    from backend.clients import cdot_sqlite as cs

    path = _versioned_db(tmp_path, ["ENST00000355349.4"])
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    plan = " | ".join(r[3] for r in con.execute(
        "EXPLAIN QUERY PLAN SELECT accession, data FROM transcripts "
        "WHERE accession >= ? AND accession < ?",
        ("ENST00000355349.", "ENST00000355349.￿")))
    con.close()
    assert "SCAN" not in plan.upper(), (
        f"the version fallback is scanning the whole table: {plan}")
    assert "INDEX" in plan.upper(), plan


def test_the_version_fallback_still_picks_the_newest(tmp_path):
    """Behaviour must be identical to the LIKE it replaced."""
    import sqlite3
    from backend.clients import cdot_sqlite as cs

    path = _versioned_db(
        tmp_path, ["ENST00000355349.2", "ENST00000355349.10", "ENST00000355349.4"])
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    blob = cs._newest_versioned(con, "ENST00000355349")
    con.close()
    assert blob is not None
    import gzip, json
    assert json.loads(gzip.decompress(blob))["id"] == "ENST00000355349.10", (
        "picked the wrong version — the range bound must not truncate the "
        "two-digit version")


def test_a_different_transcript_is_not_matched(tmp_path):
    """The upper bound must not spill into the next accession."""
    import sqlite3
    from backend.clients import cdot_sqlite as cs

    path = _versioned_db(tmp_path, ["ENST00000355349X.1", "ENST000003553490.1"])
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    assert cs._newest_versioned(con, "ENST00000355349") is None, (
        "matched an accession that merely starts with the query")
    con.close()
