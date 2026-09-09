"""The ClinVar build must not publish a database missing its query index.

WHY THIS EXISTS, and it is a process lesson rather than a code one. On
2026-08-31 `idx_variants_gene_submitters` was added to
scripts/build_clinvar_db.py, a targeted ClinVar rebuild was run, it reported
"indexes built in 38.0 s" and succeeded — and production stayed exactly as slow.
There was no way to tell from the log whether the published database actually
had the index or whether the container app had pulled a stale
`heartvar-builder:latest`. That ambiguity cost a full measure-deploy-measure
cycle.

The index is not cosmetic. The gene queries in clients/clinvar.py are
``WHERE gene_symbol = ? ORDER BY number_submitters DESC, variation_id ASC
LIMIT 100``; without it SQLite adds ``USE TEMP B-TREE FOR ORDER BY`` and reads
every row for the gene before the LIMIT applies. Measured on the real
4.46M-row database: MYH7 172.2 ms -> 4.1 ms, and over the Azure Files SMB mount
that difference presents as ~14 s per curation. So a build that omits it ships a
14-second regression with no signal.

These tests pin the DDL the build issues and the gate that refuses to publish
without it, so the next silent omission fails the build instead of the site.
"""
from __future__ import annotations

import pathlib
import re
import sqlite3

import pytest

BUILD_SRC = (pathlib.Path(__file__).resolve().parent.parent.parent
             / "scripts" / "build_clinvar_db.py")

REQUIRED_INDEX = "idx_variants_gene_covering"


def _source() -> str:
    return BUILD_SRC.read_text(encoding="utf-8")


def test_the_build_creates_the_gene_query_index():
    src = _source()
    assert REQUIRED_INDEX in src, (
        f"{REQUIRED_INDEX} is not created by the build — every gene query "
        "falls back to reading the whole gene and sorting it")


def test_the_index_leads_with_gene_symbol_then_the_sort_key():
    """Column ORDER is the whole point. (number_submitters, gene_symbol) would
    not serve `WHERE gene_symbol = ?`, and omitting DESC would not serve the
    ORDER BY without a reverse scan."""
    src = _source()
    flat = re.sub(r'["\']\s*\n\s*["\']', "", src)
    m = re.search(r"CREATE INDEX\s+" + REQUIRED_INDEX
                  + r"\s+ON variants\s*\(([^)]*)\)", flat)
    assert m, f"could not find the {REQUIRED_INDEX} DDL in the build script"
    cols = " ".join(m.group(1).replace('"', " ").split())
    assert cols.startswith("gene_symbol"), (
        f"index must lead with gene_symbol to serve the WHERE clause: {cols}")
    assert "number_submitters DESC" in cols, (
        f"index must carry the sort key, descending, to avoid the temp "
        f"B-tree: {cols}")
    for col in ("name", "start", "stop", "obj_type"):
        assert col in cols, (
            f"{col} missing from the index: the WHERE clause would be "
            f"evaluated against the table, one row fetch per candidate. {cols}")


def test_the_build_refuses_to_publish_without_it():
    """A gate that cannot fail is not a gate."""
    src = _source()
    assert "missing required index" in src, (
        "the build does not verify its own indexes before publishing")
    call = "_dbbuild.publish(tmp_path, db_path)"
    assert call in src, f"could not find the publish call ({call!r})"
    assert src.index("missing required index") < src.index(call), (
        "the index check runs AFTER publishing, so a bad database still ships")


def test_the_build_names_the_indexes_it_created():
    """"indexes built in 38.0 s" is indistinguishable between three indexes and
    four. Naming them is what makes a stale builder image visible in the log."""
    src = _source()
    assert "sqlite_master" in src and "tbl_name='variants'" in src, (
        "the build does not read back and log the indexes it created")


def test_the_index_actually_removes_the_temp_btree():
    """The claim, on a real SQLite rather than by assertion."""
    con = sqlite3.connect(":memory:")
    con.execute("CREATE TABLE variants (variation_id INTEGER, gene_symbol TEXT, "
                "number_submitters INTEGER, name TEXT)")
    con.executemany(
        "INSERT INTO variants VALUES (?,?,?,?)",
        [(i, "MYH7", i % 7, f"NM_x:c.{i}G>A") for i in range(500)])
    q = ("SELECT variation_id FROM variants WHERE gene_symbol = ? "
         "ORDER BY number_submitters DESC, variation_id ASC LIMIT 100")

    con.execute("CREATE INDEX idx_variants_gene ON variants (gene_symbol)")
    before = " ".join(r[3] for r in con.execute("EXPLAIN QUERY PLAN " + q, ("MYH7",)))
    assert "TEMP B-TREE" in before.upper(), (
        f"expected the pre-fix plan to sort; got: {before}")

    con.execute(f"CREATE INDEX {REQUIRED_INDEX} ON variants "
                "(gene_symbol, number_submitters DESC, variation_id ASC)")
    after = " ".join(r[3] for r in con.execute("EXPLAIN QUERY PLAN " + q, ("MYH7",)))
    assert "TEMP B-TREE" not in after.upper(), (
        f"the index did not remove the sort: {after}")
    assert REQUIRED_INDEX in after, after
    con.close()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


GENE_QUERY = """
SELECT variation_id, obj_type, name, gene_symbol,
       clinical_significance, review_status,
       number_submitters, phenotype_list, last_evaluated, start, stop
FROM variants
WHERE gene_symbol = ?
  AND (name LIKE ? OR name LIKE ? OR (start = ? AND stop = ?
       AND obj_type = 'single nucleotide variant'))
ORDER BY number_submitters DESC, variation_id ASC
LIMIT 100
"""


def _filter_reads(con):
    """(cols read from the table, cols read from the index) by the filter."""
    prog = list(con.execute("EXPLAIN " + GENE_QUERY,
                            ("MYH7", "%x%", "%x%", 1, 1)))
    branches = [a for a, op, p1, p2, *_ in prog if op == "Ne"]
    if not branches:
        pytest.skip("this SQLite compiled a shape this test cannot read")
    last = max(branches)
    table = index = 0
    for addr, op, p1, p2, p3, p4, p5, cmt in prog:
        if addr > last:
            break
        if op == "Column":
            if p1 == 0:
                table += 1
            elif p1 == 2:
                index += 1
    return table, index


def _fixture(tmp_path, index_cols):
    tmp_path = pathlib.Path(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "cv.db"
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE variants (variation_id INTEGER, obj_type TEXT, "
                "name TEXT, gene_symbol TEXT, clinical_significance TEXT, "
                "review_status TEXT, number_submitters INTEGER, "
                "phenotype_list TEXT, last_evaluated TEXT, start INTEGER, "
                "stop INTEGER)")
    con.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        [(i, "single nucleotide variant", f"NM_x(MYH7):c.{i}G>A", "MYH7",
          "Pathogenic", "reviewed", i % 5, "cardiomyopathy", "2026", i, i)
         for i in range(400)])
    con.execute("CREATE INDEX idx_variants_gene ON variants (gene_symbol)")
    con.execute(f"CREATE INDEX idx_variants_gene_submitters ON variants ({index_cols})")
    con.commit()
    return con


def test_the_narrow_index_still_reads_the_table(tmp_path):
    """The control. Without this, the test below could pass vacuously."""
    con = _fixture(tmp_path, "gene_symbol, number_submitters DESC, variation_id ASC")
    table, _index = _filter_reads(con)
    con.close()
    assert table > 0, (
        "expected the narrow index to force table reads for the filter; if "
        "this SQLite does not, the test below proves nothing")


def test_the_shipped_index_makes_the_filter_index_only(tmp_path):
    """The property the 13.73 s depended on."""
    con = _fixture(
        tmp_path,
        "gene_symbol, number_submitters DESC, variation_id ASC, "
        "name, start, stop, obj_type")
    table, index = _filter_reads(con)
    con.close()
    assert index > 0, "the filter reads nothing from the index at all"
    assert table == 0, (
        f"the filter still reads {table} column(s) from the TABLE, so every "
        "candidate row is fetched to evaluate it — that is one SMB round trip "
        "per row of the gene, measured at 13.73 s to return two rows on MYH7")


def test_the_wider_index_returns_the_same_rows(tmp_path):
    """A faster index that changes the answer is not a fix."""
    narrow = _fixture(tmp_path / "a", "gene_symbol, number_submitters DESC, variation_id ASC")
    wide = _fixture(tmp_path / "b",
                    "gene_symbol, number_submitters DESC, variation_id ASC, "
                    "name, start, stop, obj_type")
    args = ("MYH7", "%c.42G>A%", "%c.42G>A%", 99, 99)
    a = narrow.execute(GENE_QUERY, args).fetchall()
    b = wide.execute(GENE_QUERY, args).fetchall()
    narrow.close(); wide.close()
    assert a == b and a, f"results diverged: {a} vs {b}"


def test_every_index_is_created_before_the_gate_reads_them_back():
    src = _source()
    read_back = src.index("SELECT name FROM sqlite_master")
    creates = [m.start() for m in re.finditer(r"CREATE INDEX (\w+)", src)]
    names = re.findall(r"CREATE INDEX (\w+)", src)
    assert creates, "no CREATE INDEX statements found"
    late = [n for pos, n in zip(creates, names) if pos > read_back]
    assert not late, (
        f"index/indexes {late} are created AFTER the sqlite_master read-back, "
        "so the gate checks a stale snapshot and refuses to publish a database "
        "that is actually correct"
    )


def test_the_gate_checks_exactly_what_the_build_creates():
    """A required set that drifts from the DDL fails builds for no reason (or,
    worse, passes one that is missing an index)."""
    src = "\n".join(l for l in _source().splitlines()
                    if not l.lstrip().startswith("#"))
    created = set(re.findall(r"CREATE INDEX (\w+)", src))
    m = re.search(r"required = \{(.*?)\}", src, re.S)
    assert m, "could not find the required-index set"
    required = set(re.findall(r'"(\w+)"', m.group(1)))
    assert required == created, (
        f"the gate's required set and the CREATE INDEX statements disagree.\n"
        f"  created but not required: {sorted(created - required)}\n"
        f"  required but not created: {sorted(required - created)}"
    )


GENE_QUERIES = {
    "variant lookup": (
        "SELECT variation_id, obj_type, name, gene_symbol, clinical_significance,"
        " review_status, number_submitters, phenotype_list, last_evaluated,"
        " start, stop FROM variants WHERE gene_symbol = ?"
        " AND (name LIKE ? OR name LIKE ? OR (start = ? AND stop = ?"
        " AND obj_type = 'single nucleotide variant'))"
        " ORDER BY number_submitters DESC, variation_id ASC LIMIT 100",
        ("MYH7", "%x%", "%x%", 1, 1)),
    "gene landscape": (
        "SELECT variation_id, name, clinical_significance, review_status,"
        " number_submitters, phenotype_list, start FROM variants"
        " WHERE gene_symbol = ? AND (clinical_significance LIKE '%pathogenic%')",
        ("MYH7",)),
    "pm5 same-residue": (
        "SELECT variation_id, name, clinical_significance, review_status,"
        " number_submitters, phenotype_list, start FROM variants"
        " WHERE gene_symbol = ? AND name LIKE '%p.%'"
        " ORDER BY number_submitters DESC, variation_id ASC",
        ("MYH7",)),
    "phenotype vocabulary": (
        "SELECT DISTINCT phenotype_list FROM variants WHERE gene_symbol = ?"
        " AND phenotype_list IS NOT NULL AND phenotype_list != ''",
        ("MYH7",)),
}

COVERING_COLS = ("gene_symbol, number_submitters DESC, variation_id ASC, "
                 "name, clinical_significance, review_status, phenotype_list, "
                 "start, stop, obj_type, last_evaluated")


def _fixture(tmp_path, index_cols=None):
    """A cdot-shaped variants table. index_cols=None means only the plain
    gene_symbol index, i.e. the pre-fix state."""
    tmp_path = pathlib.Path(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(tmp_path / "cv.db")
    con.execute("CREATE TABLE variants (variation_id INTEGER, obj_type TEXT, "
                "name TEXT, gene_symbol TEXT, clinical_significance TEXT, "
                "review_status TEXT, number_submitters INTEGER, "
                "phenotype_list TEXT, last_evaluated TEXT, start INTEGER, "
                "stop INTEGER)")
    con.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        [(i, "single nucleotide variant", f"NM_x(MYH7):c.{i}G>A p.Arg{i}Cys",
          "MYH7", "Pathogenic", "reviewed", i % 5, "cardiomyopathy", "2026", i, i)
         for i in range(400)])
    con.execute("CREATE INDEX idx_variants_gene ON variants (gene_symbol)")
    if index_cols:
        con.execute(f"CREATE INDEX {REQUIRED_INDEX} ON variants ({index_cols})")
    con.commit()
    return con


def _table_reads(con, sql, args):
    """How many columns the query reads from the TABLE cursor (p1=0) rather
    than an index cursor. Zero means index-only."""
    prog = list(con.execute("EXPLAIN " + sql, args))
    return sum(1 for _a, op, p1, *_ in prog if op == "Column" and p1 == 0)


@pytest.mark.parametrize("label", sorted(GENE_QUERIES))
def test_the_pre_fix_state_does_read_the_table(tmp_path, label):
    """The CONTROL. Without it the test below could pass vacuously."""
    sql, args = GENE_QUERIES[label]
    con = _fixture(tmp_path / label.replace(" ", "_"))
    reads = _table_reads(con, sql, args)
    con.close()
    assert reads > 0, (
        f"{label}: expected the pre-fix state to fetch table rows; if this "
        "SQLite does not, the covering assertion proves nothing")


@pytest.mark.parametrize("label", sorted(GENE_QUERIES))
def test_every_gene_query_is_index_only_with_the_covering_index(tmp_path, label):
    """The property the whole fix rests on: no query may fetch a table row to
    evaluate its filter, because each such fetch is an SMB round trip."""
    sql, args = GENE_QUERIES[label]
    con = _fixture(tmp_path / label.replace(" ", "_"), COVERING_COLS)
    reads = _table_reads(con, sql, args)
    plan = " ".join(r[3] for r in con.execute("EXPLAIN QUERY PLAN " + sql, args))
    con.close()
    assert reads == 0, (
        f"{label}: still reads {reads} column(s) from the TABLE, so every "
        f"candidate row is a page fetch. Plan: {plan}")
    assert "COVERING INDEX" in plan.upper(), f"{label}: {plan}"


def test_the_covering_index_returns_the_same_rows(tmp_path):
    """A faster index that changes the answer is not a fix."""
    sql, args = GENE_QUERIES["pm5 same-residue"]
    before = _fixture(tmp_path / "a")
    after = _fixture(tmp_path / "b", COVERING_COLS)
    a, b = before.execute(sql, args).fetchall(), after.execute(sql, args).fetchall()
    before.close(); after.close()
    assert a == b and a, f"results diverged: {len(a)} vs {len(b)}"


def test_the_index_leads_with_gene_symbol_then_the_sort_key():
    """Column ORDER is load-bearing: gene_symbol first or it cannot serve the
    WHERE; number_submitters DESC next or the ORDER BY needs a temp B-tree."""
    src = _source()
    flat = re.sub(r'["\']\s*\n\s*["\']', "", src)
    m = re.search(r"CREATE INDEX\s+" + REQUIRED_INDEX
                  + r"\s+ON variants\s*\(([^)]*)\)", flat)
    assert m, f"could not find the {REQUIRED_INDEX} DDL"
    cols = " ".join(m.group(1).replace('"', " ").split())
    assert cols.startswith("gene_symbol"), cols
    assert "number_submitters DESC" in cols, cols
    for col in ("name", "clinical_significance", "review_status",
                "phenotype_list", "start", "stop", "obj_type"):
        assert col in cols, (
            f"{col} missing — its query goes back to fetching table rows. {cols}")


def test_the_table_is_not_rewritten():
    """Rewriting the table in gene order achieves the same thing and failed
    twice in the build container (OOM-killed at exit 137; then still running
    after 35 minutes). CREATE INDEX has completed on every build."""
    code = "\n".join(l for l in _source().splitlines()
                     if not l.lstrip().startswith("#"))
    assert "ORDER BY gene_symbol" not in code, (
        "the build rewrites the table in gene order again — that path OOM-killed "
        "and then hung in the build container")


def test_the_index_sorter_is_pinned_to_disk():
    """CREATE INDEX sorts every row by the index key, and the covering index
    carries name and phenotype_list — ~900 MB of sort data over 4.46 M rows.
    With temp_store at its default that sort lives in memory, and this
    container has a ceiling: the same default OOM-killed the table-rewrite
    attempt at exit 137."""
    code = "\n".join(l for l in _source().splitlines()
                     if not l.lstrip().startswith("#"))
    assert "temp_store=FILE" in code, (
        "the index sorter is not pinned to a file; a ~900 MB in-memory sort "
        "will OOM the build container")
    assert code.index("temp_store=FILE") < code.index("CREATE INDEX"), (
        "temp_store is set after the indexes are already being built")
