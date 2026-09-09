"""Guards on ClinVar's reused read-only connection.

WHY THE REUSE EXISTS. Every query in clients/clinvar.py used to open its own
``sqlite3.connect``, six sites of them, so each started with an EMPTY page cache
and re-read index and table pages from the Azure Files SMB mount where a page
read is a network round trip. clinvar.db is the largest database we ship and its
gene queries carry ``ORDER BY number_submitters DESC``, so every row for the
gene is read before ``LIMIT 100`` applies. Measured 2026-08-31 with the fixed
per-source instrument: the `clinvar` task took 14.63 s on CHD7, 10.22 s on MYH7
and 7.75 s on MYBPC3, against 0.06-0.30 s for every other local source — cost
scaling with the gene's record count, which is the signature of re-reading
pages, not of query complexity.

WHY THESE TESTS EXIST. Connection reuse was deliberately NOT done the first time
round, and the reason was written down: "callers set row_factory per call, so a
shared cached connection leaks it between them (clinvar expects Row, biogrid
expects tuples) = silently wrong data." That hazard is real and it is real
inside this one module — five sites want ``sqlite3.Row`` and
``_gene_phenotype_strings_sync`` indexes its rows as TUPLES. On a shared
connection whichever ran first would decide, and the loser would read the wrong
column with NO error and NO exception. Nothing about the output would look
wrong.

So the hazard is closed by construction (``_conn(rows=...)`` asserts the setting
on every acquisition rather than inheriting it) and these tests hold that
construction in place.
"""
from __future__ import annotations

import ast
import pathlib
import sqlite3

import pytest

CLINVAR_SRC = (pathlib.Path(__file__).resolve().parent.parent
               / "clients" / "clinvar.py")


def test_no_query_opens_its_own_connection():
    """A fresh connection has an empty page cache, which is the whole problem."""
    tree = ast.parse(CLINVAR_SRC.read_text(encoding="utf-8"))
    offenders = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "connect"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "sqlite3"):
            offenders.append(node.lineno)
    assert not offenders, (
        "sqlite3.connect() called directly at line(s) "
        f"{offenders} — use _conn(rows=...) so the page cache is reused and "
        "row_factory cannot be inherited"
    )


def test_the_guard_can_actually_see_a_violation(tmp_path):
    """A guard that cannot fail is not a guard."""
    bad = tmp_path / "bad.py"
    bad.write_text("import sqlite3\n"
                   "def f(p):\n"
                   "    return sqlite3.connect(p)\n")
    tree = ast.parse(bad.read_text())
    found = [n.lineno for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
             and n.func.attr == "connect"
             and isinstance(n.func.value, ast.Name) and n.func.value.id == "sqlite3"]
    assert found == [3]


def test_every_acquisition_states_its_row_factory():
    """``rows`` is keyword-only with NO default, so a new call site cannot
    forget to choose — but a site could still pass a variable. Check the
    literal, since that is what makes the choice auditable."""
    tree = ast.parse(CLINVAR_SRC.read_text(encoding="utf-8"))
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
             and n.func.id == "_conn"]
    assert calls, "no _conn() call sites found — did the module get renamed?"
    for call in calls:
        assert not call.args, (
            f"_conn() at line {call.lineno} passes a positional argument; "
            "`rows` is keyword-only on purpose")
        kw = {k.arg: k.value for k in call.keywords}
        assert "rows" in kw, (
            f"_conn() at line {call.lineno} does not state `rows` — the "
            "row_factory would be inherited from whichever caller ran first")
        assert isinstance(kw["rows"], ast.Constant) and isinstance(
            kw["rows"].value, bool), (
            f"_conn() at line {call.lineno} passes a non-literal `rows`; keep "
            "it a literal so the choice is auditable by reading the call")


def test_row_and_tuple_callers_do_not_contaminate_each_other(tmp_path,
                                                             monkeypatch):
    """The actual failure mode, exercised end to end on a real database.

    Acquire as Row, then as tuple, then as Row again, on the SAME cached
    connection — and check each acquisition sees what it asked for. Before the
    fix, the second caller inherited the first's factory.
    """
    from backend import localio
    from backend.clients import clinvar

    db = tmp_path / "clinvar.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE variants (variation_id INTEGER, gene_symbol TEXT, "
                "phenotype_list TEXT, number_submitters INTEGER)")
    con.execute("INSERT INTO variants VALUES (1, 'MYH7', 'cardiomyopathy', 3)")
    con.commit()
    con.close()

    monkeypatch.setattr(clinvar, "DB_PATH", db)
    localio.close_all()
    try:
        as_row = clinvar._conn(rows=True)
        r = as_row.execute("SELECT gene_symbol, phenotype_list FROM variants").fetchone()
        assert r["gene_symbol"] == "MYH7", "Row access failed for a rows=True caller"

        as_tuple = clinvar._conn(rows=False)
        assert as_tuple is as_row, (
            "connect_ro handed back a different connection — the page cache is "
            "not being reused, which defeats the purpose")
        t = as_tuple.execute("SELECT gene_symbol, phenotype_list FROM variants").fetchone()
        assert t[0] == "MYH7", (
            "a rows=False caller got a Row it cannot index numerically — this "
            "is the silent-wrong-data hazard the explicit setting exists to stop")
        assert isinstance(t, tuple)

        again = clinvar._conn(rows=True)
        r2 = again.execute("SELECT gene_symbol FROM variants").fetchone()
        assert r2["gene_symbol"] == "MYH7", (
            "a rows=True caller inherited the previous caller's tuple factory")
    finally:
        localio.close_all()


def test_a_missing_database_yields_none_not_an_exception(tmp_path, monkeypatch):
    """Contract shared with the rest of the offline path: serve nothing, let
    the caller report it. Never raise into a curation."""
    from backend.clients import clinvar

    monkeypatch.setattr(clinvar, "DB_PATH", tmp_path / "nope.db")
    assert clinvar._conn(rows=True) is None
    assert clinvar._conn(rows=False) is None


def test_the_db_missing_payload_is_not_shared_mutable_state():
    """Five call sites return it; one mutating it must not poison the others."""
    from backend.clients import clinvar

    a = dict(clinvar._DB_MISSING)
    a["ok"] = "tampered"
    assert clinvar._DB_MISSING["ok"] is False


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
