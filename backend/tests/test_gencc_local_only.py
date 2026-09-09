"""GenCC must never fetch at curation time.

WHY THIS EXISTS. gencc.py served the local snapshot only while it was under
30 days old (``CACHE_TTL_SECONDS``) and otherwise fell through to a live
thegencc.org TSV download — inside a curation. GenCC is a MONTHLY-cadence
source (build_all.sh:488, copied onto the mount at :492), so the snapshot ages
past that boundary every single month and the next curation pays for the
refresh. The switch that would have prevented it, HEARTVAR_OFFLINE_STRICT, is
not set in the deploy.

A snapshot going stale is a BUILD problem. Fixing it inside a user's request is
the wrong place: it makes one unlucky curation wear a multi-second download, and
it fails in the one situation where the network is the thing that is broken.
So the read path is local-only and staleness is REPORTED instead.

``_refresh_cache`` still exists and is still the one true parser — scripts/
build_gencc_db.py calls it directly, which is what keeps the snapshot current.
"""
from __future__ import annotations

import asyncio
import json
import time

import pytest

from backend.clients import gencc


def _snapshot(tmp_path, *, age_days: float = 0.0):
    p = tmp_path / "gencc_submissions.json"
    p.write_text(json.dumps({
        "genes": {"MYH7": [{"disease": "HCM", "classification": "Definitive",
                            "submitter": "ClinGen", "moi": "AD",
                            "gene_curie": "HGNC:7577"}]},
        "meta": {"row_count": 1},
    }))
    if age_days:
        old = time.time() - age_days * 24 * 3600
        import os
        os.utime(p, (old, old))
    return p


def _reset(monkeypatch, path):
    monkeypatch.setattr(gencc, "DATA_FILE", path, raising=False)
    monkeypatch.setattr(gencc, "_INDEX", None, raising=False)
    monkeypatch.setattr(gencc, "_META", None, raising=False)

    async def _boom():
        raise AssertionError(
            "gencc reached the live thegencc.org refresh during a curation"
        )
    monkeypatch.setattr(gencc, "_refresh_cache", _boom)


def test_a_fresh_snapshot_is_served_locally(monkeypatch, tmp_path):
    _reset(monkeypatch, _snapshot(tmp_path))
    index, meta = asyncio.run(gencc._ensure_loaded())
    assert index and "MYH7" in index
    assert not meta.get("stale")


def test_a_snapshot_PAST_THE_TTL_is_still_served_and_never_refetched(
        monkeypatch, tmp_path):
    """⚠ THE BUG. 400 days old — the old code would have gone to the network
    here, mid-curation. GenCC is monthly-cadence, so this state is reached
    every month by construction."""
    _reset(monkeypatch, _snapshot(tmp_path, age_days=400))
    index, meta = asyncio.run(gencc._ensure_loaded())
    assert index and "MYH7" in index, "a stale snapshot must still be served"
    assert meta.get("stale") is True, "staleness must be reported, not hidden"
    assert meta.get("age_days", 0) >= 399


def test_it_is_local_only_even_with_offline_strict_OFF(monkeypatch, tmp_path):
    """The guarantee must not depend on an operator flag. HEARTVAR_OFFLINE_STRICT
    is not set in the deploy, and gencc must be local-only regardless."""
    monkeypatch.delenv("HEARTVAR_OFFLINE_STRICT", raising=False)
    _reset(monkeypatch, _snapshot(tmp_path, age_days=400))
    index, _ = asyncio.run(gencc._ensure_loaded())
    assert index and "MYH7" in index


def test_a_missing_snapshot_reports_rather_than_downloading(monkeypatch, tmp_path):
    """No snapshot is a BUILD failure. Downloading one inside a curation hides
    that, and does it while wearing the user's latency."""
    _reset(monkeypatch, tmp_path / "absent.json")
    index, meta = asyncio.run(gencc._ensure_loaded())
    assert index is None
    assert "error" in meta and "snapshot" in meta["error"].lower()


def test_the_builder_still_owns_the_refresh():
    """_refresh_cache must remain the single parser — scripts/build_gencc_db.py
    calls it directly, and a second copy of the TSV schema would drift."""
    import inspect
    from pathlib import Path
    assert callable(gencc._refresh_cache)
    builder = (Path(__file__).resolve().parents[2]
               / "scripts" / "build_gencc_db.py").read_text()
    assert "_refresh_cache()" in builder
    import ast, textwrap
    tree = ast.parse(textwrap.dedent(inspect.getsource(gencc._ensure_loaded)))
    called = {
        n.func.id for n in ast.walk(tree)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
    } | {
        n.func.attr for n in ast.walk(tree)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
    }
    assert "_refresh_cache" not in called, (
        "the curation read path must never trigger a live GenCC download"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
