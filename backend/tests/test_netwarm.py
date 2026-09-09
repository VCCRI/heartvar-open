"""Tests for outbound-connectivity warming — the instrument for the last 18 s.

WHAT IS LEFT, AND WHY IT IS NOT THE MOUNT. Deployment log, 2026-08-30 10:06:

    [warmup] warmed 14/14 local files (concurrency 2) — total file-time 0.49s
    [warmup] resolver + FASTA ready in 0.34s
    [timing] DB gather 23.38s (erepo=18.00s, mgi=18.00s, vep=18.00s,
             gene_literature=18.00s, pubtator3=18.00s, medgen=18.00s,
             clinvar=18.00s, biogrid=18.00s, pubmed=18.00s, fetal_heart=18.00s,
             alphafold=18.00s, gtex=18.00s, uniprot=18.00s, panelapp=18.00s,
             ... gnomad=0.95s, gencc=0.41s, chdgene=0.00s)
    [timing] ... same_site=0.04s      <- was 15.87s before warming

So the mount is fast (0.49 s for everything) and warming worked. But FOURTEEN
sources took *exactly* 18.00 s on the first curation and 0.19 s on the second,
and every one of them makes an OUTBOUND HTTP call — either as its primary source
or as the live fallback that is still enabled (HEARTVAR_OFFLINE_STRICT unset
leaves 12 clients with one). The local-only sources were unaffected.

⚠ THIS MODULE IS DELIBERATELY AN INSTRUMENT FIRST. The identical 18.00 s across
hosts that resolve to DIFFERENT servers does not fit plain per-host DNS latency —
they release together, which points at something shared. Rather than guess again
(two hypotheses have already died here: the sync json.loads, and warming being
too shallow), this measures each layer SEPARATELY and logs it:

    getaddrinfo per host   -> is it DNS?
    TCP connect per host   -> is it reachability / IPv6 fallback?

Whichever column shows the seconds is the answer, and it lands in the deployment
log the admin page serves. Warming the successful lookups is the likely fix and
comes free with the measurement.

Contracts:
  1. Every host is probed; one failure never stops the rest or raises.
  2. It reports per-host DNS and connect timings separately.
  3. It runs off the event loop and can be disabled.
  4. No outbound HTTP REQUEST is made — DNS and TCP only, so no third-party API
     is called and no rate limit is touched.
"""
from __future__ import annotations

import asyncio

import pytest

from backend import netwarm


def test_hosts_are_real_and_deduplicated():
    hosts = netwarm.outbound_hosts()
    assert {"eutils.ncbi.nlm.nih.gov", "rest.ensembl.org"} <= set(hosts)
    assert len(hosts) == len(set(hosts)), "duplicate hosts would double the log"
    assert all("/" not in h and ":" not in h for h in hosts), (
        "hosts must be bare hostnames, not URLs"
    )


def test_reports_dns_and_connect_separately():
    report = asyncio.run(netwarm.warm_outbound(["localhost"]))
    entry = report["hosts"]["localhost"]
    assert "dns_ms" in entry and "connect_ms" in entry, (
        "the two layers must be timed separately — that separation IS the "
        "diagnostic"
    )


def test_a_bogus_host_is_recorded_not_raised():
    report = asyncio.run(netwarm.warm_outbound(
        ["this-host-does-not-exist.invalid"]))
    entry = report["hosts"]["this-host-does-not-exist.invalid"]
    assert entry["error"], "a failed lookup must be reported"
    assert isinstance(entry["dns_ms"], float), "it must still be timed"


def test_one_bad_host_does_not_stop_the_others():
    report = asyncio.run(netwarm.warm_outbound(
        ["this-host-does-not-exist.invalid", "localhost"]))
    assert report["hosts"]["localhost"]["dns_ms"] >= 0.0
    assert len(report["hosts"]) == 2


def test_makes_no_http_request(monkeypatch):
    """DNS and TCP only. A real request would hit third-party rate limits."""
    import backend.netwarm as m
    src = __import__("inspect").getsource(m)
    for forbidden in ("httpx", "urlopen", "requests.get", "AsyncClient"):
        assert forbidden not in src, (
            f"{forbidden} in netwarm would mean a real outbound request"
        )


def test_can_be_disabled(monkeypatch):
    monkeypatch.setenv("HEARTVAR_WARM_OUTBOUND", "0")
    assert netwarm.outbound_warming_enabled() is False


def test_enabled_by_default(monkeypatch):
    monkeypatch.delenv("HEARTVAR_WARM_OUTBOUND", raising=False)
    assert netwarm.outbound_warming_enabled() is True


def test_lifespan_schedules_it():
    import ast
    from pathlib import Path
    tree = ast.parse((Path(__file__).resolve().parents[1] / "app.py")
                     .read_text(encoding="utf-8"))
    src = next(ast.unparse(n) for n in ast.walk(tree)
               if isinstance(n, ast.AsyncFunctionDef) and n.name == "lifespan")
    assert "warm_outbound" in src


def test_concurrency_default_is_wide_enough_for_every_host(monkeypatch):
    """⚠ MEASURED. `[netwarm] layer totals — DNS 27.16s, TCP connect 1.63s`
    across 24 hosts, worst single host 5325 ms. At the old limit of 4, that
    27 s of resolution serialised into 9.05 s of start-up — and a curation that
    arrived meanwhile contended for the same resolver.

    These probes are pure latency: they hold no CPU and no mount. So the limit
    should be wide enough that total time is bounded by the SLOWEST host rather
    than by the queue.
    """
    monkeypatch.delenv("HEARTVAR_WARM_OUTBOUND_CONCURRENCY", raising=False)
    assert netwarm._concurrency() >= len(netwarm.outbound_hosts()), (
        "concurrency must cover every host, or DNS latency serialises"
    )


def test_netwarm_runs_before_the_mount_warm():
    """DNS is the dominant cold cost (27 s) and the mount warm is 0.5 s, so
    warming the mount first is backwards — it delays the expensive one."""
    import ast
    from pathlib import Path
    tree = ast.parse((Path(__file__).resolve().parents[1] / "app.py")
                     .read_text(encoding="utf-8"))
    src = next(ast.unparse(n) for n in ast.walk(tree)
               if isinstance(n, ast.AsyncFunctionDef) and n.name == "lifespan")
    assert src.index("warm_outbound") < src.index("warm_local_data"), (
        "warm_outbound must be scheduled before warm_local_data — DNS is 27 s "
        "of the cold cost, the mount is 0.5 s"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
