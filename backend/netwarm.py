"""Warm — and above all MEASURE — outbound connectivity at start-up.

WHAT THIS IS FOR. After the mount was warmed (PR #53/#54) the deployment log of
2026-08-30 10:06 read:

    [warmup] warmed 14/14 local files (concurrency 2) — total file-time 0.49s
    [warmup] resolver + FASTA ready in 0.34s
    [timing] DB gather 23.38s (erepo=18.00s, mgi=18.00s, vep=18.00s,
             gene_literature=18.00s, pubtator3=18.00s, medgen=18.00s,
             clinvar=18.00s, biogrid=18.00s, pubmed=18.00s, fetal_heart=18.00s,
             alphafold=18.00s, gtex=18.00s, uniprot=18.00s, panelapp=18.00s,
             ... gnomad=0.95s, gencc=0.41s, chdgene=0.00s)
    [timing] ... same_site=0.04s        <- was 15.87 s before warming

The mount is fast. Warming worked. But FOURTEEN sources took exactly 18.00 s on
the first curation and 0.19 s on the second, and every one of them makes an
outbound HTTP call — as its primary source, or as the live fallback that
``HEARTVAR_OFFLINE_STRICT`` being unset still leaves on for 12 clients. Every
unaffected source is local-only.

⚠ AN INSTRUMENT BEFORE IT IS A FIX, on purpose. Identical 18.00 s across hosts
that resolve to DIFFERENT servers does NOT fit plain per-host DNS latency; they
release together, which points at something shared. Two hypotheses have already
died on this stall (the handoff's sync ``json.load``s at 250 ms of CPU, and then
"warming is too shallow", refuted by the 0.49 s above), so this measures each
layer separately instead of assuming one:

    getaddrinfo per host  -> DNS
    TCP connect per host  -> reachability, IPv6 fallback, SNAT exhaustion

Whichever column carries the seconds is the answer, and it lands in the
deployment log that /api/admin/log-content serves — no Azure access needed.

⚠ NO HTTP REQUEST IS MADE, and a test enforces it. DNS resolution and a TCP
connect that is immediately closed touch no API, spend no rate limit, and send no
credentials. Warming the successful lookups is the likely fix and comes free with
the measurement.
"""
from __future__ import annotations

import asyncio
import logging
import os
import socket
import time

log = logging.getLogger("heartvar.netwarm")

_HOSTS: tuple[str, ...] = (
    "eutils.ncbi.nlm.nih.gov",
    "www.ncbi.nlm.nih.gov",
    "pubmed.ncbi.nlm.nih.gov",
    "pmc.ncbi.nlm.nih.gov",
    "rest.ensembl.org",
    "www.ebi.ac.uk",
    "panelapp-aus.org",
    "api.platform.opentargets.org",
    "thegencc.org",
    "search.thegencc.org",
    "www.informatics.jax.org",
    "thebiogrid.org",
    "gtexportal.org",
    "rest.uniprot.org",
    "www.uniprot.org",
    "rest.genenames.org",
    "www.omim.org",
    "www.alliancegenome.org",
    "chdgene.victorchang.edu.au",
    "spliceai-38-xwkwwwxdwq-uc.a.run.app",
    "spliceai-37-xwkwwwxdwq-uc.a.run.app",
    "login.microsoftonline.com",
    "accounts.google.com",
)

_DEFAULT_CONCURRENCY = 32
_DEFAULT_TIMEOUT = 10.0


def outbound_hosts() -> list[str]:
    """The hosts to probe, deduplicated and order-preserving."""
    seen: set[str] = set()
    out: list[str] = []
    for h in _HOSTS:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return out


def outbound_warming_enabled() -> bool:
    """ON by default. ``0``/``false``/``no``/``off`` disables it."""
    return os.environ.get("HEARTVAR_WARM_OUTBOUND", "").strip().lower() \
        not in {"0", "false", "no", "off"}


def _concurrency() -> int:
    """How many hosts to probe at once. See _DEFAULT_CONCURRENCY for why it is
    wide rather than polite: DNS probes are latency, not load."""
    raw = os.environ.get("HEARTVAR_WARM_OUTBOUND_CONCURRENCY", "").strip()
    try:
        value = int(raw)
        return value if value > 0 else _DEFAULT_CONCURRENCY
    except ValueError:
        return _DEFAULT_CONCURRENCY


async def _probe(host: str, timeout: float) -> dict:
    """Time DNS and TCP for one host. Never raises.

    The two are timed separately because that separation is the whole point: a
    slow ``dns_ms`` with a fast ``connect_ms`` means resolution, the reverse means
    reachability, and both slow means the egress path.
    """
    loop = asyncio.get_running_loop()
    entry: dict = {"dns_ms": 0.0, "connect_ms": 0.0, "error": "",
                   "family": "", "addr": ""}

    t0 = time.perf_counter()
    try:
        infos = await asyncio.wait_for(
            loop.getaddrinfo(host, 443, type=socket.SOCK_STREAM),
            timeout=timeout,
        )
    except Exception as exc:  # noqa: BLE001 — probe must never raise
        entry["dns_ms"] = (time.perf_counter() - t0) * 1000
        entry["error"] = f"dns: {type(exc).__name__}"
        return entry
    entry["dns_ms"] = (time.perf_counter() - t0) * 1000
    if not infos:
        entry["error"] = "dns: no records"
        return entry
    family, _, _, _, sockaddr = infos[0]
    entry["family"] = "IPv6" if family == socket.AF_INET6 else "IPv4"
    entry["addr"] = str(sockaddr[0])

    t1 = time.perf_counter()
    writer = None
    try:
        _, writer = await asyncio.wait_for(
            asyncio.open_connection(sockaddr[0], 443), timeout=timeout)
    except Exception as exc:  # noqa: BLE001 — probe must never raise
        entry["connect_ms"] = (time.perf_counter() - t1) * 1000
        entry["error"] = f"connect: {type(exc).__name__}"
        return entry
    finally:
        if writer is not None:
            try:
                writer.close()
            except Exception:  # noqa: BLE001
                pass
    entry["connect_ms"] = (time.perf_counter() - t1) * 1000
    return entry


async def warm_outbound(hosts: list[str] | None = None) -> dict:
    """Probe every outbound host, log the timings, return the report.

    Bounded concurrency so start-up does not open two dozen sockets at once.
    Best effort throughout: this is a background task and a probe failure must
    never affect the server.
    """
    targets = hosts if hosts is not None else outbound_hosts()
    timeout = _DEFAULT_TIMEOUT
    gate = asyncio.Semaphore(_concurrency())
    report: dict = {"hosts": {}}

    async def one(host: str) -> None:
        async with gate:
            report["hosts"][host] = await _probe(host, timeout)

    t0 = time.perf_counter()
    await asyncio.gather(*(one(h) for h in targets))
    total = time.perf_counter() - t0

    entries = report["hosts"]
    slow = sorted(entries.items(),
                  key=lambda kv: -(kv[1]["dns_ms"] + kv[1]["connect_ms"]))
    failed = [h for h, e in entries.items() if e["error"]]
    log.info(
        "[netwarm] probed %d hosts in %.2fs (%d failed) — slowest: %s",
        len(entries), total, len(failed),
        ", ".join(f"{h} dns={e['dns_ms']:.0f}ms conn={e['connect_ms']:.0f}ms"
                  f"{' ' + e['family'] if e['family'] else ''}"
                  f"{' ERR=' + e['error'] if e['error'] else ''}"
                  for h, e in slow[:8]) or "none",
    )
    if failed:
        log.warning("[netwarm] %d host(s) unreachable at start-up: %s",
                    len(failed),
                    ", ".join(f"{h} ({entries[h]['error']})" for h in failed))
    dns_total = sum(e["dns_ms"] for e in entries.values()) / 1000
    conn_total = sum(e["connect_ms"] for e in entries.values()) / 1000
    log.info("[netwarm] layer totals — DNS %.2fs, TCP connect %.2fs across %d "
             "hosts. If a first curation still costs ~18s while BOTH of these "
             "are small, the cost is NOT resolution or reachability — look at "
             "the process-wide outbound pacers (_ncbi_throttle, spliceai's "
             "global concurrency cap).", dns_total, conn_total, len(entries))
    return report
