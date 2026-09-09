"""Tests for the hardened client-IP derivation, the global daily cap, and the
request metrics (audit §11). The vulnerability fixed here: the old _client_ip
trusted the first X-Forwarded-For hop unconditionally, so any caller could set
XFF to a random IP and win a fresh rate-limit bucket per request."""

import ipaddress

import pytest
from fastapi import HTTPException
from starlette.requests import Request

import backend.app as app
import backend.ratelimit as ratelimit


def _req(peer="203.0.113.9", headers=None):
    hdrs = [(k.lower().encode(), v.encode()) for k, v in (headers or {}).items()]
    scope = {
        "type": "http",
        "method": "GET",
        "path": "/api/curate/stream",
        "headers": hdrs,
        "client": (peer, 51234),
        "scheme": "http",
        "server": ("testserver", 80),
    }
    return Request(scope)


def _configure(monkeypatch, *, nets=None, trust_all=False, header="x-azure-client-ip"):
    monkeypatch.setattr(ratelimit, "_TRUSTED_PROXY_NETS",
                        [ipaddress.ip_network(n, strict=False) for n in (nets or [])])
    monkeypatch.setattr(ratelimit, "_TRUST_ALL_PROXIES", trust_all)
    monkeypatch.setattr(ratelimit, "_PROXY_CONFIGURED", bool(nets) or trust_all)
    monkeypatch.setattr(ratelimit, "_CLIENT_IP_HEADER", header)


def test_default_ignores_spoofed_xff(monkeypatch):
    """Unconfigured (no trusted proxy): forwarded headers are ignored and we key
    on the socket peer — the spoof attempt is neutralised."""
    _configure(monkeypatch, nets=[])
    r = _req(peer="198.51.100.7", headers={
        "x-forwarded-for": "1.2.3.4",
        "x-azure-client-ip": "5.6.7.8",
    })
    assert app._client_ip(r) == "198.51.100.7"


def test_default_two_spoof_attempts_share_one_bucket(monkeypatch):
    """Two requests with different forged XFF from the same peer key identically
    (so they can't each win a fresh bucket)."""
    _configure(monkeypatch, nets=[])
    a = app._client_ip(_req(peer="198.51.100.7", headers={"x-forwarded-for": "1.1.1.1"}))
    b = app._client_ip(_req(peer="198.51.100.7", headers={"x-forwarded-for": "2.2.2.2"}))
    assert a == b == "198.51.100.7"


def test_trusted_proxy_prefers_platform_header(monkeypatch):
    _configure(monkeypatch, nets=["10.0.0.0/8"])
    r = _req(peer="10.0.0.5", headers={
        "x-azure-client-ip": "203.0.113.42",
        "x-forwarded-for": "203.0.113.42, 10.0.0.5",
    })
    assert app._client_ip(r) == "203.0.113.42"


def test_trusted_proxy_xff_rightmost_untrusted(monkeypatch):
    """XFF = client, proxyA, proxyB (peer). Trusted hops stripped right→left →
    the real client is the right-most untrusted hop."""
    _configure(monkeypatch, nets=["10.0.0.0/8"])
    r = _req(peer="10.0.0.9", headers={
        "x-forwarded-for": "203.0.113.42, 10.0.0.1, 10.0.0.9",
    })
    assert app._client_ip(r) == "203.0.113.42"


def test_trusted_proxy_strips_ipv4_port(monkeypatch):
    """Azure App Service appends the client as ip:port on XFF."""
    _configure(monkeypatch, nets=["10.0.0.0/8"])
    r = _req(peer="10.0.0.9", headers={
        "x-forwarded-for": "203.0.113.42:6010, 10.0.0.9",
    })
    assert app._client_ip(r) == "203.0.113.42"


def test_untrusted_peer_ignores_forged_headers(monkeypatch):
    """Proxy IS configured, but the request arrives from an UNtrusted peer (a
    direct hit bypassing the proxy). Forged headers must be ignored."""
    _configure(monkeypatch, nets=["10.0.0.0/8"])
    r = _req(peer="45.66.77.88", headers={
        "x-azure-client-ip": "1.1.1.1",
        "x-forwarded-for": "2.2.2.2",
    })
    assert app._client_ip(r) == "45.66.77.88"


def test_invalid_header_value_falls_back(monkeypatch):
    _configure(monkeypatch, nets=["10.0.0.0/8"])
    r = _req(peer="10.0.0.5", headers={"x-azure-client-ip": "not-an-ip"})
    assert app._client_ip(r) == "10.0.0.5"


def test_trust_all_with_azure_header(monkeypatch):
    _configure(monkeypatch, trust_all=True)
    r = _req(peer="10.0.0.5", headers={"x-azure-client-ip": "203.0.113.99"})
    assert app._client_ip(r) == "203.0.113.99"


def test_trust_all_without_header_falls_closed_to_peer(monkeypatch):
    """trust_all="*" but NO platform header → must NOT trust the leftmost
    (client-controllable) XFF hop; fail closed to the socket peer so an attacker
    can't win a fresh bucket per forged XFF. (Review finding: medium.)"""
    _configure(monkeypatch, trust_all=True)
    a = app._client_ip(_req(peer="169.254.130.1",
                            headers={"x-forwarded-for": "203.0.113.1, 1.2.3.4:5678"}))
    b = app._client_ip(_req(peer="169.254.130.1",
                            headers={"x-forwarded-for": "9.9.9.9, 1.2.3.4:5678"}))
    assert a == b == "169.254.130.1"


def test_over_broad_cidr_falls_closed_to_peer(monkeypatch):
    """An over-broad trusted CIDR (every hop parses as trusted) must not trust
    the spoofable leftmost XFF hop either. (Review finding: low.)"""
    _configure(monkeypatch, nets=["0.0.0.0/0"])
    r = _req(peer="10.0.0.5", headers={"x-forwarded-for": "203.0.113.1, 10.0.0.5"})
    assert app._client_ip(r) == "10.0.0.5"


def test_daily_cap_disabled_by_default(monkeypatch):
    monkeypatch.setattr(ratelimit, "_GLOBAL_DAILY_CAP", 0)
    monkeypatch.setattr(ratelimit, "METRICS", app._RequestMetrics())
    for _ in range(1000):
        app._check_daily_cap()


def test_daily_cap_enforced(monkeypatch):
    monkeypatch.setattr(ratelimit, "_GLOBAL_DAILY_CAP", 3)
    monkeypatch.setattr(ratelimit, "METRICS", app._RequestMetrics())
    for _ in range(3):
        app._check_daily_cap()
    with pytest.raises(HTTPException) as exc:
        app._check_daily_cap()
    assert exc.value.status_code == 429


def test_metrics_record_and_429_window():
    m = app._RequestMetrics()
    m.record(200)
    m.record(429)
    assert m.total == 2
    assert m.by_status[429] == 1
    assert m.record_429() >= 1
