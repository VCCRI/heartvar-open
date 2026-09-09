"""Inbound rate-limiting, trusted-proxy client-IP derivation, and lightweight
request metrics for the HeartVar backend.

Extracted verbatim from backend/app.py as a pure code-move. This module owns the
slowapi Limiter instance, the client-IP derivation used as its key function, the
global daily-cap backstop, and the streaming-safe metrics middleware. backend/app.py
re-imports these names (so `from backend.app import _client_ip` etc. keep working)
and performs the FastAPI wiring (app.state.limiter / exception handler / middleware /
the @limiter.limit decorator), which must stay next to the `app` object."""

from __future__ import annotations

import ipaddress
import logging
import os
from collections import deque
from datetime import datetime, timezone
from time import perf_counter

from fastapi import HTTPException, Request
from slowapi import Limiter
from slowapi.util import get_remote_address

log = logging.getLogger("heartvar.ratelimit")


def _parse_proxy_nets(raw: str):
    nets, trust_all = [], False
    for entry in (raw or "").split(","):
        entry = entry.strip()
        if not entry:
            continue
        if entry == "*":
            trust_all = True
            continue
        try:
            net = ipaddress.ip_network(entry, strict=False)
        except ValueError:
            log.warning("Ignoring invalid HEARTVAR_TRUSTED_PROXY_IPS entry: %r", entry)
            continue
        if net.prefixlen == 0:
            log.warning(
                "HEARTVAR_TRUSTED_PROXY_IPS entry %r trusts ALL addresses — this "
                "makes X-Forwarded-For spoofable; list specific proxy ranges instead.",
                entry,
            )
        nets.append(net)
    return nets, trust_all


_TRUSTED_PROXY_NETS, _TRUST_ALL_PROXIES = _parse_proxy_nets(
    os.environ.get("HEARTVAR_TRUSTED_PROXY_IPS", "")
)
_CLIENT_IP_HEADER = (
    os.environ.get("HEARTVAR_CLIENT_IP_HEADER", "x-azure-client-ip").strip().lower()
    or "x-azure-client-ip"
)
_PROXY_CONFIGURED = bool(_TRUSTED_PROXY_NETS) or _TRUST_ALL_PROXIES


def _strip_port(hop: str) -> str:
    """Drop a trailing :port from an XFF hop ('1.2.3.4:5678' or '[::1]:443')."""
    hop = hop.strip()
    if hop.startswith("["):
        return hop[1:].split("]", 1)[0]
    if hop.count(":") == 1:
        return hop.split(":", 1)[0]
    return hop


def _valid_ip(ip_str: str):
    try:
        ipaddress.ip_address(ip_str)
        return ip_str
    except ValueError:
        return None


def _is_trusted_proxy(ip_str: str) -> bool:
    if _TRUST_ALL_PROXIES:
        return True
    try:
        ip = ipaddress.ip_address(ip_str)
    except ValueError:
        return False
    return any(ip in net for net in _TRUSTED_PROXY_NETS)


def _client_ip(request: Request) -> str:
    """Best-effort real client IP for rate-limit keying. Trusts forwarded
    headers only when HEARTVAR_TRUSTED_PROXY_IPS is configured AND the immediate
    peer is itself a trusted proxy; otherwise keys on the socket peer (fail
    closed — not spoofable). See the module notes above."""
    peer = get_remote_address(request) or "unknown"
    if not _PROXY_CONFIGURED:
        return peer
    if not (_TRUST_ALL_PROXIES or _is_trusted_proxy(peer)):
        return peer
    hdr = request.headers.get(_CLIENT_IP_HEADER)
    if hdr:
        cand = _valid_ip(_strip_port(hdr.split(",")[0]))
        if cand:
            return cand
    fwd = request.headers.get("x-forwarded-for")
    if fwd:
        hops = [h for h in (p.strip() for p in fwd.split(",")) if h]
        for hop in reversed(hops):
            ip_only = _strip_port(hop)
            if not _is_trusted_proxy(ip_only):
                cand = _valid_ip(ip_only)
                if cand:
                    return cand
    return peer


def _session_account_id(request: Request) -> str | None:
    """Identity of the signed-in Microsoft Entra ID account, or None.

    Reads the session cookie directly rather than importing from backend.app,
    which would be a circular import (app imports this module). Guarded because
    request.session raises AssertionError when SessionMiddleware is absent — auth
    not being configured must mean "anonymous", never a 500.
    """
    try:
        account = request.session.get("account")
    except (AssertionError, AttributeError):
        return None
    if not account:
        return None
    for claim in ("oid", "sub", "preferred_username", "email"):
        value = account.get(claim)
        if isinstance(value, str) and value.strip():
            return f"u:{claim}:{value.strip()}"
    return "u:unidentified"


def _rate_limit_key(request: Request) -> str:
    """Rate-limit bucket: the signed-in user, else the client IP.

    Keying on identity when we have one is strictly better than keying on IP:

      * A shared institutional egress puts every curator at one hospital behind a
        single address, so an IP bucket throttles colleagues as a group.
      * It cannot be rotated. IP limits fall to anyone with a handful of
        addresses; a session is tied to a real Microsoft account.
      * It sidesteps the unresolved trusted-proxy problem for authenticated
        traffic: until HEARTVAR_TRUSTED_PROXY_IPS is set, _client_ip() correctly
        fails closed to the App Service front end, collapsing all users into one
        bucket. Signed-in users get their own regardless of that setting.

    Anonymous traffic still keys on _client_ip(), unchanged.
    """
    return _session_account_id(request) or _client_ip(request)


_RATE_PER_MIN = os.environ.get("HEARTVAR_RATE_PER_MIN", "1000").strip() or "1000"
_RATE_PER_HOUR = os.environ.get("HEARTVAR_RATE_PER_HOUR", "1000").strip() or "1000"
_RATE_PER_DAY = os.environ.get("HEARTVAR_RATE_PER_DAY", "1000").strip() or "1000"
_CURATE_RATE_LIMIT = f"{_RATE_PER_MIN}/minute;{_RATE_PER_HOUR}/hour;{_RATE_PER_DAY}/day"

_CHAT_PER_MIN = os.environ.get("HEARTVAR_CHAT_RATE_PER_MIN", "1000").strip() or "1000"
_CHAT_PER_HOUR = os.environ.get("HEARTVAR_CHAT_RATE_PER_HOUR", "1000").strip() or "1000"
_CHAT_PER_DAY = os.environ.get("HEARTVAR_CHAT_RATE_PER_DAY", "1000").strip() or "1000"
_CHAT_RATE_LIMIT = f"{_CHAT_PER_MIN}/minute;{_CHAT_PER_HOUR}/hour;{_CHAT_PER_DAY}/day"
_REDIS_URL = os.environ.get("HEARTVAR_REDIS_URL", "").strip()
if _REDIS_URL:
    limiter = Limiter(key_func=_rate_limit_key, storage_uri=_REDIS_URL)
    log.info("Inbound rate limiter using shared Redis storage")
else:
    limiter = Limiter(key_func=_rate_limit_key)


_RATE_ALERT_PER_MIN = int(os.environ.get("HEARTVAR_429_ALERT_PER_MIN", "50") or "50")
_GLOBAL_DAILY_CAP = int(os.environ.get("HEARTVAR_GLOBAL_DAILY_CAP", "1000") or "500")


class _RequestMetrics:
    def __init__(self) -> None:
        self.total = 0
        self.by_status: dict[int, int] = {}
        self._recent_429: deque[float] = deque()
        self.curate_day: str | None = None
        self.curate_day_count = 0

    def record(self, status: int) -> None:
        self.total += 1
        self.by_status[status] = self.by_status.get(status, 0) + 1

    def record_429(self) -> int:
        now = perf_counter()
        self._recent_429.append(now)
        cutoff = now - 60.0
        while self._recent_429 and self._recent_429[0] < cutoff:
            self._recent_429.popleft()
        return len(self._recent_429)


METRICS = _RequestMetrics()


def _check_daily_cap() -> None:
    """Topology-independent global backstop on the expensive endpoints. ON by
    default (HEARTVAR_GLOBAL_DAILY_CAP=500; set 0 to disable). Caps total
    AI-bearing requests (curation + chat) per UTC day across ALL callers — so
    even a perfectly-spoofed / IP-rotating flood can't exceed the day's ceiling.
    The counter is shared by /api/curate/stream and /api/chat; the user-facing
    429 message is phrased for curation but applies to both."""
    if _GLOBAL_DAILY_CAP <= 0:
        return
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    if METRICS.curate_day != today:
        METRICS.curate_day = today
        METRICS.curate_day_count = 0
    METRICS.curate_day_count += 1
    if METRICS.curate_day_count > _GLOBAL_DAILY_CAP:
        raise HTTPException(
            429, "Daily curation capacity has been reached. Please try again tomorrow."
        )


_AI_PER_USER_DAILY = int(os.environ.get("HEARTVAR_AI_PER_USER_DAILY", "1000") or "50")
_AI_CHAT_PER_USER_DAILY = int(
    os.environ.get("HEARTVAR_AI_CHAT_PER_USER_DAILY", "1000") or "60"
)
_AI_QUOTA_CAPS = {"curation": _AI_PER_USER_DAILY, "chat": _AI_CHAT_PER_USER_DAILY}
_AI_QUOTA_MESSAGES = {
    "curation": (
        "You have used your daily allowance of {cap} AI interpretations. The "
        "evidence and the rule-based classification are still available — untick "
        "\u201cInclude AI interpretation\u201d to continue, or try again tomorrow."
    ),
    "chat": (
        "You have used your daily allowance of {cap} assistant questions. The "
        "variant's evidence and classification remain available above; the "
        "assistant is back tomorrow."
    ),
}
_ai_quota: dict = {"day": None, "counts": {}}


def check_ai_user_quota(key: str, kind: str = "curation") -> None:
    """Charge one AI request to `key`, raising 429 once its daily quota is spent.

    `kind` selects the cap and the user-facing message: "curation" for an
    AI-bearing /api/curate/stream, "chat" for /api/chat. Set the corresponding env
    var to 0 to disable that quota.

    Counted at request START, so a curation that ends up degrading to evidence-only
    (budget exhausted mid-day, Anthropic error) still consumes quota. Deliberately
    conservative: metering after the fact would let a client that abandons the
    stream retry indefinitely, billing every attempt.
    """
    cap = _AI_QUOTA_CAPS.get(kind, _AI_PER_USER_DAILY)
    if cap <= 0:
        return
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    if _ai_quota["day"] != today:
        _ai_quota["day"] = today
        _ai_quota["counts"].clear()
    bucket = f"{kind}:{key}"
    spent = _ai_quota["counts"].get(bucket, 0) + 1
    _ai_quota["counts"][bucket] = spent
    if spent > cap:
        log.warning(
            "Per-user daily AI quota exhausted (kind=%s key=%s cap=%d)", kind, key, cap,
        )
        raise HTTPException(
            429,
            detail={
                "error": "ai_quota_exhausted",
                "message": _AI_QUOTA_MESSAGES[kind].format(cap=cap),
            },
        )


class _MetricsMiddleware:
    """Pure-ASGI middleware (NOT BaseHTTPMiddleware, which would buffer the SSE
    stream). Records the response status and logs a warning/alert on 429s."""

    def __init__(self, app) -> None:
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope.get("type") != "http":
            await self.app(scope, receive, send)
            return
        status = {"code": 0}

        async def _send(message):
            if message["type"] == "http.response.start":
                status["code"] = message["status"]
            await send(message)

        await self.app(scope, receive, _send)
        try:
            code = status["code"]
            METRICS.record(code)
            if code == 429:
                rate = METRICS.record_429()
                ip = _client_ip(Request(scope))
                path = scope.get("path", "")
                if rate >= _RATE_ALERT_PER_MIN:
                    log.error("429 SPIKE: %d rejections in 60s (path=%s ip=%s)", rate, path, ip)
                else:
                    log.warning("Rate-limit 429 (path=%s ip=%s, %d in last 60s)", path, ip, rate)
        except Exception:  # noqa: BLE001 — metrics must never break a response
            pass
