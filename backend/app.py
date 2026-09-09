from __future__ import annotations

import asyncio
import base64
import contextlib
import hashlib
import hmac
import json
import logging
import os
import re
import secrets
import time
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import AsyncIterator
from urllib.parse import urlencode

import httpx
import msal

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException, Request
from fastapi.exception_handlers import http_exception_handler
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    RedirectResponse,
    StreamingResponse,
)
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.middleware.sessions import SessionMiddleware
from slowapi import _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded

load_dotenv()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
log = logging.getLogger("heartvar.app")

from .logredact import (  # noqa: E402
    begin_curation as _log_begin_curation,
    install as _install_log_redaction,
)

_install_log_redaction()


from .localio import log_safe as _log_safe


from . import budget
from .claude import (
    CHAT_MODEL,
    STREAM_RESTART,
    chat_reply,
    parse_first_json_object,
    report_mode_caps,
    stream_claude,
)
from .models import (
    ChatRequest,
    CurationRequest,
    InterpretRequest,
    RescoreRequest,
)
from . import interpret_cache
from .loop_probe import maybe_enable_loop_probe
from .netwarm import outbound_warming_enabled, warm_outbound
from .warmup import warm_local_data, warm_resolver, warming_enabled
from .prewarm import maybe_start_prewarm, note_activity
from .clients.biogrid import fetch_biogrid
from .clients.vep_offline import offline_status as vep_offline_status
from .clients.clinvar import fetch_clinvar
from .clients.ensembl_vep import fetch_transcript_exons
from .clients.fetal_heart import fetch_fetal_heart
from .clients.hgnc_alias import canonicalise_gene_symbol
from .clients.panelapp import _ensure_descendant_cache
from .prompt import build_prompt

# keep importing them via ``from backend.app import X`` AND via attribute access
__all__ = [
    "app",
    "CurationRequest",
    "_TIER_POINTS",
    "_BARE_CODE_POINTS",
    "HARD_CODED_CRITERIA_CODES",
    "AI_EVALUATED_CRITERIA_CODES",
    "_CRITERION_NAMES",
    "CANONICAL_CRITERIA_ORDER",
    "_CLINVAR_STAR_STRENGTH",
    "_points_for",
    "compute_points_total",
    "classification_for",
    "apply_benign_combining_floor",
    "classification_for_criteria",
    "tier_without_clinvar_assertion",
    "apply_gene_validity_ceiling",
    "_VCEP_FREQ",
    "_any_protein_position_from_vep",
    "_bs1_threshold",
    "_clinvar_pp5_bp6_criteria",
    "_criterion_applicable",
    "_gate_criteria_applicability",
    "_gene_lof_mechanism",
    "_hgvsp_int_position",
    "_parse_pm1_hotspots",
    "_pm2_threshold",
    "_select_clinvar_assertion_record",
    "apply_cross_criterion_exclusions",
    "classify_family_history",
    "compute_hard_coded_criteria",
    "gene_mechanism",
    "mechanism_consistency_flag",
    "gene_inheritance_modes",
    "gene_validity_ceiling",
    "carrier_status",
    "merge_hard_coded_and_ai",
    "build_no_ai_criteria",
    "infer_supplementary_criteria",
    "_RequestMetrics",
    "_check_daily_cap",
    "_client_ip",
    "_gather_evidence_sse",
    "_protein_position_from_vep",
]

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _load_anthropic_key() -> None:
    """Load the owner's Anthropic key so every server-side AI call bills to one
    account. Precedence: a gitignored ``API.txt`` at the repo root (the local
    dev path — a single bare ``sk-ant-…`` line) is AUTHORITATIVE when present;
    otherwise fall back to the ``ANTHROPIC_API_KEY`` environment variable (the
    deploy path, since API.txt is gitignored AND dockerignored and never ships
    in the image). The raw key value is NEVER logged."""
    api_txt = PROJECT_ROOT / "API.txt"
    try:
        if api_txt.is_file():
            key = api_txt.read_text(encoding="utf-8").strip()
            if key:
                os.environ["ANTHROPIC_API_KEY"] = key
                log.info("Anthropic API key loaded from API.txt")
                return
    except OSError as e:
        log.warning("Could not read API.txt: %s", e)
    if os.environ.get("ANTHROPIC_API_KEY", "").strip():
        log.info("Anthropic API key loaded from the environment")
    else:
        log.warning(
            "No Anthropic API key found (API.txt missing/empty and "
            "ANTHROPIC_API_KEY unset). Server-side AI is DISABLED — the app "
            "serves evidence-only (deterministic) results until a key is set."
        )


_load_anthropic_key()
SERVER_AI_AVAILABLE = bool(os.environ.get("ANTHROPIC_API_KEY", "").strip())

try:
    RESULT_MAX_TOKENS = int(os.environ.get("HEARTVAR_RESULT_MAX_TOKENS", "") or 4000)
except ValueError:
    RESULT_MAX_TOKENS = 4000

_CHAT_INPUT_EST_TOKENS = 4000


NO_AI_SUPPLEMENTARY = os.environ.get(
    "HEARTVAR_NO_AI_SUPPLEMENTARY", "1"
).strip().lower() in ("1", "true", "yes", "on")


from .ratelimit import (
    _CHAT_RATE_LIMIT,
    _CURATE_RATE_LIMIT,
    _MetricsMiddleware,
    _RequestMetrics,
    _check_daily_cap,
    _client_ip,
    check_ai_user_quota,
    limiter,
)


from .acmg.constants import (
    _TIER_POINTS,
    _BARE_CODE_POINTS,
    HARD_CODED_CRITERIA_CODES,
    AI_EVALUATED_CRITERIA_CODES,
    _CRITERION_NAMES,
    CANONICAL_CRITERIA_ORDER,
    _CLINVAR_STAR_STRENGTH,
)
from .acmg.tiers import (
    _points_for,
    compute_points_total,
    classification_for,
    apply_benign_combining_floor,
    apply_gene_validity_ceiling,
    classification_for_criteria,
    tier_without_clinvar_assertion,
)
from .acmg.hard_coded import (
    apply_pm1_collision_rules,
    cap_all_to_spec_strength,
    _normalize_criteria,
    _select_clinvar_assertion_record,
    _clinvar_pp5_bp6_criteria,
    _VCEP_FREQ,
    _criterion_applicable,
    _parse_pm1_hotspots,
    _gate_criteria_applicability,
    _bs1_threshold,
    _pm2_threshold,
    _gene_lof_mechanism,
    gene_mechanism,
    mechanism_consistency_flag,
    gene_inheritance_modes,
    gene_validity_ceiling,
    carrier_status,
    compute_hard_coded_criteria,
    apply_cross_criterion_exclusions,
    merge_hard_coded_and_ai,
    classify_family_history,
    _hgvsp_int_position,
    _any_protein_position_from_vep,
)
from .acmg.no_ai import (
    build_no_ai_criteria,
    infer_supplementary_criteria,
    server_owned_placeholder,
)


def _variant_chromosome(evidence: dict | None) -> str:
    """Pull the variant chromosome out of the VEP block for downstream
    inference (notably hemizygosity inference on chrX). Returns the
    bare chromosome ("X", "Y", "1", …) uppercased, or "" when VEP
    didn't return a usable seq_region_name (failed lookup, malformed
    response). Strips any leading "chr" so callers can compare against
    "X" without worrying about prefix variants."""
    if not isinstance(evidence, dict):
        return ""
    vep = evidence.get("vep") if isinstance(evidence.get("vep"), dict) else {}
    chrom = (vep.get("seq_region_name") or vep.get("chromosome") or "")
    if isinstance(chrom, str):
        chrom = chrom.strip()
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        return chrom.upper()
    return ""


def _build_clinical_context(req: "CurationRequest", evidence: dict | None) -> dict:
    """Assemble the structured clinical_context dict consumed by the
    AI prompt builder. Runs three pieces of light inference over the
    raw form inputs so the prompt builder gets a single normalised
    object rather than re-deriving these flags inline:

    1. Hemizygosity inference: a male proband with an unknown
       zygosity on a chrX variant is auto-marked hemizygous and
       flagged so the prompt + UI can label the value "inferred".
    2. De novo consistency: "De novo" inheritance without parental
       data is downgraded to an explicit unconfirmed assumption +
       a note. Confirmed de novo paired with trio/duo data is
       upgraded into a denovo_confirmed flag the prompt builder
       uses to switch from PM6 to PS2 guidance.
    3. AR + het: a het variant under AR inheritance flips the
       compound_het_possible flag so the prompt nudges PM3/BP2
       handling toward phase-aware reasoning.

    The proband is by definition affected — every record always
    carries proband_affected=True so the prompt builder doesn't have
    to special-case the missing flag.
    """
    zygosity = (req.zygosity or "").strip().lower()
    inheritance_input = (req.inheritance_input or "").strip().upper()
    proband_sex = (req.proband_sex or "").strip().lower()
    trio_status = (req.trio_status or "").strip().lower()
    denovo_status = (req.denovo_status or "").strip().lower()
    notes: list[str] = []

    chrom = _variant_chromosome(evidence)
    zygosity_inferred = False
    if not zygosity and proband_sex == "male" and chrom == "X":
        zygosity = "hemi"
        zygosity_inferred = True

    if inheritance_input == "DN" and not trio_status:
        if not denovo_status:
            denovo_status = "unconfirmed"
        notes.append(
            "De novo selected without parental data — treated as unconfirmed"
        )

    denovo_confirmed = (
        denovo_status == "confirmed" and trio_status != "duo"
    )
    denovo_trio_stated = trio_status == "trio"

    compound_het_possible = inheritance_input == "AR" and zygosity == "het"

    def _nn_int(value: object) -> int:
        try:
            return max(0, int(value))  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return 0

    seg_affected_carriers = _nn_int(getattr(req, "seg_affected_carriers", 0))
    if denovo_status == "inherited_affected":
        seg_affected_carriers = max(seg_affected_carriers, 1)
    seg_affected_noncarriers = _nn_int(getattr(req, "seg_affected_noncarriers", 0))
    seg_meioses = _nn_int(getattr(req, "seg_meioses", 0))
    in_trans_pathogenic = (getattr(req, "in_trans_pathogenic", "") or "").strip().lower()
    alt_cause_present = (getattr(req, "alt_cause_present", "") or "").strip().lower()
    alt_cause_detail = (getattr(req, "alt_cause_detail", "") or "").strip()

    return {
        "zygosity": zygosity,
        "zygosity_inferred": zygosity_inferred,
        "inheritance_input": inheritance_input,
        "proband_sex": proband_sex,
        "trio_status": trio_status,
        "denovo_status": denovo_status,
        "denovo_confirmed": denovo_confirmed,
        "denovo_trio_stated": denovo_trio_stated,
        "denovo_confirmed_count": _nn_int(
            getattr(req, "denovo_confirmed_count", 0)),
        "denovo_unconfirmed_count": _nn_int(
            getattr(req, "denovo_unconfirmed_count", 0)),
        "compound_het_possible": compound_het_possible,
        "proband_affected": True,
        "chromosome": chrom,
        "family_history_summary": classify_family_history(req.family),
        "seg_affected_carriers": seg_affected_carriers,
        "seg_affected_noncarriers": seg_affected_noncarriers,
        "seg_meioses": seg_meioses,
        "in_trans_pathogenic": in_trans_pathogenic,
        "alt_cause_present": alt_cause_present,
        "alt_cause_detail": alt_cause_detail,
        "notes": notes,
    }


@contextlib.asynccontextmanager
async def lifespan(app: FastAPI):
    """Fire the JAX HPO descendant cache warm-up in the background so the
    first real curation doesn't pay the ~9 s cold-start cost. Scheduled,
    not awaited — the server is ready to accept requests immediately,
    and any request that lands before the warm-up finishes will simply
    re-enter ``_ensure_descendant_cache`` and wait on the same lock.

    Also starts the activity-gated system-prompt cache pre-warmer (C1) IFF
    ``HEARTVAR_CACHE_PREWARM`` is truthy. Default OFF: with the env unset,
    ``maybe_start_prewarm`` returns None and NO background warming task is
    created, so behaviour is byte-identical to today (no extra Anthropic
    calls). When enabled, the task is cancelled+awaited cleanly on shutdown so
    no orphaned task / exit warning remains.

    Finally, arms the event-loop probe IFF ``HEARTVAR_LOOP_PROBE`` is truthy —
    the instrument for the ~19 s production stall, where fifteen sources
    completed in the same millisecond. It turns on asyncio's own slow-callback
    reporting so the loop NAMES the callback that blocked it. Default OFF, and
    it should stay off outside a diagnostic run: debug mode times every
    callback. See backend/loop_probe.py."""
    maybe_enable_loop_probe()

    async def _go() -> None:
        try:
            await _ensure_descendant_cache()
        except Exception:  # noqa: BLE001 — non-fatal background task
            log.exception("PanelApp HPO-descendant warm-up failed")
        if outbound_warming_enabled():
            try:
                await warm_outbound()
            except Exception:  # noqa: BLE001 — non-fatal background task
                log.exception("outbound warm-up failed")
        if not warming_enabled():
            return
        try:
            await warm_local_data()
        except Exception:  # noqa: BLE001 — non-fatal background task
            log.exception("local data warm-up failed")
        try:
            await warm_resolver()
        except Exception:  # noqa: BLE001 — non-fatal background task
            log.exception("resolver warm-up failed")
    asyncio.create_task(_go())
    prewarm_task = maybe_start_prewarm()
    _record_deploy_event()
    try:
        yield
    finally:
        if prewarm_task is not None:
            prewarm_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await prewarm_task


app = FastAPI(title="HeartVar backend", lifespan=lifespan)
class _RevalidateStaticFiles(StaticFiles):
    """StaticFiles that forces browsers to revalidate every asset.

    The frontend assets (heartvar.css / heartvar.js / heartvar.landing.js) are
    referenced by fixed, UNVERSIONED URLs and the index HTML is served
    ``no-store`` (always fresh). Without a ``Cache-Control`` header on the
    assets, Starlette's StaticFiles sends only ETag/Last-Modified, so browsers
    fall back to *heuristic* caching and can serve a STALE css/js after a deploy
    — new HTML referencing new markup + old CSS renders the new elements
    unstyled (e.g. the genome-build toggle appearing raw on load). ``no-cache``
    keeps the file cacheable but forces an ETag revalidation on every load: an
    unchanged file returns a cheap 304 (no body), a changed file is re-fetched
    immediately. Avoids per-deploy URL version bumping.
    """

    async def get_response(self, path: str, scope):
        response = await super().get_response(path, scope)
        response.headers["Cache-Control"] = "no-cache"
        return response


app.mount("/static", _RevalidateStaticFiles(directory=PROJECT_ROOT / "static"), name="static")
_STRUCTURE_DIR = Path(
    os.environ.get("HEARTVAR_STRUCTURE_DIR") or (PROJECT_ROOT / "data" / "alphafold")
)
try:
    _STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
except OSError as e:
    log.warning(
        "Could not create structure dir %s (%s) — 3-D viewer disabled, app continues",
        _STRUCTURE_DIR, e,
    )
app.mount(
    "/structures",
    StaticFiles(directory=_STRUCTURE_DIR, check_dir=False),
    name="structures",
)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)
app.add_middleware(_MetricsMiddleware)
app.add_middleware(
    SessionMiddleware,
    secret_key=os.environ.get("SESSION_SECRET_KEY") or secrets.token_hex(32),
    https_only=os.environ.get("SESSION_HTTPS_ONLY", "1").strip() != "0",
    same_site="lax",
)


_msal_app: msal.ConfidentialClientApplication | None = None


def _get_msal_app() -> msal.ConfidentialClientApplication:
    global _msal_app
    if _msal_app is None:
        _msal_app = msal.ConfidentialClientApplication(
            client_id=os.environ["CLIENT_ID"],
            client_credential=os.environ["CLIENT_SECRET"],
            authority=os.environ.get(
                "AUTHORITY",
                "https://login.microsoftonline.com/common",
            ),
        )
    return _msal_app


def session_account(request: Request) -> dict | None:
    """The signed-in Microsoft Entra ID account for this request, or None.

    Wraps ``request.session`` because reading it raises AssertionError when
    SessionMiddleware is not installed — which is the case for any consumer that
    builds the app without auth configured, and for unit tests of unrelated
    endpoints. Auth being absent must degrade to "anonymous", never to a 500.
    """
    try:
        account = request.session.get("account")
    except (AssertionError, AttributeError):
        return None
    return account or None


def account_key(account: dict) -> str:
    """Stable per-user identity for quota keying.

    Prefers ``oid`` (the immutable Entra object ID) over ``preferred_username``,
    which is an email address and can be reassigned. ``sub`` is per-application so
    it is a valid fallback. Anything unrecognisable keys to a shared bucket rather
    than silently getting an unlimited one.
    """
    for claim in ("oid", "sub", "preferred_username", "email"):
        value = account.get(claim)
        if isinstance(value, str) and value.strip():
            return f"u:{claim}:{value.strip()}"
    return "u:unidentified"


def account_label(account: dict) -> str:
    """What to show in the UI and write to logs."""
    for claim in ("preferred_username", "email", "name"):
        value = account.get(claim)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return "signed-in user"


def _admin_emails() -> frozenset[str]:
    """Normalised set of email addresses that have admin-portal access.

    Configured via ``HEARTVAR_ADMIN_EMAILS`` — a comma-separated list of email
    addresses exactly as they appear in the identity provider's token claims
    (``preferred_username`` or ``email``). Comparison is case-insensitive.
    """
    raw = os.environ.get("HEARTVAR_ADMIN_EMAILS", "")
    return frozenset(e.strip().lower() for e in raw.split(",") if e.strip())


def is_admin(account: dict | None) -> bool:
    """True when the signed-in account belongs to the admin allow-list."""
    if not account:
        return False
    admins = _admin_emails()
    if not admins:
        return False
    for claim in ("preferred_username", "email"):
        value = account.get(claim)
        if isinstance(value, str) and value.strip().lower() in admins:
            return True
    return False


def ai_auth_required() -> bool:
    """Whether the AI paths demand a session.

    True only once Microsoft Entra ID is actually configured, so deploying this code without
    CLIENT_ID/CLIENT_SECRET leaves AI behaviour exactly as it was rather than
    locking the feature out of reach. Set HEARTVAR_AI_REQUIRES_SIGNIN=0 to force it
    off, or =1 to insist on it (which will 401 every AI request until the app
    registration is in place — deliberate, so a misconfiguration fails loudly
    rather than silently serving AI to anonymous callers).
    """
    override = os.environ.get("HEARTVAR_AI_REQUIRES_SIGNIN", "").strip().lower()
    if override in ("0", "false", "no", "off"):
        return False
    if override in ("1", "true", "yes", "on"):
        return True
    return bool(configured_providers())


def microsoft_configured() -> bool:
    return bool(
        os.environ.get("CLIENT_ID", "").strip()
        and os.environ.get("CLIENT_SECRET", "").strip()
    )


def google_configured() -> bool:
    return bool(
        os.environ.get("GOOGLE_CLIENT_ID", "").strip()
        and os.environ.get("GOOGLE_CLIENT_SECRET", "").strip()
    )


def configured_providers() -> list[dict]:
    """Sign-in providers this deployment can actually complete.

    Drives the chooser dialog, so an unconfigured provider is never offered as a
    button that dead-ends in a KeyError on its missing settings.
    """
    providers = []
    if microsoft_configured():
        providers.append({"key": "microsoft", "label": "Microsoft"})
    if google_configured():
        providers.append({"key": "google", "label": "Google"})
    return providers


_GOOGLE_DISCOVERY = "https://accounts.google.com/.well-known/openid-configuration"
_google_discovery_cache: dict = {}


async def _google_endpoints() -> dict:
    """Google's authorize/token endpoints, memoised (they are effectively static)."""
    if _google_discovery_cache:
        return _google_discovery_cache
    async with httpx.AsyncClient(timeout=10.0) as client:
        resp = await client.get(_GOOGLE_DISCOVERY)
    resp.raise_for_status()
    doc = resp.json()
    for field in ("authorization_endpoint", "token_endpoint", "issuer"):
        if not doc.get(field):
            raise ValueError(f"Google discovery document lacks {field}")
    _google_discovery_cache.update(doc)
    return _google_discovery_cache


def _google_redirect_uri() -> str:
    """Must match the URI registered in the Google Cloud console exactly."""
    explicit = os.environ.get("GOOGLE_REDIRECT_URI", "").strip()
    if explicit:
        return explicit
    base = os.environ.get("REDIRECT_URI", "").strip()
    if base.endswith("/auth/callback"):
        return base[: -len("/auth/callback")] + "/auth/google/callback"
    raise HTTPException(
        500,
        "Google sign-in is misconfigured: set GOOGLE_REDIRECT_URI (or REDIRECT_URI "
        "ending in /auth/callback).",
    )


def _safe_next(raw: str | None) -> str:
    """Clamp a ``?next=`` sign-in return target to a path on this site.

    Both callbacks end with ``RedirectResponse(url=next_url)``, and ``next`` is
    supplied by whoever built the link — so an unchecked value is an open
    redirect: ``/login?next=https://evil.example`` sends a user who has just
    typed their VCCRI credentials onward to an attacker's page, arriving from a
    genuine heartvar.victorchang.edu.au link. Only ever ours to hand out, but
    the client now builds it from location.pathname (heartvar.auth.js), which
    real page URLs make user-influenced, so it is checked here.

    Accepted: a single leading "/" followed by anything that is not another "/"
    or "\\". That rejects absolute URLs, protocol-relative "//host" (which a
    browser reads as a host, not a path), and the "/\\host" variant browsers
    fold to the same thing. Anything else falls back to the home page.
    """
    if not raw or not raw.startswith("/"):
        return "/"
    if raw[:2] in ("//", "/\\"):
        return "/"
    return raw


@app.get("/login/google")
async def login_google(request: Request) -> RedirectResponse:
    """Begin Google sign-in (authorization code + PKCE)."""
    if not google_configured():
        raise HTTPException(404, "Google sign-in is not configured on this server.")
    doc = await _google_endpoints()
    verifier = secrets.token_urlsafe(64)
    challenge = (
        base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest())
        .rstrip(b"=")
        .decode()
    )
    state = secrets.token_urlsafe(24)
    nonce = secrets.token_urlsafe(24)
    request.session["google_flow"] = {
        "state": state, "nonce": nonce, "verifier": verifier,
    }
    request.session["auth_next"] = _safe_next(request.query_params.get("next"))
    params = urlencode({
        "client_id": os.environ["GOOGLE_CLIENT_ID"],
        "response_type": "code",
        "redirect_uri": _google_redirect_uri(),
        "scope": "openid email profile",
        "state": state,
        "nonce": nonce,
        "code_challenge": challenge,
        "code_challenge_method": "S256",
    })
    return RedirectResponse(url=f"{doc['authorization_endpoint']}?{params}")


@app.get("/auth/google/callback")
async def auth_google_callback(request: Request) -> RedirectResponse:
    """Google OAuth callback — exchange the code and store the account."""
    flow = request.session.pop("google_flow", None)
    next_url = request.session.pop("auth_next", "/")
    if not flow:
        raise HTTPException(400, "Auth session expired — please try signing in again.")
    if request.query_params.get("error"):
        log.warning("Google auth error: %s",
                    _log_safe(request.query_params.get("error", "")))
        return RedirectResponse(url="/?auth_error=denied")
    code = request.query_params.get("code", "")
    if not code or not hmac.compare_digest(
        str(flow.get("state") or ""), request.query_params.get("state", "")
    ):
        return RedirectResponse(url="/?auth_error=state")

    doc = await _google_endpoints()
    async with httpx.AsyncClient(timeout=10.0) as client:
        token_resp = await client.post(
            doc["token_endpoint"],
            data={
                "grant_type": "authorization_code",
                "code": code,
                "redirect_uri": _google_redirect_uri(),
                "client_id": os.environ["GOOGLE_CLIENT_ID"],
                "client_secret": os.environ["GOOGLE_CLIENT_SECRET"],
                "code_verifier": flow.get("verifier", ""),
            },
            headers={"Accept": "application/json"},
        )
    if token_resp.status_code != 200:
        log.warning("Google token exchange failed: %s", token_resp.status_code)
        return RedirectResponse(url="/?auth_error=token")
    id_token = (token_resp.json() or {}).get("id_token")
    if not id_token:
        return RedirectResponse(url="/?auth_error=token")

    try:
        claims = _claims_from_google_id_token(id_token, str(flow.get("nonce") or ""))
    except ValueError as e:
        log.warning("Google ID token rejected: %s", _log_safe(str(e)))
        return RedirectResponse(url="/?auth_error=token")

    request.session["account"] = {
        "sub": claims.get("sub"),
        "preferred_username": claims.get("email"),
        "name": claims.get("name"),
        "idp": "google",
    }
    log.info("Sign-in ok: provider=google user=%s",
             _log_safe(claims.get("email") or claims.get("sub") or "?"))
    return RedirectResponse(url=next_url)


def _claims_from_google_id_token(id_token: str, nonce: str) -> dict:
    """Decode and validate Google's ID token.

    The signature is not re-checked: the token arrives in the body of a direct,
    TLS-verified, server-to-server POST to Google's token endpoint and never passes
    through the browser — the case OpenID Connect Core 3.1.3.7 item 6 explicitly
    allows TLS validation to cover. iss/aud/nonce/exp ARE validated; those are what
    bind the token to this client and this attempt. If this ever accepts a token
    from the browser instead, signature verification becomes mandatory.
    """
    parts = id_token.split(".")
    if len(parts) != 3:
        raise ValueError("malformed ID token")
    padded = parts[1] + "=" * (-len(parts[1]) % 4)
    claims = json.loads(base64.urlsafe_b64decode(padded))
    if claims.get("iss") not in ("https://accounts.google.com", "accounts.google.com"):
        raise ValueError("issuer mismatch")
    if claims.get("aud") != os.environ["GOOGLE_CLIENT_ID"]:
        raise ValueError("audience is not this client")
    if not hmac.compare_digest(str(claims.get("nonce") or ""), nonce):
        raise ValueError("nonce mismatch")
    exp = claims.get("exp")
    if not isinstance(exp, (int, float)) or exp < time.time() - 60:
        raise ValueError("expired")
    if not claims.get("sub"):
        raise ValueError("no subject")
    return claims


def require_ai_account(request: Request) -> dict | None:
    """Gate an AI-spending request. Returns the account (None when the gate is off).

    Raises 401 with a MACHINE-READABLE body rather than redirecting to /login.
    A redirect is right for a page navigation but wrong here: /api/curate/stream
    and /api/chat are called with fetch(), so a 302 would be followed
    transparently and the caller would receive the login page's HTML instead of
    an error it can act on. The frontend keys off
    ``detail.error == "signin_required"`` to open the sign-in dialog.
    """
    account = session_account(request)
    if account is not None or not ai_auth_required():
        return account
    raise HTTPException(
        401,
        detail={
            "error": "signin_required",
            "message": (
                "Sign in to use the AI interpretation. The variant evidence and the "
                "rule-based classification remain available without signing in."
            ),
        },
    )


@app.get("/login")
async def login(request: Request) -> RedirectResponse:
    """Initiate Microsoft login — redirect the browser to the Microsoft Entra ID login page."""
    next_url = _safe_next(request.query_params.get("next"))
    flow = _get_msal_app().initiate_auth_code_flow(
        scopes=[],
        redirect_uri=os.environ["REDIRECT_URI"],
    )
    request.session["auth_flow"] = flow
    request.session["auth_next"] = next_url
    return RedirectResponse(url=flow["auth_uri"])


@app.get("/auth/callback")
async def auth_callback(request: Request) -> RedirectResponse:
    """Microsoft Entra ID OAuth callback — exchange the auth code for a token and store the account."""
    flow = request.session.pop("auth_flow", None)
    next_url = request.session.pop("auth_next", "/")
    if not flow:
        raise HTTPException(400, "Auth session expired — please try logging in again.")
    result = _get_msal_app().acquire_token_by_auth_code_flow(
        flow,
        dict(request.query_params),
    )
    if "error" in result:
        log.warning(
            "MSAL auth error: %s — %s",
            _log_safe(result.get("error", "")),
            _log_safe(result.get("error_description", "")),
        )
        raise HTTPException(
            401,
            f"Authentication failed: {result.get('error_description') or result.get('error')}",
        )
    request.session["account"] = result.get("id_token_claims") or {}
    return RedirectResponse(url=next_url)


@app.get("/api/auth/status")
async def auth_status(request: Request) -> dict:
    """Whether AI needs a sign-in here, and whether anyone is signed in.

    Separate from /api/auth/me, which 401s when anonymous and so cannot
    distinguish "not signed in" from "auth is not configured on this server".
    The frontend needs that distinction: with no CLIENT_ID it must render NO
    sign-in affordance at all, because /login would raise a KeyError on the
    missing app-registration settings. Also lets a deploy without Microsoft Entra ID behave
    exactly as it did before auth existed.
    """
    account = session_account(request)
    return {
        "ai_requires_signin": ai_auth_required(),
        "authenticated": account is not None,
        "label": account_label(account) if account else "",
        "providers": configured_providers(),
        "is_admin": is_admin(account),
    }


@app.get("/api/auth/me")
async def auth_me(request: Request):
    """Return the signed-in account info, or 401 if not authenticated."""
    account = request.session.get("account")
    if not account:
        raise HTTPException(401, "Not authenticated")
    return {"name": account.get("name"), "email": account.get("preferred_username")}


@app.get("/logout")
async def logout(request: Request) -> RedirectResponse:
    """Clear the session and redirect to Microsoft logout."""
    idp = (session_account(request) or {}).get("idp")
    request.session.clear()
    _configured = os.environ.get("POST_LOGOUT_REDIRECT_URI", "").strip()
    post_logout = _configured if _configured else str(request.base_url)
    if idp == "google":
        return RedirectResponse(url=post_logout)
    authority = os.environ.get("AUTHORITY", "https://login.microsoftonline.com/common")
    return RedirectResponse(
        url=f"{authority}/oauth2/v2.0/logout?post_logout_redirect_uri={post_logout}"
    )


@app.get("/")
@app.get("/about")
@app.get("/contact")
async def serve_frontend() -> FileResponse:
    """The single-page shell, served at every addressable page path.

    /about and /contact are DECLARED EXPLICITLY rather than caught by a
    ``/{path:path}`` route. A catch-all here would be matched after the
    hand-written routes above but before nothing at all: every mistyped
    ``/api/...`` call, and any future endpoint added below this line, would
    answer 200 with an HTML body instead of a 404 — a failure mode that hides
    itself, because the browser renders the app and the missing data just looks
    like an empty result. Two literal paths cost one line each and cannot
    shadow anything. See backend/tests/test_page_routes.py.

    Section links stay client-side fragments (/about#privacy), so only the page
    paths themselves need a route. The client reads location.pathname on boot
    and shows the matching page — see showPage() in static/heartvar.js.
    """
    target = PROJECT_ROOT / "index.html"
    if not target.exists():
        raise HTTPException(404, "index.html not found")
    return FileResponse(
        target,
        headers={"Cache-Control": "no-store, no-cache, must-revalidate"},
    )


_NOT_FOUND_HTML = """<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Page not found — HeartVar</title>
<style>
  body { margin:0; min-height:100vh; display:flex; align-items:center;
         justify-content:center; background:#faf8f7;
         font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         color:#2b2b2b; }
  .box { text-align:center; padding:2rem; max-width:32rem; }
  .code { font-size:3rem; font-weight:600; color:#C8341A; margin:0 0 .5rem; }
  h1 { font-size:1.25rem; margin:0 0 .75rem; font-weight:600; }
  p { margin:0 0 1.5rem; line-height:1.5; color:#555; }
  a.btn { display:inline-block; background:#C8341A; color:#fff; text-decoration:none;
          padding:.6rem 1.2rem; border-radius:6px; font-weight:500; }
  a.btn:hover { background:#a82b16; }
  .alt { margin-top:1.25rem; font-size:.9rem; }
  .alt a { color:#C8341A; }
</style>
</head><body>
  <div class="box">
    <p class="code">404</p>
    <h1>We couldn't find that page</h1>
    <p>The address may be mistyped, or the page may have moved.</p>
    <a class="btn" href="/">Go to HeartVar</a>
    <p class="alt"><a href="/about">About</a> &nbsp;·&nbsp; <a href="/contact">Contact</a></p>
  </div>
</body></html>"""


@app.exception_handler(StarletteHTTPException)
async def _not_found_page(request: Request, exc: StarletteHTTPException):
    """Serve the styled 404 to browsers; delegate everything else unchanged.

    Registering a handler for StarletteHTTPException intercepts EVERY
    HTTPException, so the non-404 and non-browser paths must fall through to
    FastAPI's own handler rather than be reimplemented here — otherwise this
    would quietly change the body of every 4xx/5xx the app raises.
    """
    if exc.status_code == 404 and not request.url.path.startswith("/api/"):
        if "text/html" in (request.headers.get("accept") or ""):
            return HTMLResponse(_NOT_FOUND_HTML, status_code=404)
    return await http_exception_handler(request, exc)


_MIRROR_FRESHNESS_FILES = (
    "clinvar.db", "uniprot.db", "medgen.db", "mgi.db", "biogrid.db",
    "opentargets.db", "panelapp_aus_snapshot.json", "hgnc_alias_map.json",
    "erepo_all.tsv", "gencc_submissions.json", "clingen_gene_validity.json",
)


def _mirror_mtime_iso(data_dir: Path) -> str | None:
    """Newest write time across the monthly data artifacts, as a UTC ISO string.

    A fallback for the build stamp, not a replacement: a file's mtime says when it
    was last WRITTEN, which for these builders is when the source was last
    refreshed. Missing files are skipped, so a partially-provisioned mount still
    yields an answer.
    """
    newest = 0.0
    for name in _MIRROR_FRESHNESS_FILES:
        try:
            newest = max(newest, (data_dir / name).stat().st_mtime)
        except OSError:
            continue
    if newest <= 0:
        return None
    return (
        datetime.fromtimestamp(newest, timezone.utc)
        .replace(microsecond=0)
        .isoformat()
    )


_LOG_BASE_DIR = Path(
    os.environ.get("HEARTVAR_LOGS_DIR")
    or (PROJECT_ROOT / "data" / "logs")
)
_LOG_MAX_FILES = 5
_LOG_MAX_READ_BYTES = 512 * 1024
_ALLOWED_LOG_CATEGORIES: frozenset[str] = frozenset(["db_builder", "deployments"])
_LOG_CATEGORY_DIRS: dict[str, str] = {
    "db_builder": "db_builder",
    "deployments": "deployments",
}

def _log_dir(category: str) -> Path:
    return _LOG_BASE_DIR / category


def _prune_log_files(directory: Path) -> None:
    """Delete the oldest .log files, keeping only the last _LOG_MAX_FILES."""
    try:
        files = sorted(directory.glob("*.log"), key=lambda p: p.stat().st_mtime)
        for old in files[:-_LOG_MAX_FILES]:
            try:
                old.unlink()
            except OSError:
                pass
    except OSError as exc:
        log.warning("Could not prune log files in %s: %s", directory, exc)


def _record_deploy_event() -> None:
    """Create a per-deployment log file on startup and add a FileHandler so
    subsequent heartvar log records are also written to it.

    Keeps at most _LOG_MAX_FILES deployment log files; older ones are deleted.
    Failures are logged and swallowed — a broken log directory must never
    prevent the server from starting.
    """
    now = datetime.now(timezone.utc).replace(microsecond=0)
    ts_file = now.strftime("%Y-%m-%dT%H-%M-%SZ")
    log_dir = _log_dir("deployments")
    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"{ts_file}.log"
        header_lines = [
            "=== HeartVar server started ===",
            f"Time (UTC):          {now.isoformat()}",
            f"AI available:        {SERVER_AI_AVAILABLE}",
            f"Auth providers:      {', '.join(p['key'] for p in configured_providers()) or 'none'}",
            f"Admin configured:    {bool(_admin_emails())}",
            f"Offline VEP:         {vep_offline_status()}",
            "================================",
            "",
        ]
        log_path.write_text("\n".join(header_lines), encoding="utf-8")
        handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        )
        handler.setLevel(logging.INFO)
        logging.getLogger("heartvar").addHandler(handler)
        _install_log_redaction()
        _prune_log_files(log_dir)
        log.info("Deploy log: %s", log_path.name)
    except OSError as exc:
        log.warning("Could not create deploy log file: %s", exc)


def _extract_build_error(label: str, log_dir: Path) -> str | None:
    """Return the build output for *label* from the most recent db_builder log.

    Scans for the ``>>> BUILD label`` / ``>>> REFRESH label`` start marker and
    the ``!!! FAIL  label`` end marker and returns those lines (capped at 40 so
    a huge traceback doesn't swamp the UI).  Returns None if the log is absent
    or if no FAIL record for that label is found.
    """
    try:
        logs = sorted(log_dir.glob("*.log"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not logs:
            return None
        content = logs[0].read_text(encoding="utf-8", errors="replace")
        lines = content.splitlines()
        start_idx: int | None = None
        end_idx: int | None = None
        for i, line in enumerate(lines):
            if f">>> BUILD {label}" in line or f">>> REFRESH  {label}" in line:
                start_idx = i
            if start_idx is not None and f"!!! FAIL  {label}" in line:
                end_idx = i
                break
        if start_idx is None or end_idx is None:
            return None
        snippet = lines[start_idx : end_idx + 1]
        if len(snippet) > 40:
            snippet = ["... (truncated — showing last 39 lines) ..."] + snippet[-39:]
        return "\n".join(snippet)
    except OSError:
        return None


def _require_admin(request: Request) -> None:
    """Raise 403 when the caller is not an admin.  Call before touching admin data."""
    account = session_account(request)
    if not is_admin(account):
        raise HTTPException(403, detail={"error": "forbidden", "message": "Admin access required."})


def _safe_log_name(name: str) -> bool:
    """True when name contains only safe filename characters (no path traversal)."""
    return bool(re.match(r'^[\w\-+.]+\.log$', name)) and ".." not in name


@app.get("/api/admin/log-index")
async def admin_log_index(request: Request) -> dict:
    """List the last {_LOG_MAX_FILES} log files for each admin category.

    Returns ``{"db_builder": [...], "deployments": [...]}`` where each entry is
    ``{"name": "...", "ts": "...", "size": <bytes>}``.
    Accessible only to signed-in admin accounts.
    """
    _require_admin(request)
    result: dict[str, list[dict]] = {}
    for category in sorted(_ALLOWED_LOG_CATEGORIES):
        d = _log_dir(category)
        entries: list[dict] = []
        if d.exists():
            files = sorted(d.glob("*.log"), key=lambda p: p.stat().st_mtime, reverse=True)
            for p in files[:_LOG_MAX_FILES]:
                try:
                    stat = p.stat()
                    entries.append({
                        "name": p.name,
                        "ts": datetime.fromtimestamp(
                            stat.st_mtime, timezone.utc
                        ).replace(microsecond=0).isoformat(),
                        "size": stat.st_size,
                    })
                except OSError:
                    continue
        result[category] = entries
    return result


@app.get("/api/admin/log-content")
async def admin_log_content(
    request: Request, category: str, name: str
) -> dict:
    """Return the tail of a single log file.

    ``category`` must be one of the allowed values; ``name`` must match a file
    that actually exists in that category's directory.

    The path is resolved by ENUMERATION, not concatenation. Neither argument is
    ever joined onto a path: the category selects a constant directory name from
    a literal map, and the filename is compared for equality against the entries
    of that directory. So the only paths that can be opened are ones the server
    itself listed, traversal has no string to travel through, and the guarantee
    does not rest on a regex being exhaustive. (It also satisfies CodeQL's
    py/path-injection, which cannot see through a validate-then-join.)
    """
    _require_admin(request)
    directory_name = _LOG_CATEGORY_DIRS.get(category)
    if directory_name is None:
        raise HTTPException(400, "Invalid log category.")
    if not _safe_log_name(name):
        raise HTTPException(400, "Invalid log file name.")
    target = None
    try:
        for candidate in (_LOG_BASE_DIR / directory_name).iterdir():
            if (
                candidate.name == name
                and not candidate.is_symlink()
                and candidate.is_file()
            ):
                target = candidate
                break
    except OSError:
        target = None
    if target is None:
        raise HTTPException(404, "Log file not found.")
    log_path = target
    try:
        size = log_path.stat().st_size
        with log_path.open("rb") as fh:
            if size > _LOG_MAX_READ_BYTES:
                fh.seek(size - _LOG_MAX_READ_BYTES)
            raw = fh.read()
        content = raw.decode("utf-8", errors="replace")
        truncated = size > _LOG_MAX_READ_BYTES
        if truncated:
            content = content.split("\n", 1)[-1]
            content = (
                f"... (truncated \u2014 showing the last "
                f"{_LOG_MAX_READ_BYTES // 1024} KB of {size // 1024} KB) ...\n"
                + content
            )
    except OSError:
        log.exception(
            "Could not read admin log file %s/%s",
            _log_safe(category), _log_safe(name),
        )
        raise HTTPException(500, "Could not read log file.")
    return {
        "content": content, "name": name, "category": category,
        "truncated": truncated, "size": size,
    }


_DB_SOURCE_ARTIFACTS: dict[str, str] = {
    "uniprot":           "data/uniprot.db",
    "clinvar":           "data/clinvar.db",
    "gnomad_constraint": "data/gnomad_constraint.db",
    "gtex":              "data/gtex.db",
    "opentargets":       "data/opentargets.db",
    "medgen":            "data/medgen.db",
    "mgi":               "data/mgi.db",
    "biogrid":           "data/biogrid.db",
    "hpo_labels":        "backend/data/hpo_labels.json",
    "panelapp":          "data/panelapp_aus_snapshot.json",
    "gnomad_freq":       "data/gnomad_freq.db",
    "alphafold":         "data/alphafold/manifest.json",
    "spliceai":          "data/spliceai_cardiac.masked.grch38.vcf.gz",
    "fetal_heart":       "data/fetal_heart.db",
    "hgnc_alias":        "backend/data/hgnc_alias_map.json",
    "gene_id_map":       "backend/data/gene_id_map.json.gz",
    "erepo":             "backend/data/erepo_all.tsv",
    "gencc":             "backend/data/gencc_submissions.json",
    "clingen_gv":        "backend/data/clingen_gene_validity.json",
    "vep":               "data/vep/.heartvar_vep_manifest.json",
}

_DB_SOURCE_CADENCE: dict[str, str] = {
    "uniprot": "monthly", "clinvar": "monthly", "medgen": "monthly",
    "mgi": "monthly", "panelapp": "monthly", "hpo_labels": "monthly",
    "hgnc_alias": "monthly", "erepo": "monthly", "gencc": "monthly",
    "clingen_gv": "monthly", "opentargets": "monthly", "biogrid": "monthly",
    "gnomad_constraint": "versioned", "gnomad_freq": "versioned", "gtex": "versioned",
    "vep": "versioned",
    "alphafold": "static", "spliceai": "static", "fetal_heart": "static",
    "gene_id_map": "image",
}

_DB_SOURCE_ENV_OVERRIDES: dict[str, str] = {
    "uniprot":           "UNIPROT_DB_PATH",
    "clinvar":           "CLINVAR_DB_PATH",
    "gnomad_constraint": "GNOMAD_CONSTRAINT_DB_PATH",
    "gnomad_freq":       "GNOMAD_FREQ_DB_PATH",
    "gtex":              "GTEX_DB_PATH",
    "opentargets":       "OPENTARGETS_DB_PATH",
    "medgen":            "MEDGEN_DB_PATH",
    "mgi":               "MGI_DB_PATH",
    "biogrid":           "BIOGRID_DB_PATH",
    "panelapp":          "PANELAPP_SNAPSHOT_PATH",
    "spliceai":          "SPLICEAI_DB_PATH",
    "fetal_heart":       "FETAL_HEART_DB_PATH",
    "hgnc_alias":        "HGNC_ALIAS_MAP_PATH",
    "gene_id_map":       "GENE_ID_MAP_PATH",
    "erepo":             "HEARTVAR_EREPO_TSV",
    "gencc":             "GENCC_SNAPSHOT_PATH",
    "clingen_gv":        "CLINGEN_GV_PATH",
}

_DB_SOURCE_MOUNT_COPIES: frozenset[str] = frozenset(
    {"hgnc_alias", "erepo", "gencc", "clingen_gv"}
)


def _db_artifact_candidates(name: str, rel: str) -> list[Path]:
    """Every location a source's artifact can legitimately occupy.

    A set override is the ONLY candidate, because that is how the readers resolve
    it — ``os.environ.get(VAR) or <default>``, never both. Falling back to the
    declared path when an override names a missing file would report "present"
    off a copy the reader will never open, which is this bug in reverse.

    With no override there are two honest locations: the declared path, and (for
    the four artifacts build_all.sh copies onto the mount) the mount copy."""
    env_var = _DB_SOURCE_ENV_OVERRIDES.get(name, "")
    if env_var:
        override = os.environ.get(env_var, "").strip()
        if override:
            return [Path(override)]
    candidates = [PROJECT_ROOT / rel]
    if name in _DB_SOURCE_MOUNT_COPIES:
        candidates.append(PROJECT_ROOT / "data" / Path(rel).name)
    return candidates


def _db_artifact_display(path: Path) -> str:
    """Repo-relative while the artifact sits inside the deployment tree, absolute
    otherwise, so an override pointing outside it stays legible."""
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


@app.get("/api/admin/db-status")
async def admin_db_status(request: Request) -> dict:
    """Per-source data status: artifact presence, last-updated mtime, and last
    build-run timestamp from build_stamp.json.

    Returns ``{"last_run": <iso|null>, "sources": [{name, last_updated, present}, ...]}``.
    Accessible only to signed-in admin accounts.
    """
    _require_admin(request)
    stamp_path = Path(
        os.environ.get("HEARTVAR_BUILD_STAMP_PATH")
        or (PROJECT_ROOT / "data" / "build_stamp.json")
    )
    stamp: dict = {}
    try:
        loaded = json.loads(stamp_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            stamp = loaded
    except (OSError, ValueError):
        pass

    last_run: str | None = stamp.get("last_run_utc") or None
    failed_set: set[str] = set(stamp.get("failed") or [])
    db_log_dir = _log_dir("db_builder")

    sources: list[dict] = []
    for name, rel in _DB_SOURCE_ARTIFACTS.items():
        last_updated: str | None = None
        present = False
        artifact_path = rel
        for candidate in _db_artifact_candidates(name, rel):
            try:
                mtime = candidate.stat().st_mtime
            except OSError:
                continue
            present = True
            artifact_path = _db_artifact_display(candidate)
            last_updated = (
                datetime.fromtimestamp(mtime, timezone.utc)
                .replace(microsecond=0)
                .isoformat()
            )
            break

        failed_in_last_run = name in failed_set
        error_snippet: str | None = None
        if failed_in_last_run:
            error_snippet = _extract_build_error(name, db_log_dir)

        sources.append({
            "name": name,
            "last_updated": last_updated,
            "present": present,
            "artifact_path": artifact_path,
            "cadence": _DB_SOURCE_CADENCE.get(name, "monthly"),
            "failed_in_last_run": failed_in_last_run,
            "error_snippet": error_snippet,
        })

    return {"last_run": last_run, "sources": sources}


@app.get("/api/data-status")
async def data_status() -> dict:
    """When the local data mirror was last refreshed.

    Preferred source is ``data/build_stamp.json``, which scripts/build_all.sh
    writes onto the same mount the app reads. Exists so the About page can state a
    date that is true by construction: the previous hand-maintained note claimed
    "Last update: 21/07/2026" while the monthly refresh had in fact been a no-op
    since first build, so the one number users could see was the one number nobody
    was updating.

    ``last_refresh_utc`` from the stamp advances only when a builder actually
    (re)built something — a crash-retry run that skips everything must not make the
    mirror look fresher than it is.

    FALLBACK: with no stamp (a mount last built by a version that did not write
    one, or a build still in progress — the stamp is written at the end), fall back
    to the newest mtime of the monthly data files. That is a real, honest date
    rather than a blank, and it upgrades itself to the exact stamp on the next
    completed build. ``source`` records which was used. Returns ``{}`` only when
    neither is available, which the frontend renders as no date at all.
    """
    stamp_path = Path(
        os.environ.get("HEARTVAR_BUILD_STAMP_PATH")
        or (PROJECT_ROOT / "data" / "build_stamp.json")
    )
    data_dir = Path(
        os.environ.get("HEARTVAR_DATA_DIR") or stamp_path.parent
    )

    stamp: dict = {}
    try:
        loaded = json.loads(stamp_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            stamp = loaded
    except (OSError, ValueError):
        pass

    last_refresh = stamp.get("last_refresh_utc")
    if isinstance(last_refresh, str) and last_refresh:
        return {
            "last_refresh_utc": last_refresh,
            "last_run_utc": stamp.get("last_run_utc"),
            "source": "build_stamp",
            "refreshed_count": len(stamp.get("refreshed") or []),
            "built_count": len(stamp.get("built") or []),
            "failed_count": len(stamp.get("failed") or []),
        }

    mtime_iso = _mirror_mtime_iso(data_dir)
    if mtime_iso:
        return {
            "last_refresh_utc": mtime_iso,
            "last_run_utc": stamp.get("last_run_utc"),
            "source": "file_mtime",
        }
    return {}


@app.get("/health")
@app.get("/healthz")
async def healthz() -> dict:
    return {"ok": True}


@app.get("/api/clinvar")
async def clinvar_lookup(gene: str, hgvs: str = "") -> dict:
    """ClinVar lookup against the local SQLite cache.

    Rebuild the database monthly by running ``python3 build_clinvar_db.py``
    at the project root — it pulls ClinVar's latest variant_summary.txt.gz
    from NCBI's FTP and replaces ``data/clinvar.db`` in place. The same
    database powers the in-pipeline ClinVar evidence step (see
    ``backend/clients/clinvar.py``), so no separate network call is made
    for either path.
    """
    gene = (gene or "").strip()
    if not gene:
        raise HTTPException(400, "gene query parameter is required")
    result = await fetch_clinvar(gene, hgvs)
    if not result.get("ok"):
        raise HTTPException(503, result.get("error") or "ClinVar lookup failed")
    records = result.get("records") or []
    if not records:
        return {
            "status": "not_found",
            "gene": gene,
            "hgvs": hgvs,
            "matches": 0,
        }
    top = records[0]
    return {
        "status": "ok",
        "gene": gene,
        "hgvs": hgvs,
        "matches": len(records),
        "clinical_significance": top.get("clinical_significance"),
        "review_status": top.get("review_status"),
        "stars": top.get("stars"),
        "submission_count": top.get("number_submitters"),
        "phenotypes": top.get("conditions") or [],
        "variation_id": top.get("variation_id"),
        "accession": top.get("accession"),
        "last_evaluated": top.get("last_evaluated"),
    }


@app.get("/api/biogrid")
async def biogrid_lookup(gene: str) -> dict:
    """BioGRID physical-interaction lookup against the local SQLite cache.

    Rebuild the database quarterly by running ``python3 build_biogrid_db.py``
    at the project root — it pulls BioGRID's latest tab3 archive for
    Homo sapiens (freely available, no API key) and replaces
    ``data/biogrid.db`` in place. The same database powers the in-pipeline
    BioGRID evidence step (see ``backend/clients/biogrid.py``), so no
    separate network call is made for either path.
    """
    gene = (gene or "").strip()
    if not gene:
        raise HTTPException(400, "gene query parameter is required")
    result = await fetch_biogrid(gene)
    if not result.get("ok"):
        raise HTTPException(503, result.get("error") or "BioGRID lookup failed")
    return {
        "status": "ok",
        "gene": gene,
        "total_interactions": result.get("total_interactions", 0),
        "low_throughput_interactions": result.get("low_throughput_interactions", 0),
        "unique_partners": result.get("unique_partners", 0),
        "partners": result.get("partners", []),
        "url": result.get("url"),
    }


@app.get("/api/fetal_heart")
async def fetal_heart_lookup(gene: str) -> dict:
    """Fetal cardiac expression lookup against the local SQLite cache.

    Returns per-(cell_type, stage) pseudobulk mean log1p(CPM) and the
    fraction of cells in each group expressing the gene, derived from
    Farah et al. 2024 ("Heart of Cells", Nature 627:854) via the UCSC
    Cell Browser.

    Rebuild the database monthly by running ``python3 build_fetal_heart_db.py``
    at the project root — it pulls the public count matrix + metadata
    and replaces ``data/fetal_heart.db`` in place. The same database
    powers the in-pipeline fetal-heart evidence step (see
    ``backend/clients/fetal_heart.py``), so no separate network call
    is made for either path.
    """
    gene = (gene or "").strip()
    if not gene:
        raise HTTPException(400, "gene query parameter is required")
    result = await fetch_fetal_heart(gene)
    if not result.get("ok"):
        raise HTTPException(503, result.get("error") or "Fetal heart lookup failed")
    if not result.get("found"):
        return {
            "status": "not_found",
            "gene": gene,
            "dataset": result.get("dataset"),
            "url": result.get("url"),
        }
    return {
        "status": "ok",
        "gene": gene,
        "dataset": result.get("dataset"),
        "stages": result.get("stages") or [],
        "cell_types": result.get("cell_types") or [],
        "band_thresholds": result.get("band_thresholds"),
        "url": result.get("url"),
    }


@app.post("/api/acmg/rescore")
async def acmg_rescore(req: RescoreRequest) -> dict:
    """Re-score a criteria set the curator has edited in the Criteria tab.

    Backs the per-criterion on/off switch + strength picker. The adjusted tier is
    shown BESIDE HeartVar's own call, never in place of it, so this endpoint's
    job is only to answer "what would the engine say about this set".

    Server-side deliberately. The browser has no scorer and must not grow one
    (tests/test_acmg_constants_parity.test_js_has_no_second_classifier): a second
    implementation cannot be kept honest, and the last one had already drifted
    into returning a tier the backend no longer returned.

    It runs the SAME chain as the automatic call — exclusions, then points, then
    the tier with the benign combining floor — so an override cannot reach a
    classification the engine itself refuses to reach. If a curator turns on both
    PP1 and BS4, they still cancel; BA1 still forces Benign.

    One deliberate difference: ``verify_source_evidence=False``. The gates that
    ask "did the source actually cite an assay / a PMID / a proband count?" read
    fields the client does not post, so on this path they would read every
    criterion as uncited and demote it — an edit to PP1 was pulling an untouched
    PS4 down from Moderate to Supporting and switching PS3 off. Those criteria
    already passed the gates when the curate call emitted them; re-running them
    blind is not a second safety net, it is data loss.

    Cheap and stateless: no evidence gather, no AI call, no DB, no sign-in. The
    sign-in point matters — the judgement criteria that most need a manual
    override are exactly the ones marked "not assessed" in the evidence-only
    flow, which is the flow available without signing in.
    """
    criteria = [c.model_dump() for c in req.criteria]
    criteria, forced = apply_cross_criterion_exclusions(
        criteria, req.gene or None, req.inheritance_input or None,
        has_curator_segregation=False,
        verify_source_evidence=False,
    )
    points_total = compute_points_total(criteria)
    return {
        "points_total": points_total,
        "classification": classification_for_criteria(
            points_total, criteria, forced,
            validity=gene_validity_ceiling(req.gene),
        ),
        "criteria": criteria,
    }


@app.get("/api/gene/validate")
async def gene_validate(symbol: str) -> dict:
    """Validate a gene symbol against the local HGNC alias map (warn-early).

    Backs the landing form's on-blur, non-blocking gene-symbol check so the
    curator learns BEFORE the full curation runs whether they typed the
    HGNC-approved symbol, a previous/alias symbol (with the approved symbol to
    use), or something unrecognised. Mirrors the in-pipeline canonicalisation
    in ``_gather_evidence_sse`` (same ``canonicalise_gene_symbol`` source), so
    the upfront note never contradicts the post-results badge.

    Returns ``{symbol, status, approved, message}`` where ``status`` is one of:
      * ``approved``     — input is the current HGNC-approved symbol.
      * ``alias``        — input is a previous/alias symbol; ``approved`` holds
                           the symbol evidence will be gathered for.
      * ``unrecognized`` — input is not in the local HGNC map.
      * ``empty``        — no symbol supplied.

    Never raises on a bad symbol (fail-open mirrors the pipeline): an
    unrecognised symbol is a soft signal, not an error, so the curator can
    still proceed (the local map may be stale).
    """
    sym = (symbol or "").strip()
    if not sym:
        return {"symbol": sym, "status": "empty", "approved": None, "message": ""}
    canon = canonicalise_gene_symbol(sym)
    approved = canon.get("approved")
    if canon.get("is_alias") and approved:
        return {
            "symbol": sym,
            "status": "alias",
            "approved": approved,
            "message": (
                f"“{sym}” is a previous/alias symbol — the HGNC-approved symbol "
                f"is {approved}. Evidence will be gathered for {approved}."
            ),
        }
    if not canon.get("recognized"):
        return {
            "symbol": sym,
            "status": "unrecognized",
            "approved": None,
            "message": (
                f"“{sym}” is not a recognised HGNC-approved gene symbol. Check the "
                f"symbol (verify at genenames.org). You can still proceed — "
                f"gene-level evidence may be incomplete."
            ),
        }
    return {
        "symbol": sym,
        "status": "approved",
        "approved": approved or sym,
        "message": "",
    }


@app.post("/api/chat")
@limiter.limit(_CHAT_RATE_LIMIT)
async def chat(request: Request, req: ChatRequest) -> dict:
    """Server-side "Ask about this variant" chatbot.

    Runs a capped chat-model call on the OWNER's Anthropic key, scoped by a
    server-owned system prompt that refuses anything off-topic (see
    ``claude.chat_reply``). This is the highest abuse-surface endpoint —
    open-ended user text on the owner's key — so it carries its own tighter
    per-IP rate limit, is counted against the daily spend budget, and fails
    closed with a clear message (never silently spends) when AI is unavailable.
    """
    if not SERVER_AI_AVAILABLE:
        raise HTTPException(503, "The assistant is not configured on this server.")
    chat_account = require_ai_account(request)
    if chat_account is not None:
        check_ai_user_quota(account_key(chat_account), "chat")
    _check_daily_cap()
    _ctx_cap, _out_cap = report_mode_caps(req.mode)
    chat_cost_est = budget.estimate_call_cost(
        CHAT_MODEL,
        max(_CHAT_INPUT_EST_TOKENS, _ctx_cap // 4),
        _out_cap,
    )
    if not budget.reserve(chat_cost_est):
        raise HTTPException(
            429,
            "The daily AI limit has been reached — the assistant is unavailable "
            "until tomorrow. The variant's evidence and classification remain "
            "available above.",
        )
    history = [{"role": t.role, "text": t.text} for t in req.history]
    try:
        answer = await chat_reply(req.context, req.question, history, req.mode)
    except Exception as e:  # noqa: BLE001 — log detail server-side, hide from client
        log.exception("chat_reply failed for question %r", _log_safe(req.question[:80]))
        raise HTTPException(502, "The assistant is temporarily unavailable.") from e
    finally:
        budget.release(chat_cost_est)
    return {"answer": answer}


from .evidence import (
    _sse,
    _protein_position_from_vep,
    _gather_evidence_sse,
    build_same_site_evidence,
)


_PYTHON_AUTHORITATIVE_SUPPLEMENTARY = ("PM1", "PP1", "BS4")


def _python_authoritative_supplementary(
    evidence: dict,
    clinical_context: dict,
    gene: str | None,
    hard_coded_criteria: list[dict],
) -> list[dict]:
    """The PM1/PP1/BS4 entries the AI path merges in alongside the hard-coded
    set, ALWAYS all three.

    infer_supplementary_criteria emits an entry only when its derivation fires,
    which left the list at 18-21 depending on the inputs. The system prompt
    tells the model that 21 verdicts are supplied to it and forbids it from
    emitting these three, so a short list told the model a verdict existed and
    then did not show it. Every gap is filled with the same not_assessed
    placeholder the no-AI arm uses, so both arms describe an absence with the
    same words and neither claims the rule was evaluated and found absent.

    Empty when HEARTVAR_NO_AI_SUPPLEMENTARY=0, which now disarms the three on
    the AI arm too and not just the no-key one."""
    if not NO_AI_SUPPLEMENTARY:
        return []
    derived = {
        c.get("code"): c
        for c in infer_supplementary_criteria(
            evidence, clinical_context, gene, hard_coded_criteria,
        )
        if c.get("code") in _PYTHON_AUTHORITATIVE_SUPPLEMENTARY
    }
    return [
        derived.get(code) or server_owned_placeholder(code)
        for code in _PYTHON_AUTHORITATIVE_SUPPLEMENTARY
    ]


async def _interpret_sse(
    req: CurationRequest,
    state: dict,
    hard_coded_criteria: list[dict],
    t_request_start: float,
) -> AsyncIterator[str]:
    evidence = state["evidence"]
    variant_id = state.get("variant_id")
    gene = state.get("gene") or ""
    clean_hgvs_c = state.get("hgvs_c") or ""
    clinical_context = evidence.get("clinical_context") or {}
    timing = state.setdefault("timing", {})

    precomputed_criteria = _python_authoritative_supplementary(
        evidence, clinical_context, gene, hard_coded_criteria,
    ) + hard_coded_criteria

    t_prompt = perf_counter()
    user_prompt = build_prompt(
        gene, clean_hgvs_c, req.hpo, req.inheritance, req.family, evidence,
        clinical_context=clinical_context,
        hard_coded_criteria=precomputed_criteria,
        segregation_context=req.segregation_context,
    )
    timing["prompt_build"] = perf_counter() - t_prompt

    prompt_chars = len(user_prompt)
    timing["prompt_chars"] = prompt_chars
    log.info(
        "[timing] %s %s — prompt built: %d chars (~%d tokens est.)",
        _log_safe(gene), _log_safe(clean_hgvs_c), prompt_chars, prompt_chars // 4,
    )
    yield _sse("stage2_start", {})

    t_claude_start = perf_counter()
    t_first_token: float | None = None
    chunks: list[str] = []
    try:
        async for chunk in stream_claude(
            user_prompt, max_tokens=RESULT_MAX_TOKENS,
        ):
            if chunk is STREAM_RESTART:
                chunks.clear()
                t_first_token = None
                yield _sse("stage2_restart", {})
                continue
            if t_first_token is None:
                t_first_token = perf_counter()
                timing["claude_ttft"] = t_first_token - t_claude_start
                log.info(
                    "[timing] %s %s — Claude time-to-first-token %.2fs",
                    _log_safe(gene), _log_safe(clean_hgvs_c), t_first_token - t_claude_start,
                )
            chunks.append(chunk)
            yield _sse("stage2_chunk", {"chunk": chunk})
    except Exception:
        log.exception("Claude stream failed for %s %s", _log_safe(gene), _log_safe(clean_hgvs_c))
        yield _sse("error", {
            "message": "The AI interpretation could not be completed. Please try again.",
        })
        return

    t_claude_done = perf_counter()
    claude_gen = t_claude_done - (t_first_token or t_claude_start)
    claude_total = t_claude_done - t_claude_start
    timing["claude_gen"] = claude_gen
    timing["claude_total"] = claude_total

    t_parse = perf_counter()
    full_text = "".join(chunks)
    try:
        result = parse_first_json_object(full_text)
    except json.JSONDecodeError as e:
        log.error(
            "Claude returned non-JSON for %s %s: %s",
            _log_safe(gene), _log_safe(clean_hgvs_c), e,
        )
        yield _sse("error", {
            "message": "The AI interpretation could not be parsed. Please try again.",
        })
        return
    timing["json_parse"] = perf_counter() - t_parse

    ai_criteria = _normalize_criteria(
        result.get("criteria"), default_source="ai",
    )
    ai_criteria = [
        c for c in ai_criteria
        if c.get("code") not in HARD_CODED_CRITERIA_CODES
    ]
    criteria = merge_hard_coded_and_ai(precomputed_criteria, ai_criteria)
    criteria = _gate_criteria_applicability(criteria, gene, evidence, req.hpo)
    criteria = apply_pm1_collision_rules(criteria, gene)
    criteria = cap_all_to_spec_strength(criteria, gene)
    criteria, forced_classification = apply_cross_criterion_exclusions(
        criteria, gene, clinical_context.get("inheritance_input"),
        has_curator_segregation=bool(
            clinical_context.get("seg_affected_carriers")
        ),
        zygosity=clinical_context.get("zygosity"),
        in_trans_pathogenic=clinical_context.get("in_trans_pathogenic"),
    )
    points_total = compute_points_total(criteria)
    classification = classification_for_criteria(
        points_total, criteria, forced_classification,
        validity=gene_validity_ceiling(gene),
    )

    evidence["carrier_status"] = carrier_status(
        evidence, clinical_context, gene, classification, req.hpo,
    )

    total_elapsed = perf_counter() - t_request_start
    timing["total"] = total_elapsed
    log.info(
        "[timing] %s %s — Claude gen %.2fs (claude_total %.2fs) | TOTAL %.2fs | criteria=%d | points=%d | class=%s",
        _log_safe(gene), _log_safe(clean_hgvs_c), claude_gen, claude_total, total_elapsed, len(criteria),
        points_total, classification,
    )
    log.info("TIMING_REPORT: %s", json.dumps({
        "gene": _log_safe(gene), "hgvs_c": _log_safe(clean_hgvs_c), **timing,
    }, default=float))

    pts_no_pp5, tier_no_pp5 = tier_without_clinvar_assertion(
        criteria, forced_classification,
        validity=gene_validity_ceiling(gene),
    )
    yield _sse("stage2_complete", {
        "classification": classification,
        "confidence": result.get("confidence"),
        "summary": result.get("summary"),
        "criteria": criteria,
        "points_total": points_total,
        "points_total_without_pp5": pts_no_pp5,
        "classification_without_pp5": tier_no_pp5,
        "borderline_reasoning": result.get("borderline_reasoning"),
        "vus_subclassification": result.get("vus_subclassification"),
        "gene_context": result.get("gene_context") or {},
        "db_evidence": evidence,
        "variant_id": variant_id,
    })

    erepo = state.get("erepo") or {}
    if erepo.get("ok") and erepo.get("found"):
        yield _sse("erepo_verdict", {
            "source": "erepo",
            "found": True,
            "classification": erepo.get("classification"),
            "criteria": erepo.get("criteria") or [],
            "vcep": erepo.get("vcep"),
            "url": erepo.get("url"),
        })


def _mint_interpret_token(request, req, state, hard_coded_criteria) -> str | None:
    """The token that powers the "Add AI interpretation" button, or None.

    ⚠ MINTED AT PRELIMINARY TIME, NOT AFTER THE GATHER. The button is gated in
    the frontend on `payload.interpret_token` (heartvar.js: `const _tok =
    payload.interpret_token`), and the preliminary payload did not carry one. So
    in evidence-only mode the summary card appeared WITHOUT the button and only
    grew one when db_only_complete re-rendered it — i.e. after the literature
    chain, the very wait the preliminary card exists to skip. The earlier and
    more reliably the preliminary fires, the more visible that gap gets.

    Safe to mint this early because interpret_cache.store holds the payload BY
    REFERENCE and `state` is the same dict the gather keeps populating: by the
    time a curator can click, literature has landed and state["evidence"]
    contains it. `hard_coded_criteria` is refreshed on the held payload when the
    gather finishes, so a click always resumes from the final criteria even
    though the token predates them.

    Returns None when the deployment has no server AI, or when AI needs an
    identity and this caller has none — the UI then keeps the sign-in pill it
    showed before this feature existed.
    """
    if not SERVER_AI_AVAILABLE:
        return None
    account = session_account(request)
    if account is None and ai_auth_required():
        return None
    return interpret_cache.store(
        {"req": req, "state": state, "hard_coded_criteria": hard_coded_criteria},
        account_key(account) if account is not None else None,
    )


def _preliminary_payload(req, state: dict, t_request_start: float, *,
                         interpret_token: str | None = None,
                         ai_unavailable: str | None = None) -> dict | None:
    """The deterministic classification, computed before literature lands.

    Runs the SAME pipeline the evidence-only path runs — compute_hard_coded_criteria
    -> infer_supplementary_criteria -> build_no_ai_criteria ->
    apply_cross_criterion_exclusions -> _gate_criteria_applicability -> score —
    so the tier shown here is the tier that arrives in db_only_complete, not an
    approximation of it. Nothing in that chain reads pubmed / pubtator3 / pmcoa /
    gene_literature (verified: zero references in acmg/hard_coded.py and the
    supplementary inference), which is precisely why it can run early.

    Deliberately does NOT carry `db_evidence`. That payload is ~550 KB and is
    what the final event is for; this one exists to put a number and a criteria
    table on screen, and shipping the evidence twice would cost more than it
    saves.

    Returns None on ANY problem. It is an extra event on top of an unchanged
    path, so a failure here must degrade to "no preliminary banner" and never
    disturb the real result.
    """
    try:
        evidence = state.get("evidence") or {}
        if not evidence:
            return None
        gene = state.get("gene") or (req.gene or "").strip()
        clinical_context = _build_clinical_context(req, evidence)
        hard_coded_criteria = compute_hard_coded_criteria(
            evidence, clinical_context, gene,
        )
        supplementary = (
            infer_supplementary_criteria(
                evidence, clinical_context, gene, hard_coded_criteria,
            )
            if NO_AI_SUPPLEMENTARY else []
        )
        criteria = build_no_ai_criteria(hard_coded_criteria, supplementary)
        criteria, forced = apply_cross_criterion_exclusions(
            criteria, gene, clinical_context.get("inheritance_input"),
            has_curator_segregation=bool(
                clinical_context.get("seg_affected_carriers")
            ),
            zygosity=clinical_context.get("zygosity"),
            in_trans_pathogenic=clinical_context.get("in_trans_pathogenic"),
        )
        criteria = _gate_criteria_applicability(criteria, gene, evidence, req.hpo)
        criteria = apply_pm1_collision_rules(criteria, gene)
        criteria = cap_all_to_spec_strength(criteria, gene)
        points_total = compute_points_total(criteria)
        classification = classification_for_criteria(
            points_total, criteria, forced,
            validity=gene_validity_ceiling(gene),
        )
        return {
            "criteria": criteria,
            "points_total": points_total,
            "classification": classification,
            "variant_id": state.get("variant_id"),
            "gene": gene,
            "awaiting": ["pubmed", "pubtator3", "pmcoa", "gene_literature"],
            "elapsed_s": round(perf_counter() - t_request_start, 2),
            "interpret_token": interpret_token,
            "ai_unavailable": ai_unavailable,
        }
    except Exception:  # pragma: no cover - never break the real result
        log.exception("preliminary classification failed — continuing without it")
        return None


@app.post("/api/curate/stream")
@limiter.limit(_CURATE_RATE_LIMIT)
async def curate_stream(request: Request, req: CurationRequest) -> StreamingResponse:
    hgvs_c = (req.hgvs_c or "").strip()
    _log_begin_curation(
        hgvs_c, req.gene, req.amino_acid, req.hpo,
        req.segregation_context, req.alt_cause_detail,
    )
    if not hgvs_c:
        raise HTTPException(400, "hgvs_c (variant) is required")

    if req.ai_mode != "none":
        ai_account = require_ai_account(request)
        if ai_account is not None:
            check_ai_user_quota(account_key(ai_account), "curation")

    note_activity()

    _check_daily_cap()

    async def event_gen():
        t_request_start = perf_counter()
        yield _sse("stage1_start", {})

        t_stage = {}
        _t = perf_counter()
        state: dict = {}

        interpret_token: str | None = None
        async for event_str in _gather_evidence_sse(req, state):
            if '"event": "scoring_ready"' in event_str or "event: scoring_ready" in event_str:
                if req.ai_mode == "none" and "evidence" in state:
                    _early_criteria = compute_hard_coded_criteria(
                        state.get("evidence") or {},
                        (state.get("evidence") or {}).get("clinical_context")
                        or _build_clinical_context(req, state.get("evidence") or {}),
                        state.get("gene") or req.gene,
                    )
                    interpret_token = _mint_interpret_token(
                        request, req, state, _early_criteria)
                    early = _preliminary_payload(
                        req, state, t_request_start,
                        interpret_token=interpret_token,
                    )
                    if early is not None:
                        yield _sse("preliminary_classification", early)
                continue
            yield event_str
        t_stage["db_gather"] = perf_counter() - _t
        if "evidence" not in state:
            return
        evidence = state["evidence"]
        variant_id = state["variant_id"]
        clean_hgvs_c = state.get("hgvs_c") or hgvs_c
        gene = state.get("gene") or (req.gene or "").strip()

        clinical_context = _build_clinical_context(req, evidence)
        evidence["clinical_context"] = clinical_context

        _t = perf_counter()
        hard_coded_criteria = compute_hard_coded_criteria(
            evidence, clinical_context, gene,
        )
        t_stage["hard_coded"] = perf_counter() - _t

        # Deterministic, source-attributed disease mechanism (LoF / GoF / DN /
        gene_mech = gene_mechanism(evidence, gene)
        evidence["gene_mechanism"] = gene_mech
        evidence["mechanism_flag"] = mechanism_consistency_flag(
            evidence, gene, mech=gene_mech
        )

        _t = perf_counter()
        evidence["transcript_exons"] = await fetch_transcript_exons(
            (evidence.get("vep") or {}).get("transcript_id")
        )
        t_stage["transcript_exons"] = perf_counter() - _t

        _t = perf_counter()
        evidence["same_site"] = await build_same_site_evidence(
            evidence, gene, evidence.get("transcript_exons"),
        )
        t_stage["same_site"] = perf_counter() - _t

        ai_requested = req.ai_mode != "none"
        ai_unavailable_reason: str | None = None
        if ai_requested:
            if not SERVER_AI_AVAILABLE:
                ai_unavailable_reason = "unconfigured"
            elif not budget.ai_within_budget():
                ai_unavailable_reason = "budget"
                log.warning(
                    "[budget] daily AI spend cap reached — serving evidence-only "
                    "for %s %s", _log_safe(gene), _log_safe(clean_hgvs_c),
                )
        serve_evidence_only = (not ai_requested) or (ai_unavailable_reason is not None)

        if serve_evidence_only:
            supplementary = (
                infer_supplementary_criteria(
                    evidence, clinical_context, gene, hard_coded_criteria,
                )
                if NO_AI_SUPPLEMENTARY else []
            )
            no_ai_criteria = build_no_ai_criteria(
                hard_coded_criteria, supplementary,
            )
            no_ai_criteria, forced_classification = apply_cross_criterion_exclusions(
                no_ai_criteria, gene, clinical_context.get("inheritance_input"),
                has_curator_segregation=bool(
                    clinical_context.get("seg_affected_carriers")
                ),
                zygosity=clinical_context.get("zygosity"),
                in_trans_pathogenic=clinical_context.get("in_trans_pathogenic"),
            )
            no_ai_criteria = _gate_criteria_applicability(
                no_ai_criteria, gene, evidence, req.hpo,
            )
            no_ai_criteria = apply_pm1_collision_rules(
                no_ai_criteria, gene)
            no_ai_criteria = cap_all_to_spec_strength(no_ai_criteria, gene)
            points_total = compute_points_total(no_ai_criteria)
            classification = classification_for_criteria(
                points_total, no_ai_criteria, forced_classification,
                validity=gene_validity_ceiling(gene),
            )
            log.info(
                "[no-ai] %s %s — preliminary class=%s points=%d "
                "(deterministic-only)",
                _log_safe(gene), _log_safe(clean_hgvs_c), classification, points_total,
            )
            pts_no_pp5, tier_no_pp5 = tier_without_clinvar_assertion(
                no_ai_criteria, forced_classification,
                validity=gene_validity_ceiling(gene),
            )
            evidence["carrier_status"] = carrier_status(
                evidence, clinical_context, gene, classification, req.hpo,
            )
            if interpret_token is not None:
                held = interpret_cache.peek(interpret_token)
                if held is not None:
                    held["hard_coded_criteria"] = hard_coded_criteria
                else:
                    interpret_token = None
            if interpret_token is None:
                interpret_token = _mint_interpret_token(
                    request, req, state, hard_coded_criteria)

            t_stage["total"] = perf_counter() - t_request_start
            t_stage["other"] = t_stage["total"] - sum(
                v for k, v in t_stage.items() if k != "total"
            )
            log.info(
                "[timing] %s %s — evidence-only TOTAL %.2fs | %s",
                _log_safe(gene), _log_safe(clean_hgvs_c), t_stage["total"],
                " ".join(f"{k}={v:.2f}s" for k, v in t_stage.items()
                         if k != "total"),
            )
            yield _sse("db_only_complete", {
                "db_evidence": evidence,
                "variant_id": variant_id,
                "hard_coded_criteria": hard_coded_criteria,
                "criteria": no_ai_criteria,
                "points_total": points_total,
                "classification": classification,
                "points_total_without_pp5": pts_no_pp5,
                "classification_without_pp5": tier_no_pp5,
                "preliminary": True,
                "ai_unavailable": ai_unavailable_reason,
                "interpret_token": interpret_token,
            })
            return

        state["gene"] = gene
        state["hgvs_c"] = clean_hgvs_c
        async for event_str in _interpret_sse(
            req, state, hard_coded_criteria, t_request_start,
        ):
            yield event_str

    return StreamingResponse(
        event_gen(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache, no-transform",
            "X-Accel-Buffering": "no",
        },
    )


@app.post("/api/curate/interpret")
@limiter.limit(_CURATE_RATE_LIMIT)
async def curate_interpret(
    request: Request, req: InterpretRequest,
) -> StreamingResponse:
    """Add the AI interpretation to a curation whose evidence is already gathered.

    Backs the "Add AI interpretation" button on an evidence-only result. The
    curator forgot to tick the box (or signed in afterwards); without this they
    had to edit the query and re-run every reference lookup, VEP call and
    literature fetch to reach a step that needs none of them.

    Streams the SAME events as the second half of /api/curate/stream, from the
    same generator (`_interpret_sse`) over the held `state` — so the result is
    byte-identical to the one a ticked box would have produced, and there is no
    second interpretive path to keep in sync.

    The client sends only a token. Evidence never round-trips through the
    browser, so it cannot be edited on the way back to shape the prompt or the
    ACMG score.

    Gates, in this order and all of them the same as an AI-bearing curate call:
    the per-IP rate limit (decorator), the sign-in gate, the token (which is
    account-bound), AI availability + today's budget, the per-user AI quota, and
    the global daily cap. Raised BEFORE the StreamingResponse for the reason
    documented on /api/curate/stream — once the response is on the wire the
    status line is fixed and an error can only truncate a 200.
    """
    account = require_ai_account(request)
    account_id = account_key(account) if account is not None else None

    entry = interpret_cache.get(req.token, account_id)
    if entry is None:
        raise HTTPException(
            410,
            detail={
                "error": "interpret_expired",
                "message": (
                    "This result's evidence is no longer held on the server. "
                    "Re-run the variant with \u201cInclude AI interpretation\u201d "
                    "ticked."
                ),
            },
        )

    if not SERVER_AI_AVAILABLE:
        raise HTTPException(
            503,
            detail={
                "error": "unconfigured",
                "message": "AI interpretation is not configured on this server.",
            },
        )
    if not budget.ai_within_budget():
        log.warning("[budget] daily AI spend cap reached — refusing interpret")
        raise HTTPException(
            503,
            detail={
                "error": "budget",
                "message": (
                    "The shared daily AI limit has been reached. Please try "
                    "again tomorrow."
                ),
            },
        )

    if account is not None:
        check_ai_user_quota(account_id, "curation")

    note_activity()
    _check_daily_cap()

    return StreamingResponse(
        _interpret_sse(
            entry["req"], entry["state"], entry["hard_coded_criteria"],
            perf_counter(),
        ),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache, no-transform",
            "X-Accel-Buffering": "no",
        },
    )
