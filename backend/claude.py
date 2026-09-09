from __future__ import annotations

import asyncio
import json
import logging
import os
import re
from typing import AsyncIterator

from anthropic import APIConnectionError, APITimeoutError, AsyncAnthropic

from . import budget
from .prompt import SYSTEM_PROMPT

MODEL = os.environ.get("HEARTVAR_MODEL", "").strip() or "claude-sonnet-4-6"
try:
    MAX_TOKENS = int(os.environ.get("HEARTVAR_MAX_TOKENS", "") or 8000)
except ValueError:
    MAX_TOKENS = 8000

try:
    TEMPERATURE = float(os.environ.get("HEARTVAR_TEMPERATURE", "") or 0.0)
except ValueError:
    TEMPERATURE = 0.0

CHAT_MODEL = os.environ.get("HEARTVAR_CHAT_MODEL", "").strip() or "claude-sonnet-4-6"
try:
    CHAT_MAX_TOKENS = int(os.environ.get("HEARTVAR_CHAT_MAX_TOKENS", "") or 400)
except ValueError:
    CHAT_MAX_TOKENS = 400

log = logging.getLogger("heartvar.claude")

_CACHE_TTL = os.environ.get("HEARTVAR_CACHE_TTL", "").strip()

def _int_env(name: str, default: int) -> int:
    """Parse an int env var, falling back to ``default`` on unset/blank/invalid.
    Mirrors the MAX_TOKENS try/except ValueError pattern; factored out so the
    parsing is unit-testable without reloading the module."""
    try:
        return int(os.environ.get(name, "") or default)
    except ValueError:
        return default


_CLAUDE_MAX_ATTEMPTS = _int_env("HEARTVAR_CLAUDE_MAX_ATTEMPTS", 3)
_CLAUDE_BACKOFF = (0.5, 1.5, 3.0)
_CLAUDE_RETRYABLE = (APITimeoutError, APIConnectionError)

STREAM_RESTART = object()

_client: AsyncAnthropic | None = None


def _get_client() -> AsyncAnthropic:
    global _client
    if _client is None:
        _client = AsyncAnthropic()
    return _client


def strip_codefences(text: str) -> str:
    return re.sub(r"```(?:json)?", "", text).strip()


def parse_first_json_object(text: str) -> dict:
    """Extract the first valid JSON object from `text`, tolerating both
    leading code fences and trailing prose. Claude occasionally appends an
    "Here is the result:" preamble or commentary after the closing brace,
    which json.loads rejects with "Extra data". raw_decode reads exactly the
    first JSON value and stops, so trailing content is ignored."""
    stripped = strip_codefences(text)
    if stripped and stripped.rstrip()[-1:] != "}":
        tail = stripped[-100:]
        log.warning("Claude response appears truncated — last 100 chars: %r", tail)
        raise json.JSONDecodeError(
            "Response truncated — max_tokens may be too low",
            stripped,
            len(stripped),
        )
    start = stripped.find("{")
    if start < 0:
        return json.loads(stripped)
    obj, _end = json.JSONDecoder().raw_decode(stripped[start:])
    return obj


def _system_blocks() -> list[dict]:
    cache_control: dict = {"type": "ephemeral"}
    if _CACHE_TTL:
        cache_control["ttl"] = _CACHE_TTL
    return [
        {
            "type": "text",
            "text": SYSTEM_PROMPT,
            "cache_control": cache_control,
        }
    ]


async def warm_system_cache() -> None:
    """Pre-warm the prompt cache for the production system block.

    The per-curation cost is dominated by the COLD-cache WRITE of the ~18.9K-token
    system prompt (≈$0.071/call at the 1.25× write multiplier). On a WARM cache the
    SAME system block bills as a cache READ (0.10×, ≈$0.006). This sends ONE minimal
    request whose ``system`` is the IDENTICAL ``_system_blocks()`` the
    real curations use, so the cache entry the next real curation reads is the one
    written here — turning that curation's cold write into a cheap read.

    Contract for the background pre-warmer (see backend/prewarm.py):
      * NEVER raises — any failure is logged and swallowed (a warm failure must
        never disturb the server; the next real curation just pays the cold write).
      * Tallies its own spend via ``budget.record_usage`` so warming shows up in the
        daily USD budget exactly like a real call (the activity-gated loop also
        checks ``ai_within_budget`` before each warm, so warming can never push
        spend past the daily cap by more than one minimal call).
      * Logs the cache_creation / cache_read token split so warming is observable in
        the ``[budget]`` log — that read figure is the signal the deploy doc tells
        the operator to watch (warm-vs-drop-marker tradeoff).
    """
    try:
        client = _get_client()
        msg = await client.messages.create(
            model=MODEL,
            max_tokens=1,
            temperature=TEMPERATURE,
            system=_system_blocks(),
            messages=[{"role": "user", "content": "warm"}],
        )
        budget.record_usage(MODEL, msg.usage)
        usage = getattr(msg, "usage", None)
        log.info(
            "[prewarm] warmed system cache | cache_creation=%s cache_read=%s "
            "input=%s output=%s",
            getattr(usage, "cache_creation_input_tokens", "?"),
            getattr(usage, "cache_read_input_tokens", "?"),
            getattr(usage, "input_tokens", "?"),
            getattr(usage, "output_tokens", "?"),
        )
    except Exception as e:  # noqa: BLE001 — warming must never disturb the server
        log.warning("[prewarm] warm_system_cache failed (ignored): %r", e)


async def stream_claude(
    user_prompt: str, max_tokens: int = MAX_TOKENS
) -> AsyncIterator[str]:
    """Yield text chunks as Claude streams the response.

    Transient stream failures (dropped connection / read timeout) are retried
    with backoff. On a retry the response is re-emitted from scratch, so a
    ``STREAM_RESTART`` sentinel is yielded first to tell the consumer to discard
    the chunks it accumulated for the failed attempt. Non-transient errors
    propagate to the caller unchanged.
    """
    client = _get_client()
    for attempt in range(_CLAUDE_MAX_ATTEMPTS):
        try:
            if attempt:
                yield STREAM_RESTART
            async with client.messages.stream(
                model=MODEL,
                max_tokens=max_tokens,
                temperature=TEMPERATURE,
                system=_system_blocks(),
                messages=[{"role": "user", "content": user_prompt}],
            ) as stream:
                async for text in stream.text_stream:
                    yield text
                final = await stream.get_final_message()
            budget.record_usage(MODEL, final.usage)
            log.info(
                "Claude output: %s tokens (cap: %d) | Claude input: %s tokens",
                getattr(final.usage, "output_tokens", "?"),
                max_tokens,
                getattr(final.usage, "input_tokens", "?"),
            )
            if final.stop_reason == "max_tokens":
                log.warning(
                    "Claude response truncated at max_tokens=%d (output_tokens=%s) — "
                    "consider raising max_tokens or tightening prompt brevity rules",
                    max_tokens,
                    getattr(final.usage, "output_tokens", "?"),
                )
            return
        except _CLAUDE_RETRYABLE as e:
            if attempt < _CLAUDE_MAX_ATTEMPTS - 1:
                backoff = _CLAUDE_BACKOFF[min(attempt, len(_CLAUDE_BACKOFF) - 1)]
                log.warning(
                    "Claude stream attempt %d/%d failed transiently (%r) — "
                    "retrying in %.1fs",
                    attempt + 1, _CLAUDE_MAX_ATTEMPTS, e, backoff,
                )
                await asyncio.sleep(backoff)
                continue
            log.warning(
                "Claude stream: giving up after %d transient failures (%r)",
                _CLAUDE_MAX_ATTEMPTS, e,
            )
            raise


_CHAT_SYSTEM_PREAMBLE = (
    "You are HeartVar's variant-interpretation assistant. Your ONLY permitted "
    "topic is the single variant described in the EVIDENCE CONTEXT below — "
    "specifically its evidence, its ACMG/AMP criteria, and its HeartVar "
    "classification. Answer EXCLUSIVELY from that EVIDENCE CONTEXT: do not use "
    "outside knowledge and do not invent, infer, or supplement any data not "
    "present in it. If a fact is not in the context, say it is not available in "
    "the loaded evidence and recommend primary-source review. If the curator "
    "asks anything outside this variant's evidence/criteria/classification — "
    "including general genetics or medicine, other variants or genes, patient "
    "management or treatment, writing code, or any non-variant topic — you MUST "
    "refuse in a single sentence such as: 'I can only discuss this variant, its "
    "loaded evidence, and its classification.' Do not answer the off-topic "
    "part. Never give patient-management or treatment advice.\n\n"
    "ANSWER FORMAT — this is a chat box, so keep replies SHORT and skimmable:\n"
    "- Lead with the direct answer in one sentence.\n"
    "- Keep the whole reply under ~120 words. Most questions need just 2-4 "
    "sentences.\n"
    "- If you list items, use at most 3-4 short bullets ('- ' lines), one line "
    "each.\n"
    "- Do NOT use markdown headings (no '#', '##', '###') or horizontal rules. "
    "Use '**bold**' sparingly for a key term, and reference ACMG/AMP criterion "
    "codes (e.g. PM1, PS3) inline.\n"
    "- No preamble, no restating the question, no closing pleasantries. Flag "
    "when something needs primary-source review.\n\n"
    "=== EVIDENCE CONTEXT ===\n"
)

_REPORT_SYSTEM_PREAMBLE = (
    "You are HeartVar's variant-interpretation assistant, writing a CLINICAL "
    "CURATION REPORT for the single variant in the EVIDENCE CONTEXT below. Your "
    "only permitted topic is that variant. Write EXCLUSIVELY from the EVIDENCE "
    "CONTEXT: do not invent, infer from outside knowledge, or supplement any "
    "fact that is not present in it.\n\n"

    "GROUNDING — these rules outrank the output format. A curator pastes this "
    "into a patient record, so an unsupported line is worse than a missing "
    "one:\n"
    "- OMIT any line you cannot support from the EVIDENCE CONTEXT. Do not write "
    "a placeholder, do not write 'unknown', do not hedge — leave it out "
    "entirely and let the shorter report stand.\n"
    "- In particular, HeartVar loads NO penetrance or expressivity data. Never "
    "write a penetrance or expressivity statement, and never attach a citation "
    "to one. Omit those lines.\n"
    "- Identifiers may be ECHOED but never GENERATED. A PMID, MIM number, "
    "ClinVar accession (VCV/RCV/SCV) or transcript accession may appear only if "
    "that exact string appears verbatim in the EVIDENCE CONTEXT. If you would "
    "have to recall or construct one, omit the whole clause instead.\n"
    "- Use HeartVar's classification EXACTLY as given. Do NOT append or invent "
    "a sub-tier: no 'VUS-3A', '3B', 'hot/warm/cold', no numeric suffix of any "
    "kind. HeartVar computes the five ACMG/AMP tiers and nothing finer.\n"
    "- Attribute a criterion only where the context shows it was met. Do not "
    "describe a criterion HeartVar marked not met or not assessed as if it "
    "supported the call.\n\n"

    "OUTPUT FORMAT — plain text, no markdown headings, no bold, no preamble and "
    "no closing remarks.\n\n"
    "NEVER PREFIX A BULLET WITH AN ACMG CODE OR STRENGTH. Every bullet is a "
    "plain clinical sentence.\n"
    "  WRONG: '- PM1 (Moderate): The variant falls within the zinc-finger "
    "domain.'\n"
    "  WRONG: '- PM2: absent from gnomAD.'\n"
    "  RIGHT: '- The variant falls within the GATA-type 2 zinc-finger "
    "DNA-binding domain.'\n"
    "The criterion codes belong in HeartVar's Criteria tab; this document is "
    "prose for a patient record and must read as such.\n\n"
    "Emit ONLY the report, in exactly this order:\n\n"
    "1. One identity line: <transcript>(<GENE>): <c. HGVS>; <p. HGVS>\n"
    "   Include only the parts present in the context; drop the p. term if the "
    "variant has no protein consequence.\n"
    "2. One line: This variant is classified as <TIER> (ACMG/AMP <±N> points).\n"
    "3. A line reading exactly 'Evidence in support of pathogenic "
    "classification:' followed by '- ' bullets IN THIS ORDER:\n"
    "   (a) FIRST, one molecular-observation bullet built from the consequence "
    "type, whether the region is repetitive, and the conservation score — e.g. "
    "'In-frame deletion insertion in a non-repetitive region that is highly "
    "conserved.' Include this whenever the context supplies a consequence, even "
    "though it is not itself a met criterion.\n"
    "   (b) THEN one bullet per met pathogenic criterion.\n"
    "   Omit the whole section when there is nothing supporting a pathogenic "
    "reading.\n"
    "4. A line reading exactly 'Evidence in support of benign classification:' "
    "followed by '- ' bullets for met benign criteria. Omit the whole section "
    "when no benign criterion was met.\n"
    "5. A line reading exactly 'Additional information:' followed by '- ' "
    "bullets covering, ONLY where the context supports each: the gene's disease "
    "mechanism and associated conditions; the mode of inheritance; the "
    "proband's zygosity; comparable or same-residue variants and their reported "
    "classifications; whether prior evidence of pathogenicity exists; whether "
    "published segregation evidence exists; whether published functional "
    "evidence exists; and how the variant was inherited in this family.\n\n"

    "Write each bullet as one clinical sentence in a senior clinical "
    "geneticist's register. Statements are about the variant, not about "
    "databases "
    "('absent from gnomAD', not 'gnomAD shows absent'). Where the context "
    "records an ABSENCE of evidence, say so plainly ('No published segregation "
    "evidence has been identified for this variant.'): a checked-and-absent "
    "fact is worth stating, an unchecked one is not.\n\n"

    "=== EVIDENCE CONTEXT ===\n"
)

_CHAT_CONTEXT_MAX_CHARS = 8000

_REPORT_CONTEXT_MAX_CHARS = 16000

try:
    REPORT_MAX_TOKENS = int(os.environ.get("HEARTVAR_REPORT_MAX_TOKENS", "") or 1400)
except ValueError:
    REPORT_MAX_TOKENS = 1400


def report_mode_caps(mode: str | None) -> tuple[int, int]:
    """(context_char_cap, max_output_tokens) for a chat ``mode``.

    Exposed so the endpoint can size its budget RESERVATION from the same place
    the call itself is sized. Booking the chat cap for a report would under-book
    the spend it is about to consume, which is exactly what the reservation
    exists to prevent."""
    if mode == "report":
        return _REPORT_CONTEXT_MAX_CHARS, REPORT_MAX_TOKENS
    return _CHAT_CONTEXT_MAX_CHARS, CHAT_MAX_TOKENS


async def chat_reply(
    context: str,
    question: str,
    history: list[dict] | None = None,
    mode: str = "chat",
) -> str:
    """Answer one chatbot question, grounded in the client-supplied evidence
    ``context`` and scoped by the server-owned preamble. Non-streaming, cheap
    model, short cap. Records usage against the spend budget. Returns the
    answer text (may be a one-line refusal for off-topic questions).

    ``mode="report"`` swaps in the clinical-report preamble and its larger caps.
    Only the FORMAT differs — the report preamble is no less scoped, and adds an
    explicit rule to omit anything the context cannot support."""
    client = _get_client()
    ctx_cap, out_cap = report_mode_caps(mode)
    preamble = _REPORT_SYSTEM_PREAMBLE if mode == "report" else _CHAT_SYSTEM_PREAMBLE
    safe_context = (context or "")[:ctx_cap]
    system = [
        {
            "type": "text",
            "text": preamble + safe_context,
            "cache_control": {"type": "ephemeral"},
        }
    ]
    hist_lines = []
    for turn in (history or [])[-4:]:
        role = "Curator" if (turn.get("role") == "u") else "Assistant"
        text = str(turn.get("text") or "").strip()
        if text:
            hist_lines.append(f"{role}: {text}")
    hist_block = (
        "Earlier in this conversation:\n" + "\n".join(hist_lines) + "\n\n"
        if hist_lines
        else ""
    )
    user_prompt = f"{hist_block}Question: {question}"
    msg = await client.messages.create(
        model=CHAT_MODEL,
        max_tokens=out_cap,
        temperature=TEMPERATURE,
        system=system,
        messages=[{"role": "user", "content": user_prompt}],
    )
    budget.record_usage(CHAT_MODEL, msg.usage)
    return next((b.text for b in msg.content if b.type == "text"), "")
