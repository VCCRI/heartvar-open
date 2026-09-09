"""Short-lived hold on a finished evidence gather, so the AI interpretation can
be added to a result the curator already has.

The friction this removes: a curator who forgets to tick "Include AI
interpretation" before running a variant had to edit the query, tick the box and
sit through the entire annotation a second time — every reference lookup, VEP
call and literature fetch repeated to reach a step that needs none of them. The
AI half of ``/api/curate/stream`` depends only on the gathered evidence, so we
keep that evidence for half an hour and let ``/api/curate/interpret`` resume from
it (see ``_interpret_sse`` in backend/app.py — one code path, both entry points).

Design constraints, all deliberate:

* **In-process.** The app already requires a single worker for its live-source
  TTL cache (see backend/clients/_cache.py), so a module-level dict is
  consistent with the deployment rather than a new assumption. A restart drops
  every token; the UI falls back to "re-run", which is exactly today's
  behaviour.
* **Nothing is written to disk.** An evidence dict carries the curator's
  clinical context — HPO terms, family/segregation fields. Results are
  deliberately un-addressable in this app and this store keeps that property.
* **Account-bound.** ``account`` is the ``account_key`` of the session that ran
  the gather (None only where the deployment has no auth configured and AI is
  therefore open). ``get`` will not return an entry to a different account, so a
  leaked token cannot expose one curator's phenotype data to another.
* **Small and capped.** An evidence dict runs to megabytes (literature, the
  ClinVar landscape, transcript exons), so the cap is entries, not bytes, and it
  is low. Oldest goes first.

Tokens are NOT single-use: a stream that dies mid-generation, or an unparseable
model response, should be retryable inside the window. Spend is bounded by the
AI quota and daily budget at the endpoint, not by burning the token here.
"""
from __future__ import annotations

import secrets
import time
from collections import OrderedDict

TTL_SECONDS = 30 * 60

MAX_ENTRIES = 8

_store: OrderedDict[str, dict] = OrderedDict()


def _now() -> float:
    return time.monotonic()


def _purge(now: float) -> None:
    """Drop expired entries. Called on every store/get — no background task."""
    for token in [t for t, e in _store.items() if now - e["created"] > TTL_SECONDS]:
        del _store[token]


def store(payload: dict, account: str | None) -> str:
    """Hold ``payload`` and return the token that resumes it.

    ``payload`` is opaque here — app.py puts the gather's ``state``, the request
    and the deterministic criteria in it. Keeping this module ignorant of the
    shape means the resume contract lives in one place (``_interpret_sse``).
    """
    now = _now()
    _purge(now)
    token = secrets.token_urlsafe(32)
    _store[token] = {"payload": payload, "account": account, "created": now}
    while len(_store) > MAX_ENTRIES:
        _store.popitem(last=False)
    return token


def get(token: str, account: str | None) -> dict | None:
    """The held payload, or None when the token is unknown, expired, or belongs
    to someone else.

    One return value for all three cases on purpose: the caller cannot
    distinguish "expired" from "another account's token", so the endpoint cannot
    leak the existence of another curator's run. The UI treats every failure the
    same way — fall back to re-running.
    """
    now = _now()
    _purge(now)
    entry = _store.get(token)
    if entry is None or entry["account"] != account:
        return None
    return entry["payload"]


def clear() -> None:
    """Drop everything. For tests and for a clean shutdown."""
    _store.clear()


def size() -> int:
    """Live entry count, after purging. For tests."""
    _purge(_now())
    return len(_store)


def peek(token: str) -> dict | None:
    """The held payload WITHOUT an account check, or None if unknown/expired.

    Deliberately account-blind, and deliberately not a second `get`. The one
    caller is the curate endpoint refreshing the criteria on a token IT just
    minted, inside the same request — there is no cross-account question to
    answer, and threading the account through only to re-check it against
    itself would invite someone to reach for this where `get` belongs.

    `get` remains the ONLY way to resume a run from a token, so the
    account-scoping guarantee on that path is unchanged.
    """
    now = _now()
    _purge(now)
    entry = _store.get(token)
    return entry["payload"] if entry is not None else None
