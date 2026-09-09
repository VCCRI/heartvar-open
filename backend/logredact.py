"""Keep the curated gene and variant out of the log files.

WHAT THIS IS FOR. The About page tells curators that what they enter is "not
written to a database, logged with their content, or retained after the request
completes". That was not true: 60 log call sites across 18 modules interpolated
the gene symbol, the HGVS string, the gnomAD variant id or the chromosome and
position, and the deployment log keeps every record for the life of the server
process. The claim and the behaviour are now reconciled in favour of the claim.

WHY A FILTER RATHER THAN 60 EDITS. Editing every call site is a one-off that
silently rots: the next log line someone adds re-introduces the leak, and a
mis-edited format string breaks logging at runtime rather than at import. A
filter is one place, it covers call sites nobody has found yet, and it can be
tested directly by logging a gene and asserting the gene is not in the output.

WHY A RANDOM ID AND NOT A HASH. Log lines still need to be groupable — "these
five source failures belong to one curation" is the whole point of having them.
A hash of the gene would do that and is NOT acceptable: there are about 20,000
gene symbols, so a hash is brute-forced in milliseconds, and gene+HGVS is only
a little better. The correlation id here is random and unrelated to the content,
so it groups without recording anything about what was curated. It lives only in
the log line, and the mapping back to a variant exists nowhere.

WHAT IS DELIBERATELY NOT REDACTED. Single-character values. A lone "T" in a log
message tells a reader nothing they could use, and word-boundary-replacing it
would mangle unrelated text. Values of two characters or more are replaced.
"""
from __future__ import annotations

import contextvars
import logging
import re
import secrets

REDACTED = "<redacted>"

_MIN_REDACT_LEN = 2

_ctx: contextvars.ContextVar[dict | None] = contextvars.ContextVar(
    "heartvar_curation", default=None,
)


def begin_curation(*values: object) -> str:
    """Open a curation scope and return its correlation id.

    Every value passed here is scrubbed out of any log record emitted inside
    the scope. Call it once, as early in the request as the gene and variant
    are known, and pass everything that came from the curator.

    Extra values can be registered later with ``add_values`` — the resolved
    HGVS, the genomic coordinates and the gnomAD id are all derived downstream
    of the original input and each is just as identifying.
    """
    cid = secrets.token_hex(4)
    _ctx.set({"id": cid, "patterns": _compile(values)})
    return cid


def add_values(*values: object) -> None:
    """Register further identifiers against the open curation scope. No-op
    outside one, so a client can call it unconditionally."""
    cur = _ctx.get()
    if cur is None:
        return
    pats = dict(cur["patterns"])
    pats.update(_compile(values))
    _ctx.set({"id": cur["id"], "patterns": pats})


def end_curation() -> None:
    """Close the scope. The identifiers go out of scope with it."""
    _ctx.set(None)


def current_id() -> str:
    cur = _ctx.get()
    return cur["id"] if cur else "-"


def _compile(values) -> dict[str, re.Pattern]:
    out: dict[str, re.Pattern] = {}
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if len(s) < _MIN_REDACT_LEN:
            continue
        if s.isdigit() and len(s) < 4:
            continue
        left = r"\b" if s[:1].isalnum() or s[:1] == "_" else ""
        right = r"\b" if s[-1:].isalnum() or s[-1:] == "_" else ""
        out[s] = re.compile(left + re.escape(s) + right)
    return out


class RedactCurationIdentifiers(logging.Filter):
    """Strip the open curation's identifiers from every record, and stamp the
    record with the correlation id.

    ⚠ MUST BE ATTACHED TO HANDLERS, NOT TO A LOGGER. A filter on a logger runs
    only for records logged through that logger's own methods; records from a
    child logger reach an ancestor's HANDLERS but never its filters. Every
    client logs through a child ("heartvar.clinvar"), so a filter on the
    "heartvar" logger sees none of them. Measured, not assumed: with the filter
    on the logger the gene symbol appears in the output verbatim.
    """

    def filter(self, record: logging.LogRecord) -> bool:  # noqa: A003
        cur = _ctx.get()
        record.curation = cur["id"] if cur else "-"
        if cur is None or not cur["patterns"]:
            return True
        try:
            msg = record.getMessage()
        except Exception:  # noqa: BLE001 - a broken format string must not stop logging
            return True
        original = msg
        for pat in cur["patterns"].values():
            msg = pat.sub(REDACTED, msg)
        if msg != original:
            record.msg = msg
            record.args = ()
        return True


def install(*logger_names: str) -> int:
    """Attach the filter to every handler on the named loggers and on root.

    Idempotent, and safe to call again after a handler is added — the
    deployment-log handler is attached during startup, well after import, and
    it is the handler whose output is retained, so it is the one that matters
    most. Returns the number of handlers newly filtered.
    """
    names = logger_names or ("heartvar", "")
    added = 0
    for name in names:
        logger = logging.getLogger(name) if name else logging.getLogger()
        for h in list(logger.handlers):
            if any(isinstance(f, RedactCurationIdentifiers) for f in h.filters):
                continue
            h.addFilter(RedactCurationIdentifiers())
            added += 1
    return added


def unfiltered_handlers(*logger_names: str) -> list[str]:
    """Handlers that would leak, for a test to assert is empty."""
    names = logger_names or ("heartvar", "")
    bad = []
    for name in names:
        logger = logging.getLogger(name) if name else logging.getLogger()
        for h in list(logger.handlers):
            if not any(isinstance(f, RedactCurationIdentifiers) for f in h.filters):
                bad.append(f"{name or 'root'}:{type(h).__name__}")
    return bad
