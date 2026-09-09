"""HPO (Human Phenotype Ontology) descendant-set resolution.

Expands an HPO ancestor ID into the flat set of its descendant IDs via the
canonical JAX ontology ``/descendants`` endpoint. Used by the PanelApp client
to expand its handful of cardiac ancestor terms into the full HPO subtree they
cover (so a proband phenotype anywhere under a cardiac ancestor is recognised).

Failure policy: silent. Any 404 / timeout / malformed payload / network outage
returns an empty set so the caller can fall back to its curated map without
raising. Successful fetches are cached at module scope — descendant sets are
stable within a process.
"""

from __future__ import annotations

import asyncio
import logging
import re

import httpx

from ._http_retry import _with_connect_cap, make_async_client

log = logging.getLogger("heartvar.hpo")

HPO_API_BASE = "https://ontology.jax.org/api/hp/terms"
HPO_ID_RE = re.compile(r"HP:\d{7}")

_DESCENDANT_CACHE: dict[str, set[str]] = {}
_DESCENDANT_LOCK = asyncio.Lock()

DESCENDANTS_TIMEOUT = 30.0


def _normalise_id(hpo_id: str) -> str | None:
    """Strip whitespace + uppercase + validate against the canonical
    HP:NNNNNNN shape. Returns None for anything that doesn't match so
    the caller can skip malformed inputs without raising."""
    s = (hpo_id or "").strip().upper()
    return s if HPO_ID_RE.fullmatch(s) else None


class HPOClient:
    """HPO descendant-set resolver.

    Usage::

        client = HPOClient()
        kids = await client.get_descendants("HP:0030680")
    """

    def __init__(self, api_base: str = HPO_API_BASE) -> None:
        self._api_base = api_base.rstrip("/")

    async def get_descendants(self, hpo_id: str) -> set[str]:
        """Fetch every descendant of ``hpo_id`` in a single API call via
        JAX's ``/descendants`` endpoint and return a flat set of the
        descendant IDs (excluding the root itself).

        Failure-silent: 404, timeout, malformed payload, network outage
        all return an empty set so the caller can fall back to its
        curated map without raising. Cached at module scope on first
        successful fetch — descendant sets are stable within a session.

        Implementation note: JAX exposes both
        ``/api/hp/terms/{id}/children`` (one level) and
        ``/api/hp/terms/{id}/descendants`` (whole subtree). The
        descendants endpoint is a single call that returns the same
        data a recursive children-walk would produce, but avoids the
        ~100× request amplification — especially relevant for broad
        ancestor terms (HP:0030680 has ~800 descendants).
        """
        normalised = _normalise_id(hpo_id)
        if not normalised:
            return set()
        if normalised in _DESCENDANT_CACHE:
            return _DESCENDANT_CACHE[normalised]

        try:
            async with make_async_client() as client:
                r = await client.get(
                    f"{self._api_base}/{normalised}/descendants",
                    headers={"Accept": "application/json"},
                    timeout=_with_connect_cap(DESCENDANTS_TIMEOUT),
                )
        except (httpx.HTTPError, asyncio.TimeoutError) as exc:
            log.info("HPO descendants %s failed: %r", normalised, exc)
            return set()
        if r.status_code != 200:
            log.info("HPO descendants %s returned HTTP %s", normalised, r.status_code)
            return set()
        try:
            payload = r.json()
        except ValueError:
            log.info("HPO descendants %s returned non-JSON body", normalised)
            return set()
        if not isinstance(payload, list):
            return set()
        descendants: set[str] = set()
        for entry in payload:
            if not isinstance(entry, dict):
                continue
            d_id = (entry.get("id") or "").strip().upper()
            if HPO_ID_RE.fullmatch(d_id):
                descendants.add(d_id)
        async with _DESCENDANT_LOCK:
            _DESCENDANT_CACHE[normalised] = descendants
        return descendants


_DEFAULT_CLIENT = HPOClient()


async def fetch_hpo_descendants(hpo_id: str) -> set[str]:
    """Convenience wrapper for the descendant traversal used by the
    panelapp client to expand its cardiac ancestor terms into the
    full HPO subtree they cover."""
    return await _DEFAULT_CLIENT.get_descendants(hpo_id)
