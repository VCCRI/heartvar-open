"""Offline-strict mode toggle.

When ``HEARTVAR_OFFLINE_STRICT`` is set (1/true/yes/on), the offline-primary
data clients must serve ONLY their local build and MUST NOT fall back to a live
network call on a local miss. The design goal is that, with every local mirror
provisioned, the ONLY outbound network access at curation time is the literature
triplet (PubMed / PubTator3 / PMC).

Default OFF: every client keeps its existing live fallback and behaviour is
byte-for-byte unchanged — so the flag is safe to ship dark and is flipped by the
operator only once IT has provisioned all the local data files (see
``scripts/README.md``).

Read at CALL time (not import) so a test or a deploy can flip it via the
environment without re-importing the module. Cheap enough to call per lookup.
"""

from __future__ import annotations

import os

_TRUE = {"1", "true", "yes", "on"}


def offline_strict() -> bool:
    """True when offline-strict mode is enabled — live fallbacks must be skipped."""
    return os.environ.get("HEARTVAR_OFFLINE_STRICT", "").strip().lower() in _TRUE


def vep_offline_enabled() -> bool:
    """True when ``HEARTVAR_VEP_OFFLINE`` is set — the annotation backbone is
    local and must not reach ``rest.ensembl.org``.

    DELIBERATELY NARROWER THAN ``offline_strict``, and the distinction matters.
    ``offline_strict`` answers "may a client fall back to the network when its
    LOCAL MIRROR misses?", and it is off because twelve clients still want that
    fallback. This answers a different question: "has the operator declared the
    Ensembl annotation path local?" A client can honour the second while still
    keeping its own live fallback under the first.

    It lives here rather than being imported from ``vep_offline`` so that
    callers like ``gnomad`` can ask without importing the VEP module, which
    would pull ``ensembl_vep`` and create an import cycle. Same env var, read at
    CALL time for the same reason as above.
    """
    return os.environ.get("HEARTVAR_VEP_OFFLINE", "").strip().lower() in _TRUE
