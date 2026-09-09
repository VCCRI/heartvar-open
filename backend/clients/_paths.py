"""Shared filesystem paths for the local-DB clients.

``clinvar.py`` / ``biogrid.py`` / ``fetal_heart.py`` read the large *mounted*
databases from the repo-root ``data/`` directory, resolved three parents up from
this file (``backend/clients/_paths.py`` → repo root). They previously each
repeated the same ``Path(__file__).resolve().parent.parent.parent`` line; this is
the single source of truth.

Note the two data roots: small *shipped* JSON caches live in ``backend/data/``
(one directory shallower) and are resolved by their own clients (chdgene, gencc,
hpo_labels, panelapp). Do NOT route those through here — the depth differs.
"""
from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"


def data_path(name: str) -> Path:
    """Path to a file under the repo-root mounted ``data/`` directory."""
    return DATA_DIR / name


def db_path(name: str, env_var: str) -> Path:
    """Resolve a mounted-DB path, honouring an env override for deploy
    flexibility (mirrors the ``ALPHAMISSENSE_PATH`` pattern in
    ``alphamissense.py``). Falls back to ``data/<name>`` when ``env_var`` is
    unset or empty.

    NOTE: this is the *read* side only. The matching build scripts
    (``scripts/build_<name>_db.py``) still write to the default ``data/``
    location; if you point a reader at a custom path, regenerate the DB there.
    """
    return Path(os.environ.get(env_var) or (DATA_DIR / name))
