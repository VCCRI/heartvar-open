#!/usr/bin/env python3
"""build_gencc_db.py — refresh the offline GenCC submissions snapshot.

GenCC publishes no API — only a bulk submissions TSV. The client
(backend/clients/gencc.py) lazily downloads + caches it to
``backend/data/gencc_submissions.json`` on first use and refreshes when the cache
is older than 30 days. This builder performs that refresh explicitly so the
monthly data job keeps the snapshot current, instead of a live thegencc.org
download firing mid-service when the cache ages out.

It reuses the client's own ``_refresh_cache()`` — identical parse + output — so
there is no second copy of the TSV schema to drift. Writes to ``gencc.DATA_FILE``
(``backend/data`` by default, or ``GENCC_SNAPSHOT_PATH``); build_all.sh then
copies it onto the data mount for the ``GENCC_SNAPSHOT_PATH`` override.

Cadence: MONTHLY.
"""

from __future__ import annotations

import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from backend.clients import gencc  # noqa: E402


def main() -> int:
    index, meta = asyncio.run(gencc._refresh_cache())
    n_genes = len(index or {})
    print(
        f"[gencc] wrote {gencc.DATA_FILE} — genes={n_genes}, "
        f"rows={meta.get('row_count')}"
    )
    return 0 if n_genes else 1


if __name__ == "__main__":
    raise SystemExit(main())
