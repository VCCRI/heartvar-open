"""A cdot data provider that keeps its transcripts on disk instead of in memory.

WHY. ``cdot.hgvs.dataproviders.JSONDataProvider`` parses its transcript files
fully into RAM. Measured 2026-08-28 on the two GRCh38 files (906,754
transcripts): **4480 MB RSS**, against a 4 GiB container. The resolver was not
deployable. Resolution itself is cheap — 0.2 ms per variant — so the cost is
entirely the resident transcript set.

Restricting the set to the cardiac panel got it to 498 MB but drops any
variant whose gene is off-panel to REST. That trades correctness for memory when
the problem is the storage layout. A gene
being outside our panel is not a reason to annotate it worse.

cdot is built for this: every lookup routes through ``_get_transcript(tx_ac)``,
with ``_get_transcript_ids_for_gene(gene)`` for the gene walk. Implementing those
two against an indexed SQLite table gives full coverage at flat memory.

Everything above those hooks — the c.→g. arithmetic, exon handling, MANE tag
parsing, versioned-accession fallback — is cdot's own unmodified code. That is
deliberate: reshaping the records here would mean re-deriving cdot's semantics,
which is precisely the class of mistake that has cost this project a week.
"""
from __future__ import annotations

import functools
import gzip
import json
import logging
import os
import sqlite3
import threading
from typing import NamedTuple

log = logging.getLogger("heartvar.cdot_sqlite")

RSS_CEILING_MB = int(os.environ.get("HEARTVAR_RESOLVER_RSS_CEILING_MB", "512"))


def current_rss_mb() -> float:
    """Peak resident set of this process, in MB.

    ru_maxrss is BYTES on Darwin and KILOBYTES on Linux — getting that backwards
    makes the guard either never fire or always fire, so it is handled here once
    rather than at each call site.
    """
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    divisor = 1024 * 1024 if os.uname().sysname == "Darwin" else 1024
    return peak / divisor


def _rows_to_record(blob: bytes) -> dict:
    return json.loads(gzip.decompress(blob).decode("utf-8"))


def _newest_versioned(conn: sqlite3.Connection, accession: str) -> bytes | None:
    """The blob for ``accession``, falling back to its newest version.

    The VEP CLI reports transcript ids WITHOUT a version (`ENST00000355349`)
    while cdot keys on the versioned accession (`ENST00000355349.4`). Missing
    that would silently reinstate the exon-number-only NMD behaviour this
    function exists to replace — a fallback that looks like a working result.
    An EXPLICIT version is never second-guessed.
    """
    row = conn.execute(
        "SELECT data FROM transcripts WHERE accession = ?", (accession,)
    ).fetchone()
    if row:
        return row[0]
    if "." in accession:
        return None
    rows = conn.execute(
        "SELECT accession, data FROM transcripts "
        "WHERE accession >= ? AND accession < ?",
        (f"{accession}.", f"{accession}.￿"),
    ).fetchall()
    best, best_v = None, -1
    for acc, blob in rows:
        try:
            v = int(acc.rsplit(".", 1)[1])
        except (IndexError, ValueError):
            continue
        if v > best_v:
            best, best_v = blob, v
    return best


class ManeSelect(NamedTuple):
    """A gene's MANE Select transcript in both namespaces. Either may be None."""

    ensembl: str | None
    refseq: str | None


def _is_mane_select(tag: str | None) -> bool:
    """True for a cdot build-block ``tag`` naming MANE Select.

    The two namespaces spell it differently and BOTH have to match. RefSeq rows
    carry a human-readable ``"MANE Select"``; Ensembl rows carry a comma list,
    ``"CCDS,gencode_basic,gencode_primary,MANE_Select,Ensembl_canonical"``.
    Normalising space to underscore and casefolding covers both without
    per-format parsing. ``MANE Plus Clinical`` deliberately does NOT match: it is
    a second transcript for genes needing one, not the primary, and picking it
    would change which coding change the classification runs on.
    """
    return "mane_select" in (tag or "").replace(" ", "_").casefold()


@functools.lru_cache(maxsize=4096)
def mane_select_for_gene(
    db_path: str | None, gene: str | None, assembly: str = "GRCh38",
) -> ManeSelect:
    """The gene's MANE Select transcripts, read locally. Never raises.

    WHY THIS EXISTS. A curator typing ``CHD7`` + ``c.5058del`` — no transcript,
    which is what the web form invites — had nothing to anchor to, so
    hgvs_resolver returned None and the curation fell through to
    ``rest.ensembl.org``. On 2026-08-31 Ensembl was returning ReadTimeout, then
    HTTP 500, then 503, and that one variant cost 89 s. The whole point of the
    offline cache is not to depend on that host, and bare-gene HGVS was the hole
    in it.

    This is a LOOKUP, not a guess. Verified against the failing production log:
    the local answer for CHD7 is ``ENST00000423902.7``, which is byte-identical
    to the transcript REST itself eventually selected. Only MANE Select is
    returned — no "first transcript" or "longest CDS" heuristic — and a gene
    without one yields (None, None) so the caller keeps its REST fallback rather
    than running the classification on a transcript nobody designated.

    Cached because a gene walk decompresses every transcript record the gene has
    (CHD7 has 53), and one curation can ask more than once.
    """
    if not db_path or not gene:
        return ManeSelect(None, None)
    from ..localio import connect_ro

    conn = connect_ro(db_path)
    if conn is None:
        return ManeSelect(None, None)
    try:
        rows = conn.execute(
            "SELECT accession, data FROM transcripts WHERE gene = ?", (gene,),
        ).fetchall()
    except sqlite3.Error as exc:
        log.warning("cdot: MANE lookup failed for gene %r (%r)", gene, exc)
        return ManeSelect(None, None)

    ensembl = refseq = None
    for accession, blob in rows:
        try:
            record = _rows_to_record(blob)
        except (OSError, ValueError, TypeError):
            continue
        build = ((record.get("genome_builds") or {}).get(assembly)) or {}
        if not _is_mane_select(build.get("tag")):
            continue
        if _is_refseq_accession(accession):
            if refseq is None:
                refseq = accession
            elif refseq != accession:
                log.warning("cdot: gene %r has multiple RefSeq MANE Select rows "
                            "(%s, %s) — keeping %s", gene, refseq, accession, refseq)
        else:
            if ensembl is None:
                ensembl = accession
            elif ensembl != accession:
                log.warning("cdot: gene %r has multiple Ensembl MANE Select rows "
                            "(%s, %s) — keeping %s", gene, ensembl, accession, ensembl)
    return ManeSelect(ensembl, refseq)


def _is_refseq_accession(accession: str | None) -> bool:
    return bool(accession) and str(accession).upper().startswith(
        ("NM_", "NR_", "XM_", "XR_"))


def exon_cdna_lengths(db_path: str | None, accession: str | None,
                      assembly: str = "GRCh38") -> list[int] | None:
    """Per-exon cDNA lengths in transcript (5'->3') order, or None.

    WHAT THIS IS FOR. It is the offline replacement for
    ``ensembl_vep._fetch_exon_lengths``, which asks Ensembl
    ``/lookup/id?expand=1``. Without it clients/vep_offline had to leave
    ``exon_lengths`` None and compute ``nmd_escape`` from the exon "N/total"
    string alone, so the third ClinGen PVS1 rule — a PTC in the LAST 50 nt of
    the PENULTIMATE exon escapes NMD — could never fire offline. Offline then
    returned nmd_escape False where REST returned True, and the direction is the
    unsafe one: False keeps PVS1 at full strength.

    ⚠⚠ THE OFF-BY-ONE, AND IT WOULD HAVE BEEN SILENT. cdot stores each exon as
        [genomic_start, genomic_end, transcript_index, cdna_start, cdna_end, gap]
    with genomic starts 0-BASED HALF-OPEN and cDNA 1-based inclusive. Measured
    on NM_000257.4: every exon's ``end - start + 1`` is exactly ONE MORE than
    its cDNA span, and the two sums differ by exactly the exon count — 6067
    against 6027 over 40 exons. Using the genomic form would have shifted the
    cumulative cDNA offset by 40 bases on MYH7, which is most of the 50 nt
    window the rule turns on.

    So the cDNA columns are used, and they reproduce Ensembl's own numbers
    exactly: REST for ENST00000355349 gives [41, 56, 209, 144, 157, ...] summing
    to 6027, and so does this.

    ⚠ ORDER COMES FROM THE INDEX, not the list. cdot stores exons in ascending
    GENOMIC order, which on the minus strand is the reverse of transcript order
    — and both ``_predict_nmd_escape`` and ``hard_coded``'s splice grading index
    this list by transcript rank.

    STDLIB ONLY, deliberately: ``make_provider`` drags in cdot -> hgvs ->
    psycopg2 -> ipython, and a blob read must not pay for that. It also means
    exon lengths still work on a deployment where cdot is not importable.
    """
    build, exons = _build_for(db_path, accession, assembly)
    if not exons:
        return None
    try:
        lengths = [int(e[4]) - int(e[3]) + 1 for e in exons]
    except (IndexError, TypeError, ValueError) as exc:
        log.warning("cdot exon lengths: unusable row for %r (%r)", accession, exc)
        return None
    if not lengths or any(n <= 0 for n in lengths):
        return None
    return lengths


def _build_for(db_path: str | None, accession: str | None, assembly: str):
    """(record's build block, exons sorted into transcript order) or (None, None)."""
    if not db_path or not accession or not os.path.isfile(db_path):
        return None, None
    from ..localio import connect_ro

    conn = connect_ro(db_path)
    if conn is None:
        return None, None
    conn.row_factory = None
    try:
        blob = _newest_versioned(conn, accession)
    except sqlite3.Error as exc:
        log.warning("cdot: query failed for %r (%r)", accession, exc)
        return None, None
    if not blob:
        return None, None
    try:
        record = _rows_to_record(blob)
        build = ((record.get("genome_builds") or {}).get(assembly)) or {}
        exons = build.get("exons") or []
        if not exons:
            return None, None
        return build, sorted(exons, key=lambda e: e[2])
    except (AttributeError, IndexError, KeyError, TypeError, ValueError,
            OSError, json.JSONDecodeError) as exc:
        log.warning("cdot: unusable record for %r (%r)", accession, exc)
        return None, None


def transcript_exons(db_path: str | None, accession: str | None,
                     assembly: str = "GRCh38") -> dict | None:
    """The exon ladder in ``fetch_transcript_exons``'s shape, from cdot.

    ``{"ok": True, "transcript_id": str, "strand": ±1,
       "exons": [{"start", "end", "rank"}, ...]}`` with rank in 5'->3' order.

    WHY: ensembl_vep.fetch_transcript_exons asks Ensembl /lookup/id, and that
    call is BOTH slow and, for the transcripts HeartVar mostly picks, doomed.
    Measured 2026-08-29 against rest.ensembl.org:
        /lookup/id/NM_000257.4     -> HTTP 400 in 1.3-1.8 s  ("ID not found")
        /lookup/id/ENST00000355349 -> HTTP 200 in 1.3-1.4 s
    Ensembl resolves ENSEMBL ids only, so on the RefSeq path — which the HGVS
    resolver now makes the norm — the app spent a round trip on a call that
    CANNOT succeed and then a second one on the MANE fallback. On a warm-cache
    curation that pair was ~1.4 s of a ~3.0 s total, i.e. nearly half.

    It also gets a BETTER answer: cdot holds RefSeq, so the ladder is the PICKED
    transcript's rather than a MANE substitute, which is what
    build_same_site_evidence had to flag as "proven contiguous against the MANE
    transcript rather than against the picked one".

    ⚠⚠ THE SECOND OFF-BY-ONE IN THIS RECORD, and it is not the same one as in
    exon_cdna_lengths. cdot's genomic starts are 0-BASED HALF-OPEN; Ensembl's
    Exon array is 1-BASED INCLUSIVE, and codon_genomic_positions compares
    genomic positions against these spans directly (``host[0] <= pos <=
    host[1]``). So ``start + 1``, END UNCHANGED. Verified against live REST for
    ENST00000355349: the shifted form reproduces Ensembl's 40 exons EXACTLY
    (first 23435620-23435660, last 23412740-23412871); the unshifted form does
    not match at all.

    STDLIB ONLY, like exon_cdna_lengths — no cdot/hgvs import.
    """
    build, exons = _build_for(db_path, accession, assembly)
    if not exons:
        return None
    strand = -1 if str(build.get("strand", "")).strip() in ("-", "-1") else 1
    try:
        rows = [
            {"start": int(e[0]) + 1, "end": int(e[1]), "rank": i + 1}
            for i, e in enumerate(exons)
        ]
    except (IndexError, TypeError, ValueError) as exc:
        log.warning("cdot: unusable exon row for %r (%r)", accession, exc)
        return None
    if any(r["end"] < r["start"] for r in rows):
        return None
    return {"ok": True, "transcript_id": accession, "strand": strand,
            "exons": rows}


def make_provider(db_path: str):
    """Build a cdot ``LocalDataProvider`` backed by ``db_path``, or None.

    The class is defined INSIDE the function because its base class comes from
    ``cdot``, which drags in ``hgvs`` → ``psycopg2`` → ``ipython``. A deployment
    with the resolver disabled must not pay that import cost, and
    ``hgvs_resolver`` defers its imports for the same reason.
    """
    if not db_path or not os.path.isfile(db_path):
        log.warning("cdot database not found at %r — HGVS input will use REST",
                    db_path)
        return None
    try:
        from cdot.hgvs.dataproviders.json_data_provider import LocalDataProvider
        from hgvs.exceptions import HGVSDataNotAvailableError
    except Exception as exc:  # pragma: no cover - import environment
        log.warning("cdot is not importable (%r) — HGVS input will use REST", exc)
        return None

    class _OfflineSeqFetcher:
        """A seqfetcher that refuses, so ``hgvs`` never reaches the network.

        ⚠ THIS IS AN OFFLINE HOLE THAT WAS INVISIBLE. ``del``, ``dup``, ``inv``
        and ``delins`` state no reference allele, so ``hgvs`` fills it in by
        FETCHING the transcript sequence — and its default
        ``bioutils.seqfetcher`` does that over HTTP against
        ``rest.ensembl.org:80`` and NCBI. Three things are wrong with that, all
        observed in production on 2026-08-31:

          * it is a live Ensembl dependency, which is precisely what
            HEARTVAR_VEP_OFFLINE exists to remove — and it never appeared in any
            audit of our Ensembl call sites because bioutils does not use httpx,
            so none of our request logging or retry accounting sees it;
          * it fetches VERSION-STRIPPED (``/sequence/id/ENST00000423902``), so
            it can return a sequence from a different Ensembl release than the
            pinned 113 cache. That is a silently WRONG reference base, not a
            slow one. Measured: CHD7 ``c.5058del`` came back stating ref ``T``
            where GRCh38 holds ``A`` at 8:60845257 (Ensembl's own VEP reports
            ``A/-`` at that position, so the locus was right and only the base
            was wrong). Nothing but the FASTA cross-check in
            hgvs_resolver.resolve caught it, and its reward for being right was
            to send the curation to REST;
          * it is a network round trip on the critical path.

        Raising makes ``c_to_g`` fail fast for exactly the edit classes that omit
        their reference, which is what hands them to the FASTA interval path in
        ``hgvs_resolver`` — the path built for them, and the one whose docstring
        records being verified with sockets hard-blocked. The fetcher was never
        meant to be reached; the network simply being available was enough for it
        to pre-empt the offline path.
        """

        source = "HeartVar offline (network fetching disabled)"

        def fetch_seq(self, ac, start_i=None, end_i=None):
            raise HGVSDataNotAvailableError(
                f"sequence fetch for {ac} is disabled: HeartVar resolves an "
                "omitted reference allele from the local FASTA, never over the "
                "network (see cdot_sqlite._OfflineSeqFetcher)"
            )

    class _SqliteProvider(LocalDataProvider):
        """cdot transcripts served a row at a time from SQLite."""

        def __init__(self, path: str, assemblies=None, mode=None, cache=None,
                     seqfetcher=None):
            super().__init__(assemblies=assemblies or ["GRCh38"], mode=mode,
                             cache=cache, seqfetcher=seqfetcher)
            self._path = path
            self._local = threading.local()


        @property
        def _conn(self) -> sqlite3.Connection:
            conn = getattr(self._local, "conn", None)
            if conn is None:
                conn = sqlite3.connect(f"file:{self._path}?mode=ro", uri=True,
                                       check_same_thread=False)
                from ..localio import _tune_for_slow_storage

                _tune_for_slow_storage(conn, self._path)
                self._local.conn = conn
            return conn

        def _get_transcript(self, tx_ac):
            if not tx_ac:
                return None
            row = self._conn.execute(
                "SELECT data FROM transcripts WHERE accession = ?", (tx_ac,)
            ).fetchone()
            return _rows_to_record(row[0]) if row else None

        def _get_transcripts(self, tx_acs) -> dict:
            """Retrieve-many hook. Overridden so a gene walk is ONE query rather
            than one per transcript — cdot's default calls _get_transcript in a
            loop, which is fine for an in-memory dict and needless round trips
            here."""
            acs = [a for a in dict.fromkeys(tx_acs or []) if a]
            if not acs:
                return {}
            out: dict = {}
            for i in range(0, len(acs), 500):
                chunk = acs[i:i + 500]
                placeholders = ",".join("?" * len(chunk))
                for accession, blob in self._conn.execute(
                    f"SELECT accession, data FROM transcripts "
                    f"WHERE accession IN ({placeholders})", chunk
                ):
                    out[accession] = _rows_to_record(blob)
            return out

        def _get_transcript_ids_for_gene(self, gene):
            if not gene:
                return []
            return [r[0] for r in self._conn.execute(
                "SELECT accession FROM transcripts WHERE gene = ?", (gene,)
            )]

        def _get_gene(self, gene):
            """Exact-key lookup, matching JSONDataProvider's bare
            ``self.genes.get(gene)``.

            cdot keys genes by Entrez id (RefSeq) or ENSG (Ensembl), NOT by
            symbol — so a symbol lookup misses here exactly as it does upstream.
            Re-keying by symbol would be a behaviour change disguised as a fix,
            and nothing in HeartVar calls get_gene_info anyway; this exists
            because the base class requires it.
            """
            if not gene:
                return None
            row = self._conn.execute(
                "SELECT data FROM genes WHERE key = ?", (str(gene),)
            ).fetchone()
            return _rows_to_record(row[0]) if row else None

        def _get_contig_interval_tree(self, alt_ac):
            """Only used by ``get_tx_for_region``, which nothing in HeartVar
            calls. Raising is honest: building an interval tree over every
            transcript would reintroduce exactly the resident-memory cost this
            class exists to remove, and silently returning an empty tree would
            report "no transcripts in this region" as a fact."""
            raise NotImplementedError(
                "get_tx_for_region is not supported by the SQLite cdot provider "
                "— it would require an in-memory interval tree over the whole "
                "transcript set, which is what this provider avoids"
            )


        def transcript_count(self) -> int:
            row = self._conn.execute(
                "SELECT value FROM meta WHERE key = 'transcript_count'"
            ).fetchone()
            if row:
                return int(row[0])
            return int(self._conn.execute(
                "SELECT count(*) FROM transcripts").fetchone()[0])

        def meta(self) -> dict:
            try:
                return {k: v for k, v in
                        self._conn.execute("SELECT key, value FROM meta")}
            except sqlite3.Error:
                return {}

    try:
        provider = _SqliteProvider(db_path, seqfetcher=_OfflineSeqFetcher())
        count = provider.transcript_count()
        if count <= 0:
            log.warning("cdot database at %r has no transcripts — HGVS input "
                        "will use REST", db_path)
            return None
        log.info("cdot SQLite provider ready: %s transcripts, meta=%s",
                 f"{count:,}", provider.meta())
        return provider
    except Exception as exc:
        log.warning("could not open the cdot database at %r (%r) — HGVS input "
                    "will use REST", db_path, exc)
        return None
