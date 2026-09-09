"""Local HGVS → genomic coordinate resolution, so offline VEP can serve HGVS input.

WHY THIS EXISTS. ``vep --offline`` refuses ``--format hgvs`` outright — VEP 113
``Parser/HGVS.pm:104`` throws ``ERROR: Cannot use HGVS format in offline mode``
unconditionally, verified in ``ensemblorg/ensembl-vep:release_113.0``, because
resolving a transcript-anchored HGVS needs a database connection. So the offline
cache could only ever serve COORDINATE input, and curators type HGVS. Every
HGVS-entered curation went to ``rest.ensembl.org``, which is the shared public
quota whose exhaustion or IP-block would take HeartVar down for every user
(``clients/_throttle.py:21``).

This module closes that gap: it converts ``NM_000257.4:c.1208G>A`` to
``14:23429278 C>T`` locally, and ``clients/vep_offline.fetch_hgvs`` then annotates
those coordinates from the local cache. Together with the merged cache, that
takes both input shapes off Ensembl for annotation.

WHAT IT IS NOT. It does not annotate, does not call VEP, and knows nothing about
ACMG. HGVS string in, genomic locus out. That keeps it a pure function, testable
without VEP, Docker or a network — which is what makes the parity gate cheap.

DEPENDENCIES, and the reason they are lazy. ``hgvs`` (biocommons) is the reference
implementation of c.→g. arithmetic; reimplementing intronic-offset, UTR and
strand handling ourselves is exactly the off-by-one class that silently moves
ACMG calls, so we do not. It costs us ``psycopg2`` (required BY NAME — ``psycopg2
-binary`` does not satisfy it) and ``ipython`` as runtime dependencies even though
UTA is never touched. Both imports are deferred to first use so a deployment with
the resolver disabled pays nothing for them.

REPRESENTATION IS THE SUBTLE PART, and it is why ``Resolved`` carries two tuples.

  * ``vep_*`` is what goes to ``vep --format vcf``. VEP re-derives HGVS from the
    genomic change, so an un-normalised indel still round-trips EXACTLY —
    measured: ``NM_000257.4:c.1207_1209del`` → ``14:23429277 CCG>''`` → VEP →
    ``NM_000257.4:c.1207_1209del`` / ``NP_000248.2:p.Arg403del``. Annotation
    parity therefore needs no normalisation from us.
  * ``vcf_*`` is the LEFT-ANCHORED form, and it is NOT optional. gnomAD stores
    anchored left-aligned keys, and ``ensembl_vep.gnomad_variant_id_with_provenance``
    prefers VEP's own ``vcf_string`` for exactly this reason — but the VEP CLI
    does NOT emit ``vcf_string`` (verified: with and without ``--var_synonyms``
    the top-level keys are only allele_string, assembly_name, end, id, input,
    most_severe_consequence, seq_region_name, start, strand,
    transcript_consequences). It is a REST-only field. So the offline path has to
    produce the anchored key itself or every indel gets a gnomAD key that "can
    never match" — silently costing PM2/BA1/BS1 on indels.
"""
from __future__ import annotations

import functools
import gzip
import json
import logging
import os
import re
from dataclasses import dataclass

log = logging.getLogger("heartvar.hgvs_resolver")

_NC_SPECIAL = {23: "X", 24: "Y", 12920: "MT"}

_REFSEQ_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "NP_", "NC_")


@dataclass(frozen=True)
class Resolved:
    """A resolved genomic locus, in the two representations that have different
    consumers. See the module docstring on why both exist."""

    chrom: str
    vep_pos: int
    vep_ref: str
    vep_alt: str
    vcf_pos: int
    vcf_ref: str
    vcf_alt: str
    accession: str
    note: str | None = None

    def vep_row(self) -> str:
        """The stdin row for ``vep --format vcf``.

        ⚠ TABS. ``Parser/VCF.pm`` splits on TAB; a space-separated row gives exit
        0, zero bytes of stdout and seven "uninitialized value $ref" warnings —
        a completely silent failure. This is the same trap that made
        ``fetch_region`` inert (see ``vep_offline.fetch_region``), so the row is
        built here once rather than by each caller.
        """
        return "\t".join((
            self.chrom, str(self.vep_pos), ".",
            self.vep_ref or "-", self.vep_alt or "-", ".", ".", ".",
        ))

    def gnomad_id(self) -> str:
        """The anchored left-aligned ``chrom-pos-ref-alt`` key gnomAD stores."""
        return f"{self.chrom}-{self.vcf_pos}-{self.vcf_ref}-{self.vcf_alt}"


def resolver_enabled() -> bool:
    """Off unless HEARTVAR_HGVS_RESOLVER is truthy.

    Default OFF for the same reason HEARTVAR_VEP_OFFLINE is: this changes the
    coordinates the annotation backbone runs on, and the parity gate has to
    pass first.
    """
    return os.environ.get("HEARTVAR_HGVS_RESOLVER", "").strip().lower() in {
        "1", "true", "yes", "on",
    }


def _is_refseq(accession: str | None) -> bool:
    return bool(accession) and str(accession).upper().startswith(_REFSEQ_PREFIXES)


def _nc_to_chrom(accession: str) -> str | None:
    """``NC_000014.9`` → ``"14"``. Returns None for anything that is not an
    NC_ chromosome accession, so a non-chromosomal alignment is a resolution
    failure rather than a silently wrong chromosome."""
    m = re.fullmatch(r"NC_0+(\d+)\.\d+", str(accession or ""))
    if not m:
        return None
    n = int(m.group(1))
    if n in _NC_SPECIAL:
        return _NC_SPECIAL[n]
    return str(n) if 1 <= n <= 22 else None


def _db_path() -> str:
    """The cdot SQLite database. Preferred over the JSON files, by a lot.

    HEARTVAR_CDOT_FILES (the JSON route) is still honoured, because a build that
    staged the JSONs before the .db existed should not simply stop working — but
    it costs 4480 MB of RSS against a 4 GiB container, so it warns loudly. See
    clients/cdot_sqlite.
    """
    return os.environ.get("HEARTVAR_CDOT_DB", "").strip()


def _data_files() -> list[str]:
    """The cdot transcript JSONs, RefSeq first.

    Two files, deliberately. The RefSeq set serves ``NM_``/``NR_`` input, which is
    the clinical norm; the Ensembl set serves ``ENST`` input, which curators can
    also supply (``_is_refseq`` exists in the REST client for that reason). With
    only the RefSeq file, ENST-anchored HGVS raises
    ``HGVSDataNotAvailableError: No alignments`` — measured.

    The Ensembl file must be RELEASE-MATCHED to the VEP cache (113), because a
    transcript present in one and absent in the other is a resolution that VEP
    then cannot annotate.
    """
    raw = os.environ.get("HEARTVAR_CDOT_FILES", "").strip()
    if raw:
        return [p for p in (s.strip() for s in raw.split(",")) if p]
    return []


@functools.lru_cache(maxsize=1)
def _providers():
    """(data_provider, parser, mapper) or None. Cached: building the provider
    parses ~104 MB of gzipped JSON, so it must happen once per process."""
    db = _db_path()
    files = _data_files()
    if not db and not files:
        log.warning("HEARTVAR_HGVS_RESOLVER is set but neither HEARTVAR_CDOT_DB "
                    "nor HEARTVAR_CDOT_FILES is set — HGVS input will use REST")
        return None
    if not db:
        missing = [p for p in files if not os.path.isfile(p)]
        if missing:
            log.warning("cdot transcript file(s) missing (%s) — HGVS input will "
                        "use REST", ", ".join(missing))
            return None
    try:
        import hgvs.assemblymapper
        import hgvs.parser

        if db:
            from . import cdot_sqlite

            hdp = cdot_sqlite.make_provider(db)
            if hdp is None:
                return None
        else:
            from cdot.hgvs.dataproviders import JSONDataProvider

            log.warning("using HEARTVAR_CDOT_FILES (JSON) — this loads the whole "
                        "transcript set into memory (~4.5 GB measured). Set "
                        "HEARTVAR_CDOT_DB to the cdot SQLite instead.")
            hdp = JSONDataProvider(files)
        parser = hgvs.parser.Parser()
        mapper = hgvs.assemblymapper.AssemblyMapper(
            hdp, assembly_name="GRCh38", alt_aln_method="splign",
            normalize=False, replace_reference=False,
            prevalidation_level="NONE",
        )
        return hdp, parser, mapper
    except Exception as exc:  # pragma: no cover - import/build failure
        log.warning("could not build the HGVS resolver (%r) — HGVS input will "
                    "use REST", exc)
        return None


@functools.lru_cache(maxsize=1)
def _versions_by_accession() -> dict[str, list[str]]:
    """base accession → its versions present in the cdot files, ascending.

    WHY: a curator may type ``NM_000257`` with no version, and the REST path
    DELIBERATELY strips versions before querying ("RefSeq versions are frequently
    stale" — ``vep_offline.fetch_hgvs``). cdot needs an explicit version:
    ``NM_000257:c.1208G>A`` raises ``HGVSDataNotAvailableError: No alignments for
    NM_000257``, while ``NM_000257.3`` (a STALE version) resolves correctly —
    both measured. So a bare accession is filled in from the newest version
    actually present, and the substitution is reported back in ``Resolved.note``
    rather than made silently.
    """
    out: dict[str, list[str]] = {}
    db = _db_path()
    if db:
        import sqlite3

        try:
            conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
            for (acc,) in conn.execute("SELECT accession FROM transcripts"):
                base, _, ver = str(acc).rpartition(".")
                if base and ver.isdigit():
                    out.setdefault(base, []).append(str(acc))
            conn.close()
        except sqlite3.Error as exc:
            log.warning("could not index accession versions in %s (%r)", db, exc)
        for base in out:
            out[base].sort(key=lambda a: int(a.rpartition(".")[2]))
        return out
    for path in _data_files():
        try:
            opener = gzip.open if path.endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8") as fh:
                blob = json.load(fh)
        except (OSError, ValueError) as exc:
            log.warning("could not index versions in %s (%r)", path, exc)
            continue
        for acc in (blob.get("transcripts") or {}):
            base, _, ver = str(acc).rpartition(".")
            if base and ver.isdigit():
                out.setdefault(base, []).append(str(acc))
    for base in out:
        out[base].sort(key=lambda a: int(a.rpartition(".")[2]))
    return out


def _mane_select_for_gene(gene: str | None) -> tuple[str | None, str | None]:
    """(accession, note) for the gene's MANE Select transcript, or (None, None).

    ENSEMBL IS PREFERRED, and that is a compatibility choice rather than a
    clinical one. For bare-gene HGVS the REST cascade queried
    ``/vep/human/hgvs/{gene}:{hgvs}`` with ``VEP_PARAMS`` — NOT
    ``_REFSEQ_PARAMS`` — so it returned ENST/ENSP accessions. vep_offline keys
    ``refseq=`` off the accession this returns, so handing back the RefSeq
    NM_ here would silently flip the transcript set on every transcript-less
    curation. That exact swap was measured once before: transcript_id 100/101
    mismatched, hgvsp 89/101 (ENSP vs NP), and because hard_coded.py's MANE
    lookup reads the first MANE row's hgvsp, PS1/PM5 would residue-match an
    ENSP accession against ClinVar's NP_ one. A deterministic rule change, not
    a display difference.

    So this preserves the namespace the curator was already getting. Whether
    bare-gene input SHOULD instead resolve to RefSeq — which is what ClinVar
    names on, and arguably the better clinical answer — is a separate decision
    that needs its own parity run.
    """
    if not gene:
        return None, None
    db = _db_path()
    if not db:
        return None, None
    try:
        from . import cdot_sqlite

        mane = cdot_sqlite.mane_select_for_gene(db, str(gene).strip())
    except Exception as exc:
        log.warning("resolver: MANE lookup failed for gene %r (%r) — using REST",
                    gene, type(exc).__name__)
        return None, None
    accession = mane.ensembl or mane.refseq
    if not accession:
        log.info("resolver: no MANE Select transcript for gene %r — using REST",
                 gene)
        return None, None
    return accession, (
        f"no transcript was supplied; resolved against {accession}, the MANE "
        f"Select transcript for {gene}"
    )


def _pin_version(accession: str) -> tuple[str, str | None]:
    """Return (accession_to_use, note). Unchanged when already versioned."""
    if "." in accession:
        return accession, None
    versions = _versions_by_accession().get(accession)
    if not versions:
        return accession, None
    chosen = versions[-1]
    return chosen, (
        f"{accession} was supplied without a version; resolved against "
        f"{chosen} (newest available locally)"
    )


def _ref_alt_from_edit(edit) -> tuple[str, str] | None:
    """Genomic ref/alt from an hgvs posedit edit, for every edit class.

    ⚠ THE PROTOTYPE ONLY HANDLED ONE CLASS. It read ``edit.ref``/``edit.alt``
    directly, which exists on ``NARefAlt`` and NOT on ``Dup`` or ``Inv``:
    ``c.1208dup`` and ``c.1207_1209inv`` both raise
    ``AttributeError: 'Dup' object has no attribute 'alt'`` — measured. It went
    unnoticed because the inputs it had been exercised on are all NARefAlt, so
    the two crashing classes were never reached. Substitutions, deletions,
    delins, insertions, intronic offsets (``c.1208+5``) and 3'UTR (``c.*34``) are
    all NARefAlt and were fine.

    Returning None rather than raising keeps the module's contract: it either
    resolves or says it cannot.
    """
    kind = type(edit).__name__
    ref = getattr(edit, "ref", None)
    if kind == "NARefAlt":
        return str(ref or ""), str(getattr(edit, "alt", None) or "")
    if kind == "Dup":
        if not ref:
            return None
        return str(ref), str(ref) * 2
    if kind == "Inv":
        if not ref:
            return None
        return str(ref), _revcomp(str(ref))
    return None


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _left_anchor(fasta, chrom: str, pos: int, ref: str, alt: str
                 ) -> tuple[int, str, str] | None:
    """Convert a possibly-unanchored change to the anchored, left-aligned VCF
    form gnomAD stores.

    Needed because the VEP CLI does not emit ``vcf_string`` (REST-only), and
    without an anchored key an indel's gnomAD lookup "can never match" — see
    ``ensembl_vep.gnomad_variant_id_with_provenance``. A meaningless miss is worse
    than no lookup, because PM2 reads absence as evidence.

    Standard VCF normalisation: trim a shared suffix, then a shared prefix, then
    while either allele is empty roll one base to the left and re-prepend the
    reference base. ``fasta`` is any object with pysam's
    ``fetch(chrom, start0, end0)``.
    """
    if fasta is None:
        return None
    ref, alt = ref or "", alt or ""
    if ref == alt:
        return None
    try:
        while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        while not ref or not alt:
            if pos <= 1:
                return None
            base = fasta.fetch(chrom, pos - 2, pos - 1).upper()
            if not base:
                return None
            ref, alt, pos = base + ref, base + alt, pos - 1
            while (len(ref) > 1 and len(alt) > 1
                   and ref[-1] == alt[-1] and pos > 1):
                ref, alt = ref[:-1], alt[:-1]
                if not ref[1:] and not alt[1:]:
                    break
                base = fasta.fetch(chrom, pos - 2, pos - 1).upper()
                if not base:
                    break
                ref, alt, pos = base + ref, base + alt, pos - 1
            break
    except (OSError, ValueError) as exc:
        log.warning("left-anchoring failed at %s:%s (%r)", chrom, pos, exc)
        return None
    if not ref or not alt:
        return None
    return pos, ref, alt


@functools.lru_cache(maxsize=1)
def _fasta():
    """The reference, for left-anchoring only. Reuses HEARTVAR_VEP_FASTA — the
    same file ``--fasta`` already needs, so this adds no new artefact."""
    path = os.environ.get("HEARTVAR_VEP_FASTA", "").strip()
    if not path or not os.path.isfile(path):
        return None
    try:
        import pysam
        return pysam.FastaFile(path)
    except Exception as exc:  # pragma: no cover - depends on the deployed file
        log.warning("could not open HEARTVAR_VEP_FASTA=%r for anchoring (%r) — "
                    "indel gnomAD keys will be unavailable", path, exc)
        return None


def _nc_for(accession: str, hdp) -> str | None:
    """The NC_ chromosome accession a transcript aligns to, from cdot alone.

    Needed by the offline path, which has to name the target accession itself
    rather than letting AssemblyMapper pick it. Returns None when the transcript
    has no chromosomal alignment (a scaffold-only transcript), which is a
    resolution failure rather than a guess.
    """
    try:
        options = hdp.get_tx_mapping_options(accession) or []
    except Exception:
        return None
    for row in options:
        acc = row.get("alt_ac") if isinstance(row, dict) else None
        if acc and _nc_to_chrom(acc):
            return acc
    return None


def _offline_ref_alt(hdp, parsed, accession: str, nc_ac: str, chrom: str,
                     fasta) -> tuple[int, str, str] | None:
    """Genomic (pos, ref, alt) for the edit classes whose HGVS OMITS the reference
    allele, computed without any sequence fetch.

    WHY A SECOND PATH EXISTS. ``del``, ``dup``, ``inv`` and ``delins`` do not state
    their reference bases, so ``hgvs``'s own ``c_to_g`` goes and fetches the
    TRANSCRIPT sequence to fill them — measured, it asks for
    ``get_seq('NM_000257.4', 1311, 1314)`` and that lands on
    ``bioutils.seqfetcher`` → NCBI. That is a per-variant network call, which is
    the thing this module exists to remove.

    ``AlignmentMapper`` maps the INTERVAL with no sequence access at all
    (verified with sockets hard-blocked: RefSeq and Ensembl, del/dup/inv, zero
    network), and the reference then comes from the FASTA — the same genome VEP
    annotates against, which is a better authority for a genomic ref than a
    transcript record anyway.

    Only these four classes are handled here. Everything ``hgvs`` can do offline
    (substitutions in CDS/intron/UTR, insertions) stays with ``hgvs``, because
    those need strand-aware allele handling and re-implementing that is the
    off-by-one class we are avoiding.
    """
    if fasta is None:
        return None
    try:
        from hgvs.alignmentmapper import AlignmentMapper

        am = AlignmentMapper(hdp, accession, nc_ac, "splign")
        gi = am.c_to_g(parsed.posedit.pos)
        start, end = int(gi.start.base), int(gi.end.base)
        if end < start:
            return None
        gref = (fasta.fetch(chrom, start - 1, end) or "").upper()
        if not gref:
            return None
        strand = getattr(am, "strand", 1)
    except Exception as exc:
        log.info("offline interval mapping failed for %s (%r)", accession, exc)
        return None

    edit = parsed.posedit.edit
    kind = type(edit).__name__
    stated_alt = getattr(edit, "alt", None)

    if kind == "Dup":
        return start, gref, gref * 2
    if kind == "Inv":
        return start, gref, _revcomp(gref)
    if kind == "NARefAlt":
        if not stated_alt:
            return start, gref, ""
        alt = _revcomp(str(stated_alt)) if strand == -1 else str(stated_alt)
        return start, gref, alt
    return None


def resolve(hgvs_c: str, transcript: str | None, gene: str | None = None
            ) -> Resolved | None:
    """``NM_000257.4``, ``c.1208G>A`` → a genomic locus, or None.

    None means "use REST" and is never an error: a resolver that cannot resolve
    this one variant must not take the curation down with it. Every failure path
    logs at most a warning.
    """
    if not resolver_enabled():
        return None
    if not hgvs_c:
        return None
    mane_note: str | None = None
    if not transcript:
        transcript, mane_note = _mane_select_for_gene(gene)
        if not transcript:
            return None
    built = _providers()
    if built is None:
        return None
    hdp, parser, mapper = built
    accession, note = _pin_version(str(transcript).strip())
    note = "; ".join(n for n in (mane_note, note) if n) or None
    var = f"{accession}:{str(hgvs_c).strip().split(':', 1)[-1]}"
    try:
        parsed = parser.parse_hgvs_variant(var)
    except Exception as exc:
        log.info("resolver could not parse %s (%s) — using REST",
                 var, type(exc).__name__)
        return None

    pos = ref = alt = chrom = None
    try:
        g = mapper.c_to_g(parsed)
        chrom = _nc_to_chrom(getattr(g, "ac", ""))
        ra = _ref_alt_from_edit(g.posedit.edit)
        if chrom is not None and ra is not None:
            ref, alt = ra
            pos = int(g.posedit.pos.start.base)
    except Exception as exc:
        log.debug("hgvs c_to_g needed sequence for %s (%s) — trying the offline "
                  "interval path", var, type(exc).__name__)

    if pos is not None and ref:
        fa = _fasta()
        if fa is not None:
            try:
                actual = (fa.fetch(chrom, pos - 1, pos - 1 + len(ref)) or "").upper()
            except (OSError, ValueError):
                actual = ""
            if actual and actual != ref.upper():
                log.info(
                    "resolver: %s states a reference that disagrees with GRCh38 "
                    "at %s:%s (stated %r, genome %r) — using REST so the input "
                    "error is reported", var, chrom, pos, ref, actual)
                return None

    if pos is None:
        nc_ac = _nc_for(accession, hdp)
        chrom = _nc_to_chrom(nc_ac or "")
        if chrom is None:
            log.info("resolver found no chromosomal alignment for %s — using REST",
                     accession)
            return None
        got = _offline_ref_alt(hdp, parsed, accession, nc_ac, chrom, _fasta())
        if got is None:
            log.info("resolver cannot express %s offline — using REST", var)
            return None
        pos, ref, alt = got

    anchored = _left_anchor(_fasta(), chrom, pos, ref, alt)
    if anchored is None:
        if len(ref) == 1 and len(alt) == 1:
            anchored = (pos, ref, alt)
        else:
            log.info("resolver could not anchor %s — using REST so the gnomAD "
                     "key stays meaningful", var)
            return None

    return Resolved(
        chrom=chrom, vep_pos=pos, vep_ref=ref, vep_alt=alt,
        vcf_pos=anchored[0], vcf_ref=anchored[1], vcf_alt=anchored[2],
        accession=accession, note=note,
    )
