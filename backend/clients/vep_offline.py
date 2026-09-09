"""Offline Ensembl VEP via the local ``vep`` CLI + indexed cache.

The single biggest source of live traffic in HeartVar is Ensembl VEP: every
curation hits ``rest.ensembl.org`` for consequence / HGVS / MANE / REVEL /
phyloP. ``scripts/setup_offline_vep.sh`` installs the ~27 GB indexed GRCh38
cache (measured: ~25 GB cache + ~0.9 GB FASTA + ~1.5 GB REVEL) as the final step
of ``scripts/build_all.sh``; this module is the runtime half that *uses* it, so
annotation can run entirely on-box.

Design — safe by construction:
  * Gated behind ``HEARTVAR_VEP_OFFLINE`` (default OFF). With it off this module
    is never called and production behaviour is byte-for-byte unchanged.
  * Every entry point returns ``None`` on ANY problem (binary missing, cache
    missing, non-zero exit, timeout, unparseable output, empty result). The
    caller in ``ensembl_vep.py`` then falls through to the live REST path. So
    turning the flag on can only add a local fast-path — it can never break
    annotation.
  * The result dict is the SAME shape ``ensembl_vep`` returns, built from the
    same pure helpers (``_pick_transcript_consequence`` / ``_build_transcript_table``
    / ``_max_dbnsfp`` …), so downstream code (evidence.py, ACMG) can't tell the
    two apart.

Runtime requirement: the ``vep`` binary must be on PATH inside the app
container. It IS — the runtime image is built ``FROM ensemblorg/ensembl-vep``
(see the Dockerfile), and the data-builder image uses the same release via one
``VEP_RELEASE`` build arg. ``vep --offline`` refuses a cache built by a different
release, and it refuses it by returning nothing, which this module reads as
"fall back to REST" — so a release skew costs the offline path silently.

PLUGIN FIELD NAMES (was the "validation gap"): the exact keys in VEP-CLI
``--json`` output (REVEL / phyloP / CADD casing) can differ from the REST
response's lowercased ones. ``_first`` below looks up BOTH spellings so either
works, and ``setup_offline_vep.sh`` now proves it during the build: it refuses to
record REVEL as installed unless a REVEL score actually comes through
``vep --json``, printing ``REVEL OK`` when it does. So that prerequisite is
settled by the build log, not by hand.

``HEARTVAR_VEP_OFFLINE`` IS ON in the deployment as of 2026-08-26, to test this
path. STILL OUTSTANDING, and now measurable rather than blocking: the
offline-vs-REST parity comparison. Field names being right is not the same as the
two paths agreeing on ACMG-affecting values, and REST fallback hides a
disagreement as easily as it rescues a failure — so parity has to be compared on
real curations, not inferred from the absence of errors. Config knobs:

  HEARTVAR_VEP_OFFLINE    "1"/"true"/"yes"/"on" to enable (default off)
  HEARTVAR_VEP_BINARY     vep executable (default "vep")
  HEARTVAR_VEP_DATA       --dir_cache (required when enabled)
  HEARTVAR_VEP_FASTA      --fasta reference (optional but recommended for HGVS)
  HEARTVAR_VEP_ASSEMBLY   --assembly (default "GRCh38")
  HEARTVAR_VEP_REVEL      path to the tabixed REVEL data file. When set, the
                          REVEL plugin is added (--plugin REVEL,file=<path>) so PP3/BP4
                          keeps its calibrated-REVEL signal offline. The
                          release-113 plugin emits a **lowercase** `revel` key
                          (verified against the live cache 2026-08-26); _assemble
                          accepts that plus `REVEL` and dbNSFP's `REVEL_score` /
                          `revel_score`. The claim that it emits `REVEL` is what
                          setup_offline_vep.sh's verification check was written
                          against, and that check then failed for five runs on a
                          working install. ~2-3 GB — far smaller than full
                          dbNSFP; see scripts/setup_offline_vep.sh.
  HEARTVAR_VEP_PLUGINS_DIR  --dir_plugins, where REVEL.pm lives (the setup script
                          installs it under <cache>/Plugins). Optional; VEP's
                          built-in plugin path is used when unset.
  HEARTVAR_VEP_EXTRA_ARGS extra argv appended verbatim, shlex-split (e.g. the
                          --plugin CADD,… / --plugin dbNSFP,… flags once those
                          data files are staged)
  HEARTVAR_VEP_TIMEOUT    per-call subprocess timeout seconds (default 45)
"""
from __future__ import annotations

import asyncio
import glob
import json
import logging
import os
import shlex
import shutil

from ..localio import run_local
from .ensembl_vep import (
    _HIGH_TERMS,
    _complement,
    _MODERATE_TERMS,
    _build_transcript_table,
    _consequence_impact_rank,
    _is_refseq,
    _max_dbnsfp,
    _pick_transcript_consequence,
    _strip_version,
    nmd_escape_from_lengths,
    nmd_escape_needs_exon_lengths,
)

log = logging.getLogger("heartvar.vep.offline")

_TRUE = {"1", "true", "yes", "on"}


def offline_enabled() -> bool:
    """True when HEARTVAR_VEP_OFFLINE is set to a truthy value."""
    return os.environ.get("HEARTVAR_VEP_OFFLINE", "").strip().lower() in _TRUE


def offline_status() -> str:
    """One line describing whether offline VEP will actually be used, for the
    startup banner.

    WHY: turning HEARTVAR_VEP_OFFLINE on changes the annotation backbone, and
    every way it can fail to take effect is SILENT — a missing binary, a cache dir
    that is not mounted, or a release skew all just fall back to live REST. The
    banner reported AI, auth and admin state but nothing about VEP, so "is it
    actually annotating offline?" had no answer short of reading curation output.

    Reports the manifest's release too, because ``vep --offline`` refuses a cache
    built by a different release by returning NOTHING, which reads here as a
    fallback rather than an error. Cheap: two stats and a small JSON read.
    """
    if not offline_enabled():
        return "off (live Ensembl REST)"
    cfg = _config()
    if cfg is None:
        return "SET BUT UNUSABLE — falling back to live REST (see warnings above)"
    data_dir = os.environ.get("HEARTVAR_VEP_DATA", "").strip()

    if data_dir:
        trees = [
            t for t in glob.glob(os.path.join(data_dir, "homo_sapiens*", "*_*"))
            if os.path.isdir(t)
            and any(os.path.isdir(os.path.join(t, e)) for e in os.listdir(t))
        ]
        if trees and not any(os.path.isfile(os.path.join(t, "info.txt")) for t in trees):
            log.warning(
                "offline VEP: no info.txt in any of %s — every call will exit 2 "
                "with 'ERROR: SIFT not available'. The cache stream almost "
                "certainly did not finish (info.txt is the last member of the "
                "tarball). Re-run the data build with FORCE_VEP=1.",
                ", ".join(sorted(trees)),
            )
            return ("SET BUT THE CACHE HAS NO info.txt — every call will exit 2 "
                    "('SIFT not available'); re-run the build with FORCE_VEP=1")

    needy = extra_args_needing_variation_cache(cfg)
    if needy:
        log.warning(
            "offline VEP: HEARTVAR_VEP_EXTRA_ARGS asks for %s, which read the "
            "variation cache. The data build excludes it (17.46 GB of a 21.5 GB "
            "tree, never read otherwise). These flags will return no data rather "
            "than erroring — rebuild with VEP_CACHE_KEEP_VARIATION=1 if they are "
            "genuinely wanted.", " ".join(needy),
        )

    manifest = os.path.join(data_dir, ".heartvar_vep_manifest.json")
    try:
        with open(manifest, encoding="utf-8") as fh:
            rec = json.load(fh)
        release = rec.get("vep_release") or "?"
        revel = "with REVEL" if rec.get("revel") else "no REVEL"
        return f"on for COORDINATE input only (cache release {release}, {revel})"
    except (OSError, ValueError):
        return "on, but no build manifest at " + manifest


def _config() -> dict | None:
    """Resolve the offline VEP invocation config, or ``None`` when it is not
    usable (disabled, or the binary / cache dir is missing). A ``None`` here
    means "silently fall back to REST" — never an error."""
    if not offline_enabled():
        return None
    binary = os.environ.get("HEARTVAR_VEP_BINARY", "vep").strip() or "vep"
    if shutil.which(binary) is None:
        log.warning("HEARTVAR_VEP_OFFLINE set but '%s' not on PATH — using REST", binary)
        return None
    data_dir = os.environ.get("HEARTVAR_VEP_DATA", "").strip()
    if not data_dir or not os.path.isdir(data_dir):
        log.warning("HEARTVAR_VEP_OFFLINE set but HEARTVAR_VEP_DATA=%r is not a "
                    "directory — using REST", data_dir)
        return None
    fasta = os.environ.get("HEARTVAR_VEP_FASTA", "").strip() or None
    if not fasta:
        log.warning("HEARTVAR_VEP_OFFLINE set but HEARTVAR_VEP_FASTA is unset — "
                    "`vep --offline --hgvs` exits 2 without it, so offline VEP "
                    "would fail on every curation. Using REST.")
        return None
    if not os.path.isfile(fasta):
        log.warning("HEARTVAR_VEP_OFFLINE set but HEARTVAR_VEP_FASTA=%r does not "
                    "exist — `vep --offline --hgvs` exits 2 without a readable "
                    "reference. Using REST. (This path embeds the VEP release, "
                    "so check it after a release bump.)", fasta)
        return None
    extra = os.environ.get("HEARTVAR_VEP_EXTRA_ARGS", "").strip()
    try:
        extra_args = shlex.split(extra) if extra else []
    except ValueError:
        log.warning("HEARTVAR_VEP_EXTRA_ARGS is not shell-parseable — ignoring it")
        extra_args = []
    try:
        timeout = float(os.environ.get("HEARTVAR_VEP_TIMEOUT", "45"))
    except ValueError:
        timeout = 45.0
    return {
        "binary": binary,
        "data_dir": data_dir,
        "fasta": fasta,
        "assembly": os.environ.get("HEARTVAR_VEP_ASSEMBLY", "GRCh38").strip() or "GRCh38",
        "revel": os.environ.get("HEARTVAR_VEP_REVEL", "").strip() or None,
        "plugins_dir": os.environ.get("HEARTVAR_VEP_PLUGINS_DIR", "").strip() or None,
        "extra_args": extra_args,
        "timeout": timeout,
    }


_VARIATION_CACHE_FLAGS = frozenset({
    "--check_existing", "--af", "--af_1kg", "--af_gnomad", "--af_gnomade",
    "--af_gnomadg", "--max_af", "--pubmed", "--var_synonyms", "--everything",
    "--check_svs", "--failed",
})


def extra_args_needing_variation_cache(cfg: dict) -> list[str]:
    """Any HEARTVAR_VEP_EXTRA_ARGS entries that need the variation cache.

    Reported rather than stripped: the operator asked for them, and silently
    dropping a flag is the same class of mistake as silently ignoring a missing
    file. Set VEP_CACHE_KEEP_VARIATION=1 on the data build to install the full
    tree if one of these is genuinely wanted.
    """
    return [a for a in cfg.get("extra_args", []) if a.split("=", 1)[0] in _VARIATION_CACHE_FLAGS]


def _argv(cfg: dict, fmt: str, *, refseq: bool) -> list[str]:
    """Build the ``vep`` argv for a single STDIN→STDOUT JSON annotation run.

    ``--json`` emits one JSON object per input variant on stdout;
    ``--no_stats`` keeps the stats summary off stdout so the output is clean
    JSONL. The annotation flags mirror the REST ``VEP_PARAMS`` (hgvs / symbol /
    mane / numbers) so the parsed fields line up with the REST path. Plugin
    flags (CADD / dbNSFP) are supplied by the operator via
    HEARTVAR_VEP_EXTRA_ARGS once their data files are staged."""
    argv = [
        cfg["binary"], "--offline", "--cache",
        "--dir_cache", cfg["data_dir"],
        "--assembly", cfg["assembly"],
        "--json", "--no_stats", "--force_overwrite",
        "-o", "STDOUT", "-i", "STDIN",
        "--format", fmt,
        "--hgvs", "--symbol", "--mane", "--numbers",
        "--canonical", "--biotype",
        "--sift", "b", "--polyphen", "b",
    ]
    if cfg["fasta"]:
        argv += ["--fasta", cfg["fasta"]]
    argv.append("--merged")
    if cfg.get("plugins_dir"):
        argv += ["--dir_plugins", cfg["plugins_dir"]]
    if cfg.get("revel"):
        argv += ["--plugin", f"REVEL,file={cfg['revel']},no_match=1"]
    argv += cfg["extra_args"]
    return argv


def _most_severe_within(reported: str | None, kept: list[dict]) -> str | None:
    """``most_severe_consequence`` restricted to the transcripts we KEPT.

    THE ONE FIELD THE NARROWING FORGOT. _to_rest_transcript_set rebuilt each
    entry as ``{**entry, "transcript_consequences": kept}``, so this field came
    through untouched — and VEP computed it across the whole MERGED
    RefSeq+Ensembl set. When the severe call belonged to a transcript the filter
    then dropped, the offline path reported a consequence REST could not have
    returned, in a field that is read by acmg/hard_coded.py (four sites) and
    evidence.py (four sites), and that gates clients/protvar.py entirely
    (``most_severe_consequence == "missense_variant"``).

    DELIBERATELY CONSERVATIVE. VEP's severity ordering is finer than anything
    worth reimplementing here, and a hand-maintained ranking would be a new
    source of truth that could be silently wrong in the same field. So the
    reported value is replaced ONLY when it is provably unattainable — when no
    kept row claims that term at all. Otherwise it stands exactly as VEP
    computed it.

    When a replacement is needed, the term comes from the highest-ranked KEPT
    row using ``_consequence_impact_rank`` — the same ranking
    ``_pick_transcript_consequence`` already uses — so the answer is always a
    consequence some retained transcript actually carries.
    """
    kept_terms = {t for tc in kept for t in (tc.get("consequence_terms") or [])}
    if not reported or reported in kept_terms:
        return reported
    best = max(kept, key=_consequence_impact_rank, default=None)
    if best is None:
        return reported
    terms = list(best.get("consequence_terms") or [])
    if not terms:
        return reported
    for tier in (_HIGH_TERMS, _MODERATE_TERMS):
        for t in terms:
            if t in tier:
                return t
    return terms[0]


def _to_rest_transcript_set(entries: list[dict], *, refseq: bool) -> list[dict]:
    """Restrict each entry's ``transcript_consequences`` to the transcript set
    the REST path would have returned, so the merged cache does not silently
    widen the annotation.

    WHY THIS EXISTS. The REST client sends two different parameter sets
    (``ensembl_vep.py``): ``_REFSEQ_PARAMS`` (``refseq=1``) for RefSeq input and
    ``VEP_PARAMS`` for everything else, including the whole coordinate path. Each
    returns ONE transcript set — measured on MYH7 c.1208G>A, ``refseq=1`` returns
    2 RefSeq transcripts and no ENST; the default returns 62 ENST and no NM_.
    The MERGED cache returns BOTH sets at once, so without this filter the
    offline path would hand downstream code ~64 rows where REST gave it 2.

    That is not cosmetic. ``_build_transcript_table`` marks a row
    ``is_mane_select`` from a truthy ``mane_select``, and in a merged response
    BOTH members of a MANE pair carry one — NM_000257.4 names
    ENST00000355349.4 and ENST00000355349 names NM_000257.4. Two MANE-Select
    rows then tie through the sort down to ``str(transcript_id)``, where
    "ENST…" sorts before "NM_…", so the Ensembl row wins. ``hard_coded.py``'s
    MANE lookup takes the first MANE row's ``hgvsp`` and PS1/PM5 would
    residue-match ``ENSP00000347507.3:p.Arg403Gln`` against ClinVar instead of
    ``NP_000248.2:p.Arg403Gln`` — a DETERMINISTIC rule change, decided by an
    alphabetical tiebreak. ``transcript_set`` would also read 'mixed' rather
    than 'RefSeq'.

    So the merged cache buys one 25.6 GB tree instead of two 23.2 GB ones, and
    this function is the price: it puts the response back to the shape the REST
    path produced, which is the shape the published numbers were measured on.
    """
    out: list[dict] = []
    for entry in entries:
        tcs = entry.get("transcript_consequences") or []
        if not tcs:
            out.append(entry)
            continue
        kept = [
            tc for tc in tcs
            if _is_refseq(tc.get("transcript_id")) == refseq
            and (refseq or str(tc.get("transcript_id") or "").upper().startswith("ENST"))
        ]
        if not kept:
            continue
        out.append({
            **entry,
            "transcript_consequences": kept,
            "most_severe_consequence": _most_severe_within(
                entry.get("most_severe_consequence"), kept),
        })
    return out


async def _run(cfg: dict, fmt: str, stdin_text: str, *, refseq: bool) -> list[dict] | None:
    """Run ``vep`` once, feeding ``stdin_text`` and parsing the JSONL stdout
    into a list of entry dicts. Returns ``None`` on any failure so the caller
    falls back to REST."""
    argv = _argv(cfg, fmt, refseq=refseq)
    try:
        proc = await asyncio.create_subprocess_exec(
            *argv,
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
    except (OSError, ValueError) as exc:
        log.warning("offline VEP spawn failed (%r) — using REST", exc)
        return None
    try:
        out, err = await asyncio.wait_for(
            proc.communicate(stdin_text.encode()), timeout=cfg["timeout"]
        )
    except asyncio.TimeoutError:
        log.warning("offline VEP timed out after %ss — using REST", cfg["timeout"])
        try:
            proc.kill()
        except ProcessLookupError:
            pass
        return None
    if proc.returncode != 0:
        _err = (err or b"").decode(errors="replace").strip()
        head, tail = _err[:700], _err[-200:]
        log.warning("offline VEP exited %s — using REST. stderr HEAD: %s%s",
                    proc.returncode, head,
                    f" ... TAIL: {tail}" if len(_err) > 900 else "")
        return None
    entries: list[dict] = []
    for line in (out or b"").decode(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            entries.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    if not entries:
        _err = (err or b"").decode(errors="replace").strip()
        log.warning(
            "offline VEP exited 0 but returned no usable annotation (%s bytes "
            "stdout) — using REST. This is a silent-failure shape: check the "
            "input format before the cache. stderr HEAD: %s",
            len(out or b""), _err[:700] or "<empty>",
        )
        return None
    entries = _to_rest_transcript_set(entries, refseq=refseq)
    return entries or None


def _first(d: dict, *keys):
    """Return the first present, non-None value among ``keys``. VEP REST
    lowercases plugin output (revel_score, cadd_phred) while the CLI plugins
    emit their native casing (REVEL_score, CADD_PHRED); checking both makes the
    parser indifferent to which produced the JSON."""
    for k in keys:
        v = d.get(k)
        if v is not None:
            return v
    return None


def _cdot_exon_lengths(tc: dict) -> list[int] | None:
    """Exon lengths for the picked transcript, from the cdot database.

    THE OFFLINE REPLACEMENT FOR ensembl_vep._fetch_exon_lengths, which asks
    Ensembl ``/lookup/id?expand=1``. Without it this module had to leave
    ``exon_lengths`` None, which cost TWO things:

      * the third ClinGen PVS1 NMD rule (last 50 nt of the penultimate exon)
        could never fire, so ``nmd_escape`` was False where REST said True —
        and False keeps PVS1 at FULL STRENGTH, so offline OVER-CALLED
        pathogenicity for that class;
      * acmg/hard_coded.py's splice-skip grading reads ``exon_lengths`` and
        returned None without one, dropping splice PVS1 to base PVS1.

    ⚠ AND IT IS BETTER THAN REST HERE, WHICH IS A BEHAVIOUR CHANGE WORTH
    NAMING. Ensembl's /lookup/id REJECTS RefSeq accessions — measured:
    `{"error":"ID 'NM_000257.4' not found"}` — so on the RefSeq path, which is
    the clinical norm, REST cannot apply rule 3 either and silently degrades.
    cdot holds both transcript sets, so offline now applies the full rule where
    REST never could. The difference is in the SAFE direction: it can only turn
    nmd_escape from False to True, which DOWNGRADES PVS1.

    VERIFIED AGAINST REST on 12 cardiac transcripts: 9 byte-identical
    (including TTN's 363 exons and RYR2's 105). The other 3 differ only because
    cdot is release-113-matched to the VEP cache while REST is on 116 — the
    correlation is exact, every transcript whose VERSION matches gives an
    identical ladder (ENST00000355349 .4/v4 identical; ENST00000334701 .11/v12,
    ENST00000333535 .9/v10, ENST00000070846 .11/v12 all differ). Offline is
    self-consistent at 113; the skew is the one already documented.

    ⚠ KEYED ON THE TRANSCRIPT **VEP** ANNOTATED AGAINST, not the one the curator
    typed. The two genuinely differ: for input NM_002880.3 the resolver honours
    the explicit version for coordinates, while the release-113 cache carries
    NM_002880.4 and reports that — so ``exon`` and ``cdna_start`` are .4-based
    and the ladder must be too. Measured: .3 puts the penultimate 3' end at cDNA
    2218 (last 50 nt 2169..2218) and .4 at 2134 (2085..2134), and the PTC is at
    2086 — so the wrong version flips the verdict. tc["transcript_id"] is the
    only field guaranteed to describe the annotation the numbers came from.
    """
    db = os.environ.get("HEARTVAR_CDOT_DB", "").strip()
    if not db:
        return None
    from . import cdot_sqlite

    return cdot_sqlite.exon_cdna_lengths(db, tc.get("transcript_id"))


def _apply_transcript_strand_repr(result: dict, tc: dict) -> None:
    """Re-express ``strand``/``allele_string`` the way ``/vep/human/hgvs`` does.

    ⚠ THE TWO REST ENDPOINTS DISAGREE, AND BOTH ANSWERS ARE "RIGHT". Measured
    2026-08-29 for MYH7 R403Q — MYH7 is on the minus strand:
        /vep/human/region : strand= 1  allele_string=C/T   (the VCF row)
        /vep/human/hgvs   : strand=-1  allele_string=G/A   (the transcript)
    `vep --format vcf` necessarily gives the first, because that is the input it
    was handed; the per-transcript ``strand`` in its output carries the second.

    So the offline coordinate path already matches REST, and the offline HGVS
    path did not — it is replacing the *hgvs* endpoint, and reported +1 / C/T.

    WHY THAT IS NOT COSMETIC. Three clients read strand and allele_string as a
    PAIR and complement the allele when strand == -1 to get back to the forward
    strand (protvar.py:290, alphamissense.py:153), so either
    convention is safe for them *as long as the pair is consistent*. But
    ensembl_vep.codon_genomic_positions uses ``strand`` ALONE, as the direction
    the CDS runs in:
        first = pos - offset * strand
        walk  = [first + i * strand for i in range(3)]
    Reproduced for MYH7 at 14:23429278 — TWO OF THREE codon positions differ:
        c.1207 (offset 0): strand=-1 -> 277,278,279   strand=1 -> 279,280,281
        c.1208 (offset 1): identical
        c.1209 (offset 2): strand=-1 -> 277,278,279   strand=1 -> 275,276,277
    Turning the flag on would therefore have MOVED THE CODON for minus-strand
    genes on HGVS input, in the field feeding the same-codon / same-residue
    evidence. Silent, and in the direction of a plausible-looking wrong answer.

    ⚠ SEPARATE, PRE-EXISTING, NOT FIXED HERE: the REST *coordinate* path has the
    same problem for the same reason — it stores /vep/human/region's top-level
    +1, so codon_genomic_positions already walks the wrong way for minus-strand
    genes on coordinate input, offline or not. The docstring's own verification
    used GATA4, a PLUS-strand gene, which cannot expose it. Fixing that means
    deciding that codon_genomic_positions should read the transcript's strand
    rather than the result's, which changes REST behaviour and published
    numbers — a deliberate change, not a side effect of this one.
    """
    if tc.get("strand") != -1:
        return
    result["strand"] = -1
    allele = result.get("allele_string") or ""
    if "/" in allele:
        result["allele_string"] = "/".join(
            _complement(part) for part in allele.split("/")
        )


def _assemble(entry: dict, supplied_transcript: str | None, *,
              source_label: str, queried_as: str, assembly_fallback: str,
              coord_extra: dict | None = None,
              transcript_strand_repr: bool = False,
              gene: str | None = None) -> dict | None:
    """Map one VEP entry dict to the HeartVar VEP result shape, reusing the
    same pure helpers the REST path uses. Returns ``None`` when the entry has
    no usable transcript consequence (caller falls back to REST).

    ``transcript_strand_repr`` selects which of REST's two conventions for
    ``strand``/``allele_string`` to reproduce — see _transcript_strand_repr.
    The HGVS path sets it; the coordinate path does not."""
    consequences = entry.get("transcript_consequences") or []
    if not consequences:
        return None
    tc = _pick_transcript_consequence(consequences, supplied_transcript, gene)
    if not tc:
        return None

    mane_select_acc = tc.get("mane_select")
    mane_clinical_acc = tc.get("mane_plus_clinical")
    transcript_table, conseq_differs, transcript_set = _build_transcript_table(
        consequences, tc
    )
    requested_honored = None
    if supplied_transcript:
        _sb = _strip_version(supplied_transcript)
        _picked_mane = mane_select_acc or mane_clinical_acc
        requested_honored = bool(
            _strip_version(tc.get("transcript_id")) == _sb
            or (_picked_mane and _strip_version(_picked_mane) == _sb)
        )

    _exon_lengths = None
    if nmd_escape_needs_exon_lengths(tc) or tc.get("intron"):
        _exon_lengths = _cdot_exon_lengths(tc)

    result = {
        "ok": True,
        "queried_as": queried_as,
        "assembly_name": entry.get("assembly_name") or assembly_fallback,
        "seq_region_name": entry.get("seq_region_name"),
        "start": entry.get("start"),
        "end": entry.get("end"),
        "allele_string": entry.get("allele_string"),
        "strand": entry.get("strand"),
        "transcript_strand": tc.get("strand"),
        "most_severe_consequence": entry.get("most_severe_consequence"),
        "transcript_id": tc.get("transcript_id"),
        "gene_id": tc.get("gene_id"),
        "gene_symbol": tc.get("gene_symbol"),
        "biotype": tc.get("biotype"),
        "hgvsc": tc.get("hgvsc"),
        "hgvsp": tc.get("hgvsp"),
        "consequence_terms": tc.get("consequence_terms"),
        "exon": tc.get("exon"),
        "intron": tc.get("intron"),
        "exon_lengths": _exon_lengths,
        "cdna_start": tc.get("cdna_start"),
        "cdna_end": tc.get("cdna_end"),
        "cds_start": tc.get("cds_start"),
        "cds_end": tc.get("cds_end"),
        "protein_start": tc.get("protein_start"),
        "protein_end": tc.get("protein_end"),
        "nmd_escape": nmd_escape_from_lengths(tc, _exon_lengths),
        "impact": tc.get("impact"),
        "sift_prediction": tc.get("sift_prediction"),
        "sift_score": tc.get("sift_score"),
        "polyphen_prediction": tc.get("polyphen_prediction"),
        "polyphen_score": tc.get("polyphen_score"),
        "cadd_phred": _first(tc, "cadd_phred", "CADD_PHRED"),
        "cadd_raw": _first(tc, "cadd_raw", "CADD_RAW"),
        "revel_score": _max_dbnsfp(
            _first(tc, "revel_score", "REVEL_score", "REVEL", "revel")
        ),
        "phylop100way": _max_dbnsfp(
            _first(tc, "phylop100way_vertebrate", "phyloP100way_vertebrate")
        ),
        "selected_transcript_id": (
            supplied_transcript
            if (supplied_transcript and source_label == "user-supplied")
            else tc.get("transcript_id")
        ),
        "selected_transcript_source": source_label,
        "is_mane_select": bool(mane_select_acc),
        "is_mane_clinical": bool(mane_clinical_acc),
        "mane_select_accession": mane_select_acc,
        "mane_clinical_accession": mane_clinical_acc,
        "requested_transcript": supplied_transcript,
        "requested_transcript_honored": requested_honored,
        "transcript_consequences_all": transcript_table,
        "consequence_differs_significantly": conseq_differs,
        "transcript_set": transcript_set,
        "annotation_source": "offline",
    }
    if transcript_strand_repr:
        _apply_transcript_strand_repr(result, tc)
    if coord_extra:
        result.update(coord_extra)
    return result


async def fetch_hgvs(gene: str, hgvs_c: str,
                     supplied_transcript: str | None) -> dict | None:
    """Offline equivalent of ``_fetch_vep_uncached`` for HGVS input.

    Only attempts when a transcript accession is available (supplied by the
    curator), because the VEP CLI ``--format hgvs`` parser needs a
    transcript-anchored HGVS (``NM_…:c.…`` / ``ENST…:c.…``) — it has no
    gene-symbol resolution like the REST ``/vep/human/hgvs/{gene}:{hgvs}``
    convenience endpoint. Bare ``gene:c.…`` input returns ``None`` here so the
    caller uses REST (which can resolve the symbol)."""
    cfg = _config()
    if cfg is None:
        return None
    from . import hgvs_resolver

    resolved = await run_local(
        hgvs_resolver.resolve, hgvs_c, supplied_transcript, gene)
    if resolved is None:
        return None
    entries = await _run(cfg, "vcf", resolved.vep_row(),
                         refseq=_is_refseq(resolved.accession))
    if not entries:
        return None
    result = _assemble(
        entries[0], supplied_transcript,
        source_label="user-supplied" if supplied_transcript else "VEP canonical",
        queried_as=f"{resolved.accession}:{hgvs_c}",
        assembly_fallback=cfg["assembly"],
        gene=gene,
        transcript_strand_repr=True,
        coord_extra={
            "input_format": "hgvs",
            "input_build": cfg["assembly"],
            "forward_variant_id": resolved.gnomad_id(),
            "resolved_locally": True,
        },
    )
    if result is not None and resolved.note:
        result["resolver_note"] = resolved.note
    return result
    cfg = _config()  # pragma: no cover
    if cfg is None or not supplied_transcript:
        return None
    refseq = _is_refseq(supplied_transcript)
    tx = _strip_version(supplied_transcript) if refseq else supplied_transcript
    queried = f"{tx}:{hgvs_c}"
    entries = await _run(cfg, "hgvs", queried, refseq=refseq)
    if not entries:
        return None
    return _assemble(
        entries[0], supplied_transcript,
        source_label="user-supplied", queried_as=queried,
        assembly_fallback=cfg["assembly"],
    )


async def fetch_region(chrom: str, pos: int, ref: str, alt: str,
                       build: str, supplied_transcript: str | None) -> dict | None:
    """Offline equivalent of ``_fetch_vep_by_coordinates_uncached``.

    The GRCh37 cache is only available offline when the operator installed it
    (VEP_ASSEMBLY=GRCh37 in the setup script); when the request build doesn't
    match the installed assembly this returns ``None`` and the REST path (which
    routes GRCh37 to grch37.rest.ensembl.org) handles it."""
    cfg = _config()
    if cfg is None:
        return None
    if build != cfg["assembly"]:
        return None
    vcf_line = "\t".join((str(chrom), str(pos), ".", ref, alt, ".", ".", "."))
    input_coords = f"{chrom}-{pos}-{ref}-{alt}"
    entries = await _run(cfg, "vcf", vcf_line, refseq=False)
    if not entries:
        return None
    entry = entries[0]
    tc = _pick_transcript_consequence(
        entry.get("transcript_consequences") or [], supplied_transcript
    )
    hgvsc = (tc.get("hgvsc") or "") if tc else ""
    coord_extra = {
        "input_format": "coordinates",
        "input_build": build,
        "input_coords": input_coords,
        "derived_hgvs": hgvsc.split(":")[-1] if hgvsc else None,
        "derived_gene_symbol": tc.get("gene_symbol") if tc else None,
        "forward_variant_id": input_coords,
    }
    return _assemble(
        entry, supplied_transcript,
        source_label="user-supplied" if supplied_transcript else "VEP canonical",
        queried_as=" ".join(vcf_line.split("\t")),
        assembly_fallback=cfg["assembly"],
        coord_extra=coord_extra,
    )
