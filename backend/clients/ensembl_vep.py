from __future__ import annotations

import logging
import math
import os
import re
import urllib.parse

import httpx

from ._cache import EXTERNAL_CACHE, TTL_VEP, _ok
from ._http_retry import _with_connect_cap, request_with_retry
from .conservation import fetch_phylop100way

log = logging.getLogger("heartvar.vep")

ENSEMBL_BASE = "https://rest.ensembl.org"
ENSEMBL_GRCH37_BASE = "https://grch37.rest.ensembl.org"
HEADERS = {
    "Accept": "application/json",
    "User-Agent": "HeartVar/1.0 (variant-curation; heartvar@victorchang.edu.au)",
}


async def _ensembl_get(
    client: httpx.AsyncClient,
    url: str,
    params: dict,
    timeout: float = 30.0,
) -> httpx.Response | None:
    """GET an Ensembl REST URL, retrying transient failures (429, 5xx, and
    transport/read timeouts) with backoff.

    Returns the final ``Response`` on a non-retryable outcome (2xx, or a 4xx
    other than 429 — both of which the caller interprets), or ``None`` when
    every attempt failed transiently.
    """
    return await request_with_retry(
        client, "GET", url, timeout=timeout,
        headers=HEADERS, params=params,
        name=f"Ensembl {url.rsplit('/', 1)[-1]}", logger=log,
    )

_TRANSCRIPT_PREFIX_RE = re.compile(
    r"^(?P<prefix>(?:NM|NR|XM|XR)_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?)"
    r":?(?P<rest>[cgnpr]\..+)$"
)

_GENE_SYMBOL_PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9-]*:(?P<rest>[cnr]\..+)$")

_COORD_INPUT_RE = re.compile(
    r"^(?:chr)?(?P<chrom>\d+|X|Y|MT|M)[:_\-]"
    r"(?P<pos>\d+)[:_\-]"
    r"(?P<ref>[ACGT]+)[:_\-]"
    r"(?P<alt>[ACGT]+)$",
    re.IGNORECASE,
)

_UNICODE_DASH_TABLE = str.maketrans(
    {"‐": "-", "‑": "-", "‒": "-", "–": "-",
     "—": "-", "―": "-", "−": "-"}
)

_RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)


def parse_variant_input(raw_input: str) -> dict:
    """Classify the curator's variant string as HGVS or genomic coordinates.

    Two supported forms:
      - HGVS, optionally with a transcript prefix:
          NM_000257.4:c.1208G>A, ENST00000355349:c.1208G>A, c.1208G>A
      - Genomic coordinates in chr:pos:ref:alt form (also -/_ separators):
          chr7:117548628:C:T, 7-117548628-C-T, 7_117548628_C_T

    Returns:
      - ``{"format": "coordinates", "chrom": "7", "pos": 117548628,
           "ref": "C", "alt": "T", "build": None}`` for coord input. The
        ``build`` field is left ``None`` here — the caller pairs the
        parsed coords with the user-selected genome build from the
        request body.
      - ``{"format": "hgvs", "hgvs": "c.1208G>A",
           "transcript": "NM_000257.4"}`` for HGVS input (``transcript``
        may be ``None`` when no NM_/ENST prefix was supplied).
      - ``{"format": "unknown"}`` when neither form matches.

    Coordinate parsing strips any "chr" prefix, normalises the chromosome
    to bare 1-22 / X / Y / MT (uppercase), and converts the position to
    int. Lowercase ref/alt are upper-cased so downstream services
    (gnomAD, Ensembl VEP /vep/human/region) get the form they expect.
    """
    if not raw_input:
        return {"format": "unknown"}
    s = raw_input.strip().translate(_UNICODE_DASH_TABLE).replace(",", "")
    if _RSID_RE.match(s):
        return {"format": "rsid", "rsid": s.lower()}
    m = _COORD_INPUT_RE.match(s) or _COORD_INPUT_RE.match(re.sub(r"\s+", "", s))
    if m:
        chrom = m.group("chrom").upper()
        if chrom == "M":
            chrom = "MT"
        return {
            "format": "coordinates",
            "chrom": chrom,
            "pos": int(m.group("pos")),
            "ref": m.group("ref").upper(),
            "alt": m.group("alt").upper(),
            "build": None,
        }
    if "c." in s or "g." in s or "n." in s:
        stripped, transcript = strip_transcript_prefix(s)
        if transcript is None:
            gm = _GENE_SYMBOL_PREFIX_RE.match(stripped)
            if gm:
                log.debug("Stripped gene-symbol prefix from HGVS input %r", stripped)
                stripped = gm.group("rest")
        return {
            "format": "hgvs",
            "hgvs": stripped,
            "transcript": transcript,
        }
    if "p." in s:
        return {"format": "protein", "hgvs_p": s}
    return {"format": "unknown"}


def strip_transcript_prefix(hgvs: str) -> tuple[str, str | None]:
    """Detach a leading transcript accession from an HGVS string.

    Returns ``(stripped, prefix)``. When no prefix is present, ``prefix`` is
    ``None`` and ``stripped`` equals the input. A debug log records each
    strip so the curator can confirm via the server console which transcript
    the user supplied.
    """
    if not hgvs:
        return hgvs, None
    m = _TRANSCRIPT_PREFIX_RE.match(hgvs.strip())
    if not m:
        return hgvs, None
    prefix = m.group("prefix")
    rest = m.group("rest")
    log.debug("Stripped transcript prefix %s from HGVS input", prefix)
    return rest, prefix

VEP_PARAMS = {
    "hgvs": "1", "CADD": "1",
    "dbNSFP": "REVEL_score,phyloP100way_vertebrate",
    "mane": "1", "numbers": "1",
    "vcf_string": "1",
}

_REFSEQ_PARAMS = {**VEP_PARAMS, "refseq": "1"}

_REFSEQ_PREFIXES = frozenset({"NM", "NR", "XM", "XR"})


def _dbsnp_rsid_from_entry(entry: dict) -> str | None:
    """Pull the dbSNP rsID from a VEP entry's ``colocated_variants`` (Ensembl VEP
    REST returns these by default for known variants). The rsID is the single
    highest-value tmVar3 entity term for variant→literature retrieval, so it is
    surfaced on the VEP result and threaded into the PubTator3 fetch
    (``evidence._fetch_literature_with_vep_aa``). Returns the first ``rs…`` id,
    or None. Never raises — a missing rsID just falls back to variant_recoder."""
    try:
        for cv in entry.get("colocated_variants") or []:
            if not isinstance(cv, dict):
                continue
            cid = cv.get("id")
            if isinstance(cid, str) and cid.lower().startswith("rs"):
                return cid
    except Exception:
        log.debug("colocated_variants rsID parse failed", exc_info=True)
    return None


def _is_refseq(accession: str | None) -> bool:
    """True if ``accession`` is a RefSeq transcript ID (NM_/NR_/XM_/XR_).
    ENST… and bare/None inputs return False."""
    if not accession:
        return False
    head = accession.split("_", 1)[0].upper()
    return "_" in accession and head in _REFSEQ_PREFIXES


def _max_dbnsfp(score_field) -> float | None:
    """dbNSFP columns come back per-transcript as a comma-separated string like
    '.,0.651,0.651,.,.'. Take the maximum numeric value, ignoring the
    '.' placeholders. Some responses return a plain float — handle both.

    The per-position dbNSFP scores (REVEL, BayesDel, phyloP/phastCons/GERP)
    are identical across a variant's transcripts, so max() simply recovers the
    single real value; it also picks the most-pathogenic / most-conserved value
    in the rare case transcripts disagree, which is the conservative choice."""
    if score_field is None:
        return None
    if isinstance(score_field, (int, float)):
        f = float(score_field)
        return f if math.isfinite(f) else None
    if isinstance(score_field, str):
        vals = []
        for tok in score_field.split(","):
            tok = tok.strip()
            if not tok or tok == ".":
                continue
            try:
                f = float(tok)
            except ValueError:
                continue
            if math.isfinite(f):
                vals.append(f)
        return max(vals) if vals else None
    return None


async def _lookup_canonical_transcript(client: httpx.AsyncClient, gene: str) -> str | None:
    r = await _ensembl_get(
        client,
        f"{ENSEMBL_BASE}/lookup/symbol/homo_sapiens/{urllib.parse.quote(gene, safe='')}",
        {"expand": "1"},
        timeout=25.0,
    )
    if r is None or r.status_code != 200:
        return None
    data = r.json()
    transcripts = data.get("Transcript", []) or []
    for t in transcripts:
        if t.get("is_canonical"):
            return t.get("id")
    return transcripts[0].get("id") if transcripts else None


def _strip_version(accession: str | None) -> str | None:
    """Drop the trailing ".N" version suffix on an accession (NM_000257.4 →
    NM_000257). Returns None unchanged."""
    if not accession:
        return None
    return accession.split(".", 1)[0]


def codon_genomic_positions(vep: dict, exons: list[dict] | None = None) -> list[int] | None:
    """Genomic positions (GRCh38, ascending) of the three bases of the codon
    the variant sits in, or None when they can't be determined.

    Derived entirely from fields the VEP result already carries:
      * ``cds_start``  — 1-based CDS coordinate of the variant, so the
        variant's offset within its codon is ``(cds_start - 1) % 3``.
      * ``transcript_strand`` — which genomic direction the CDS runs in.
      * ``start``      — the variant's genomic position.
    Verified against GATA4 NM_002052.5:c.907G>T: cds_start 907 → offset 0,
    strand +1, start 11750234 → [11750234, 11750235, 11750236], which is
    exactly the span of MANE codon 304 (c.910-912).

    ⚠⚠ IT READS ``transcript_strand``, NOT ``strand``, AND THAT WAS A REAL BUG.
    ``strand`` is the ALLELE-REPRESENTATION strand and means different things
    per endpoint — measured for MYH7 R403Q, a minus-strand gene:
        /vep/human/region : strand= 1  allele_string=C/T   (the VCF row)
        /vep/human/hgvs   : strand=-1  allele_string=G/A   (the transcript)
    So on the COORDINATE path this function received +1 for every gene
    regardless of orientation, and walked the codon the wrong way. Reproduced
    exactly for MYH7 at 14:23429278:
        c.1207 (offset 0): correct 277,278,279   got 279,280,281
        c.1208 (offset 1): identical
        c.1209 (offset 2): correct 277,278,279   got 275,276,277
    TWO OF THREE codon positions wrong, for roughly half of all genes, feeding
    the same-codon / same-residue evidence. It survived because the docstring's
    own verification used GATA4 — a PLUS-strand gene, where the two conventions
    coincide — and because the three OTHER consumers of ``strand``
    (protvar.py:290, alphamissense.py:153) read it PAIRED with
    allele_string and complement together, so they are correct either way.

    ``transcript_strand`` comes straight from the transcript consequence and is
    unambiguous. It is REQUIRED rather than falling back to ``strand``: a field
    that means two different things cannot be trusted as the CDS direction, and
    declining costs only the codon view (the caller degrades to the variant's
    own position) whereas guessing produces a plausible, wrong codon.

    ``exons`` (from :func:`fetch_transcript_exons`) is used only to detect a
    codon that STRADDLES an intron: if a naive walk leaves the exon containing
    the variant, the codon is split across an exon junction and the remaining
    bases are taken from the start of the next exon (or the end of the
    previous one, on the reverse strand). When no exon ladder is supplied we
    return the naive walk but only if it stays clear of needing one — i.e.
    ``None`` when we cannot prove the codon is contiguous. Callers must treat
    None as "cannot offer a codon view", not as "no other alleles".
    """
    if not isinstance(vep, dict) or not vep.get("ok"):
        return None
    cds_start = vep.get("cds_start")
    pos = vep.get("start")
    strand = vep.get("transcript_strand")
    if not isinstance(cds_start, int) or not isinstance(pos, int):
        return None
    if strand not in (1, -1):
        return None
    if vep.get("end") != pos:
        return None
    offset = (cds_start - 1) % 3
    first = pos - offset * strand
    walk = [first + i * strand for i in range(3)]
    ordered = sorted(walk)
    if exons:
        spans = [
            (min(e["start"], e["end"]), max(e["start"], e["end"]))
            for e in exons
            if isinstance(e, dict) and e.get("start") is not None
            and e.get("end") is not None
        ]
        host = next((sp for sp in spans if sp[0] <= pos <= sp[1]), None)
        if host is None:
            return None
        if all(host[0] <= g <= host[1] for g in ordered):
            return ordered
        spans.sort()
        idx = spans.index(host)
        inside = [g for g in walk if host[0] <= g <= host[1]]
        need = 3 - len(inside)
        nxt = spans[idx + 1] if strand == 1 else spans[idx - 1] if idx else None
        if need <= 0 or nxt is None:
            return None
        extra = (
            [nxt[0] + i for i in range(need)] if strand == 1
            else [nxt[1] - i for i in range(need)]
        )
        return sorted(inside + extra)
    return None


def _pick_transcript_consequence(
    consequences: list[dict],
    supplied_transcript: str | None,
    gene: str | None = None,
) -> dict:
    """Pick the most relevant transcript_consequence from VEP's list.

    Preference order:
      1. The transcript matching ``supplied_transcript`` (compare both
         versioned and base accession so NM_000257.4 matches NM_000257.3
         when only the base accession is supplied or vice versa).
      2. The transcript flagged as MANE Select.
      3. The first canonical Ensembl transcript.
      4. The first entry in the list.

    ``gene`` scopes steps 2-4 to the gene the CALLER ASKED FOR. VEP returns a
    block for every transcript overlapping the locus, neighbouring genes
    included, and each gene carries its own MANE Select flag — so the MANE loop
    would otherwise return whichever gene's MANE block happens to be ordered
    first. Measured on production 2026-09-07: NRAS c.34G>A resolved against
    CSDE1's MANE Select (chr1:114716127 is downstream of CSDE1), which has no
    hgvsp / hgvsc / protein_start because the position is outside its CDS. PM5
    and PS1 both read the residue out of hgvsp, so both went "not evaluated"
    and the variant fell from LP to VUS — silently, because
    ``most_severe_consequence`` is computed across ALL blocks and still said
    missense_variant.

    The filter FAILS OPEN: when no block names ``gene`` the full list is used,
    so a symbol alias can never leave the caller with nothing. Step 1 is
    deliberately outside the filter — an explicit accession is the strongest
    signal a caller can give and still wins.
    """
    if not consequences:
        return {}
    if supplied_transcript:
        target_base = _strip_version(supplied_transcript)
        for tc in consequences:
            tx_id = tc.get("transcript_id")
            if tx_id == supplied_transcript:
                return tc
            if _strip_version(tx_id) == target_base:
                return tc
            mane = tc.get("mane_select") or tc.get("mane_plus_clinical")
            if mane and (_strip_version(mane) == target_base or mane == supplied_transcript):
                return tc
    candidates = consequences
    if gene:
        want = gene.strip().upper()
        same_gene = [
            tc for tc in consequences
            if (tc.get("gene_symbol") or "").strip().upper() == want
        ]
        if same_gene:
            candidates = same_gene
    for tc in candidates:
        if tc.get("mane_select"):
            return tc
    for tc in candidates:
        if tc.get("canonical"):
            return tc
    return candidates[0]


_IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
_HIGH_TERMS = frozenset({
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification",
})
_MODERATE_TERMS = frozenset({
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant",
})


def _consequence_impact_rank(tc: dict) -> int:
    """Severity rank (3=HIGH … 0=MODIFIER) for one transcript consequence.
    Prefers VEP's own ``impact``; falls back to the consequence-term map when
    impact is absent so a missense-vs-synonymous split still registers."""
    impact = (tc.get("impact") or "").upper()
    if impact in _IMPACT_RANK:
        return _IMPACT_RANK[impact]
    terms = set(tc.get("consequence_terms") or [])
    if terms & _HIGH_TERMS:
        return 3
    if terms & _MODERATE_TERMS:
        return 2
    return 1 if terms else 0


def _build_transcript_table(
    consequences: list[dict], picked: dict
) -> tuple[list[dict], bool, str | None]:
    """Build the per-transcript consequence table for the UI from VEP's full
    ``transcript_consequences`` list, and decide whether the predicted
    consequence differs SIGNIFICANTLY across transcripts.

    Returns ``(rows, differs_significantly, transcript_set)``:
      * ``rows`` — one dict per transcript with the fields the frontend table
        renders (id, source, biotype, consequence, HGVSc/HGVSp, impact, the
        MANE flags + accessions, canonical, and whether it is the scored
        ("picked") transcript), sorted MANE Select → MANE Plus Clinical →
        canonical → protein_coding → impact severity → accession.
      * ``differs_significantly`` — True when the protein-coding transcripts
        span the HIGH boundary (LoF on some isoforms but not others) or a
        ≥2-tier impact gap (e.g. MODERATE missense vs MODIFIER intronic). A
        single transcript / uniform impact is never flagged.
      * ``transcript_set`` — 'RefSeq' / 'GENCODE' / 'mixed' / None: which
        transcript database the rows came from. Ensembl REST returns one set
        per query (RefSeq only when a RefSeq accession was supplied), so the UI
        can note the table is scoped to that set.
    """
    rows: list[dict] = []
    picked_id = picked.get("transcript_id") if picked else None
    for tc in consequences or []:
        tx_id = tc.get("transcript_id")
        if not tx_id:
            continue
        mane_select_acc = tc.get("mane_select")
        mane_clinical_acc = tc.get("mane_plus_clinical")
        source = (
            "RefSeq" if _is_refseq(tx_id)
            else "GENCODE" if str(tx_id).upper().startswith("ENST")
            else ""
        )
        terms = tc.get("consequence_terms") or []
        rows.append({
            "transcript_id": tx_id,
            "source": source,
            "biotype": tc.get("biotype"),
            "consequence_terms": terms,
            "consequence": terms[0] if terms else None,
            "hgvsc": tc.get("hgvsc"),
            "hgvsp": tc.get("hgvsp"),
            "impact": tc.get("impact"),
            "is_mane_select": bool(mane_select_acc),
            "is_mane_plus_clinical": bool(mane_clinical_acc),
            "mane_select_accession": mane_select_acc,
            "mane_plus_clinical_accession": mane_clinical_acc,
            "canonical": bool(tc.get("canonical")),
            "is_picked": bool(tx_id == picked_id),
            "_rank": _consequence_impact_rank(tc),
        })

    rows.sort(key=lambda r: (
        0 if r["is_mane_select"] else 1,
        0 if r["is_mane_plus_clinical"] else 1,
        0 if r["canonical"] else 1,
        0 if r.get("biotype") == "protein_coding" else 1,
        -r["_rank"],
        str(r["transcript_id"]),
    ))

    coding = [r for r in rows if r.get("biotype") == "protein_coding"] or rows
    differs = False
    if len(coding) > 1:
        ranks = {r["_rank"] for r in coding}
        if len(ranks) > 1:
            differs = (3 in ranks and min(ranks) < 3) or (max(ranks) - min(ranks) >= 2)

    sources = {r["source"] for r in rows if r["source"]}
    transcript_set = (
        "RefSeq" if sources == {"RefSeq"}
        else "GENCODE" if sources == {"GENCODE"}
        else "mixed" if sources
        else None
    )

    for r in rows:
        r.pop("_rank", None)
    return rows, differs, transcript_set


async def _lookup_mane_status(
    client: httpx.AsyncClient, transcript_id: str | None
) -> tuple[bool, bool]:
    """Fallback MANE check via Ensembl /lookup/id when the VEP response
    omits the MANE flags (older NM_ versions sometimes lack them).

    Returns (is_mane_select, is_mane_clinical). Failures are swallowed —
    both flags default to False.
    """
    if not transcript_id:
        return False, False
    try:
        r = await client.get(
            f"{ENSEMBL_BASE}/lookup/id/{transcript_id}",
            headers=HEADERS,
            timeout=_with_connect_cap(10.0),
        )
        if r.status_code != 200:
            return False, False
        data = r.json() or {}
        return bool(data.get("is_mane_select")), bool(data.get("is_mane_plus_clinical"))
    except Exception:
        return False, False


_EXON_LENGTHS_CACHE: dict[tuple[str, str], list[int] | None] = {}


async def _fetch_exon_lengths(
    client: httpx.AsyncClient, transcript_id: str | None, base_url: str
) -> list[int] | None:
    """Return the transcript's exon lengths in transcript (rank, 5'→3')
    order, or ``None`` when the lookup fails or is empty.

    Ensembl's ``/lookup/id?expand=1`` returns the ``Exon`` array; we sort it
    into rank order defensively (ascending genomic start on the + strand,
    descending on the - strand) rather than trusting the response order, so
    the cumulative-length arithmetic is correct for reverse-strand genes
    (MYH7, RAF1, …). Failures are swallowed — the caller degrades to the
    exon-number-only NMD rules."""
    if not transcript_id:
        return None
    key = (base_url, transcript_id)
    if key in _EXON_LENGTHS_CACHE:
        return _EXON_LENGTHS_CACHE[key]
    lengths: list[int] | None = None
    try:
        r = await client.get(
            f"{base_url}/lookup/id/{transcript_id}",
            params={"expand": "1"},
            headers=HEADERS,
            timeout=_with_connect_cap(25.0),
        )
        if r.status_code == 200:
            data = r.json() or {}
            exons = data.get("Exon") or []
            if exons:
                strand = data.get("strand")
                ordered = sorted(
                    exons, key=lambda e: e["start"], reverse=(strand == -1)
                )
                lengths = [e["end"] - e["start"] + 1 for e in ordered]
    except Exception:
        lengths = None
    _EXON_LENGTHS_CACHE[key] = lengths
    return lengths


_TRANSCRIPT_EXONS_CACHE: dict[str, dict] = {}


async def fetch_transcript_exons(transcript_id: str | None) -> dict | None:
    """Exon genomic coordinates for the gene-model track drawn under the
    ClinVar-landscape lollipop (so the curator can see which exon each variant
    — and the proband — falls in, like ClinVar's own viewer).

    Returns ``{"ok": True, "transcript_id": str, "strand": ±1, "exons":
    [{"start", "end", "rank"}, ...]}`` with ``rank`` in 5'→3' transcription
    order, or ``None`` when the lookup fails / is empty (the frontend then draws
    the lollipop without the ladder). One memoised ``/lookup/id?expand=1`` call
    per transcript.

    ``transcript_id`` is echoed back so the chart can caption itself with the
    model it was ACTUALLY built from. The caller passes the picked transcript
    from the VEP block, and the frontend could read that same field — but then
    the caption would be an assumption rather than a fact, and would keep
    claiming the right transcript if the two ever diverged. Exon ranks and
    residue numbers only mean anything relative to one model, so the label has
    to come from the same place the coordinates did."""
    if not transcript_id:
        return None
    if transcript_id in _TRANSCRIPT_EXONS_CACHE:
        return _TRANSCRIPT_EXONS_CACHE[transcript_id]

    from . import vep_offline

    if vep_offline.offline_enabled():
        db = os.environ.get("HEARTVAR_CDOT_DB", "").strip()
        if db:
            from . import cdot_sqlite

            local = cdot_sqlite.transcript_exons(db, transcript_id)
            if local is not None:
                _TRANSCRIPT_EXONS_CACHE[transcript_id] = local
                return local

        _TRANSCRIPT_EXONS_CACHE[transcript_id] = None
        return None

    result: dict | None = None
    try:
        async with httpx.AsyncClient() as client:
            r = await client.get(
                f"{ENSEMBL_BASE}/lookup/id/{transcript_id}",
                params={"expand": "1"},
                headers=HEADERS,
                timeout=_with_connect_cap(25.0),
            )
        if r.status_code == 200:
            data = r.json() or {}
            exons = data.get("Exon") or []
            if exons:
                strand = data.get("strand") or 1
                ordered = sorted(
                    exons, key=lambda e: e["start"], reverse=(strand == -1)
                )
                result = {
                    "ok": True,
                    "transcript_id": transcript_id,
                    "strand": strand,
                    "exons": [
                        {"start": e["start"], "end": e["end"], "rank": i + 1}
                        for i, e in enumerate(ordered)
                    ],
                }
    except Exception:
        result = None
    if result is not None:
        _TRANSCRIPT_EXONS_CACHE[transcript_id] = result
    return result


def nmd_escape_needs_exon_lengths(tc: dict) -> bool:
    """True only when the PENULTIMATE-exon rule is the one that decides.

    Exists so the REST path keeps firing at most one ``/lookup/id`` per
    curation, and so the offline path only reads cdot when it must."""
    exon = tc.get("exon")
    if not exon or "/" not in str(exon):
        return False
    try:
        n, total = (int(x) for x in str(exon).split("/", 1))
    except ValueError:
        return False
    return total > 1 and n == total - 1


def nmd_escape_from_lengths(tc: dict, lengths: list[int] | None) -> bool:
    """The ClinGen PVS1 NMD-escape decision, given exon lengths (or None).

    A PTC escapes NMD when it lies in a single-exon transcript, in the last
    (3'-most) exon, or in the last 50 nt of the penultimate exon
    (Abou Tayoun 2018). The first two rules need only VEP's ``exon`` "N/total"
    string; the third needs ``lengths`` in transcript rank order, and returns
    False without them rather than guessing.

    ⚠ ONE COPY, SHARED BY BOTH PATHS. This used to be inlined in
    _predict_nmd_escape and reimplemented in vep_offline._nmd_escape_from_exon,
    and the two drifted: the offline copy omitted the third rule entirely, so
    offline said False where REST said True — and that direction keeps PVS1 at
    FULL STRENGTH, i.e. offline over-called
    pathogenicity for a PTC in the last 50 nt of a penultimate exon. Two copies
    of an ACMG rule is one too many.

    ``cdna_start`` is used as the PTC position — the standard simplification;
    frameshift PTCs may sit a few codons 3' of it.
    """
    exon = tc.get("exon")
    if not exon or "/" not in str(exon):
        return False
    try:
        n, total = (int(x) for x in str(exon).split("/", 1))
    except ValueError:
        return False
    if total <= 1 or n == total:
        return True
    if n != total - 1:
        return False
    cdna_start = tc.get("cdna_start")
    if cdna_start is None or not lengths or len(lengths) < total:
        return False
    penultimate_3p = sum(lengths[: total - 1])
    return cdna_start > penultimate_3p - 50


async def _predict_nmd_escape(
    client: httpx.AsyncClient, tc: dict, base_url: str
) -> bool:
    """REST-side NMD escape: fetch exon lengths only when rule 3 decides it,
    then apply the shared rule. See :func:`nmd_escape_from_lengths`."""
    if not nmd_escape_needs_exon_lengths(tc):
        return nmd_escape_from_lengths(tc, None)
    lengths = await _fetch_exon_lengths(client, tc.get("transcript_id"), base_url)
    return nmd_escape_from_lengths(tc, lengths)


async def _splice_exon_lengths(
    client: httpx.AsyncClient, tc: dict, base_url: str
) -> list[int] | None:
    """Exon-length list (rank order) surfaced for canonical ±1,2 splice
    variants only, so the PVS1 splice arm in ``compute_hard_coded_criteria``
    can score the affected exon's frame and fraction-of-transcript removed.

    Returns ``None`` for non-splice variants and on any lookup failure — the
    splice arm then falls back to base PVS1. Reuses the cached exon-structure
    fetch, so a canonical-splice variant whose transcript was already seen
    costs nothing."""
    conseq = tc.get("consequence_terms") or []
    if not any(
        c in ("splice_acceptor_variant", "splice_donor_variant") for c in conseq
    ):
        return None
    return await _fetch_exon_lengths(client, tc.get("transcript_id"), base_url)


_BP7_PHYLOP_BACKFILL_TOKENS = (
    "synonymous", "intron_variant", "splice_region", "splice_polypyrimidine",
    "3_prime_utr_variant", "5_prime_utr_variant", "non_coding_transcript",
)


def _needs_phylop_backfill(consequence: str) -> bool:
    c = (consequence or "").lower()
    return any(t in c for t in _BP7_PHYLOP_BACKFILL_TOKENS)


async def _enrich_phylop_conservation(result: dict) -> dict:
    """Backfill ``result["phylop100way"]`` from the conservation bigWig for
    BP7-eligible synonymous / non-coding variants where dbNSFP returned nothing.
    No-op (and never raises) for everything else; the value is then cached with
    the VEP result so BP7 reads it transparently via ``vep.get("phylop100way")``."""
    if not isinstance(result, dict) or not result.get("ok"):
        return result
    if result.get("phylop100way") is not None:
        return result
    if not _needs_phylop_backfill(result.get("most_severe_consequence") or ""):
        return result
    parts = (result.get("allele_string") or "").split("/")
    if not (len(parts) == 2 and all(len(p) == 1 and p != "-" for p in parts)):
        return result
    val = await fetch_phylop100way(result.get("seq_region_name"), result.get("start"))
    if val is not None:
        result["phylop100way"] = val
    return result


async def _fetch_vep_uncached_enriched(
    gene: str, hgvs_c: str, supplied_transcript: str | None,
) -> dict:
    return await _enrich_phylop_conservation(
        await _fetch_vep_uncached(gene, hgvs_c, supplied_transcript))


async def _fetch_vep_by_coordinates_uncached_enriched(
    chrom: str, pos: int, ref: str, alt: str, build: str,
    supplied_transcript: str | None,
) -> dict:
    return await _enrich_phylop_conservation(
        await _fetch_vep_by_coordinates_uncached(
            chrom, pos, ref, alt, build, supplied_transcript))


async def fetch_vep(
    gene: str,
    hgvs_c: str,
    supplied_transcript: str | None = None,
) -> dict:
    """Cached wrapper around :func:`_fetch_vep_uncached`.

    VEP is the first call on every curation and the single point of failure
    (gnomAD / SpliceAI ids, the literature AA token, and PVS1 exon structure
    all chain off it), so its result is cached for ``TTL_VEP`` keyed on
    (gene, HGVSc, supplied transcript). Single-flight coalescing collapses
    concurrent identical lookups to one upstream call. Only successful results
    are cached — a transient failure re-resolves on the next attempt.
    """
    key = ("vep", (gene or "").upper(), hgvs_c or "", supplied_transcript or "")
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _fetch_vep_uncached_enriched(gene, hgvs_c, supplied_transcript),
        ttl=TTL_VEP,
        should_cache=_ok,
    )


async def _fetch_vep_uncached(
    gene: str,
    hgvs_c: str,
    supplied_transcript: str | None = None,
) -> dict:
    """Resolve gene+HGVSc to genomic coords + consequence via Ensembl VEP.

    When ``supplied_transcript`` is provided (the user's NM_ / ENST prefix
    detached upstream by ``strip_transcript_prefix``), VEP is queried
    against that transcript first so the per-transcript annotation matches
    the curator's intended reading frame. Otherwise VEP picks its canonical.

    Returns a dict with at minimum {"queried_as": ..., "ok": bool}. On success,
    includes assembly, chrom, pos, ref, alt, consequence, transcript fields,
    and the MANE-tracking fields selected_transcript_id /
    selected_transcript_source / is_mane_select / is_mane_clinical.
    On failure, includes "error".
    """
    from . import vep_offline
    if vep_offline.offline_enabled():
        offline = await vep_offline.fetch_hgvs(gene, hgvs_c, supplied_transcript)
        if offline is not None:
            return offline

    async with httpx.AsyncClient() as client:
        attempts: list[tuple[str, str, dict]] = []
        if supplied_transcript:
            if _is_refseq(supplied_transcript):
                base = _strip_version(supplied_transcript)
                attempts.append((f"{base}:{hgvs_c}", "user-supplied", _REFSEQ_PARAMS))
                if base != supplied_transcript:
                    attempts.append(
                        (f"{supplied_transcript}:{hgvs_c}",
                         "user-supplied (versioned)", _REFSEQ_PARAMS))
            else:
                attempts.append((f"{supplied_transcript}:{hgvs_c}", "user-supplied", VEP_PARAMS))
        attempts.append((f"{gene}:{hgvs_c}", "VEP canonical", VEP_PARAMS))
        canonical_transcript = await _lookup_canonical_transcript(client, gene)
        if canonical_transcript:
            attempts.append((f"{canonical_transcript}:{hgvs_c}", "VEP canonical", VEP_PARAMS))

        last_error: str | None = None
        last_kind: str | None = None
        saw_transient = False
        for hgvs, source_label, params in attempts:
            r = await _ensembl_get(
                client,
                f"{ENSEMBL_BASE}/vep/human/hgvs/{hgvs}",
                params,
                timeout=30.0,
            )
            if r is None:
                last_error = "transient Ensembl failure (429/5xx/timeout) after retries"
                last_kind = "transient"
                saw_transient = True
                continue
            if r.status_code != 200:
                last_error = f"{r.status_code}: {r.text[:200]}"
                if 400 <= r.status_code < 500:
                    last_kind = "not_found"
                else:
                    last_kind = "transient"
                    saw_transient = True
                continue
            payload = r.json()
            if not payload:
                last_error = "empty response"
                last_kind = "empty"
                continue
            entry = payload[0]
            most_severe = entry.get("most_severe_consequence")
            consequences = entry.get("transcript_consequences") or []
            tc = _pick_transcript_consequence(consequences, supplied_transcript, gene)

            mane_select_acc = tc.get("mane_select")
            mane_clinical_acc = tc.get("mane_plus_clinical")
            is_mane_select = bool(mane_select_acc)
            is_mane_clinical = bool(mane_clinical_acc)
            if not is_mane_select and not is_mane_clinical:
                fb_select, fb_clinical = await _lookup_mane_status(
                    client, tc.get("transcript_id")
                )
                is_mane_select = fb_select
                is_mane_clinical = fb_clinical

            nmd_escape = await _predict_nmd_escape(client, tc, ENSEMBL_BASE)
            exon_lengths = await _splice_exon_lengths(client, tc, ENSEMBL_BASE)
            transcript_table, conseq_differs, transcript_set = _build_transcript_table(
                consequences, tc
            )
            requested_honored = None
            if supplied_transcript:
                _sb = _strip_version(supplied_transcript)
                _picked_mane = tc.get("mane_select") or tc.get("mane_plus_clinical")
                requested_honored = bool(
                    _strip_version(tc.get("transcript_id")) == _sb
                    or (_picked_mane and _strip_version(_picked_mane) == _sb)
                )
            return {
                "ok": True,
                "queried_as": hgvs,
                "assembly_name": entry.get("assembly_name"),
                "seq_region_name": entry.get("seq_region_name"),
                "start": entry.get("start"),
                "end": entry.get("end"),
                "allele_string": entry.get("allele_string"),
                "strand": entry.get("strand"),
                "transcript_strand": tc.get("strand"),
                "vcf_string": entry.get("vcf_string"),
                "most_severe_consequence": most_severe,
                "transcript_id": tc.get("transcript_id"),
                "gene_id": tc.get("gene_id"),
                "gene_symbol": tc.get("gene_symbol"),
                "biotype": tc.get("biotype"),
                "hgvsc": tc.get("hgvsc"),
                "hgvsp": tc.get("hgvsp"),
                "consequence_terms": tc.get("consequence_terms"),
                "exon": tc.get("exon"),
                "intron": tc.get("intron"),
                "exon_lengths": exon_lengths,
                "cdna_start": tc.get("cdna_start"),
                "cdna_end": tc.get("cdna_end"),
                "cds_start": tc.get("cds_start"),
                "cds_end": tc.get("cds_end"),
                "protein_start": tc.get("protein_start"),
                "protein_end": tc.get("protein_end"),
                "nmd_escape": nmd_escape,
                "impact": tc.get("impact"),
                "sift_prediction": tc.get("sift_prediction"),
                "sift_score": tc.get("sift_score"),
                "polyphen_prediction": tc.get("polyphen_prediction"),
                "polyphen_score": tc.get("polyphen_score"),
                "cadd_phred": tc.get("cadd_phred"),
                "cadd_raw": tc.get("cadd_raw"),
                "revel_score": _max_dbnsfp(tc.get("revel_score")),
                "phylop100way": _max_dbnsfp(tc.get("phylop100way_vertebrate")),
                "dbsnp_rsid": _dbsnp_rsid_from_entry(entry),
                "selected_transcript_id": (
                    supplied_transcript
                    if (
                        supplied_transcript
                        and source_label == "user-supplied"
                    )
                    else tc.get("transcript_id")
                ),
                "selected_transcript_source": source_label,
                "is_mane_select": is_mane_select,
                "is_mane_clinical": is_mane_clinical,
                "mane_select_accession": mane_select_acc,
                "mane_clinical_accession": mane_clinical_acc,
                "requested_transcript": supplied_transcript,
                "requested_transcript_honored": requested_honored,
                "transcript_consequences_all": transcript_table,
                "consequence_differs_significantly": conseq_differs,
                "transcript_set": transcript_set,
            }

        return {
            "ok": False,
            "error": last_error or "VEP lookup failed",
            "failure_kind": "transient" if saw_transient else (last_kind or "not_found"),
        }


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _complement(allele: str) -> str:
    """Reverse-complement an allele string. Single nucleotides and short
    indel sequences both flow through `translate`; the order doesn't matter
    for SNV alleles, and indel orientation flips correctly because the
    sequence is read 5'→3' on the opposite strand."""
    return allele.translate(_COMPLEMENT)[::-1] if len(allele) > 1 else allele.translate(_COMPLEMENT)


async def liftover_grch37_to_grch38(
    chrom: str, pos: int, client: httpx.AsyncClient | None = None
) -> tuple[str, int] | None:
    """Lift a single GRCh37 coordinate to GRCh38 via Ensembl's REST
    ``/map`` endpoint.

    Ensembl exposes a coordinate-mapping service at
    ``https://rest.ensembl.org/map/human/GRCh37/{chrom}:{start}..{end}:1/GRCh38``.
    It returns a JSON object with ``mappings`` listing the lifted
    region(s). Returns ``(chrom, pos)`` on success, ``None`` on any
    failure (network, empty mapping, multi-mapping ambiguity).

    (An earlier draft hit a UCSC REST endpoint; that URL now redirects
    to a CGI form. The Ensembl service is the canonical REST liftover
    surface and is already in our dependency surface.)

    A pre-opened ``httpx.AsyncClient`` may be passed in so the lookup
    rides the caller's connection pool; otherwise a temporary client is
    spun up.
    """
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient()
    try:
        bare_chrom = chrom.removeprefix("chr") if chrom.lower().startswith("chr") else chrom
        url = (
            f"{ENSEMBL_BASE}/map/human/GRCh37/"
            f"{bare_chrom}:{pos}..{pos}:1/GRCh38"
        )
        r = await client.get(url, headers=HEADERS, timeout=_with_connect_cap(25.0))
        if r.status_code != 200:
            log.warning(
                "Ensembl liftover returned %s for %s:%s — %s",
                r.status_code, bare_chrom, pos, r.text[:160],
            )
            return None
        data = r.json() or {}
        mappings = data.get("mappings") or []
        if len(mappings) != 1:
            log.warning(
                "Ensembl liftover returned %d mappings for %s:%s",
                len(mappings), bare_chrom, pos,
            )
            return None
        mapped = mappings[0].get("mapped") or {}
        new_chrom = mapped.get("seq_region_name")
        new_start = mapped.get("start")
        if not new_chrom or new_start is None:
            return None
        return str(new_chrom), int(new_start)
    except Exception as exc:
        log.warning("Ensembl liftover failed for %s:%s — %r", chrom, pos, exc)
        return None
    finally:
        if own_client:
            await client.aclose()


def _vcf_input_for_region(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build the VCF-like input row that Ensembl VEP's /vep/human/region
    POST endpoint expects: 'chrom pos . ref alt . . .'. VEP REST docs
    state the trailing three fields (QUAL/FILTER/INFO) are optional but
    must be present as placeholders when supplied — '.' satisfies the
    parser without imposing a quality threshold."""
    return f"{chrom} {pos} . {ref} {alt} . . ."


async def fetch_vep_by_coordinates(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    build: str = "GRCh38",
    supplied_transcript: str | None = None,
) -> dict:
    """Cached wrapper around :func:`_fetch_vep_by_coordinates_uncached` — same
    TTL/single-flight policy as :func:`fetch_vep`, keyed on the input
    coordinates and build."""
    key = ("vep_coord", chrom, pos, ref, alt, build, supplied_transcript or "")
    return await EXTERNAL_CACHE.get_or_set(
        key,
        lambda: _fetch_vep_by_coordinates_uncached_enriched(
            chrom, pos, ref, alt, build, supplied_transcript),
        ttl=TTL_VEP,
        should_cache=_ok,
    )


async def _fetch_vep_by_coordinates_uncached(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    build: str = "GRCh38",
    supplied_transcript: str | None = None,
) -> dict:
    """Resolve a chr/pos/ref/alt variant through Ensembl VEP's
    ``/vep/human/region`` POST endpoint.

    For ``build == "GRCh37"`` the call is routed at
    ``grch37.rest.ensembl.org`` so VEP's consequence prediction lines up
    with the curator's input build. Downstream services (gnomAD,
    SpliceAI, ClinVar) work in GRCh38 only — the caller is expected to
    have lifted the GRCh37 input to GRCh38 before any of those run.

    Returns the same shape as ``fetch_vep`` with the extra fields:
      - ``input_format``      → "coordinates"
      - ``input_build``       → "GRCh38" / "GRCh37"
      - ``input_coords``      → "chrom-pos-ref-alt" string of the
                                INPUT coordinates (build-as-supplied).
      - ``derived_hgvs``      → the gene-symbol-prefixed HGVS c. notation
                                pulled from VEP for display + downstream
                                HGVS-based lookups.
      - ``derived_gene_symbol``→ gene symbol pulled from VEP's
                                transcript_consequences (the curator may
                                not have supplied one for coord input).

    On failure returns ``{"ok": False, "error": ...}``.
    """
    from . import vep_offline
    if vep_offline.offline_enabled():
        offline = await vep_offline.fetch_region(
            chrom, pos, ref, alt, build, supplied_transcript
        )
        if offline is not None:
            return offline

    base = ENSEMBL_GRCH37_BASE if build == "GRCh37" else ENSEMBL_BASE
    vcf_line = _vcf_input_for_region(chrom, pos, ref, alt)
    input_coords = f"{chrom}-{pos}-{ref}-{alt}"
    async with httpx.AsyncClient() as client:
        try:
            r = await client.post(
                f"{base}/vep/human/region",
                headers={**HEADERS, "Content-Type": "application/json"},
                params=VEP_PARAMS,
                json={"variants": [vcf_line]},
                timeout=_with_connect_cap(30.0),
            )
        except Exception as exc:
            log.warning("VEP region POST failed for %s: %r", input_coords, exc)
            return {"ok": False, "error": "VEP region POST failed",
                    "failure_kind": "transient"}
        if r.status_code != 200:
            return {
                "ok": False,
                "error": f"{r.status_code}: {r.text[:200]}",
                "failure_kind": "not_found" if 400 <= r.status_code < 500 else "transient",
                "input_format": "coordinates",
                "input_build": build,
                "input_coords": input_coords,
            }
        payload = r.json()
        if not payload:
            return {
                "ok": False,
                "error": "empty VEP region response",
                "failure_kind": "empty",
                "input_format": "coordinates",
                "input_build": build,
                "input_coords": input_coords,
            }
        entry = payload[0]
        consequences = entry.get("transcript_consequences") or []
        tc = _pick_transcript_consequence(consequences, supplied_transcript)

        mane_select_acc = tc.get("mane_select")
        mane_clinical_acc = tc.get("mane_plus_clinical")
        is_mane_select = bool(mane_select_acc)
        is_mane_clinical = bool(mane_clinical_acc)
        if not is_mane_select and not is_mane_clinical:
            fb_select, fb_clinical = await _lookup_mane_status(
                client, tc.get("transcript_id")
            )
            is_mane_select = fb_select
            is_mane_clinical = fb_clinical

        gene_symbol = tc.get("gene_symbol")
        hgvsc = tc.get("hgvsc") or ""
        derived_hgvs = hgvsc.split(":")[-1] if hgvsc else None

        nmd_escape = await _predict_nmd_escape(client, tc, base)
        exon_lengths = await _splice_exon_lengths(client, tc, base)
        transcript_table, conseq_differs, transcript_set = _build_transcript_table(
            consequences, tc
        )
        return {
            "ok": True,
            "queried_as": vcf_line,
            "assembly_name": entry.get("assembly_name"),
            "seq_region_name": entry.get("seq_region_name"),
            "start": entry.get("start"),
            "end": entry.get("end"),
            "allele_string": entry.get("allele_string"),
            "strand": entry.get("strand"),
            "transcript_strand": tc.get("strand"),
            "vcf_string": entry.get("vcf_string"),
            "most_severe_consequence": entry.get("most_severe_consequence"),
            "transcript_id": tc.get("transcript_id"),
            "gene_id": tc.get("gene_id"),
            "gene_symbol": gene_symbol,
            "biotype": tc.get("biotype"),
            "hgvsc": tc.get("hgvsc"),
            "hgvsp": tc.get("hgvsp"),
            "consequence_terms": tc.get("consequence_terms"),
            "exon": tc.get("exon"),
            "intron": tc.get("intron"),
            "exon_lengths": exon_lengths,
            "cdna_start": tc.get("cdna_start"),
            "cdna_end": tc.get("cdna_end"),
            "cds_start": tc.get("cds_start"),
            "cds_end": tc.get("cds_end"),
            "protein_start": tc.get("protein_start"),
            "protein_end": tc.get("protein_end"),
            "nmd_escape": nmd_escape,
            "impact": tc.get("impact"),
            "sift_prediction": tc.get("sift_prediction"),
            "sift_score": tc.get("sift_score"),
            "polyphen_prediction": tc.get("polyphen_prediction"),
            "polyphen_score": tc.get("polyphen_score"),
            "cadd_phred": tc.get("cadd_phred"),
            "cadd_raw": tc.get("cadd_raw"),
            "revel_score": _max_dbnsfp(tc.get("revel_score")),
            "phylop100way": _max_dbnsfp(tc.get("phylop100way_vertebrate")),
            "selected_transcript_id": (
                supplied_transcript
                if supplied_transcript else tc.get("transcript_id")
            ),
            "selected_transcript_source": (
                "user-supplied" if supplied_transcript else "VEP canonical"
            ),
            "is_mane_select": is_mane_select,
            "is_mane_clinical": is_mane_clinical,
            "mane_select_accession": mane_select_acc,
            "mane_clinical_accession": mane_clinical_acc,
            "input_format": "coordinates",
            "input_build": build,
            "input_coords": input_coords,
            "derived_hgvs": derived_hgvs,
            "derived_gene_symbol": gene_symbol,
            "transcript_consequences_all": transcript_table,
            "consequence_differs_significantly": conseq_differs,
            "transcript_set": transcript_set,
            "forward_variant_id": f"{chrom}-{pos}-{ref}-{alt}",
        }


def _id_from_vcf_string(vcf_string) -> str | None:
    """Turn VEP's ``vcf_string`` into a forward-strand "chrom-pos-ref-alt" id.

    REST returns it as a LIST (e.g. ``["15-48425438-TG-T"]``); the recoder
    path can hand it over as a plain string. It is already left-aligned and on
    the reference forward strand — the same normalization gnomAD v4.1 applies —
    so it is returned verbatim (leading ``chr`` stripped). Returns None unless
    it is a well-formed 4-field record with ACGT-only ref/alt, so a malformed
    or symbolic representation falls through to the allele_string path.
    """
    vs = vcf_string
    if isinstance(vs, (list, tuple)):
        vs = vs[0] if vs else None
    if not isinstance(vs, str) or not vs.strip():
        return None
    vs = vs.strip()
    if vs[:3].lower() == "chr":
        vs = vs[3:]
    parts = vs.split("-")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    if not (chrom and pos.isdigit() and ref and alt):
        return None
    if not (re.fullmatch(r"[ACGTacgt]+", ref) and re.fullmatch(r"[ACGTacgt]+", alt)):
        return None
    return f"{chrom}-{pos}-{ref}-{alt}"


def gnomad_variant_id_with_provenance(vep: dict) -> tuple[str | None, bool]:
    """Build a forward-strand "chrom-pos-ref-alt" variant ID from a VEP result,
    and signal whether it came from VEP's CANONICAL left-aligned ``vcf_string``.

    Order of preference:
      1. ``vcf_string`` — VEP's own left-aligned VCF representation. This is the
         only form that produces a gnomAD-matching anchored key for indels/dups
         (gnomAD stores anchored left-aligned keys; the hand-built
         allele_string/_complement key for an indel is non-anchored garbage and
         can never match). Returns ``canonical=True``.
      2. ``forward_variant_id`` — precomputed for coord input (the VCF-style
         /vep/human/region echoes its input on the forward strand). Only reached
         when vcf_string is absent. ``canonical=False``.
      3. ``allele_string`` (+ complement for reverse-strand genes like MYH7,
         RAF1) — Ensembl REST returns allele_string on the TRANSCRIPT strand for
         HGVS-c input, so reverse-strand alleles are complemented to the forward
         strand. Correct for SNVs; non-anchored for indels. ``canonical=False``.

    ``canonical`` MUST be checked by any caller that gates a "confirmed absent"
    decision on a gnomAD lookup miss — a miss on a non-canonical indel key is
    meaningless (the key never could have matched).
    """
    if not vep.get("ok"):
        return None, False
    vid = _id_from_vcf_string(vep.get("vcf_string"))
    if vid:
        return vid, True
    forward = vep.get("forward_variant_id")
    if forward:
        return forward, False
    chrom = vep.get("seq_region_name")
    pos = vep.get("start")
    allele = vep.get("allele_string") or ""
    if not (chrom and pos and "/" in allele):
        return None, False
    ref, alt = allele.split("/", 1)
    if vep.get("strand") == -1:
        ref, alt = _complement(ref), _complement(alt)
    return f"{chrom}-{pos}-{ref}-{alt}", False


def gnomad_variant_id(vep: dict) -> str | None:
    """Back-compat wrapper returning only the id (drops the canonical flag).
    Callers that gate 'confirmed absent' on a lookup miss must instead use
    ``gnomad_variant_id_with_provenance`` and honour the canonical flag."""
    return gnomad_variant_id_with_provenance(vep)[0]


async def recode_rsid_to_coords(rsid: str) -> dict | None:
    """Resolve a dbSNP rsID to GRCh38 coordinates via Ensembl ``variant_recoder``.

    Returns ``{"chrom","pos","ref","alt"}`` on the PRIMARY assembly, or ``None``
    when the rsID is unknown / multi-mapping / non-SNV-on-alt-contig. The caller
    runs the coordinate pipeline on the result. Only a representation whose SPDI
    is on a RefSeq ``NC_0000…`` chromosome is trusted (alt-contig mappings are
    skipped, mirroring fetch_variant_recoder_rsid). Never raises."""
    if not rsid:
        return None
    quoted = urllib.parse.quote(rsid, safe="")
    url = f"{ENSEMBL_BASE}/variant_recoder/human/{quoted}"
    try:
        async with httpx.AsyncClient() as client:
            r = await _ensembl_get(client, url, {"fields": "spdi,vcf_string"}, timeout=25.0)
    except httpx.HTTPError:
        return None
    if r is None or r.status_code != 200:
        return None
    try:
        data = r.json()
    except ValueError:
        return None
    for entry in (data if isinstance(data, list) else [data]):
        if not isinstance(entry, dict):
            continue
        for info in entry.values():
            if not isinstance(info, dict):
                continue
            primary = any(
                isinstance(s, str) and s.startswith("NC_0000")
                for s in (info.get("spdi") or [])
            )
            if not primary:
                continue
            for vs in (info.get("vcf_string") or []):
                if not isinstance(vs, str):
                    continue
                parts = vs.split("-")
                if len(parts) != 4:
                    continue
                chrom, pos_s, ref, alt = parts
                ref, alt = ref.upper(), alt.upper()
                if not (set(ref) <= set("ACGT") and set(alt) <= set("ACGT")):
                    continue
                if chrom.upper() in ("M", "MT"):
                    chrom = "MT"
                try:
                    pos = int(pos_s)
                except ValueError:
                    continue
                return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
    return None


async def fetch_variant_recoder_rsid(hgvs: str) -> dict:
    """Resolve a c.HGVS to its dbSNP rsID (+ SPDI) via Ensembl ``variant_recoder``.

    Used as the gnomAD-frequency fallback for INDELS / dups, where a hand-built
    ``chrom-pos-ref-alt`` id is unreliable: gnomAD keys on the left-aligned,
    anchored VCF representation, and for indels in homopolymers / repeats the
    correct left-alignment cannot be derived without the reference sequence (and
    a wrong id silently matches the WRONG variant). gnomAD can instead be queried
    by rsID, and ``variant_recoder`` returns the dbSNP rsID for known (i.e. the
    common, BA1/BS1-relevant) variants — which is exactly the set we were missing.

    Returns ``{"ok": bool, "rsid": str|None, "spdi": str|None}``. Only trusts a
    SPDI/rsID on the PRIMARY assembly (RefSeq ``NC_0000…`` chromosome) so alt-
    contig mappings don't leak in. Never raises (frequency fallback must not
    abort the gather)."""
    if not hgvs:
        return {"ok": False, "rsid": None, "spdi": None}
    quoted = urllib.parse.quote(hgvs, safe="")
    url = f"{ENSEMBL_BASE}/variant_recoder/human/{quoted}"
    try:
        async with httpx.AsyncClient() as client:
            r = await _ensembl_get(client, url, {"fields": "id,spdi"}, timeout=30.0)
    except httpx.HTTPError:
        return {"ok": False, "rsid": None, "spdi": None}
    if r is None or r.status_code != 200:
        return {"ok": False, "rsid": None, "spdi": None}
    try:
        data = r.json()
    except ValueError:
        return {"ok": False, "rsid": None, "spdi": None}
    rsid = spdi = None
    for entry in (data if isinstance(data, list) else [data]):
        if not isinstance(entry, dict):
            continue
        for info in entry.values():
            if not isinstance(info, dict):
                continue
            primary = [
                s for s in (info.get("spdi") or [])
                if isinstance(s, str) and s.startswith("NC_0000")
            ]
            if not primary:
                continue
            if spdi is None:
                spdi = primary[0]
            for i in (info.get("id") or []):
                if isinstance(i, str) and i.lower().startswith("rs"):
                    rsid = i
                    break
            if rsid:
                break
        if rsid:
            break
    return {"ok": True, "rsid": rsid, "spdi": spdi}
