"""Unit tests for backend.clients.vep_offline — the offline VEP fast-path.

Fully offline, no pytest, no real ``vep`` binary: the subprocess is faked by
monkeypatching ``asyncio.create_subprocess_exec`` with a stub process, and the
higher-level entry points patch ``vep_offline._run`` directly. Runnable with
``python -m backend.tests.test_vep_offline``.

Coverage:
  * env gating: offline_enabled / _config (disabled, missing binary/dir)
  * _argv flag composition
  * _first tolerant key lookup + _nmd_escape_from_exon
  * _assemble: CLI-cased plugin fields (REVEL_score/CADD_PHRED) map through;
    result shape + annotation_source
  * _run: JSONL parse, non-zero exit → None, timeout → None
  * fetch_hgvs / fetch_region: gating, transcript requirement, build mismatch,
    happy path
"""
from __future__ import annotations

import asyncio
import os
import tempfile

import backend.clients.vep_offline as vo


_FASTA_STANDIN = tempfile.NamedTemporaryFile(  # noqa: SIM115 - module lifetime
    suffix=".fa.gz", delete=False)
_FASTA_STANDIN.close()


class _env:
    """Temporarily set/clear env vars, restoring exact prior state."""
    def __init__(self, **kv):
        if kv.get("HEARTVAR_VEP_OFFLINE") and "HEARTVAR_VEP_FASTA" not in kv:
            kv = {**kv, "HEARTVAR_VEP_FASTA": _FASTA_STANDIN.name}
        self.kv = kv
        self._saved: dict[str, str | None] = {}

    def __enter__(self):
        for k, v in self.kv.items():
            self._saved[k] = os.environ.get(k)
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return self

    def __exit__(self, *exc):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_offline_disabled_by_default():
    with _env(HEARTVAR_VEP_OFFLINE=None):
        assert vo.offline_enabled() is False
        assert vo._config() is None


def test_offline_enabled_truthy_values():
    for val in ("1", "true", "YES", "On"):
        with _env(HEARTVAR_VEP_OFFLINE=val):
            assert vo.offline_enabled() is True


def test_config_none_when_binary_missing():
    with _env(HEARTVAR_VEP_OFFLINE="1",
              HEARTVAR_VEP_BINARY="definitely-not-a-real-binary-xyz",
              HEARTVAR_VEP_DATA="/tmp"):
        assert vo._config() is None


def test_config_none_when_data_dir_missing():
    with _env(HEARTVAR_VEP_OFFLINE="1",
              HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/no/such/dir/really"):
        assert vo._config() is None


def test_config_ok_and_argv_flags():
    with tempfile.NamedTemporaryFile(suffix=".fa.gz") as fa, \
         _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38",
              HEARTVAR_VEP_FASTA=fa.name,
              HEARTVAR_VEP_EXTRA_ARGS="--plugin CADD,x,y"):
        cfg = vo._config()
        assert cfg is not None
        assert cfg["assembly"] == "GRCh38"
        assert cfg["extra_args"] == ["--plugin", "CADD,x,y"]
        argv = vo._argv(cfg, "vcf", refseq=True)
        for flag in ("--offline", "--cache", "--json", "--no_stats",
                     "--mane", "--numbers", "--hgvs", "--merged"):
            assert flag in argv, flag
        assert "--refseq" not in argv, (
            "--refseq asks for homo_sapiens_refseq, which the build does not "
            "install; the RefSeq transcript set now comes from the merged cache"
        )
        assert "--merged" in vo._argv(cfg, "vcf", refseq=False), (
            "--merged must be unconditional — the merged tree replaces the "
            "unsuffixed species dir, so every call needs it"
        )
        assert argv[argv.index("--format") + 1] == "vcf"
        assert argv[argv.index("--dir_cache") + 1] == "/tmp"
        assert argv[-2:] == ["--plugin", "CADD,x,y"]


def test_config_and_argv_wire_revel_plugin():
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp",
              HEARTVAR_VEP_REVEL="/tmp/new_tabbed_revel_grch38.tsv.gz",
              HEARTVAR_VEP_PLUGINS_DIR="/tmp/Plugins",
              HEARTVAR_VEP_EXTRA_ARGS=None):
        cfg = vo._config()
        assert cfg["revel"] == "/tmp/new_tabbed_revel_grch38.tsv.gz"
        assert cfg["plugins_dir"] == "/tmp/Plugins"
        argv = vo._argv(cfg, "vcf", refseq=False)
        assert argv[argv.index("--dir_plugins") + 1] == "/tmp/Plugins"
        assert argv[argv.index("--plugin") + 1] == \
            "REVEL,file=/tmp/new_tabbed_revel_grch38.tsv.gz,no_match=1"


def test_config_requires_a_fasta_because_argv_always_passes_hgvs():
    """⚠ MEASURED IN THE REAL IMAGE — this is not a theoretical precondition.

        vep --offline --cache ... --format vcf --hgvs ... --merged   (no --fasta)
        -> EXIT 2, 0 bytes
        MSG: ERROR: Cannot generate HGVS coordinates (--hgvs and --hgvsg) in
             offline mode without a FASTA file (see --fasta)

    _argv passes --hgvs on EVERY call and cannot stop: hgvsc/hgvsp feed PVS1,
    PM4 and the literature amino-acid token. So a missing or mistyped
    HEARTVAR_VEP_FASTA is not a degraded offline path, it is a dead one — and
    _config used to return a usable config anyway, so every curation paid for a
    doomed subprocess, fell back to REST, and the banner still read "on".

    That is not hypothetical drift: the path in restart-webapp.yml EMBEDS the
    VEP release (.../113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz), so a
    release bump silently invalidates it — the same skew the cache is protected
    against by deriving its release from the binary."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_FASTA=None):
        assert vo._config() is None, (
            "offline VEP must report itself unusable without a FASTA rather "
            "than spawning a vep that exits 2 on every curation"
        )


def test_config_requires_the_fasta_to_actually_exist():
    """A path that is set but wrong fails identically to one that is unset, and
    is the likelier of the two — see the release-embedding note above."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp",
              HEARTVAR_VEP_FASTA="/no/such/reference/really.fa.gz"):
        assert vo._config() is None


def test_revel_plugin_disables_transcript_matching():
    """⚠ MEASURED: without this, PP3/BP4 silently loses REVEL on newer transcripts.

    REVEL.pm's default requires the REVEL row to LIST the transcript being
    annotated *and* the alt amino acid to match:

        if ($_->{transcript_ids} && !$self->{no_match}) {
          foreach my $tr_id (split /;/, $_->{transcript_ids}) {
            return $_->{result} if ($tr_id eq $tr_stable_id && $_->{altaa} eq $tva->peptide);

    REVEL v1.3 predates a lot of current MANE Select accessions, so for those the
    transcript id is simply absent from its rows and NO score comes back. REST
    does not have this problem: it reads REVEL out of dbNSFP, which is not
    transcript-keyed.

    Offline vs REST: TNNT2 variants whose MANE Select is ENST00000656932
    returned revel=None offline where REST returned a score. With no_match=1 the
    offline path returns the same values REST does.

    no_match=1 drops ONLY the transcript-id requirement. The nucleotide match
    (get_matched_variant_alleles: ref, alts, pos, strand) and the alt-amino-acid
    match ($_->{altaa} eq $tva->peptide) both still apply — which is the correct
    semantics anyway, because REVEL scores an amino-acid substitution, not a
    transcript accession."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp",
              HEARTVAR_VEP_REVEL="/tmp/revel.tsv.gz",
              HEARTVAR_VEP_EXTRA_ARGS=None):
        cfg = vo._config()
        argv = vo._argv(cfg, "vcf", refseq=False)
        spec = argv[argv.index("--plugin") + 1]
        assert spec == "REVEL,file=/tmp/revel.tsv.gz,no_match=1", (
            f"REVEL plugin spec must disable transcript matching, got {spec!r}"
        )


def _entry(most_severe, tcs):
    return {"most_severe_consequence": most_severe, "transcript_consequences": tcs}


def test_most_severe_consequence_is_dropped_when_only_a_filtered_transcript_had_it():
    """⚠ THE ONE FIELD THE NARROWING FORGOT.

    _to_rest_transcript_set exists so the merged cache "does not silently widen
    the annotation", but it rebuilt the entry as
        {**entry, "transcript_consequences": kept}
    — so most_severe_consequence came through UNTOUCHED, computed by VEP across
    the merged RefSeq+Ensembl set. If the severe call belonged to a transcript
    the filter then removed, the offline path reports a consequence REST could
    not have returned.

    That is not a cosmetic field. It is read by acmg/hard_coded.py (four sites)
    and evidence.py (four sites), and clients/protvar.py gates its ENTIRE
    invocation on `most_severe_consequence == "missense_variant"`.

    Here stop_gained belongs only to the RefSeq row, which the default
    (Ensembl) transcript set drops."""
    entries = [_entry("stop_gained", [
        {"transcript_id": "ENST00000355349", "consequence_terms": ["missense_variant"],
         "impact": "MODERATE"},
        {"transcript_id": "NM_000257.4", "consequence_terms": ["stop_gained"],
         "impact": "HIGH"},
    ])]
    out = vo._to_rest_transcript_set(entries, refseq=False)
    assert len(out) == 1
    assert [tc["transcript_id"] for tc in out[0]["transcript_consequences"]] \
        == ["ENST00000355349"]
    assert out[0]["most_severe_consequence"] == "missense_variant", (
        "most_severe_consequence still reports the dropped RefSeq transcript's "
        "stop_gained; ProtVar and four ACMG branches read this field"
    )


def test_most_severe_consequence_is_left_alone_when_the_kept_set_supports_it():
    """The other direction, and the more important one: this must NOT rewrite a
    value that the retained transcripts genuinely justify. VEP's own severity
    ordering is finer than anything reimplemented here, so the value is only
    replaced when it is PROVABLY unattainable from the kept rows."""
    entries = [_entry("stop_gained", [
        {"transcript_id": "ENST00000355349", "consequence_terms": ["stop_gained"],
         "impact": "HIGH"},
        {"transcript_id": "NM_000257.4", "consequence_terms": ["stop_gained"],
         "impact": "HIGH"},
    ])]
    out = vo._to_rest_transcript_set(entries, refseq=False)
    assert out[0]["most_severe_consequence"] == "stop_gained"


def test_most_severe_consequence_prefers_the_highest_impact_kept_row():
    """When it does have to choose, it picks from the KEPT rows by the same rank
    _pick_transcript_consequence uses — never from a dropped one."""
    entries = [_entry("transcript_ablation", [
        {"transcript_id": "ENST1", "consequence_terms": ["synonymous_variant"],
         "impact": "LOW"},
        {"transcript_id": "ENST2", "consequence_terms": ["frameshift_variant"],
         "impact": "HIGH"},
        {"transcript_id": "XM_1", "consequence_terms": ["transcript_ablation"],
         "impact": "HIGH"},
    ])]
    out = vo._to_rest_transcript_set(entries, refseq=False)
    assert out[0]["most_severe_consequence"] == "frameshift_variant"


_FLAGS_NEEDING_THE_VARIATION_CACHE = (
    "--check_existing", "--af", "--af_1kg", "--af_gnomad", "--af_gnomade",
    "--af_gnomadg", "--max_af", "--pubmed", "--var_synonyms", "--everything",
    "--check_svs", "--failed",
)


def test_the_runtime_argv_never_needs_the_variation_cache():
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_EXTRA_ARGS=None):
        cfg = vo._config()
        for fmt in ("vcf", "hgvs"):
            for refseq in (True, False):
                argv = vo._argv(cfg, fmt, refseq=refseq)
                for flag in _FLAGS_NEEDING_THE_VARIATION_CACHE:
                    assert flag not in argv, (
                        f"{flag} reads the variation cache, which "
                        "setup_offline_vep.sh no longer installs (17.46 GB of "
                        "the 21.5 GB tree). Either drop the flag or set "
                        "VEP_CACHE_KEEP_VARIATION=1 on the build."
                    )


def test_extra_args_asking_for_variation_data_is_reported():
    """HEARTVAR_VEP_EXTRA_ARGS is passed through verbatim, so it is the one way
    a variation-cache flag can reach vep without touching this repo. It must not
    do so silently: without the co-located data VEP returns the annotation minus
    the field, which reads downstream as "this variant is not in dbSNP"."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp",
              HEARTVAR_VEP_EXTRA_ARGS="--check_existing"):
        cfg = vo._config()
        assert cfg is not None
        assert vo.extra_args_needing_variation_cache(cfg) == ["--check_existing"], (
            "vep_offline must be able to name a variation-cache flag that "
            "arrived via HEARTVAR_VEP_EXTRA_ARGS"
        )


def test_first_prefers_present_nonnull():
    assert vo._first({"a": None, "b": 5}, "a", "b") == 5
    assert vo._first({"REVEL_score": "0.9"}, "revel_score", "REVEL_score") == "0.9"
    assert vo._first({}, "x", "y") is None


def test_nmd_escape_rules_are_the_shared_ones():
    """⚠ ONE COPY OF THE RULE, NOT TWO. This module used to reimplement the
    ClinGen PVS1 NMD decision in _nmd_escape_from_exon, and the copy OMITTED
    the third rule — a PTC in the last 50 nt of the penultimate exon escapes
    NMD. Offline returned False where REST returned True, and False keeps PVS1
    at FULL STRENGTH, so offline over-called pathogenicity for that class. Both paths now call
    ensembl_vep.nmd_escape_from_lengths."""
    from backend.clients.ensembl_vep import (
        nmd_escape_from_lengths as rule,
        nmd_escape_needs_exon_lengths as needs,
    )
    assert rule({"exon": "1/1"}, None) is True
    assert rule({"exon": "40/40"}, None) is True
    assert rule({"exon": "13/40"}, None) is False
    assert rule({}, None) is False
    assert rule({"exon": "35-36/40"}, None) is False

    assert needs({"exon": "39/40"}) is True
    for e in ("1/1", "40/40", "13/40", "35-36/40"):
        assert needs({"exon": e}) is False, e


def test_the_penultimate_exon_rule_now_fires_offline():
    """THE GAP THIS CLOSES. Ten exons of 100 nt each: the penultimate exon's 3'
    end is at cDNA 900, so the last 50 nt is 851..900. Without lengths the
    answer is False whatever the position — which is what offline used to
    return for every penultimate-exon PTC."""
    from backend.clients.ensembl_vep import nmd_escape_from_lengths as rule
    lengths = [100] * 10
    inside = {"exon": "9/10", "cdna_start": 875}
    outside = {"exon": "9/10", "cdna_start": 700}
    assert rule(inside, lengths) is True
    assert rule(outside, lengths) is False
    assert rule({"exon": "9/10", "cdna_start": 850}, lengths) is False
    assert rule({"exon": "9/10", "cdna_start": 851}, lengths) is True
    assert rule(inside, None) is False
    assert rule(inside, [100] * 8) is False


def _cli_entry():
    """A synthetic VEP-CLI --json entry using CLI/plugin key casing
    (REVEL_score / CADD_PHRED / phyloP100way_vertebrate)."""
    return {
        "assembly_name": "GRCh38",
        "seq_region_name": "14",
        "start": 23424081,
        "end": 23424081,
        "allele_string": "G/A",
        "strand": -1,
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [{
            "transcript_id": "NM_000257.4",
            "gene_symbol": "MYH7",
            "biotype": "protein_coding",
            "mane_select": "NM_000257.4",
            "hgvsc": "NM_000257.4:c.1208G>A",
            "hgvsp": "NP_000248.2:p.Arg403Gln",
            "consequence_terms": ["missense_variant"],
            "exon": "13/40",
            "impact": "MODERATE",
            "REVEL_score": "0.9",
            "CADD_PHRED": 28.5,
            "phyloP100way_vertebrate": "7.8",
        }],
    }


def test_assemble_maps_cli_cased_plugin_fields():
    res = vo._assemble(_cli_entry(), "NM_000257.4",
                       source_label="user-supplied", queried_as="NM_000257:c.1208G>A",
                       assembly_fallback="GRCh38")
    assert res is not None
    assert res["ok"] is True
    assert res["revel_score"] == 0.9
    assert res["cadd_phred"] == 28.5
    assert res["phylop100way"] == 7.8
    assert res["most_severe_consequence"] == "missense_variant"
    assert res["is_mane_select"] is True
    assert res["nmd_escape"] is False
    assert res["annotation_source"] == "offline"
    assert res["requested_transcript_honored"] is True


def test_assemble_reads_standalone_revel_plugin_key():
    entry = _cli_entry()
    tc = entry["transcript_consequences"][0]
    del tc["REVEL_score"]
    tc["REVEL"] = "0.812"
    res = vo._assemble(entry, "NM_000257.4", source_label="user-supplied",
                       queried_as="q", assembly_fallback="GRCh38")
    assert res is not None
    assert res["revel_score"] == 0.812


def test_assemble_none_when_no_consequences():
    entry = {"seq_region_name": "1", "transcript_consequences": []}
    assert vo._assemble(entry, None, source_label="x", queried_as="y",
                        assembly_fallback="GRCh38") is None


class _FakeProc:
    def __init__(self, out=b"", err=b"", rc=0, hang=False):
        self._out, self._err, self.returncode, self._hang = out, err, rc, hang

    async def communicate(self, _input=None):
        if self._hang:
            await asyncio.sleep(10)
        return self._out, self._err

    def kill(self):
        pass


def _patch_exec(monkey_proc):
    async def _fake_exec(*a, **k):
        return monkey_proc
    return _fake_exec


def _run_async(coro):
    return asyncio.run(coro)


def test_run_parses_jsonl():
    cfg = {"binary": "sh", "data_dir": "/tmp", "fasta": None,
           "assembly": "GRCh38", "extra_args": [], "timeout": 5.0}
    orig = vo.asyncio.create_subprocess_exec
    vo.asyncio.create_subprocess_exec = _patch_exec(
        _FakeProc(out=b'{"a":1}\n\n{"b":2}\n', rc=0)
    )
    try:
        entries = _run_async(vo._run(cfg, "vcf", "x", refseq=False))
    finally:
        vo.asyncio.create_subprocess_exec = orig
    assert entries == [{"a": 1}, {"b": 2}]


def test_run_nonzero_exit_returns_none():
    cfg = {"binary": "sh", "data_dir": "/tmp", "fasta": None,
           "assembly": "GRCh38", "extra_args": [], "timeout": 5.0}
    orig = vo.asyncio.create_subprocess_exec
    vo.asyncio.create_subprocess_exec = _patch_exec(_FakeProc(err=b"boom", rc=2))
    try:
        assert _run_async(vo._run(cfg, "vcf", "x", refseq=False)) is None
    finally:
        vo.asyncio.create_subprocess_exec = orig


def test_run_timeout_returns_none():
    cfg = {"binary": "sh", "data_dir": "/tmp", "fasta": None,
           "assembly": "GRCh38", "extra_args": [], "timeout": 0.05}
    orig = vo.asyncio.create_subprocess_exec
    vo.asyncio.create_subprocess_exec = _patch_exec(_FakeProc(hang=True))
    try:
        assert _run_async(vo._run(cfg, "vcf", "x", refseq=False)) is None
    finally:
        vo.asyncio.create_subprocess_exec = orig


def test_fetch_hgvs_none_when_disabled():
    with _env(HEARTVAR_VEP_OFFLINE=None):
        assert _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257.4")) is None


def test_fetch_hgvs_none_without_transcript():
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp"):
        assert _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", None)) is None


def test_fetch_hgvs_never_spawns_when_the_resolver_is_off():
    """`vep --format hgvs` needs a DB connection, so VEP refuses it outright:

        Parser/HGVS.pm:104
        throw("ERROR: Cannot use HGVS format in offline mode") if $self->param('offline');

    Unconditional — not gated on --cache, --fasta or the cache flavour. Verified
    in ensemblorg/ensembl-vep:release_113.0 on 2026-08-27 with --fasta supplied
    and the cache directory resolving.

    This was the SECOND fatal error behind the 2026-08-26 exit-2, hidden behind
    the first: --refseq died on a missing cache directory before the parser was
    ever built, so fixing the cache flavour alone would have moved the failure
    rather than removed it.

    The assertion that matters is NOT SPAWNING. Attempting the call cost every
    curation a wasted subprocess plus a full REST annotation from scratch —
    median 23.2s -> 63.8s, error rate 5.7% -> 11.8%. Falling back is correct;
    paying for a subprocess that cannot succeed first is not.
    """
    spawned = False

    async def _must_not_run(*a, **k):
        nonlocal spawned
        spawned = True
        raise AssertionError("offline HGVS must not spawn vep — VEP refuses the format")

    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_HGVS_RESOLVER=None):
        orig = vo._run
        vo._run = _must_not_run
        try:
            res = _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257.4"))
        finally:
            vo._run = orig
    assert res is None
    assert not spawned, "a subprocess was spawned for a format VEP cannot accept"


def test_a_local_resolution_reaches_vep_as_vcf_and_NEVER_as_hgvs(monkeypatch):
    """With the resolver ON, HGVS input DOES reach vep — as coordinates.

    This is the shape of the fix: VEP still cannot parse HGVS and never will, so
    the HGVS is resolved locally and handed over as `--format vcf`, the one input
    shape offline VEP serves. The invariant worth pinning is the format: if any
    future change routes HGVS input to `--format hgvs`, every curation exits 2
    again (that was the 2026-08-26 outage).
    """
    from backend.clients import hgvs_resolver as hr

    seen = {}

    async def _fake_run(cfg, fmt, stdin, *, refseq):
        seen["fmt"], seen["stdin"] = fmt, stdin
        return [_cli_entry()]

    resolved = hr.Resolved(
        chrom="14", vep_pos=23429278, vep_ref="C", vep_alt="T",
        vcf_pos=23429278, vcf_ref="C", vcf_alt="T", accession="NM_000257.4")
    monkeypatch.setattr(hr, "resolve", lambda *a, **k: resolved)

    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        vo._run = _fake_run
        try:
            res = _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257.4"))
        finally:
            vo._run = orig

    assert seen["fmt"] == "vcf", (
        "HGVS input reached vep as --format "
        f"{seen.get('fmt')!r}; VEP refuses hgvs offline (Parser/HGVS.pm:104)")
    assert "\t" in seen["stdin"], "the VCF row must be tab-delimited"
    assert res is not None and res["ok"] is True
    assert res["input_format"] == "hgvs"


def test_the_gnomad_key_for_a_resolved_variant_is_the_ANCHORED_form(monkeypatch):
    """gnomAD stores anchored left-aligned keys and the VEP CLI does not emit
    vcf_string, so forward_variant_id has to carry the resolver's anchored form.
    Using the VEP row's position instead would give indels a key that can never
    match — and PM2 reads absence as evidence, so the miss would look like
    rarity."""
    from backend.clients import hgvs_resolver as hr

    async def _fake_run(cfg, fmt, stdin, *, refseq):
        return [_cli_entry()]

    resolved = hr.Resolved(
        chrom="14", vep_pos=23429277, vep_ref="CCG", vep_alt="",
        vcf_pos=23429276, vcf_ref="CCCG", vcf_alt="C", accession="NM_000257.4")
    monkeypatch.setattr(hr, "resolve", lambda *a, **k: resolved)

    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        vo._run = _fake_run
        try:
            res = _run_async(vo.fetch_hgvs("MYH7", "c.1207_1209del", "NM_000257.4"))
        finally:
            vo._run = orig
    assert res["forward_variant_id"] == "14-23429276-CCCG-C"


def test_the_resolver_path_keeps_the_transcript_set_the_INPUT_asked_for(monkeypatch):
    """⚠⚠ THE HGVS PATH WAS SWAPPING THE TRANSCRIPT NAMESPACE ON EVERY REFSEQ INPUT.

    fetch_hgvs hard-coded `refseq=False`, so _to_rest_transcript_set narrowed the
    merged response to the ENSEMBL set even when the curator typed an NM_
    accession — which is the clinical norm. REST does the opposite: RefSeq input
    goes out with _REFSEQ_PARAMS (refseq=1) and comes back as NM_/NP_ rows.

    HGVS path, offline vs REST — the accession namespace mismatched on
    essentially every variant:
        transcript_id                       ENST00000290378 vs NM_005159.5
        hgvsc          100/101              ENST00000290378.6:c.166G>A
                                            vs NM_005159.5:c.166G>A
        hgvsp           89/101              ENSP00000290378.4:p.Val56Ile
                                            vs NP_005150.1:p.Val56Ile

    NOT COSMETIC, and _to_rest_transcript_set's own docstring says why:
    hard_coded.py's MANE lookup takes the first MANE row's hgvsp, so PS1/PM5
    would residue-match an ENSP accession against ClinVar's NP_ — a
    deterministic rule change. The dead `--format hgvs` branch further down this
    same function had it right all along (`refseq = _is_refseq(...)`); only the
    live resolver path did not.
    """
    from backend.clients import hgvs_resolver as hr

    seen: dict[str, bool] = {}

    async def _fake_run(cfg, fmt, stdin, *, refseq):
        seen["refseq"] = refseq
        return [_cli_entry()]

    for accession, expected in (("NM_000257.4", True),
                                ("NR_123456.1", True),
                                ("ENST00000355349.4", False)):
        resolved = hr.Resolved(
            chrom="14", vep_pos=23429278, vep_ref="C", vep_alt="T",
            vcf_pos=23429278, vcf_ref="C", vcf_alt="T", accession=accession)
        monkeypatch.setattr(hr, "resolve", lambda *a, **k: resolved)
        with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
                  HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
            orig = vo._run
            vo._run = _fake_run
            try:
                _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", accession))
            finally:
                vo._run = orig
        assert seen["refseq"] is expected, (
            f"{accession} must be annotated against "
            f"{'the RefSeq' if expected else 'the Ensembl'} transcript set, "
            f"matching what REST returns for the same input"
        )


def _minus_strand_cli_entry():
    """What `vep --format vcf` really returns for a MINUS-strand gene.

    Measured in the release-113 image, MYH7 R403Q (14:23429278 C>T):
        top-level: strand=1  allele_string=C/T
          tc ENST00000355349  strand=-1
          tc NM_000257.4      strand=-1
    The top level describes the VCF row (forward strand); the transcript
    consequences carry the gene's actual strand."""
    e = _cli_entry()
    e["start"] = e["end"] = 23429278
    e["allele_string"] = "C/T"
    e["strand"] = 1
    e["transcript_consequences"][0]["strand"] = -1
    return e


def test_the_hgvs_path_reports_the_transcript_strand_representation(monkeypatch):
    """⚠ FLIPPING THE FLAG WOULD OTHERWISE CHANGE codon_genomic_positions.

    The two REST endpoints disagree on what `strand`/`allele_string` mean, and
    both were measured on 2026-08-29 for MYH7 R403Q, a MINUS-strand gene:
        /vep/human/region : strand= 1, allele_string=C/T   (forward)
        /vep/human/hgvs   : strand=-1, allele_string=G/A   (transcript strand)

    Three clients rely on that PAIR being consistent — protvar.py:290,
    protvar.py:290 and alphamissense.py:153 both complement the allele when
    strand == -1 to get back to the forward strand — so either convention is
    safe for them.

    ensembl_vep.codon_genomic_positions is NOT: it uses `strand` ALONE as the
    direction the CDS runs in (`first = pos - offset * strand`). Reproduced
    exactly, for MYH7 at 14:23429278:
        c.1207 (offset 0): strand=-1 -> [..277,278,279]  strand=1 -> [..279,280,281]
        c.1208 (offset 1): identical
        c.1209 (offset 2): strand=-1 -> [..277,278,279]  strand=1 -> [..275,276,277]
    i.e. TWO OF THREE codon positions differ. So an offline HGVS path reporting
    the VCF row's +1 would silently move the codon for minus-strand genes —
    a regression against the REST path it is replacing, in the field feeding
    the same-codon/same-residue evidence.

    fetch_region keeps the forward representation, because that is what
    /vep/human/region returns and what the coordinate path is replacing."""
    from backend.clients import hgvs_resolver as hr

    async def _fake_run(cfg, fmt, stdin, *, refseq):
        return [_minus_strand_cli_entry()]

    resolved = hr.Resolved(
        chrom="14", vep_pos=23429278, vep_ref="C", vep_alt="T",
        vcf_pos=23429278, vcf_ref="C", vcf_alt="T", accession="NM_000257.4")
    monkeypatch.setattr(hr, "resolve", lambda *a, **k: resolved)

    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        vo._run = _fake_run
        try:
            hgvs_res = _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257.4"))
            coord_res = _run_async(
                vo.fetch_region("14", 23429278, "C", "T", "GRCh38", "NM_000257.4"))
        finally:
            vo._run = orig

    assert hgvs_res["strand"] == -1, (
        "the HGVS path must report the transcript strand, as /vep/human/hgvs "
        "does — codon_genomic_positions reads this as the CDS direction"
    )
    assert hgvs_res["allele_string"] == "G/A", (
        "strand and allele_string are a PAIR: protvar/alphamissense "
        "complement the allele when strand == -1, so reporting strand=-1 with a "
        "forward allele_string would make them complement a forward allele"
    )
    assert coord_res["strand"] == 1
    assert coord_res["allele_string"] == "C/T"


def test_a_version_substitution_is_surfaced_on_the_result(monkeypatch):
    """Resolving a bare accession picks a version, and that choice decides which
    transcript the entire ACMG evaluation ran on. It must not be silent."""
    from backend.clients import hgvs_resolver as hr

    async def _fake_run(cfg, fmt, stdin, *, refseq):
        return [_cli_entry()]

    resolved = hr.Resolved(
        chrom="14", vep_pos=23429278, vep_ref="C", vep_alt="T",
        vcf_pos=23429278, vcf_ref="C", vcf_alt="T", accession="NM_000257.4",
        note="NM_000257 was supplied without a version; resolved against NM_000257.4")
    monkeypatch.setattr(hr, "resolve", lambda *a, **k: resolved)

    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        vo._run = _fake_run
        try:
            res = _run_async(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257"))
        finally:
            vo._run = orig
    assert "without a version" in (res.get("resolver_note") or "")


def test_the_banner_does_not_claim_more_than_the_offline_path_delivers(tmp_path):
    """A banner reading plain "on" invited the wrong conclusion on 2026-08-26 —
    the cache was complete and the annotation path still was not working. It has
    to say which input it actually covers."""
    (tmp_path / ".heartvar_vep_manifest.json").write_text(
        '{"vep_release": "113", "revel": "/x/revel.tsv.gz"}'
    )
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA=str(tmp_path)):
        status = vo.offline_status()
    assert "COORDINATE input only" in status, status
    assert "113" in status and "with REVEL" in status, status


def test_fetch_region_build_mismatch_returns_none():
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        assert _run_async(
            vo.fetch_region("7", 117548628, "C", "T", "GRCh37", None)
        ) is None


def test_fetch_region_happy_path_sets_coord_fields():
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        async def _fake_run(cfg, fmt, stdin, *, refseq):
            assert fmt == "vcf"
            assert stdin == "7\t117548628\t.\tC\tT\t.\t.\t."
            return [_cli_entry()]
        vo._run = _fake_run
        try:
            res = _run_async(
                vo.fetch_region("7", 117548628, "C", "T", "GRCh38", None)
            )
        finally:
            vo._run = orig
        assert res is not None and res["ok"] is True
        assert res["input_format"] == "coordinates"
        assert res["input_coords"] == "7-117548628-C-T"
        assert res["forward_variant_id"] == "7-117548628-C-T"
        assert res["derived_hgvs"] == "c.1208G>A"


def test_the_cli_vcf_row_is_tab_delimited_because_the_parser_splits_on_TAB():
    """THE SILENT ONE, measured 2026-08-28 in ensemblorg/ensembl-vep:release_113.0.

    `vep --format vcf` runs the real VCF parser, which splits on TAB. Feeding it
    the space-separated row (what this code sent, copied from the REST helper)
    gives, on an identical cache with identical flags:

        space-separated : exit 0, 0 bytes stdout, 7 x
                          "Use of uninitialized value $ref ... Parser/VCF.pm"
        tab-separated   : exit 0, 15 transcript_consequences

    EXIT 0 IS THE POINT. _run only logged stderr on a non-zero status, so this
    produced NO log line: entries came back empty, fetch_region returned None,
    the caller fell back to REST. Offline VEP was permanently inert on the only
    input path it can serve, while offline_status() still reported "on" and every
    curation paid for a wasted subprocess. Strictly more silent than the
    2026-08-26 exit 2, which at least had a status to notice.
    """
    captured = {}
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp", HEARTVAR_VEP_ASSEMBLY="GRCh38"):
        orig = vo._run
        async def _fake_run(cfg, fmt, stdin, *, refseq):
            captured["stdin"] = stdin
            return [_cli_entry()]
        vo._run = _fake_run
        try:
            _run_async(vo.fetch_region("14", 23429278, "C", "T", "GRCh38", None))
        finally:
            vo._run = orig

    row = captured["stdin"]
    assert "\t" in row, (
        "the CLI VCF row is not tab-delimited; VEP's Parser/VCF.pm splits on TAB "
        "and a space-separated row yields exit 0 with zero output"
    )
    assert row.split("\t") == ["14", "23429278", ".", "C", "T", ".", ".", "."], row
    assert " " not in row, f"a space survives in the CLI VCF row: {row!r}"


def test_the_rest_helper_keeps_spaces_and_the_cli_does_not():
    """The two builders must DIFFER, deliberately. ensembl_vep's REST region
    endpoint accepts a whitespace-separated row; the CLI parser does not. Anyone
    "unifying" these two would reintroduce the silent failure above, so the
    difference is asserted rather than left to a comment."""
    from backend.clients.ensembl_vep import _vcf_input_for_region

    rest_row = _vcf_input_for_region("14", 23429278, "C", "T")
    assert "\t" not in rest_row, (
        "the REST row gained tabs; that endpoint takes whitespace-separated "
        "input and this change was probably meant for the CLI path"
    )
    assert rest_row == "14 23429278 . C T . . ."


def test_a_clean_exit_that_annotated_nothing_is_logged(caplog):
    """exit 0 + empty stdout was the one failure shape with no log line at all.
    _run logged stderr only when returncode != 0, so the space-delimiter bug was
    invisible for as long as it existed."""
    import logging
    cfg = {"binary": "sh", "data_dir": "/tmp", "fasta": None,
           "assembly": "GRCh38", "extra_args": [], "timeout": 5.0,
           "revel": None, "plugins_dir": None}
    orig = vo.asyncio.create_subprocess_exec
    vo.asyncio.create_subprocess_exec = _patch_exec(
        _FakeProc(out=b"", err=b"Use of uninitialized value $ref", rc=0)
    )
    try:
        with caplog.at_level(logging.WARNING):
            assert _run_async(vo._run(cfg, "vcf", "x", refseq=False)) is None
    finally:
        vo.asyncio.create_subprocess_exec = orig

    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "exited 0 but returned no usable annotation" in msgs, (
        "a vep that ran cleanly and annotated nothing is still silent:\n" + msgs
    )
    assert "uninitialized value" in msgs, "the stderr HEAD is not being logged"


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()


def test_module_docstring_does_not_claim_a_slim_runtime_image():
    """The docstring told readers the vep binary was absent and that a multi-stage
    build had to be added first. The runtime image has been FROM
    ensemblorg/ensembl-vep for a while, so that prerequisite read as unmet when it
    was met — and this project has twice been bitten by a comment that quietly
    stopped being true (the "~3 GB" REVEL figure became a too-low disk floor).
    """
    from backend.clients import vep_offline
    from .conftest import require_unstripped_prose
    doc = vep_offline.__doc__ or ""
    require_unstripped_prose(doc, "offline-vs-REST parity")
    assert "python:3.11-slim" not in doc
    assert "~30 GB" not in doc, "measured size is ~27 GB"
    assert "offline-vs-REST parity" in doc, (
        "the docstring must still name the one gate that IS outstanding"
    )


def test_revel_plugin_uses_the_named_file_parameter():
    """The release-113 REVEL plugin's SYNOPSIS is `--plugin REVEL,file=<path>`.
    A bare positional path is accepted on the command line and then yields no
    score at all — so PP3/BP4 loses its calibrated signal with nothing erroring.
    Caught by the build-time check in setup_offline_vep.sh on the first real run.
    """
    from backend.clients.vep_offline import _argv
    cfg = {
        "binary": "vep", "data_dir": "/d", "fasta": None, "assembly": "GRCh38",
        "revel": "/r.tsv.gz", "plugins_dir": "/p", "extra_args": [], "timeout": 45.0,
    }
    argv = _argv(cfg, "vcf", refseq=False)
    joined = " ".join(argv)
    assert "REVEL,file=/r.tsv.gz" in joined, joined
    assert "REVEL,/r.tsv.gz" not in joined, joined


def test_status_says_off_when_the_flag_is_unset():
    with _env(HEARTVAR_VEP_OFFLINE=None):
        assert "off" in vo.offline_status()


def test_status_flags_set_but_unusable_rather_than_claiming_on():
    """The dangerous state: configured, silently inert, annotating via REST. It
    must NOT read as 'on'."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/nonexistent-cache-dir"):
        status = vo.offline_status()
        assert "UNUSABLE" in status, status
        assert not status.startswith("on"), status


def test_status_reports_the_cache_release_and_revel(tmp_path):
    """`vep --offline` refuses a cache from a different release by returning
    NOTHING, which reads here as a REST fallback. So the release the cache was
    built for is the single most useful thing to print."""
    (tmp_path / ".heartvar_vep_manifest.json").write_text(
        '{"vep_release": "113", "revel": "/data/vep/revel/new_tabbed_revel_grch38.tsv.gz"}'
    )
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA=str(tmp_path)):
        status = vo.offline_status()
        assert status.startswith("on"), status
        assert "113" in status, status
        assert "with REVEL" in status, status


def test_status_still_says_on_without_a_manifest(tmp_path):
    """The manifest is the BUILD's record, not a runtime requirement — a cache
    can annotate perfectly well without one, so this is not 'unusable'."""
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA=str(tmp_path)):
        status = vo.offline_status()
        assert status.startswith("on"), status
        assert "no build manifest" in status, status


_MERGED_ENTRY = {
    "input": "NM_000257:c.1208G>A",
    "transcript_consequences": [
        {"transcript_id": "ENST00000355349", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "canonical": 1,
         "mane_select": "NM_000257.4", "amino_acids": "R/Q",
         "hgvsp": "ENSP00000347507.3:p.Arg403Gln",
         "hgvsc": "ENST00000355349.4:c.1208G>A"},
        {"transcript_id": "NM_000257.4", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "canonical": 1,
         "mane_select": "ENST00000355349.4", "amino_acids": "R/Q",
         "hgvsp": "NP_000248.2:p.Arg403Gln",
         "hgvsc": "NM_000257.4:c.1208G>A"},
        {"transcript_id": "ENST00000000001", "biotype": "protein_coding",
         "consequence_terms": ["intron_variant"]},
    ],
}


def test_merged_response_is_narrowed_to_the_refseq_set_for_refseq_input():
    out = vo._to_rest_transcript_set([_MERGED_ENTRY], refseq=True)
    ids = [tc["transcript_id"] for tc in out[0]["transcript_consequences"]]
    assert ids == ["NM_000257.4"], ids


def test_merged_response_is_narrowed_to_the_ensembl_set_otherwise():
    """Mirrors VEP_PARAMS, which the client uses for ENST input AND for every
    coordinate lookup — REST returned no NM_ rows for either."""
    out = vo._to_rest_transcript_set([_MERGED_ENTRY], refseq=False)
    ids = [tc["transcript_id"] for tc in out[0]["transcript_consequences"]]
    assert ids == ["ENST00000355349", "ENST00000000001"], ids


def test_the_MANE_hgvsp_stays_RefSeq_so_PS1_PM5_match_the_right_residue():
    """THE REASON THE FILTER EXISTS, and it is a deterministic-rule bug.

    ``_build_transcript_table`` sets ``is_mane_select`` from a truthy
    ``mane_select``, and in a MERGED response BOTH members of a MANE pair carry
    one — NM_000257.4 names ENST00000355349.4 and ENST00000355349 names
    NM_000257.4. Two MANE-Select rows then tie all the way down the sort key to
    ``str(transcript_id)``, where "ENST…" sorts before "NM_…". So the Ensembl row
    wins, and ``hard_coded.py``'s MANE lookup — which takes the FIRST MANE row's
    hgvsp and uses it to residue-match ClinVar for PS1/PM5 — would compare
    ENSP numbering against a ClinVar NP_ record.

    An alphabetical tiebreak silently changing which residue PS1/PM5 matches on
    is exactly the class of bug that cost 2026-08-26. Unfiltered is asserted
    first here, because a test that only checks the fixed path would pass just as
    well if the filter were deleted.
    """
    from backend.clients.ensembl_vep import _build_transcript_table

    def _mane_hgvsp(consequences):
        rows, _differs, tset = _build_transcript_table(consequences, consequences[0])
        for row in rows:
            if row.get("is_mane_select") and row.get("hgvsp"):
                return row["hgvsp"], tset
        return None, tset

    unfiltered = _MERGED_ENTRY["transcript_consequences"]
    bad, bad_set = _mane_hgvsp(unfiltered)
    assert bad == "ENSP00000347507.3:p.Arg403Gln", (
        "if this no longer holds, the merged-cache hazard has changed shape — "
        "re-derive the filter rather than deleting this test"
    )
    assert bad_set == "mixed", bad_set

    filtered = vo._to_rest_transcript_set(
        [_MERGED_ENTRY], refseq=True)[0]["transcript_consequences"]
    good, good_set = _mane_hgvsp(filtered)
    assert good == "NP_000248.2:p.Arg403Gln", good
    assert good_set == "RefSeq", good_set


def test_the_filter_drops_an_entry_whose_whole_set_is_wrong():
    """Wrong transcript set is worse than no answer: dropping the entry makes the
    caller fall back to REST, which is the safe direction. Keeping it would
    annotate the RefSeq slice against Ensembl transcripts."""
    ensembl_only = {"transcript_consequences": [
        {"transcript_id": "ENST00000355349", "consequence_terms": ["missense_variant"]},
    ]}
    assert vo._to_rest_transcript_set([ensembl_only], refseq=True) == []


def test_the_filter_passes_through_an_entry_with_no_consequences():
    """Not a transcript-set problem — _assemble owns the "no usable
    consequence" decision, and this function must not pre-empt it."""
    bare = {"input": "x", "most_severe_consequence": "intergenic_variant"}
    assert vo._to_rest_transcript_set([bare], refseq=True) == [bare]

def test_a_failing_offline_run_logs_the_CAUSE_not_the_stack(monkeypatch, caplog):
    """VEP is a Perl program: the cause is on the FIRST stderr line and the last
    few hundred characters are the stack trace and the "Ensembl API version"
    banner. The log took the TAIL, so when offline VEP began exiting 2 on every
    curation in production the log said only "...Runner.pm:128 / STACK toplevel /
    Date / Ensembl API version = 113" — which names nothing. Combined with the
    REST fallback swallowing the failure, a completely broken offline path looked
    like mere slowness for as long as nobody measured it."""
    import logging
    src = ("ERROR: Cannot detect format of input\n"
           + "filler line\n" * 40
           + "STACK Bio::EnsEMBL::VEP::Runner::run Runner.pm:200\n"
           "STACK toplevel /opt/vep/src/ensembl-vep/vep:46\n"
           "Ensembl API version = 113\n")

    class _Proc:
        returncode = 2
        async def communicate(self, _input=None):
            return b"", src.encode()

    async def _fake_exec(*a, **k):
        return _Proc()

    monkeypatch.setattr(asyncio, "create_subprocess_exec", _fake_exec)
    with _env(HEARTVAR_VEP_OFFLINE="1", HEARTVAR_VEP_BINARY="sh",
              HEARTVAR_VEP_DATA="/tmp"):
        with caplog.at_level(logging.WARNING, logger="heartvar.vep.offline"):
            cfg = {"binary": "sh", "data_dir": "/tmp", "fasta": None,
                   "assembly": "GRCh38", "revel": None, "plugins_dir": None,
                   "extra_args": [], "timeout": 30.0}
            asyncio.run(vo._run(cfg, "vcf", "in", refseq=False))
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "Cannot detect format of input" in msg, msg
