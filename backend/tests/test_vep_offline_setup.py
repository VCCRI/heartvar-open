"""Guards on scripts/setup_offline_vep.sh — the resume/completion state machine.

The bug being fixed: the old completion check was "does $VEP_DATA/homo_sapiens
exist and is it non-empty?". The data-builder sidecar restarts on EVERY deploy,
so an unpack interrupted halfway leaves a non-empty directory that every later
run reads as finished — a permanently half-installed cache that nothing reports.
That is the same silent-staleness class as the 2026-08-01 refresh, which claimed
success and rebuilt nothing.

The replacement is a FUNCTIONAL check: a component counts as complete only when
`vep` actually annotates MYH7 R403Q and returns the expected field. Requiring a
REVEL score in that output doubles as the plugin field-name validation that gates
HEARTVAR_VEP_OFFLINE (see backend/clients/vep_offline.py, "VALIDATION GAP").

Functions are extracted and evaluated in a throwaway bash with a STUB vep on
PATH — never sourced, so nothing here can start a ~25 GB download.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

from .conftest import (  # noqa: E402
    require_file, require_unstripped_prose,
)
SETUP = ROOT / "scripts" / "setup_offline_vep.sh"
SCRIPT = SETUP.read_text()

VEP_JSON_PLAIN = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"gene_symbol":"MYH7",'
    '"amino_acids":"R/Q","protein_start":403}]}'
)
VEP_JSON_HGVS = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"amino_acids":"R/Q",'
    '"hgvsc":"ENST00000355349.4:c.1208G>A"}]}'
)
VEP_JSON_REVEL = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"amino_acids":"R/Q",'
    '"transcript_id":"ENST00000355349","revel":0.886}]}'
)
VEP_JSON_RUNTIME = (
    '{"input":"14 23429278 . C T . . .","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"gene_symbol":"MYH7",'
    '"transcript_id":"ENST00000355349.4","amino_acids":"R/Q",'
    '"hgvsp":"ENSP00000347507.3:p.Arg403Gln","mane_select":"NM_000257.4",'
    '"sift_prediction":"deleterious_low_confidence","sift_score":0,'
    '"polyphen_prediction":"probably_damaging","polyphen_score":1,'
    '"revel":0.886},'
    '{"consequence_terms":["missense_variant"],"gene_symbol":"MYH7",'
    '"transcript_id":"NM_000257.4","amino_acids":"R/Q",'
    '"hgvsp":"NP_000248.2:p.Arg403Gln","mane_select":"ENST00000355349.4",'
    '"sift_prediction":"deleterious_low_confidence","sift_score":0,'
    '"polyphen_prediction":"probably_damaging","polyphen_score":1,'
    '"revel":0.886}]}'
)
VEP_JSON_REVEL_UPPER = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"amino_acids":"R/Q","REVEL":"0.913"}]}'
)
VEP_JSON_NO_REVEL = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"amino_acids":"R/Q",'
    '"transcript_id":"ENST00000355349"}]}'
)
VEP_JSON_SIFT = (
    '{"input":"14 23429278 . C T","transcript_consequences":'
    '[{"consequence_terms":["missense_variant"],"gene_symbol":"MYH7",'
    '"amino_acids":"R/Q","protein_start":403,'
    '"sift_prediction":"deleterious_low_confidence","sift_score":0,'
    '"polyphen_prediction":"probably_damaging","polyphen_score":1}]}'
)
VEP_JSON_WRONG_VARIANT = (
    '{"input":"14 23424081 . G A","transcript_consequences":'
    '[{"consequence_terms":["intron_variant"],"gene_symbol":"MYH7",'
    '"hgvsc":"ENST00000355349.4:c.1223-45G>A"}]}'
)


def _extract_function(name: str) -> str:
    match = re.search(rf"^{name}\(\)\s*\{{.*?^\}}", SCRIPT, re.MULTILINE | re.DOTALL)
    assert match, f"{name}() not found in setup_offline_vep.sh"
    return match.group(0)


def _stub_vep(tmp_path: Path, *, version: str = "113", emit: str = "") -> Path:
    """A fake `vep` on PATH. Prints a version banner for --help, else `emit`."""
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    vep = bindir / "vep"
    vep.write_text(
        "#!/usr/bin/env bash\n"
        'for a in "$@"; do [[ "$a" == "--help" ]] && '
        f'{{ echo "ensembl-vep : {version}"; exit 0; }}; done\n'
        "cat >/dev/null\n"
        f"printf '%s\\n' {json.dumps(emit)}\n"
    )
    vep.chmod(0o755)
    return bindir


def _bash(tmp_path: Path, functions: list[str], body: str, *,
          bindir: Path | None = None, env: str = "") -> subprocess.CompletedProcess:
    preamble = "\n".join([
        f'VEP_DATA="{tmp_path}/vep"',
        "VEP_ASSEMBLY=GRCh38",
        "VEP_CACHE_SPECIES=homo_sapiens_merged",
        "VEP_FASTA_SPECIES=homo_sapiens",
        "VEP_CACHE_FLAVOUR_FLAG=--merged",
        "VEP_BINARY=vep",
        f'VEP_STATE_DIR="{tmp_path}/vep/.heartvar_vep_state"',
        f'VEP_MANIFEST="{tmp_path}/vep/.heartvar_vep_manifest.json"',
        f'PLUGINS_DIR="{tmp_path}/vep/Plugins"',
        r"PROBE_VCF='14\t23429278\t.\tC\tT\t.\t.\t.'",
        'PROBE_AA="R/Q"',
        'PROBE_HGVSC="c.1208G>A"',
        "PROBE_RUNTIME_VCF='14 23429278 . C T . . .'",
        'PROBE_RUNTIME_RESIDUE="p.Arg403Gln"',
        'PROBE_RUNTIME_ENSP="ENSP00000347507"',
        'PROBE_RUNTIME_MERGED_TX="NM_000257"',
        f'CACHE_DIR="{tmp_path}/vep/homo_sapiens_merged/113_GRCh38"',
        f'FASTA_DIR="{tmp_path}/vep/homo_sapiens/113_GRCh38"',
        f'FASTA_PATH="{tmp_path}/vep/homo_sapiens/113_GRCh38/ref.fa.gz"',
        f'REVEL_TSV="{tmp_path}/vep/revel/r.tsv.gz"',
        env,
    ])
    path_line = f'export PATH="{bindir}:$PATH"\n' if bindir else ""
    script = path_line + preamble + "\n" + "\n".join(functions) + "\n" + body
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          timeout=60)


def test_release_comes_from_the_installed_binary(tmp_path):
    """Derived, never hardcoded: a cache built for a release the local vep cannot
    read is silently useless — vep --offline returns nothing, and vep_offline.py
    reads that as "fall back to REST"."""
    bindir = _stub_vep(tmp_path, version="114")
    out = _bash(tmp_path, [_extract_function("_vep_release")],
                "_vep_release", bindir=bindir)
    assert out.stdout.strip() == "114", out.stderr


def test_release_prefers_constants_pm_over_the_banner(tmp_path):
    """Constants.pm is the value VEP itself uses, it is instant to read, and its
    format is stable. `vep --help` loads BioPerl — seconds on real hardware,
    minutes under emulation — and this runs on every build."""
    bindir = _stub_vep(tmp_path, version="999")
    constants = bindir / "modules" / "Bio" / "EnsEMBL" / "VEP" / "Constants.pm"
    constants.parent.mkdir(parents=True)
    constants.write_text("use constant {\n  VEP_VERSION     = 113,\n};\n")
    out = _bash(tmp_path, [_extract_function("_vep_release")], "_vep_release",
                bindir=bindir)
    assert out.stdout.strip() == "113", f"{out.stdout!r} {out.stderr}"


def test_release_falls_back_to_the_banner_without_constants_pm(tmp_path):
    bindir = _stub_vep(tmp_path, version="114")
    out = _bash(tmp_path, [_extract_function("_vep_release")], "_vep_release",
                bindir=bindir)
    assert out.stdout.strip() == "114", f"{out.stdout!r} {out.stderr}"


def test_release_falls_back_to_env_when_binary_is_absent(tmp_path):
    out = _bash(tmp_path, [_extract_function("_vep_release")], "_vep_release",
                env="VEP_RELEASE=112; VEP_BINARY=/nonexistent/vep")
    assert out.stdout.strip() == "112", out.stderr


VALID_INFO_TXT = (
    "# comment line, skipped\n"
    "assembly\tGRCh38\n"
    "polyphen\tb\n"
    "sift\tb\n"
    "species\thomo_sapiens\n"
    "cache_region_size\t1000000\n"
    "source_sift\tsift5.2.2\n"
    "source_polyphen\t2.2.3\n"
    "var_type\ttabix\n"
)


def _write_info(tmp_path: Path, text: str | None = VALID_INFO_TXT) -> Path:
    cache_dir = tmp_path / "vep" / "homo_sapiens_merged" / "113_GRCh38"
    cache_dir.mkdir(parents=True, exist_ok=True)
    info = cache_dir / "info.txt"
    if text is None:
        if info.exists():
            info.unlink()
    else:
        info.write_text(text)
    return info


def _check(tmp_path, component: str, emit: str, *, info: str | None = VALID_INFO_TXT):
    bindir = _stub_vep(tmp_path, emit=emit)
    (tmp_path / "vep").mkdir(parents=True, exist_ok=True)
    _write_info(tmp_path, info)
    funcs = [_extract_function("_vep_smoke"),
             _extract_function("_check_cache_info"),
             _extract_function(f"_check_{component}")]
    fa = _write_fasta(tmp_path)
    return _bash(tmp_path, funcs,
                 f"_check_{component} && echo COMPLETE || echo INCOMPLETE",
                 bindir=bindir, env=f'FASTA_PATH="{fa}"; REVEL_TSV=/r.tsv.gz')


def test_revel_check_accepts_the_plugins_actual_lowercase_key(tmp_path):
    """THE BUG. release-113 REVEL.pm emits `"revel":0.886`; the check grepped for
    `"REVEL"` case-sensitively, so it declared a perfectly working offline install
    broken five runs in a row while the cache, FASTA and REVEL data were all
    fine."""
    assert "COMPLETE" in _check(tmp_path, "revel", VEP_JSON_REVEL).stdout


def test_revel_check_still_accepts_the_upper_case_key(tmp_path):
    """Both spellings are ones clients/vep_offline.py::_first reads, so both must
    count as installed."""
    assert "COMPLETE" in _check(tmp_path, "revel", VEP_JSON_REVEL_UPPER).stdout


def test_revel_check_fails_when_no_score_comes_back(tmp_path):
    """The check must still catch the real failure: a missense call with no REVEL
    key, which would silently cost PP3/BP4 its calibrated signal."""
    assert "INCOMPLETE" in _check(tmp_path, "revel", VEP_JSON_NO_REVEL).stdout


def test_the_accepted_revel_spellings_match_what_the_client_reads():
    """These two lists drifting apart IS the bug. The build must not certify a
    key the runtime cannot find, nor reject one it can."""
    client = (ROOT / "backend" / "clients" / "vep_offline.py").read_text()
    accepted = re.search(r"grep -qE '\\?\"\(([^)]*)\)", SCRIPT)
    assert accepted, "could not find the REVEL key alternation in _check_revel"
    for key in accepted.group(1).split("|"):
        assert f'"{key}"' in client, (
            f"_check_revel accepts {key!r} but clients/vep_offline.py never reads it"
        )


def test_the_probe_is_the_real_r403q_coordinate():
    """THE BUG, pinned. The probe was `14 23424081 G A` until 2026-08-26. MYH7
    p.Arg403Gln is at 14:23,429,278 (GRCh38; ClinVar 14087,
    NM_000257.4:c.1208G>A), and MYH7 is on the minus strand so the genomic change
    is C>T. The old line declared reference G where the genome has C — it was not
    a variant, so no REVEL score could ever exist for it and offline VEP could
    never complete."""
    code = "\n".join(l for l in SCRIPT.splitlines() if not l.lstrip().startswith("#"))
    assert "23429278" in code
    assert "23424081" not in code, (
        "23424081 is 5.2 kb from R403Q and its declared reference base is wrong"
    )


def test_cache_check_rejects_a_variant_that_is_not_the_probe(tmp_path):
    """`grep -q '"consequence_terms"'` passed on ANY position overlapping a
    transcript, which is how a non-variant got recorded as a working cache. The
    check must identify the residue, not just that VEP said something."""
    out = _check(tmp_path, "cache", VEP_JSON_WRONG_VARIANT).stdout
    assert "INCOMPLETE" in out, out


def test_fasta_check_rejects_a_different_variants_hgvs(tmp_path):
    """Same failure mode one component later: any hgvsc satisfied the old check.
    The expected HGVS is what proves the reference is both readable and right."""
    out = _check(tmp_path, "fasta", VEP_JSON_WRONG_VARIANT).stdout
    assert "INCOMPLETE" in out, out


def test_cache_complete_when_vep_returns_a_consequence(tmp_path):
    assert "COMPLETE" in _check(tmp_path, "cache", VEP_JSON_PLAIN).stdout


def test_cache_incomplete_when_vep_returns_nothing(tmp_path):
    """The half-unpacked-cache case: files present, annotation does not work."""
    assert "INCOMPLETE" in _check(tmp_path, "cache", "").stdout


def test_fasta_complete_only_with_an_hgvs_string(tmp_path):
    assert "COMPLETE" in _check(tmp_path, "fasta", VEP_JSON_HGVS).stdout
    assert "INCOMPLETE" in _check(tmp_path, "fasta", VEP_JSON_PLAIN).stdout


def test_revel_complete_only_with_a_revel_score(tmp_path):
    """Doubles as the plugin field-name validation gating HEARTVAR_VEP_OFFLINE: if
    the plugin key casing differs from what the client expects, it shows up here,
    in the build log, instead of silently dropping PP3/BP4's calibrated signal."""
    assert "COMPLETE" in _check(tmp_path, "revel", VEP_JSON_REVEL).stdout
    assert "INCOMPLETE" in _check(tmp_path, "revel", VEP_JSON_PLAIN).stdout


def _argv_recording_vep(tmp_path: Path):
    """A stub vep that records the argv it was called with."""
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    log = tmp_path / "argv.txt"
    vep = bindir / "vep"
    vep.write_text(
        "#!/usr/bin/env bash\n"
        f'printf "%s\\n" "$*" >> "{log}"\n'
        "cat >/dev/null\n"
        "echo '{}'\n"
    )
    vep.chmod(0o755)
    return bindir, log


def test_fasta_check_requests_the_field_it_greps_for(tmp_path):
    """Regression, found in production. _check_fasta grepped for "hgvsc" but never
    passed --hgvs, so vep never emitted the field and the FASTA was reported
    unusable on every build — with INSTALL.pl having placed it correctly. A failed
    component means no manifest, so the step could never finish."""
    bindir, log = _argv_recording_vep(tmp_path)
    (tmp_path / "vep").mkdir(parents=True, exist_ok=True)
    _bash(tmp_path, [_extract_function("_vep_smoke"), _extract_function("_check_fasta")],
          "_check_fasta || true", bindir=bindir, env="FASTA_PATH=/f.fa.gz")
    argv = log.read_text()
    assert "--hgvs" in argv, f"--hgvs not requested; vep cannot emit hgvsc: {argv}"
    assert "--fasta /f.fa.gz" in argv, argv


def test_cache_check_does_not_request_hgvs(tmp_path):
    """--hgvs in offline mode requires --fasta, and the cache check runs BEFORE the
    FASTA is installed. Adding it to _vep_smoke globally would make the cache
    check fail on a fresh mount.

    info.txt is written here because _check_cache now gates on it first — with
    it absent the function returns before ever invoking vep, and this test would
    pass for the wrong reason (an empty argv log contains no --hgvs either)."""
    bindir, log = _argv_recording_vep(tmp_path)
    (tmp_path / "vep").mkdir(parents=True, exist_ok=True)
    _write_info(tmp_path)
    _bash(tmp_path, [_extract_function("_vep_smoke"),
                     _extract_function("_check_cache_info"),
                     _extract_function("_check_cache")],
          "_check_cache || true", bindir=bindir)
    recorded = log.read_text()
    assert recorded.strip(), (
        "vep was never invoked, so this test proves nothing — the info.txt gate "
        "returned first"
    )
    assert "--hgvs" not in recorded


def test_revel_check_requests_the_plugin(tmp_path):
    bindir, log = _argv_recording_vep(tmp_path)
    (tmp_path / "vep").mkdir(parents=True, exist_ok=True)
    _bash(tmp_path, [_extract_function("_vep_smoke"), _extract_function("_check_revel")],
          "_check_revel || true", bindir=bindir, env="REVEL_TSV=/r.tsv.gz")
    argv = log.read_text()
    assert "--plugin REVEL,file=/r.tsv.gz" in argv, (
        "the release-113 plugin takes NAMED parameters; a bare path yields no "
        f"score and silently costs PP3/BP4 its signal: {argv}"
    )


STATE_FUNCS = ("_state_ok", "_state_write")


def _state(tmp_path, body, env=""):
    return _bash(tmp_path, [_extract_function(f) for f in STATE_FUNCS], body,
                 env=f"REL=113; {env}")


def test_state_write_then_ok_roundtrips(tmp_path):
    out = _state(tmp_path, "_state_write cache; _state_ok cache && echo OK || echo NO")
    assert "OK" in out.stdout, out.stderr


def test_state_absent_is_not_ok(tmp_path):
    assert "NO" in _state(tmp_path, "_state_ok cache && echo OK || echo NO").stdout


def test_state_from_a_different_release_is_not_ok(tmp_path):
    """A VEP release bump must invalidate the cache. The old non-empty-directory
    check never noticed one."""
    _state(tmp_path, "_state_write cache")
    out = _state(tmp_path, "_state_ok cache && echo OK || echo NO", env="REL=114")
    assert "NO" in out.stdout


def test_state_from_a_different_assembly_is_not_ok(tmp_path):
    _state(tmp_path, "_state_write cache")
    out = _bash(tmp_path, [_extract_function(f) for f in STATE_FUNCS],
                "_state_ok cache && echo OK || echo NO",
                env="REL=113; VEP_ASSEMBLY=GRCh37")
    assert "NO" in out.stdout


def test_partial_state_leaves_revel_incomplete(tmp_path):
    """THE resume case: cache done, REVEL not. The next run must pick up at REVEL
    rather than either re-downloading 25 GB or declaring itself finished."""
    out = _state(
        tmp_path,
        "_state_write cache\n"
        "_state_ok cache && echo CACHE_OK || echo CACHE_NO\n"
        "_state_ok revel && echo REVEL_OK || echo REVEL_NO",
    )
    assert "CACHE_OK" in out.stdout
    assert "REVEL_NO" in out.stdout


def test_manifest_is_not_written_on_partial_state():
    """build_all.sh's `versioned` KEEP branch skips a step whose output artifact
    merely EXISTS. If the manifest appeared after a partial install the step would
    be skipped forever with REVEL missing — the original bug in a new place. So
    the manifest is written only once every component verifies."""
    body = _extract_function("_manifest_write")
    assert "vep_release" in body and "assembly" in body
    assert re.search(r"_state_ok\s+cache", SCRIPT)
    assert re.search(r"_state_ok\s+revel", SCRIPT)
    assert re.search(r"^\s*_manifest_write\s*$", SCRIPT, re.MULTILINE), (
        "_manifest_write is never called"
    )


SPACE_FUNCS = ("_free_gb", "_require_gb")


def test_require_gb_passes_when_space_is_ample(tmp_path):
    out = _bash(tmp_path, [_extract_function(f) for f in SPACE_FUNCS],
                f'_require_gb "{tmp_path}" 0 "cache" && echo PASS || echo BLOCK')
    assert "PASS" in out.stdout, out.stderr


def test_require_gb_blocks_when_space_is_short(tmp_path):
    """A hard stop, so "ran out of space mid-unpack and left a broken cache"
    becomes a line in the log instead of a silent half-install."""
    out = _bash(tmp_path, [_extract_function(f) for f in SPACE_FUNCS],
                f'_require_gb "{tmp_path}" 999999 "cache" && echo PASS || echo BLOCK')
    assert "BLOCK" in out.stdout
    assert "999999" in out.stderr


def test_require_gb_resolves_a_directory_that_does_not_exist_yet(tmp_path):
    """On a FRESH mount $VEP_DATA/revel does not exist when the guard runs — the
    mkdir comes after it. df fails on a missing path, which made _free_gb return
    empty and the guard report "only 0 GB free", so REVEL was skipped on every
    first run and could never be built. Free space is a property of the volume,
    not of the leaf directory, so walk up to the nearest existing ancestor.
    """
    missing = tmp_path / "vep" / "revel"
    assert not missing.exists()
    out = _bash(tmp_path, [_extract_function(f) for f in SPACE_FUNCS],
                f'_require_gb "{missing}" 0 "REVEL" && echo PASS || echo BLOCK')
    assert "PASS" in out.stdout, (
        f"guard reported no space for a not-yet-created dir: {out.stdout}{out.stderr}"
    )


def test_free_gb_reports_the_same_volume_for_a_missing_child(tmp_path):
    """The walked-up answer must be the real volume figure, not a fallback 0."""
    funcs = [_extract_function("_free_gb")]
    here = _bash(tmp_path, funcs, f'_free_gb "{tmp_path}"').stdout.strip()
    below = _bash(tmp_path, funcs, f'_free_gb "{tmp_path}/a/b/c"').stdout.strip()
    assert here and here == below, f"{here!r} != {below!r}"
    assert int(here) > 0


def test_release_detection_survives_a_failing_binary(tmp_path):
    """A vep that exits non-zero must not kill the script through set -e/pipefail
    before the friendlier INSTALL.pl check can report anything."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "vep").write_text(
        "#!/usr/bin/env bash\necho boom >&2\nexit 3\n"
    )
    (bindir / "vep").chmod(0o755)
    out = _bash(tmp_path, [_extract_function("_vep_release")],
                'set -euo pipefail\nREL="$(_vep_release)"\necho "REACHED REL=$REL"',
                bindir=bindir, env="VEP_RELEASE=113")
    assert "REACHED REL=113" in out.stdout, f"{out.stdout}{out.stderr}"


def test_space_floors_match_the_spec():
    """The cache floor tracks whether the variation cache is being installed.

    28 GB was measured against the FULL merged tree (25.6 GB Content-Length,
    ftp.ensembl.org release-113); the floor has to clear the extracted tree, not
    the tarball, because the stream is piped straight into tar. The default
    install now excludes all_vars.gz — measured 17.46 GB of a 21.5 GB tree on
    2026-08-29 — so ~5 GB lands and 8 GB is the floor that matters.

    Both numbers live behind ONE switch on purpose. A flat 28 GB left here after
    the strip is the stale-comment-becomes-stale-floor bug the script's DISK
    block exists to prevent, and an 8 GB floor left here if the strip is ever
    reverted is the same bug pointing the other way."""
    assert 'CACHE_FLOOR_GB' in SCRIPT
    assert re.search(r'_require_gb "\$VEP_DATA" "\$CACHE_FLOOR_GB"', SCRIPT), (
        "the cache floor must come from CACHE_FLOOR_GB, not a literal"
    )
    assert re.search(r'CACHE_FLOOR_GB="\$\{CACHE_FLOOR_GB:-28\}"', SCRIPT), (
        "the KEEP_VARIATION branch must still require 28 GB"
    )
    assert re.search(r'CACHE_FLOOR_GB="\$\{CACHE_FLOOR_GB:-8\}"', SCRIPT), (
        "the stripped branch must require 8 GB"
    )
    assert re.search(r"_require_gb[^\n]*\s10\s", SCRIPT), "REVEL floor must be 10 GB"


def _write_fasta(tmp_path: Path) -> Path:
    """A readable stand-in reference. It has to EXIST: _vep_smoke attaches
    --fasta only when it can read one, and _check_cache_info defers its
    SIFT probe entirely when it cannot."""
    fa = tmp_path / "vep" / "ref.fa.gz"
    fa.parent.mkdir(parents=True, exist_ok=True)
    fa.write_text("")
    return fa


def _info_check(tmp_path, info: str | None, *, emit: str = VEP_JSON_SIFT,
                with_fasta: bool = True):
    """Run _check_cache_info against a given info.txt (None = absent).

    `emit` is what the stub vep returns, because the check is FUNCTIONAL: it
    runs the runtime's --sift/--polyphen flags rather than grepping info.txt for
    a field name."""
    _write_info(tmp_path, info)
    fa = _write_fasta(tmp_path) if with_fasta else (tmp_path / "vep" / "absent.fa.gz")
    bindir = _stub_vep(tmp_path, emit=emit)
    return _bash(tmp_path, [_extract_function("_vep_smoke"),
                            _extract_function("_check_cache_info")],
                 "_check_cache_info && echo OK || echo BAD", bindir=bindir,
                 env=f'FASTA_PATH="{fa}"')


def test_cache_check_fails_when_info_txt_is_absent(tmp_path):
    """⚠⚠ THE GAP THAT LET A BROKEN CACHE BE RECORDED COMPLETE.

    Measured 2026-08-29 on a real merged tree whose stream stopped in chr7:
    every chromosome directory present, every transcript block readable, and
    `vep` annotated MYH7 R403Q correctly — while the app's own invocation exited
    2 with `ERROR: SIFT not available` and zero bytes. Adding ONLY info.txt made
    the identical call return 3076 bytes with
    sift=deleterious_low_confidence / polyphen=probably_damaging.

    _vep_smoke passes neither --sift nor --polyphen, so _check_cache used to pass
    on exactly that tree and _state_write recorded the cache as done. The stub
    below returns a PERFECT answer, so the only thing that can fail this test is
    the info.txt gate."""
    out = _check(tmp_path, "cache", VEP_JSON_PLAIN, info=None)
    assert "INCOMPLETE" in out.stdout, (
        "a cache with no info.txt was recorded complete; that is the tree whose "
        "every curation exits 2 on 'SIFT not available'"
    )
    assert "info.txt" in out.stderr, "the failure must name the file, not SIFT"


def test_cache_check_passes_with_a_valid_info_txt(tmp_path):
    """The other direction — the gate must not reject a good cache."""
    out = _check(tmp_path, "cache", VEP_JSON_SIFT)
    assert "COMPLETE" in out.stdout, out.stderr


def test_info_check_runs_the_runtimes_sift_flags(tmp_path):
    """The check must ASK VEP, not read a field name.

    An earlier version of this check grepped info.txt for `source_sift`, a name
    derived from the Perl and then confirmed only against a stub info.txt
    written for that test. The real Ensembl file turned out to carry bare
    `sift`/`polyphen` keys too. Asserting on a schema seen only in one's own
    fixture is the same mistake as the `--format hgvs` probe and the
    case-sensitive "REVEL" grep, so this pins the behaviour that cannot be
    wrong about the schema: it runs the flags."""
    bindir, log = _argv_recording_vep(tmp_path)
    _write_info(tmp_path)
    fa = _write_fasta(tmp_path)
    _bash(tmp_path, [_extract_function("_vep_smoke"),
                     _extract_function("_check_cache_info")],
          "_check_cache_info || true", bindir=bindir, env=f'FASTA_PATH="{fa}"')
    argv = log.read_text()
    assert "--fasta" in argv, (
        "the merged cache is BAM-edited, so VEP enables --use_transcript_ref "
        "and THROWS without a FASTA (CacheDir.pm:400, Runner.pm:741). Every "
        f"probe must carry one once it exists: {argv}"
    )
    assert "--sift b" in argv and "--polyphen b" in argv, (
        f"the info check must exercise the runtime's own SIFT/PolyPhen flags: {argv}"
    )


def test_info_check_fails_when_sift_produces_nothing(tmp_path):
    """A cache that declares no SIFT — which is how Ensembl ships species
    without the data, e.g. ciona's empty `sift` line — must fail here rather
    than at the first curation. VEP throws instead of degrading unless
    --everything or REST 'safe' mode is set."""
    out = _info_check(tmp_path, VALID_INFO_TXT, emit=VEP_JSON_PLAIN)
    assert "BAD" in out.stdout, (
        "no sift_prediction came back and the check still passed"
    )
    assert "sift_prediction" in out.stderr


def test_the_probe_falls_back_to_use_given_ref_without_a_fasta(tmp_path):
    """⚠ THE INSTALL DEADLOCKS WITHOUT THIS, AND THE DEADLOCK IS INVISIBLE.

    The merged cache is BAM-edited, so CacheDir.pm:400 enables
    --use_transcript_ref and Runner.pm:741 throws unless a FASTA is supplied.
    Combined with the install order that is a cycle:

        step 1 installs the cache, _check_cache throws (no FASTA yet)
               -> the cache is NOT recorded
        step 2 refuses to install the FASTA — "Skipping FASTA: the cache is
               not usable yet" — because it gates on _state_ok cache

    so the cache can never be recorded and the FASTA is never installed. The
    live mount escaped it only because its FASTA predates the 2026-08-27 switch
    to the merged flavour; a fresh share would hang forever.

    --use_given_ref is VEP's own override and is right for a probe: PROBE_VCF
    states ref C at 14:23429278, which IS the GRCh38 reference, and SIFT and
    PolyPhen come from the cache rather than the reference."""
    bindir, log = _argv_recording_vep(tmp_path)
    _write_info(tmp_path)
    absent = tmp_path / "vep" / "not-installed-yet.fa.gz"
    _bash(tmp_path, [_extract_function("_vep_smoke"),
                     _extract_function("_check_cache_info")],
          "_check_cache_info || true", bindir=bindir, env=f'FASTA_PATH="{absent}"')
    argv = log.read_text()
    assert "--use_given_ref" in argv, (
        f"without a FASTA the probe must override use_transcript_ref: {argv}"
    )
    assert "--fasta" not in argv, "there is no FASTA to pass"
    assert "--sift b" in argv, "the probe must still exercise the runtime's flags"


def test_info_check_rejects_an_assembly_mismatch(tmp_path):
    """VEP cross-checks the assembly ONLY when info.txt supplies it
    (CacheDir.pm:397), so an absent or wrong value is how a GRCh37 tree gets
    annotated as GRCh38 — silently wrong coordinates, not an error."""
    out = _info_check(tmp_path, VALID_INFO_TXT.replace("GRCh38", "GRCh37"))
    assert "BAD" in out.stdout
    assert "GRCh37" in out.stderr


def test_the_install_strips_the_variation_cache_behind_one_switch():
    """all_vars.gz measured at 17.46 GB of a 21.5 GB tree (2026-08-29), and
    CacheDir.pm:167 only builds the Variation annotation source inside
    `if($self->param('check_existing') && $info->{variation_cols})`. So it is
    ~81% of the install that is never opened."""
    assert "all_vars.gz*" in SCRIPT, "the tar exclude for the variation cache is gone"
    assert "VEP_CACHE_KEEP_VARIATION" in SCRIPT, (
        "there must be a documented way back to the full tree, for the day a "
        "--check_existing-class flag is added"
    )


def test_non_empty_directory_is_no_longer_a_completion_check():
    """Regression guard on the bug this whole change exists to fix."""
    assert 'ls -A "$VEP_DATA/$VEP_SPECIES"' not in SCRIPT, (
        "A non-empty species directory must not count as a finished install — an "
        "unpack interrupted by a deploy restart leaves exactly that."
    )


def test_builder_image_has_the_vep_toolchain():
    """The whole reason setup_offline_vep.sh could not live in build_all.sh."""
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert "ensemblorg/ensembl-vep" in dockerfile, (
        "the builder image must carry INSTALL.pl and vep; python:3.11-slim has "
        "neither, which is why the VEP install needed a separate manual job."
    )
    assert not re.search(r"^FROM python:3\.11-slim", dockerfile, re.MULTILINE)


def test_builder_installs_python_311():
    """The VEP base is Ubuntu 22.04 (Python 3.10) and pinned numpy needs >= 3.11."""
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert "deadsnakes" in dockerfile
    assert "/opt/venv" in dockerfile


def test_builder_has_the_revel_build_tools():
    """_build_revel_data.sh needs unzip; htslib comes with the VEP image."""
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert "unzip" in dockerfile


def test_builder_runs_as_root_to_write_the_mount():
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert not re.search(r"^USER (?!root)", dockerfile, re.MULTILINE), (
        "the builder writes to the /app/data mount; the runtime image drops "
        "privileges because it only reads, but this one must not."
    )


def test_builder_points_build_all_at_the_venv_python():
    """build_all.sh defaults PY to python3, which on this base is 3.10."""
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert "PYTHON=/opt/venv/bin/python" in dockerfile


def test_builder_clears_the_vep_entrypoint():
    """The official VEP image's ENTRYPOINT launches vep; the builder must run
    entrypoint_builder.sh instead."""
    dockerfile = (ROOT / "Dockerfile.builder").read_text()
    assert 'ENTRYPOINT ["bash", "scripts/entrypoint_builder.sh"]' in dockerfile


def test_both_images_take_the_vep_release_as_a_build_arg():
    for name in ("Dockerfile", "Dockerfile.builder"):
        text = (ROOT / name).read_text()
        assert re.search(r"^ARG VEP_RELEASE=", text, re.MULTILINE), name
        assert re.search(r"^FROM ensemblorg/ensembl-vep:\$\{VEP_RELEASE\}", text,
                         re.MULTILINE), name


def test_workflow_passes_one_release_to_both_builds():
    wf = require_file(".github/workflows/push-image.yml").read_text()
    assert wf.count("VEP_RELEASE=${{ env.VEP_RELEASE }}") == 2, (
        "both image builds must take the SAME release from one env value"
    )
    assert re.search(r"^\s*VEP_RELEASE:\s*release_\d+", wf, re.MULTILINE)


WORKFLOW = ROOT / ".github/workflows/restart-webapp.yml"
DATA_BUILD_WORKFLOW = ROOT / ".github/workflows/data-build.yml"


def test_vep_paths_are_set_as_app_settings():
    """Setting the paths while the flag is off is inert (_config returns early),
    and it makes the eventual switch-over a one-setting change instead of five."""
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    for var in ("HEARTVAR_VEP_DATA", "HEARTVAR_VEP_ASSEMBLY", "HEARTVAR_VEP_FASTA",
                "HEARTVAR_VEP_REVEL", "HEARTVAR_VEP_PLUGINS_DIR"):
        assert re.search(rf"^\s+{var}=", wf, re.MULTILINE), f"{var} is not set"


def test_the_offline_flag_is_on_and_the_resolver_with_it():
    """⚠ THIS TEST WAS THE GATE, AND THE GATE IS NOW SATISFIED.

    It previously asserted HEARTVAR_VEP_OFFLINE was UNSET, with the condition
    for flipping written into its own docstring: "the flag stays unset until the
    runtime call itself is demonstrated". Keeping the history, because it is the
    reason the condition was worded that way:

    2026-08-26 — the flag was turned on and turned OFF the same day. `vep
    --offline` exited 2 on EVERY curation in production. Not slow, broken. The
    REST fallback covered for it, so the only visible symptom was curations
    taking ~64 s instead of ~23 s and the error rate doubling. Every build-time
    check was GREEN throughout, and none was lying: CHECK_ONLY proves the CACHE
    is annotatable by the vep binary in the BUILDER image. It does not prove the
    WEB APP's runtime invocation works — different process, different image,
    different environment — and that gap is exactly where it hid.

    WHAT CLOSED IT, 2026-08-29:
      * setup_offline_vep.sh::_check_runtime_call reconstructs vep_offline._argv
        flag for flag and PASSES ON THE AZURE MOUNT — cache / fasta / revel /
        runtime all functional=yes, manifest written, FAILED (0).
      * Full HTTP curations through the running app returned IDENTICAL
        classification, points, criteria met and every VEP field to REST for
        MYH7 c.1208G>A, MYBPC3 c.1504C>T and ACTC1 c.301G>A — all minus-strand,
        all confirming resolved_locally=True, i.e. genuinely annotated locally.
      * The parity gate ran BOTH paths with fetch_region/fetch_hgvs
        called directly so a fallback counts as a fallback and never as
        agreement. Coordinate path: 2 residuals, both attributed — LZTR1 REVEL
        where OFFLINE IS BETTER (REST's dbNSFP returns no data), and MAP2K1
        most_severe from ENST00001072373, an NMD-biotype transcript that exists
        only in release 116.
      * The nmd_escape gap is closed (exon lengths from cdot), so offline no
        longer over-calls PVS1 for a PTC in the last 50 nt of a penultimate exon.

    THE RESOLVER GOES ON WITH IT, and that pairing is asserted rather than left
    to judgement: `vep --offline` refuses --format hgvs outright
    (Parser/HGVS.pm:104, unconditional), so with the resolver off every
    HGVS-entered curation — the clinical norm — would still go to REST and the
    flag would buy almost nothing.

    REST remains the fallback on every failure path, so the worst case is the
    behaviour that preceded this. Reverting is one line in the workflow.
    """
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    assert re.search(r'^\s+HEARTVAR_VEP_OFFLINE="1"', wf, re.MULTILINE), (
        "offline VEP is verified on the mount and in-app; the deploy must set "
        "the flag or the cache is inert"
    )
    assert re.search(r'^\s+HEARTVAR_HGVS_RESOLVER="1"', wf, re.MULTILINE), (
        "`vep --offline` cannot parse HGVS, so without the resolver the offline "
        "flag leaves the clinical-norm input path on REST"
    )
    assert re.search(r"if _check_revel; then\s*\n\s*_state_write revel", SCRIPT), (
        "REVEL must be recorded ONLY after a score verifiably comes through"
    )


def test_the_resolver_is_useless_without_the_cdot_database():
    """The resolver reads HEARTVAR_CDOT_DB; without it hgvs_resolver falls back
    to the 4.5 GB in-memory provider or returns None. Turning the resolver on
    while that path is unset would be a silent REST fallback on every HGVS
    curation — the failure this whole sequence exists to stop."""
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    assert 'HEARTVAR_CDOT_DB="/app/data/cdot/cdot_transcripts.db"' in wf


def test_settings_paths_agree_with_the_setup_script():
    """A path set here that the builder never writes is a silent REST fallback."""
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    assert 'HEARTVAR_VEP_DATA="/app/data/vep"' in wf
    assert "/app/data/vep/revel/new_tabbed_revel_grch38.tsv.gz" in wf
    assert 'HEARTVAR_VEP_PLUGINS_DIR="/app/data/vep/Plugins"' in wf


def test_fasta_setting_matches_the_release_in_push_image():
    """The FASTA path embeds the VEP release, so it must track VEP_RELEASE. A
    stale release here points --fasta at a file that does not exist, which costs
    HGVS silently rather than erroring."""
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    push = require_file(".github/workflows/push-image.yml").read_text()
    rel = re.search(r"^\s*VEP_RELEASE:\s*release_(\d+)", push, re.MULTILINE)
    assert rel, "push-image.yml has no VEP_RELEASE"
    assert f"/{rel.group(1)}_GRCh38/" in wf, (
        f"HEARTVAR_VEP_FASTA must use release {rel.group(1)}"
    )


def test_scripts_readme_no_longer_says_vep_is_unwired():
    text = (ROOT / "scripts/README.md").read_text()
    assert "**NOT wired into the backend**" not in text


def test_setup_script_header_documents_the_measured_sizes():
    """The old header's "~3 GB" REVEL figure is what made a single 33 GB floor look
    sufficient. Wrong numbers in a comment become wrong numbers in a guard."""
    require_unstripped_prose(SCRIPT, "25.6 GB tarball")
    header = SCRIPT[:SCRIPT.index("set -euo pipefail")]
    assert "25.6 GB tarball" in header, "the header must spec the MERGED cache size"
    assert "24.9 GB tarball" not in header, (
        "24.9 GB is the superseded Ensembl-only cache; the spec line must be 25.6"
    )
    assert "667 MB" in header


def test_the_single_global_disk_floor_is_gone():
    assert "MIN_FREE_GB" not in SCRIPT, (
        "one global floor is what allowed 5 GB immediately before a step needing "
        "~25 GB; the floors are per-step now"
    )


def _run_script(tmp_path: Path, *, cache_url: str, extra_env: dict | None = None):
    bindir = _stub_vep(tmp_path, emit="")
    (bindir / "INSTALL.pl").write_text("#!/usr/bin/env bash\nexit 0\n")
    (bindir / "INSTALL.pl").chmod(0o755)
    env = {
        "PATH": f"{bindir}:{os.environ['PATH']}",
        "VEP_DATA": str(tmp_path / "vep"),
        "CACHE_URL": cache_url,
        "SKIP_REVEL": "1",
    }
    env.update(extra_env or {})
    return subprocess.run(["bash", str(SETUP)], capture_output=True, text=True,
                          timeout=120, env=env)


def test_a_failed_download_exits_one_with_a_summary_not_a_raw_tar_error(tmp_path):
    res = _run_script(tmp_path, cache_url="file:///nonexistent/cache.tar.gz")
    assert res.returncode == 1, (
        f"expected a clean exit 1, got {res.returncode}\n{res.stdout}\n{res.stderr}"
    )
    combined = res.stdout + res.stderr
    assert "No manifest written" in combined, combined
    assert "failed component" in combined, combined


def test_no_manifest_is_written_when_verification_fails(tmp_path):
    """The half-install guard: files may exist, but if vep cannot annotate then
    build_all.sh must NOT see a manifest, or its `versioned` KEEP would skip the
    step forever in a broken state."""
    _run_script(tmp_path, cache_url="file:///nonexistent/cache.tar.gz")
    assert not (tmp_path / "vep" / ".heartvar_vep_manifest.json").exists()
    state = tmp_path / "vep" / ".heartvar_vep_state"
    assert not (state / "cache.json").exists(), (
        "cache must not be recorded complete when vep cannot annotate"
    )


def test_missing_toolchain_is_reported_clearly(tmp_path):
    """Without INSTALL.pl the script must say so rather than fail obscurely."""
    res = subprocess.run(
        ["bash", str(SETUP)], capture_output=True, text=True, timeout=60,
        env={"PATH": "/usr/bin:/bin", "VEP_DATA": str(tmp_path / "vep")},
    )
    assert res.returncode == 1
    assert "INSTALL.pl is not on PATH" in res.stderr, res.stderr


def _sabotage(bindir: Path) -> None:
    """Replace every installer with one that fails loudly if it is ever called."""
    for tool in ("curl", "INSTALL.pl"):
        path = bindir / tool
        path.write_text(
            "#!/usr/bin/env bash\n"
            f'echo "SABOTAGE: {tool} was invoked" >&2\nexit 99\n'
        )
        path.chmod(0o755)


def _stub_vep_answering(tmp_path: Path) -> Path:
    """A vep that answers each probe correctly, so every check passes —
    including the runtime-call probe, which is the only one that asks for the
    RefSeq protein HGVS and so the only one that would have caught the
    2026-08-26 flavour bug."""
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    vep = bindir / "vep"
    vep.write_text(
        "#!/usr/bin/env bash\n"
        "hgvs=0; revel=0; runtime=0\n"
        'for a in "$@"; do\n'
        '  [[ "$a" == "--help" ]] && { echo "ensembl-vep : 113"; exit 0; }\n'
        '  [[ "$a" == "--hgvs" ]] && hgvs=1\n'
        '  [[ "$a" == REVEL,* ]] && revel=1\n'
        '  [[ "$a" == "--sift" ]] && runtime=1\n'
        "done\n"
        "cat >/dev/null\n"
        f"if (( runtime )); then printf '%s\\n' {json.dumps(VEP_JSON_RUNTIME)}\n"
        f"elif (( revel )); then printf '%s\\n' {json.dumps(VEP_JSON_REVEL)}\n"
        f"elif (( hgvs )); then printf '%s\\n' {json.dumps(VEP_JSON_HGVS)}\n"
        f"else printf '%s\\n' {json.dumps(VEP_JSON_PLAIN)}; fi\n"
    )
    vep.chmod(0o755)
    return bindir


def _write_state(tmp_path: Path, *components: str) -> None:
    _write_info(tmp_path)

    state = tmp_path / "vep" / ".heartvar_vep_state"
    state.mkdir(parents=True, exist_ok=True)
    for comp in components:
        (state / f"{comp}.json").write_text(json.dumps({
            "component": comp, "vep_release": "113", "assembly": "GRCh38",
        }, indent=2))


def _check_only(tmp_path: Path, bindir: Path):
    return subprocess.run(
        ["bash", str(SETUP)], capture_output=True, text=True, timeout=120,
        env={
            "PATH": f"{bindir}:{os.environ['PATH']}",
            "VEP_DATA": str(tmp_path / "vep"),
            "CHECK_ONLY": "1",
        },
    )


def test_check_only_downloads_nothing(tmp_path):
    """The whole point: a probe that could start a 25 GB stream is not a probe."""
    bindir = _stub_vep(tmp_path, emit="")
    _sabotage(bindir)
    res = _check_only(tmp_path, bindir)
    combined = res.stdout + res.stderr
    assert "SABOTAGE" not in combined, combined
    assert "CHECK ONLY" in combined, combined
    assert res.returncode == 1, combined


def test_check_only_names_every_component_and_its_state(tmp_path):
    """'state=ABSENT' and 'functional=NO' are different failures — a component
    recorded complete but no longer working means the mount changed underneath
    it, which is not the same problem as one that never finished."""
    bindir = _stub_vep(tmp_path, emit="")
    _sabotage(bindir)
    out = _check_only(tmp_path, bindir).stdout
    for comp in ("cache", "fasta", "revel"):
        assert comp in out, out
    assert "state=ABSENT" in out, out
    assert "functional=NO" in out, out


def test_check_only_exits_zero_when_every_component_verifies(tmp_path):
    """Exit status matches a real run, so build_all.sh's per-source summary does
    not report a working mirror as failed."""
    bindir = _stub_vep_answering(tmp_path)
    _sabotage(bindir)
    _write_state(tmp_path, "cache", "fasta", "revel")
    res = _check_only(tmp_path, bindir)
    combined = res.stdout + res.stderr
    assert res.returncode == 0, combined
    assert "SABOTAGE" not in combined, combined
    assert "functional=NO" not in combined, combined


def test_check_only_reports_verified_components_with_no_manifest(tmp_path):
    """The exact state the mount was left in on 2026-08-25: components verify
    individually and the all-or-nothing manifest is absent, which is what
    build_all.sh reports as `FAILED vep`. Saying so is the diagnosis."""
    bindir = _stub_vep_answering(tmp_path)
    _sabotage(bindir)
    _write_state(tmp_path, "cache", "fasta", "revel")
    out = _check_only(tmp_path, bindir).stdout
    assert "NO MANIFEST" in out, out


def test_check_only_does_not_need_install_pl(tmp_path):
    """A diagnostic must not be blocked by the toolchain it does not use — that
    guard is for the installer."""
    bindir = _stub_vep(tmp_path, emit="")
    res = subprocess.run(
        ["bash", str(SETUP)], capture_output=True, text=True, timeout=120,
        env={"PATH": f"{bindir}:/usr/bin:/bin", "VEP_DATA": str(tmp_path / "vep"),
             "CHECK_ONLY": "1"},
    )
    assert "INSTALL.pl is not on PATH" not in res.stderr, res.stderr
    assert "CHECK ONLY" in res.stdout, res.stdout


def test_the_runtime_never_asks_for_a_cache_FLAVOUR_the_build_does_not_install():
    """THE 2026-08-26 EXIT-2 BUG, pinned at the boundary it crossed.

    `vep --offline` exited 2 on EVERY curation in production while every
    build-time check stayed green. Neither was lying, and the cache was never the
    problem:

      * The build installed `homo_sapiens_vep_<REL>_<ASM>.tar.gz` — the
        ENSEMBL-only cache — into `$VEP_DATA/homo_sapiens/<REL>_<ASM>`.
      * The runtime appended `--refseq` whenever the curator's transcript was a
        RefSeq accession (`_is_refseq("NM_...")`), which for this corpus is EVERY
        variant.
      * VEP 113 `CacheDir.pm` does
            $species_dir_name .= '_'.$_ for grep { $self->param($_) } qw(refseq merged);
        then `throw("ERROR: Cache directory $dir not found\\n")`.

    Reproduced locally in `ensemblorg/ensembl-vep:release_113.0`, exit 2:
        MSG: ERROR: Cache directory /app/data/vep/homo_sapiens_refseq not found
    Note there is no release subdirectory in that path — VEP fails at the species
    directory, before it ever looks at a cache version.

    Every build-time probe ran WITHOUT the flavour flag, so not one exercised the
    thing the runtime added. That is the whole gap: a verified cache is not a
    verified path, because what fails is THE CALL.

    This test is the cross-file assertion that closes it — the flavour the
    runtime asks for must be the flavour the build installs.
    """
    import backend.clients.vep_offline as vo

    m = re.search(r'VEP_CACHE_SPECIES="\$\{VEP_CACHE_SPECIES:-([a-z_]+)\}"', SCRIPT)
    assert m, "could not determine which cache species the build installs"
    installed = m.group(1)
    assert "${VEP_CACHE_SPECIES}_vep_${REL}_${VEP_ASSEMBLY}.tar.gz" in SCRIPT, (
        "the cache tarball URL no longer derives from VEP_CACHE_SPECIES — "
        "re-derive this test against whatever now decides the flavour"
    )
    flag = re.search(r'VEP_CACHE_FLAVOUR_FLAG="\$\{VEP_CACHE_FLAVOUR_FLAG:-(--[a-z]+)\}"',
                     SCRIPT)
    assert flag, "the build no longer pins the flavour flag in one variable"

    cfg = {"binary": "vep", "data_dir": "/app/data/vep", "fasta": None,
           "assembly": "GRCh38", "revel": None, "plugins_dir": None,
           "extra_args": [], "timeout": 45.0}

    for refseq in (True, False):
        argv = vo._argv(cfg, "vcf", refseq=refseq)
        asked = [a for a in argv if a in ("--refseq", "--merged")]
        assert asked, (
            "the runtime passes no cache-flavour flag, so VEP will look for the "
            f"unsuffixed species dir, but the build installs {installed!r} — "
            "which is a missing-directory exit 2 on every call"
        )
        assert len(asked) == 1, f"--refseq and --merged are mutually exclusive: {asked}"
        assert asked[0] == flag.group(1), (
            f"runtime passes {asked[0]} but the build installs the "
            f"{flag.group(1)} flavour; VEP derives the cache directory from the "
            "flag, so these cannot differ"
        )
        wanted = "homo_sapiens" + asked[0].removeprefix("-" * 2).join(("_", ""))[:0] + \
                 "_" + asked[0].lstrip("-")
        assert wanted == installed, (
            f"runtime passes {asked[0]}, so VEP looks for {wanted!r} under "
            f"--dir_cache, but the build installs {installed!r}"
        )


def test_the_build_probes_the_RUNTIME_call_not_just_the_cache():
    """The lesson of 2026-08-26, as a guard.

    CHECK_ONLY proved the CACHE was annotatable — by the vep binary in the
    BUILDER image, from the mount, at build time — and said nothing about the web
    app's invocation. 8 of the runtime's 20 flags and its input format appeared in
    no check at all, and the one that mattered (the cache flavour) took a day and
    five attempts to find. So the build now runs the app's own call, and the
    manifest is gated on it.
    """
    assert "_check_runtime_call()" in SCRIPT, (
        "the build no longer probes the runtime call; without it a broken "
        "offline path passes every component check, which is the 2026-08-26 bug"
    )
    for flag in ("--sift b", "--polyphen b", "--mane",
                 "--symbol", "--numbers", "--canonical", "--biotype"):
        assert flag in SCRIPT, f"the runtime probe does not exercise {flag}"
    assert re.search(r"if _check_runtime_call; then", SCRIPT), \
        "the runtime check must run in the installer, not only under CHECK_ONLY"
    assert "THE RUNTIME CALL FAILS even though every component verified" in SCRIPT, \
        "a runtime-call failure must be loud — it is the invisible one"


def test_the_runtime_probe_uses_a_format_vep_can_actually_serve_offline():
    """THE 2026-08-28 BUG: the manifest gate could never pass.

    `_check_runtime_call` gates the manifest, and between 8e3f6d1 and 2026-08-28
    it ran

        vep --offline ... --format hgvs   (input NM_000257:c.1208G>A)

    and required `NP_000248.2:p.Arg403Gln` back. VEP 113 refuses that
    combination unconditionally:

        Parser/HGVS.pm:104
        throw("ERROR: Cannot use HGVS format in offline mode")
            if $self->param('offline');

    Not gated on --cache, --fasta, --dir_cache or the cache flavour. So a
    flawless ~26 GB merged-cache install would still have ended `FAILED vep`
    with no manifest — and because build_all.sh's `versioned` KEEP is
    conditional on that manifest, every later run would re-enter the step.

    WHY NOTHING CAUGHT IT, which is the part worth keeping. 8e3f6d1 found this
    same VEP limit and fixed the CLIENT for it — vep_offline.fetch_hgvs returns
    None BEFORE spawning — but left the probe reconstructing the call the client
    had just stopped making. And the probe cannot be caught by the suite's stub
    `vep`, which happily answers any argv, so the test that asserted
    `"--format hgvs" in SCRIPT` was actively pinning the bug in place.

    This test asserts the property that actually matters: the probe must send the
    format the runtime sends, and must not send the one VEP refuses.
    """
    import re as _re
    import backend.clients.vep_offline as vo

    body = _re.search(r"_check_runtime_call\(\) \{(.*?)\n\}", SCRIPT, _re.S)
    assert body, "could not isolate _check_runtime_call"
    body = body.group(1)

    assert "--format vcf" in body, (
        "the runtime probe must use --format vcf: fetch_region is the ONLY "
        "runtime entry point offline VEP can serve, and it sends vcf"
    )
    assert "--format hgvs" not in body, (
        "the runtime probe asks for --format hgvs, which VEP refuses in offline "
        "mode (Parser/HGVS.pm:104). This check gates the manifest, so it would "
        "fail on every install no matter how complete the cache is"
    )

    import asyncio
    assert asyncio.run(vo.fetch_hgvs("MYH7", "c.1208G>A", "NM_000257.4")) is None, (
        "fetch_hgvs now attempts offline VEP; if HGVS input became servable "
        "offline, revisit the runtime probe's format too"
    )


def test_the_runtime_probe_proves_the_cache_is_MERGED_not_just_present():
    """The flavour flag resolving a directory proves the NAME is right; it does
    not prove the CONTENTS are. A merged tree carries both transcript sets, so
    the probe requires a RefSeq transcript in the raw response alongside the
    Ensembl protein the coordinate path keeps. Without the RefSeq assertion the
    probe would pass against a plain cache that had merely been renamed."""
    body = re.search(r"_check_runtime_call\(\) \{(.*?)\n\}", SCRIPT, re.S)
    assert body, "could not isolate _check_runtime_call"
    body = body.group(1)
    assert "PROBE_RUNTIME_MERGED_TX_RE" in body, (
        "the runtime probe no longer requires a RefSeq transcript, so it cannot "
        "tell a merged cache from a plain one"
    )
    m = re.search(r"PROBE_RUNTIME_MERGED_TX_RE='([^']+)'", SCRIPT)
    assert m, "the merged-cache proof is no longer a single pinned pattern"
    pattern = m.group(1)
    assert "transcript_id" in pattern, (
        "the merged-cache proof is not anchored on the transcript_id field, so "
        '"mane_select":"NM_..." on an Ensembl row satisfies it and a plain '
        "cache passes: " + pattern
    )
    merged = ('{"transcript_consequences":[{"transcript_id":"ENST00000355349.4",'
              '"mane_select":"NM_000257.4"},{"transcript_id":"NM_000257.4"}]}')
    plain = ('{"transcript_consequences":[{"transcript_id":"ENST00000355349.4",'
             '"mane_select":"NM_000257.4"}]}')
    assert subprocess.run(["grep", "-qE", pattern], input=merged,
                          text=True).returncode == 0, \
        "the pattern does not match a genuinely merged response"
    assert subprocess.run(["grep", "-qE", pattern], input=plain,
                          text=True).returncode != 0, \
        ("the pattern matches an ENSEMBL-ONLY response — it is being satisfied "
         "by the mane_select cross-reference, not by a RefSeq transcript")
    assert "PROBE_RUNTIME_ENSP" in body, (
        "the runtime probe no longer requires the Ensembl protein, which is the "
        "set _to_rest_transcript_set keeps for coordinate input"
    )


def test_the_fasta_path_does_not_follow_the_cache_flavour():
    """INSTALL.pl --AUTO f publishes the reference under the UNSUFFIXED species
    dir, and HEARTVAR_VEP_FASTA in restart-webapp.yml points there. When the
    cache moved to the merged flavour, deriving FASTA_PATH from CACHE_DIR would
    have sent the FASTA check hunting under homo_sapiens_merged/ and failed a
    correctly-installed reference — the same false-negative shape as the
    case-sensitive REVEL grep."""
    assert 'FASTA_PATH="$FASTA_DIR/' in SCRIPT, \
        "FASTA_PATH must derive from FASTA_DIR, not from the cache dir"
    assert 'FASTA_DIR="$VEP_DATA/$VEP_FASTA_SPECIES/' in SCRIPT
    assert '--SPECIES "$VEP_FASTA_SPECIES"' in SCRIPT, \
        "INSTALL.pl --AUTO f must be given the unsuffixed species"
    wf = require_file(".github/workflows/restart-webapp.yml").read_text()
    assert "/app/data/vep/homo_sapiens/" in wf, (
        "HEARTVAR_VEP_FASTA must still point at the unsuffixed species dir"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_the_data_build_workflow_can_override_the_once_a_month_guard():
    """entrypoint_builder.sh refreshes at most once per calendar month, decided
    from build_stamp.json. A MANUAL dispatch is by definition an off-cycle
    request, and on 2026-08-28 the guard swallowed one whole and reported
    success. FORCE_BUILD is how the workflow says "I meant it", so it has to be
    reachable from the dispatch form — it was not."""
    wf = require_file(".github/workflows/data-build.yml").read_text(encoding="utf-8")
    assert "force_build:" in wf, (
        "the data-build workflow cannot set FORCE_BUILD, so any dispatch in a "
        "month that already refreshed is a silent no-op reported as success"
    )
    assert "FORCE_BUILD=1" in wf, "force_build is an input but never reaches the container"


def test_the_data_build_workflow_can_force_a_vep_reinstall():
    """build_all.sh's `versioned` KEEP skips the whole vep step when
    .heartvar_vep_manifest.json merely EXISTS, and setup_offline_vep.sh skips a
    component its state file calls complete. The manifest written 2026-08-26
    records the WRONG cache flavour (Ensembl-only, not merged), so without a way
    to force it the corrected install can never start."""
    wf = require_file(".github/workflows/data-build.yml").read_text(encoding="utf-8")
    assert "force_vep:" in wf, (
        "no way to reinstall the VEP cache over a stale manifest, so a wrong "
        "install is permanent"
    )
    assert "FORCE_VEP=1" in wf, (
        "force_vep is an input but never reaches the container as FORCE_VEP; "
        "without it build_all.sh KEEPs the vep step and setup_offline_vep.sh "
        "never runs at all"
    )
