#!/usr/bin/env python3
"""Harvest per-gene ACMG criteria specifications from the ClinGen CSpec Registry.

WHY THIS EXISTS. ``backend/data/vcep_criteria_applicability.json`` was
hand-transcribed. Reading two specifications directly on 2026-08-25 found three
classes of defect in it, each capable of producing a wrong clinical call:

  * MERGED ROWS. FBN1's PM1 arrived as one string with a "|" in it, silently
    joining two DIFFERENT strengths — PM1_strong ("Cysteine residues in
    cbEGF-like domains") and PM1 Moderate (the Cys-in-EGF-like/TB/hybrid list).
    Implementing from that copy would have scored the Strong rule at Moderate.
  * MISSING ROWS. FBN1's PM4 ("None (at any strength level)") and PM5 ("None")
    are absent entirely, so `_criterion_applicable` defaults them OPEN and both
    can fire on FBN1 today.
  * NO STRENGTH DIMENSION. The file records applicability only, so a spec that
    says "PM1 at Strong here, Moderate there" cannot be represented at all.

The registry serves each rule set as server-rendered HTML at
``/cspec/ui/svi/doc/<GN>``; there is no JSON API (every plausible /api path
returns 400). So this parses HTML — which is exactly the "prose parsing is where
a silent mis-parse becomes a clinical over-call" risk the audits flagged. Two
mitigations, both mandatory:

  1. ``--verify`` checks the output against values READ BY HAND from the source
     PDFs. If the registry changes its markup, the assertions fail loudly rather
     than emitting a plausible-looking but wrong table.
  2. Nothing is inferred. A criterion with no strength row and a "Not Applicable"
     marker is recorded not_applicable; a criterion the page does not describe is
     recorded "unknown", never "applicable".

Usage:
    python3 scripts/build_vcep_criteria.py --verify          # fetch + self-check
    python3 scripts/build_vcep_criteria.py -o out.json       # write the table
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time

BASE = "https://cspec.genome.network/cspec/ui/svi/doc/"

CODES = (
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)

STRENGTHS = ("Stand Alone", "Very Strong", "Strong", "Moderate", "Supporting")

GN_GENES: dict[str, tuple[str, ...]] = {
    "GN002": ("MYH7",),     "GN022": ("FBN1",),     "GN038": ("SHOC2",),
    "GN039": ("NRAS",),     "GN040": ("RAF1",),     "GN041": ("SOS1",),
    "GN042": ("SOS2",),     "GN043": ("PTPN11",),   "GN044": ("KRAS",),
    "GN045": ("MAP2K1",),   "GN046": ("HRAS",),     "GN047": ("RIT1",),
    "GN048": ("MAP2K2",),   "GN049": ("BRAF",),     "GN087": ("MRAS",),
    "GN094": ("LZTR1",),    "GN095": ("MYBPC3",),   "GN098": ("TNNI3",),
    "GN099": ("TNNT2",),    "GN100": ("TPM1",),     "GN101": ("ACTC1",),
    "GN102": ("MYL2",),     "GN103": ("MYL3",),     "GN112": ("KCNQ1",),
    "GN127": ("RRAS2",),
    "GN128": ("PPP1CB",),
    "GN125": ("BMPR2",),
}


def _clean(html: str) -> str:
    txt = re.sub(r"<[^>]+>", " ", html)
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
                 ("&nbsp;", " "), ("&quot;", '"'), ("&#39;", "'")):
        txt = txt.replace(a, b)
    return re.sub(r"\s+", " ", txt).strip()


_WITHDRAWN_RE = re.compile(r"^\W*not\s+applicable\b", re.I)


def _row_withdraws(row_html: str) -> bool:
    """True when this strength row's FIRST bullet withdraws the criterion."""
    first_li = re.search(r"<li[^>]*>(.*?)</li>", row_html, re.S)
    return bool(first_li and _WITHDRAWN_RE.match(_clean(first_li.group(1))))


def fetch(gn: str, retries: int = 3) -> str:
    """Fetch one rule set via curl.

    curl, not urllib: on this network the TLS chain carries a self-signed
    certificate (corporate proxy) and urllib fails CERTIFICATE_VERIFY_FAILED,
    while curl succeeds through the macOS keychain. Verification stays ON — this
    is not an --insecure shortcut.
    """
    url = f"{BASE}{gn}"
    for attempt in range(retries):
        proc = subprocess.run(
            ["curl", "-sSL", "--fail", "--max-time", "90",
             "-A", "HeartVar/1.0", url],
            capture_output=True, text=True,
        )
        if proc.returncode == 0 and proc.stdout:
            return proc.stdout
        if attempt == retries - 1:
            raise SystemExit(
                f"{gn}: fetch failed (curl exit {proc.returncode}): "
                f"{proc.stderr.strip()[:300]}"
            )
        print(f"  retry {gn} after curl exit {proc.returncode}", file=sys.stderr)
        time.sleep(2 * (attempt + 1))
    raise AssertionError("unreachable")


def _criterion_blocks(html: str) -> dict[str, str]:
    """Slice the page into one HTML block per criterion code.

    Each code appears MORE THAN ONCE in the document (parent row plus the
    goto-dropdown), so anchor on the parent <tr> specifically and bound each
    block by the next parent in document order. Slicing between two occurrences
    of the SAME code yields the wrong span and silently empty results.
    """
    starts: list[tuple[int, str]] = []
    for code in CODES:
        m = re.search(
            rf'<tr class="criteria-code parent[^"]*"[^>]*data-cspec-cc-label="{code}"',
            html,
        )
        if m:
            starts.append((m.start(), code))
    starts.sort()
    out: dict[str, str] = {}
    for i, (pos, code) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(html)
        out[code] = html[pos:end]
    return out


def parse_criterion(block: str, gene_symbol: str | None = None) -> dict:
    """One criterion -> {applicability, strengths{...}, acmg_summary}.

    Parsed per ROW so each strength keeps its own text. Collecting all labels
    and all texts separately and zipping them misaligns as soon as one strength
    has no text, which is the common case.
    """
    strengths: dict[str, str] = {}
    withdrawn: dict[str, bool] = {}
    for row in re.findall(r'<div class="row strength[^"]*">(.*?)(?=<div class="row |\Z)',
                          block, re.S):
        lab = re.search(r'strength-label[^>]*>(.*?)</div>', row, re.S)
        txt = re.search(r'strength-text[^>]*>(.*?)</div>', row, re.S)
        if not lab:
            continue
        label = _clean(lab.group(1))
        raw = txt.group(1) if txt else ""
        body = _clean(raw)
        if label in STRENGTHS and body:
            strengths[label] = body
            withdrawn[label] = _row_withdraws(raw)

    summary = ""
    for row in re.findall(r'<div class="row summary">(.*?)(?=<div class="row |\Z)',
                          block, re.S):
        t = re.search(r'strength-text[^>]*>(.*?)</div>', row, re.S)
        if t:
            summary = _clean(t.group(1))
            break

    specs = comments = ""
    for _label, _key in (("VCEP Specifications", "specs"), ("Comments", "comments")):
        m = re.search(
            _label + r'\s*:\s*</?span[^>]*>\s*<span class="col[^"]*value">(.*?)</span>',
            block, re.S,
        )
        if not m:
            m = re.search(
                _label + r'\s*:.{0,200}?class="col[^"]*value">(.*?)</span>',
                block, re.S,
            )
        if m:
            if _key == "specs":
                specs = _clean(m.group(1))
            else:
                comments = _clean(m.group(1))

    not_applicable = bool(
        re.search(r'class="[^"]*\bna\b[^"]*">\s*Not Applicable', block)
    )
    _live_row = any(
        lab for lab in (strengths or {}) if not withdrawn.get(lab)
    )
    _cw = ""
    for _txt in (comments, specs):
        _t = (_txt or "").strip()
        if not re.match(r"(?i)^not applicable\b", _t):
            continue
        if re.match(r"(?i)^not applicable\s+because\b", _t):
            if _live_row:
                continue
            _cw = _t
            break
        _m = re.match(r"(?i)^not applicable\s+for\s+([A-Z0-9]+)\s*\.?\s*$", _t)
        if _m and _m.group(1).upper() == (gene_symbol or "").upper():
            _cw = _t
            break
    if _cw:
        return {
            "applicability": "not_applicable",
            "strengths": strengths,
            "acmg_summary": summary,
            **({"vcep_specifications": specs} if specs else {}),
            **({"comments": comments} if comments else {}),
        }

    if strengths and all(withdrawn.get(lab) for lab in strengths):
        applicability = "not_applicable"
    elif strengths:
        applicability = "applicable"
    elif not_applicable:
        applicability = "not_applicable"
    else:
        applicability = "unknown"
    out = {
        "applicability": applicability,
        "strengths": strengths,
        "acmg_summary": summary,
    }
    if specs:
        out["vcep_specifications"] = specs
    if comments:
        out["comments"] = comments
    return out


def build(gns: list[str]) -> dict:
    genes: dict[str, dict] = {}
    for gn in gns:
        print(f"fetching {gn} ...", file=sys.stderr)
        html = fetch(gn)
        blocks = _criterion_blocks(html)
        if len(blocks) < 20:
            raise SystemExit(
                f"{gn}: only {len(blocks)} criteria parsed — the registry markup "
                "has probably changed. Refusing to emit a partial table."
            )
        parsed_by_gene = {
            gene: {code: parse_criterion(b, gene)
                   for code, b in blocks.items()}
            for gene in GN_GENES[gn]
        }
        for gene in GN_GENES[gn]:
            entry: dict = {"gn": gn}
            entry.update(parsed_by_gene[gene])
            genes[gene] = entry
    return {
        "_README": (
            "Per-gene ACMG criteria specifications harvested from the ClinGen "
            "CSpec Registry by scripts/build_vcep_criteria.py. Each criterion "
            "carries applicability AND per-strength rule text. applicability is "
            "'applicable' when the registry publishes at least one strength row "
            "whose first bullet is not a withdrawal, 'not_applicable' when it "
            "renders Not Applicable OR when every published strength row opens "
            "with a 'Not applicable ...' bullet (KCNQ1 GN112 withdraws "
            "PP2/BP1/BP2/BP3/BS2 that way), and 'unknown' when "
            "the page describes neither — 'unknown' must never be read as "
            "permission."
        ),
        "_source": BASE,
        "genes": genes,
    }


def verify(table: dict) -> int:
    g = table["genes"]
    failures: list[str] = []

    def check(cond: bool, msg: str) -> None:
        if not cond:
            failures.append(msg)

    fbn1_pm1 = g["FBN1"]["PM1"]["strengths"]
    check("Strong" in fbn1_pm1, "FBN1 PM1 should have a Strong row (cbEGF Cys)")
    check("cbEGF" in fbn1_pm1.get("Strong", ""),
          "FBN1 PM1_Strong should name cbEGF-like domains")
    check("Moderate" in fbn1_pm1, "FBN1 PM1 should have a Moderate row")
    check("TB domain" in fbn1_pm1.get("Moderate", ""),
          "FBN1 PM1 Moderate should name the TB domain")

    check(g["FBN1"]["PM4"]["applicability"] == "applicable",
          f"FBN1 PM4 should be applicable, got {g['FBN1']['PM4']['applicability']}")
    check("PVS1" in g["FBN1"]["PM4"]["strengths"].get("Moderate", ""),
          "FBN1 PM4's caveat should bar co-application with PVS1 — the "
          "gene-specific echo of Abou Tayoun 2018's global PM4/PVS1 exclusion")
    check(g["FBN1"]["PM5"]["applicability"] == "applicable",
          f"FBN1 PM5 should be applicable, got {g['FBN1']['PM5']['applicability']}")
    check("caution" in g["FBN1"]["PM5"]["strengths"].get("Moderate", "").lower(),
          "FBN1 PM5 should carry the cysteine-creating caution caveat")

    check(g["MYH7"]["PVS1"]["applicability"] == "not_applicable",
          f"MYH7 PVS1 should be not_applicable, got {g['MYH7']['PVS1']['applicability']}")

    ps3 = g["MYH7"]["PS3"]["strengths"]
    check("splicing" in ps3.get("Strong", "").lower(),
          "MYH7 PS3 Strong should be the in-vitro splicing / RNA route")
    check("knock-in" in ps3.get("Moderate", "").lower(),
          "MYH7 PS3 Moderate should be the knock-in animal route")
    check("myofilament" in ps3.get("Supporting", "").lower()
          or "ipsc" in ps3.get("Supporting", "").lower(),
          "MYH7 PS3 Supporting should be the in-vitro myofilament / iPSC-CM route")

    check("167-931" in g["MYH7"]["PM1"]["strengths"].get("Moderate", "")
          .replace("–", "-"),
          "MYH7 PM1 Moderate should carry the codon 167-931 hotspot")

    bm = g["BMPR2"]
    check("1%" in bm["BA1"]["strengths"].get("Stand Alone", ""),
          "BMPR2 BA1 Stand Alone should be the 1% gnomAD threshold")
    check("0.1%" in bm["BS1"]["strengths"].get("Strong", ""),
          "BMPR2 BS1 Strong should be the 0.1% threshold")
    check("0.01%" in bm["PM2"]["strengths"].get("Supporting", ""),
          "BMPR2 PM2 Supporting should be the 0.01% threshold")
    for code in ("BA1", "BS1", "PM2"):
        check("1,000 allele counts" in " ".join(
                  bm[code]["strengths"].values()),
              f"BMPR2 {code} should carry the 1,000-allele minimum")
    for code in ("PM3", "PP4", "PP5", "BP6"):
        check(bm[code]["applicability"] == "not_applicable",
              f"BMPR2 {code} should be not_applicable, got "
              f"{bm[code]['applicability']}")
    check(bm["BP5"]["applicability"] == "applicable",
          f"BMPR2 BP5 should be applicable, got {bm['BP5']['applicability']}")
    check(bm["PVS1"]["applicability"] == "applicable",
          "BMPR2 PVS1 should be applicable (LOF is the PAH mechanism)")
    check("critical amino acid" in bm["PM1"]["strengths"].get("Strong", ""),
          "BMPR2 PM1 Strong should be the critical-residue enumeration")
    check("33-131" in bm["PM1"]["strengths"].get("Moderate", "")
          .replace("–", "-"),
          "BMPR2 PM1 Moderate should carry the aa 33-131 extracellular domain")

    if "KCNQ1" in g:
        kq = g["KCNQ1"]
        for code, quote in (
            ("PP2", "presence of benign variation throughout the"),
            ("BP1", "not limited to truncating variants"),
            ("BP2", "biallelic cases"),
            ("BP3", "Not applicable to"),
            ("BS2", "incomplete penetrance"),
        ):
            check(kq[code]["applicability"] == "not_applicable",
                  f"KCNQ1 {code} should be not_applicable (GN112 withdraws it "
                  f"in the first bullet), got {kq[code]['applicability']}")
            check(quote in " ".join(kq[code]["strengths"].values()),
                  f"KCNQ1 {code}'s withdrawal text should still be recorded "
                  f"(expected {quote!r})")
        check(kq["BP5"]["applicability"] == "applicable",
              "KCNQ1 BP5 is LIMITED to phenotypes matching another LQTS gene, "
              f"not withdrawn, got {kq['BP5']['applicability']}")
        for code in ("PM1", "PP4", "PVS1", "PS4"):
            check(kq[code]["applicability"] == "applicable",
                  f"KCNQ1 {code} should stay applicable, got "
                  f"{kq[code]['applicability']}")
    for gene, code in (("MYH7", "PM4"), ("MYH7", "PM5"), ("PTPN11", "PM1"),
                       ("BMPR2", "BP7")):
        if gene in g:
            check(g[gene][code]["applicability"] == "applicable",
                  f"{gene} {code} carries a 'not applicable' CAVEAT, not a "
                  f"withdrawal, got {g[gene][code]['applicability']}")

    for f in failures:
        print(f"FAIL  {f}", file=sys.stderr)
    print(
        f"\n{len(failures)} failure(s); {len(g)} genes, "
        f"{sum(len(v) - 1 for v in g.values())} criterion entries",
        file=sys.stderr,
    )
    return 1 if failures else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", help="write the table here")
    ap.add_argument("--verify", action="store_true",
                    help="check the parse against hand-read spec values")
    ap.add_argument("--gn", nargs="*", help="limit to these GN ids")
    args = ap.parse_args()

    gns = args.gn or sorted(GN_GENES)
    table = build(gns)
    rc = verify(table) if args.verify else 0
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(table, fh, indent=1, ensure_ascii=False, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.out}", file=sys.stderr)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
