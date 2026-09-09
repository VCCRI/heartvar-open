
'use strict';
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const JS_PATH = process.env.HEARTVAR_JS_PATH
  || path.resolve(__dirname, '..', '..', '..', 'static', 'heartvar.js');

function makeElement(tag) {
  const el = {
    tagName: String(tag || 'div').toUpperCase(),
    style: {}, dataset: {}, classList: {
      add() {}, remove() {}, toggle() {}, contains() { return false; },
    },
    children: [], attributes: {},
    innerHTML: '', textContent: '', value: '', checked: false,
    scrollHeight: 0, offsetWidth: 0, clientWidth: 0,
    appendChild(c) { this.children.push(c); return c; },
    removeChild(c) { return c; },
    insertAdjacentHTML() {},
    setAttribute(k, v) { this.attributes[k] = v; },
    getAttribute(k) { return Object.prototype.hasOwnProperty.call(this.attributes, k) ? this.attributes[k] : null; },
    removeAttribute(k) { delete this.attributes[k]; },
    addEventListener() {}, removeEventListener() {},
    querySelector() { return makeElement('div'); },
    querySelectorAll() { return []; },
    closest() { return null; },
    focus() {}, blur() {}, click() {}, remove() {},
    getBoundingClientRect() { return { top: 0, left: 0, width: 0, height: 0, bottom: 0, right: 0 }; },
    scrollIntoView() {},
  };
  el.parentElement = null;
  return el;
}

const documentStub = {
  readyState: 'complete',
  body: makeElement('body'),
  documentElement: makeElement('html'),
  getElementById() { return makeElement('div'); },
  querySelector() { return makeElement('div'); },
  querySelectorAll() { return []; },
  createElement(t) { return makeElement(t); },
  createElementNS(_ns, t) { return makeElement(t); },
  createTextNode(t) { return { textContent: t }; },
  addEventListener() {}, removeEventListener() {},
  head: makeElement('head'),
  cookie: '',
};

const sandbox = {
  console,
  setTimeout, clearTimeout, setInterval, clearInterval,
  queueMicrotask,
  Promise, Math, JSON, Date, RegExp, Error, Object, Array, String, Number, Boolean, Map, Set, WeakMap,
  encodeURIComponent, decodeURIComponent, encodeURI, decodeURI,
  parseInt, parseFloat, isNaN, isFinite,
  document: documentStub,
  navigator: { userAgent: 'node-render-smoke', clipboard: { writeText: () => Promise.resolve() } },
  location: { href: 'http://localhost/', search: '', hash: '', pathname: '/', origin: 'http://localhost' },
  history: { pushState() {}, replaceState() {} },
  localStorage: {
    _d: {},
    getItem(k) { return Object.prototype.hasOwnProperty.call(this._d, k) ? this._d[k] : null; },
    setItem(k, v) { this._d[k] = String(v); },
    removeItem(k) { delete this._d[k]; },
    clear() { this._d = {}; },
  },
  sessionStorage: { getItem: () => null, setItem() {}, removeItem() {}, clear() {} },
  fetch: () => Promise.resolve({ ok: true, json: () => Promise.resolve({}), text: () => Promise.resolve('') }),
  EventSource: function EventSource() { return { addEventListener() {}, close() {} }; },
  requestAnimationFrame: (fn) => setTimeout(fn, 0),
  cancelAnimationFrame: () => {},
  getComputedStyle: () => ({ getPropertyValue: () => '' }),
  matchMedia: () => ({ matches: false, addEventListener() {}, addListener() {} }),
  alert() {}, confirm: () => true, prompt: () => null,
  URL, URLSearchParams, TextEncoder, TextDecoder,
  performance: { now: () => 0 },
  addEventListener() {}, removeEventListener() {},
  scrollTo() {}, innerWidth: 1280, innerHeight: 900, devicePixelRatio: 1,
};

const context = vm.createContext(sandbox);

vm.runInContext('globalThis.window = globalThis; globalThis.self = globalThis;', context);

let source = fs.readFileSync(JS_PATH, 'utf8');
if (source.charCodeAt(0) === 0xfeff) source = source.slice(1);
vm.runInContext(source, context, { filename: 'heartvar.js' });

function gata4Evidence() {
  return {
    vep: {
      ok: true, most_severe_consequence: 'missense_variant',
      seq_region_name: '8', start: 11750234, end: 11750234, strand: 1,
      allele_string: 'G/T', vcf_string: '8-11750234-G-T',
      transcript_id: 'NM_002052.5', selected_transcript_id: 'NM_002052.5',
      selected_transcript_source: 'user-supplied',
      user_transcript: 'NM_002052.5',
      hgvsc: 'NM_002052.5:c.907G>T', hgvsp: 'NP_002043.2:p.Gly303Trp',
      cds_start: 907, protein_start: 303, impact: 'MODERATE',
      biotype: 'protein_coding', gene_symbol: 'GATA4', exon: '6/7',
      is_mane_select: false, is_mane_clinical: false,
      mane_select_accession: null, mane_clinical_accession: null,
      consequence_terms: ['missense_variant', 'splice_region_variant'],
      transcript_consequences_all: [
        { transcript_id: 'NM_001308093.3', source: 'RefSeq', biotype: 'protein_coding',
          consequence_terms: ['missense_variant'], consequence: 'missense_variant',
          hgvsc: 'NM_001308093.3:c.910G>T', hgvsp: 'NP_001295022.1:p.Gly304Trp',
          impact: 'MODERATE', is_mane_select: true, is_mane_plus_clinical: false,
          mane_select_accession: 'ENST00000532059.6', mane_plus_clinical_accession: null,
          canonical: false, is_picked: false },
        { transcript_id: 'NM_002052.5', source: 'RefSeq', biotype: 'protein_coding',
          consequence_terms: ['missense_variant'], consequence: 'missense_variant',
          hgvsc: 'NM_002052.5:c.907G>T', hgvsp: 'NP_002043.2:p.Gly303Trp',
          impact: 'MODERATE', is_mane_select: false, is_mane_plus_clinical: false,
          mane_select_accession: null, mane_plus_clinical_accession: null,
          canonical: false, is_picked: true },
      ],
    },
    gnomad: {
      ok: true,
      variant: {
        variantId: '8-11750234-G-T', rsid: null,
        exome: { ac: 0, an: 1460238, af: 0, ac_hom: 0, ac_hemi: null, faf95: { popmax: null }, populations: [] },
        genome: null,
      },
      gene: { gnomad_constraint: { pli: 0.86, oe_lof_upper: 0.42 } },
    },
    same_site: {
      ok: true, available: true,
      chrom: '8', variant_position: 11750234,
      codon_positions: [11750234, 11750235, 11750236],
      codon_span_available: true,
      codon_span_source: 'MANE ENST00000532059',
      codon_span_note: 'Ensembl serves an exon model only for an unversioned ENST accession, so the codon span was proven contiguous against the MANE transcript ENST00000532059 rather than against NM_002052.5.',
      codon_span_confirmed_on_picked_transcript: false,
      gnomad_available: true, gnomad_reason: null, gnomad_note: null,
      same_nucleotide: [
        { variant_id: '8-11750234-G-A', position: 11750234, ref: 'G', alt: 'A',
          rsid: 'rs1205549216', ac: 15, an: 1612518, af: 9.3e-6, faf95_popmax: 5.75e-6 },
      ],
      same_codon: [
        { variant_id: '8-11750236-G-A', position: 11750236, ref: 'G', alt: 'A',
          rsid: 'rs773684507', ac: 2, an: 1460216, af: 1.37e-6, faf95_popmax: 7.41e-6 },
      ],
      display_residue: 303, matched_residue: 304, numbering_differs: true,
      clinvar_same_residue: {
        ok: true, gene: 'GATA4', display_position: 303, matched_position: 304,
        numbering_differs: true, count: 2, pm5_eligible_count: 0,
        records: [
          { variation_id: 3897161, accession: 'VCV003897161',
            name: 'NM_001308093.3(GATA4):c.910G>A (p.Gly304Arg)',
            clinical_significance: 'Uncertain significance', tier: 'VUS',
            review_status: 'criteria provided, multiple submitters, no conflicts',
            stars: 2, consequence: 'missense', alt_aa: 'Arg', position: 11750234,
            is_proband_own_record: false, pm5_eligible: false,
            pm5_ineligible_reason: "classified 'Uncertain significance' — PM5 requires an established Pathogenic/Likely-pathogenic comparator" },
          { variation_id: 1519604, accession: 'VCV001519604',
            name: 'NM_001308093.3(GATA4):c.912G>A (p.Gly304=)',
            clinical_significance: 'Uncertain significance', tier: 'VUS',
            review_status: 'criteria provided, single submitter',
            stars: 1, consequence: 'synonymous', alt_aa: null, position: 11750236,
            is_proband_own_record: false, pm5_eligible: false,
            pm5_ineligible_reason: 'synonymous (no amino-acid change) — PM5 requires a DIFFERENT MISSENSE change at the residue' },
        ],
      },
    },
    clinvar: {
      ok: true, found: false, records: [], total_records: 0, total_submissions: 0,
      all_conditions: [], matched_on_mane: false, position_matching_used: true,
      same_position_count: 1, same_position_allele_confirmed: false,
      same_position_records: [
        { uid: '3897161', accession: 'VCV003897161',
          title: 'NM_001308093.3(GATA4):c.910G>A (p.Gly304Arg)',
          variation_id: 3897161, clinical_significance: 'Uncertain significance',
          review_status: 'criteria provided, multiple submitters, no conflicts',
          stars: 2, last_evaluated: 'Jan 21, 2026', position: 11750234,
          consequence: 'missense', allele_confirmed: false },
      ],
    },
    clinvar_pm5_candidates: {
      ok: true, gene: 'GATA4', protein_position: 303,
      matched_protein_position: 304, numbering_differs: true,
      candidates: [], count: 0, count_two_star: 0,
      ps1_candidates: [], ps1_count: 0, ps1_count_two_star: 0,
    },
    chdgene: { ok: true, listed: true, gene: 'GATA4' },
    hpo_resolved: {},
  };
}

function offPanelEvidence() {
  const ev = gata4Evidence();
  ev.same_site = {
    ok: true, available: false,
    chrom: '8', variant_position: 11750234,
    codon_positions: null, codon_span_available: false,
    codon_span_source: null, codon_span_confirmed_on_picked_transcript: false,
    gnomad_available: false, gnomad_reason: 'off-panel',
    gnomad_note: 'same-site data not available off-panel — the local gnomAD frequency DB covers the cardiac panel intervals only',
    same_nucleotide: [], same_codon: [],
    display_residue: 303, matched_residue: 304, numbering_differs: true,
    clinvar_same_residue: { ok: true, records: [], count: 0, pm5_eligible_count: 0,
                            display_position: 303, matched_position: 304, numbering_differs: true },
  };
  return ev;
}

function maneInputEvidence() {
  const ev = gata4Evidence();
  ev.vep.transcript_id = 'NM_001308093.3';
  ev.vep.selected_transcript_id = 'NM_001308093.3';
  ev.vep.user_transcript = 'NM_001308093.3';
  ev.vep.hgvsc = 'NM_001308093.3:c.910G>T';
  ev.vep.hgvsp = 'NP_001295022.1:p.Gly304Trp';
  ev.vep.is_mane_select = true;
  ev.vep.mane_select_accession = 'ENST00000532059.6';
  ev.same_site.numbering_differs = false;
  ev.same_site.clinvar_same_residue.numbering_differs = false;
  ev.same_site.clinvar_same_residue.display_position = 304;
  ev.clinvar_pm5_candidates.numbering_differs = false;
  ev.clinvar_pm5_candidates.protein_position = 304;
  return ev;
}

function myh7Evidence() {
  const ev = gata4Evidence();
  ev.vep.gene_symbol = 'MYH7';
  ev.vep.seq_region_name = '14';
  ev.vep.start = 23424001; ev.vep.end = 23424001; ev.vep.strand = -1;
  ev.vep.vcf_string = '14-23424001-C-T';
  ev.vep.transcript_id = 'NM_000257.4';
  ev.vep.selected_transcript_id = 'NM_000257.4';
  ev.vep.user_transcript = 'NM_000257.4';
  ev.vep.hgvsc = 'NM_000257.4:c.1208G>A';
  ev.vep.hgvsp = 'NP_000248.2:p.Arg403Gln';

  ev.gnomad.variant.variantId = '14-23424001-C-T';
  ev.gnomad.variant.rsid = 'rs121913627';
  ev.gnomad.variant.exome = { ac: 2, an: 1400000, af: 1.4e-6, ac_hom: 0,
                              ac_hemi: null, faf95: { popmax: null }, populations: [] };
  ev.same_site.chrom = '14';
  ev.same_site.variant_position = 23424001;
  ev.same_site.codon_positions = [23424000, 23424001, 23424002];
  ev.same_site.same_nucleotide = [
    { variant_id: '14-23424001-C-A', position: 23424001, ref: 'C', alt: 'A',
      rsid: null, ac: 3, an: 1400000, af: 2.1e-6, faf95_popmax: null },
  ];
  ev.same_site.same_codon = [];
  return ev;
}

const cases = [
  ['GATA4 c.907G>T (the reported variant)', gata4Evidence()],
  ['GATA4 off-panel / no codon span', offPanelEvidence()],
  ['GATA4 on the MANE transcript', maneInputEvidence()],
  ['MYH7 c.1208G>A (different gene, same-base allele)', myh7Evidence()],
];

const failures = [];
const results = {};

if (typeof context._renderEvidence !== 'function') {
  failures.push('_renderEvidence is not reachable on the script global — '
    + 'heartvar.js may have been wrapped in an IIFE; update this harness.');
}

for (const [label, ev] of cases) {
  if (failures.length) break;
  let html;
  try {
    html = context._renderEvidence(ev);
  } catch (e) {
    failures.push(`${label}: _renderEvidence THREW ${e && e.name}: ${e && e.message}`);
    continue;
  }
  if (typeof html !== 'string' || html.length < 200) {
    failures.push(`${label}: returned ${typeof html} of length `
      + `${html && html.length} — expected markup`);
    continue;
  }
  results[label] = html;
}

function want(label, needle) {
  const html = results[label];
  if (html === undefined) return;
  if (!html.includes(needle)) {
    failures.push(`${label}: rendered markup is missing ${JSON.stringify(needle)}`);
  }
}
function reject(label, needle) {
  const html = results[label];
  if (html === undefined) return;
  if (html.includes(needle)) {
    failures.push(`${label}: rendered markup unexpectedly contains ${JSON.stringify(needle)}`);
  }
}

const G = 'GATA4 c.907G>T (the reported variant)';

want(G, 'Other alleles at this base');
want(G, 'rs1205549216');
want(G, 'Elsewhere in this codon');
want(G, 'rs773684507');

want(G, 'Not MANE Select');
want(G, 'ENST00000532059');

reject('GATA4 on the MANE transcript', 'Not MANE Select');

want('GATA4 off-panel / no codon span', 'Not available');
want('GATA4 off-panel / no codon span', 'off-panel');
reject('GATA4 off-panel / no codon span', 'Elsewhere in this codon');

want('MYH7 c.1208G>A (different gene, same-base allele)', 'Other alleles at this base');

want(G, '<span title="');
want(G, 'codon span from MANE ENST00000532059');
reject(G, 'Elsewhere in this codon (codon span from');
reject(G, 'has no exon model))');

want(G, 'chr8:11750234');
reject(G, '(8:11750234)');

want(G, 'rsID note');
want('MYH7 c.1208G>A (different gene, same-base allele)', 'rsID note');

{
  const evAll = gata4Evidence();
  evAll.gnomad.variant.rsid = 'rs999999999';
  let h;
  try { h = context._renderEvidence(evAll); } catch (e) {
    failures.push(`rsID-note control THREW ${e && e.name}: ${e && e.message}`);
  }
  if (typeof h === 'string' && h.includes('rsID note')) {
    failures.push('rsID note fired although every allele at the base has an rsID');
  }
}

if (typeof context._buildLandscapeLollipop === 'function') {
  const lsData = {
    positions: [

      { gpos: 11750234, aa: 304, tier: 'VUS', csq: 'missense', stars: 2, count: 1 },

      { gpos: 11750236, aa: null, tier: 'VUS', csq: 'synonymous', stars: 1, count: 1 },
      { gpos: 11750213, aa: 297, tier: 'P', csq: 'missense', stars: 2, count: 4 },
    ],
    tierCounts: { P: 1, LP: 0, VUS: 2, LB: 0, B: 0 },
    tierByCsq: { missense: { P: 1, VUS: 1 }, synonymous: { VUS: 1 } },
    gMin: 11750213, gMax: 11750236, strand: 1,
    probandGpos: 11750234, probandExon: '6/7', exons: null,
    gene: 'GATA4', truncated: false, modelLabel: '',

    clinvarTxId: 'NM_001308093.3',
    probandResidueClinvar: 304,
    probandResidueModel: 303,
    codonPositions: [11750234, 11750235, 11750236],
  };
  let svg;
  try {
    svg = context._buildLandscapeLollipop(lsData, 'all', true);
  } catch (e) {
    failures.push(`_buildLandscapeLollipop THREW ${e && e.name}: ${e && e.message}`);
  }
  if (typeof svg === 'string' && svg.length) {

    const marks = [...svg.matchAll(
      /<circle[^>]*opacity="([\d.]+)"[^>]*>\s*<title>([^<]*?)([0-9]+)★<\/title>/g)]
      .map(m => ({ opacity: parseFloat(m[1]), stars: parseInt(m[3], 10),
                   tip: m[2], isVus: /^VUS/.test(m[2]) }));
    if (marks.length < 3) {
      failures.push(`D6: expected 3 marks with tooltips, parsed ${marks.length}`);
    }
    for (const m of marks) {
      const wanted = m.stars >= 2 ? 0.95 : 0.5;
      if (m.opacity !== wanted) {
        failures.push(`D6: a ${m.stars}★ mark rendered at opacity ${m.opacity}, `
          + `expected ${wanted} — the "faded = <2★" caption is false for it`);
      }
    }
    if (!marks.some(m => m.isVus && m.stars === 2 && m.opacity === 0.95)) {
      failures.push('D6: the 2★ VUS mark is not rendered at full opacity');
    }
    if (!marks.some(m => m.isVus && m.stars === 1 && m.opacity === 0.5)) {
      failures.push('D6: the 1★ VUS mark should still be faded');
    }
    if (svg.includes('size ∝ records') && !svg.includes('needle size ∝ records')) {
      failures.push('D6: caption still claims size scales for EVERY mark, '
        + 'but VUS dots are fixed-radius');
    }

    if (!svg.includes('residue 303 on this model (304 in ClinVar · NM_001308093.3)')) {
      failures.push('D2: the mark at the proband residue does not show both '
        + 'residue numbers');
    }

    if (/<title>VUS · residue 304 ·/.test(svg)) {
      failures.push('D2: a mark still reports a bare ClinVar residue number '
        + 'with no indication of which numbering it is in');
    }

    if (!svg.includes('residue 297 (ClinVar numbering · NM_001308093.3)')) {
      failures.push('D2: an off-residue mark does not name its numbering');
    }
    if (svg.includes('residue 296')) {
      failures.push('D2: a distant mark was translated into the model\'s '
        + 'numbering — that offset is not known to hold there');
    }

    if (!svg.includes('residue from position; the ClinVar name has none')) {
      failures.push('D2: the aa:null mark inside the proband codon carries no '
        + 'residue at all');
    }

    if (!svg.includes('Residue numbers on the marks are ClinVar')) {
      failures.push('D2: the chart caption does not explain the numbering split');
    }

    const same = Object.assign({}, lsData, {
      probandResidueClinvar: 304, probandResidueModel: 304,
    });
    let svgSame;
    try {
      svgSame = context._buildLandscapeLollipop(same, 'all', true);
    } catch (e) {
      failures.push(`_buildLandscapeLollipop (same-numbering) THREW ${e && e.name}: ${e && e.message}`);
    }
    if (typeof svgSame === 'string') {
      if (svgSame.includes('in ClinVar · NM_001308093.3')
          || svgSame.includes('(ClinVar numbering')) {
        failures.push('D2: the second residue number is shown even though the '
          + 'two numberings agree');
      }
      if (svgSame.includes('Residue numbers on the marks are ClinVar')) {
        failures.push('D2: the numbering caption shows when nothing differs');
      }
      if (!/<title>VUS · residue 304 ·/.test(svgSame)) {
        failures.push('D2: the plain single-number tooltip was lost');
      }
    }
  } else if (!failures.length) {
    failures.push('_buildLandscapeLollipop returned no markup');
  }
} else {
  failures.push('_buildLandscapeLollipop is not reachable on the script global');
}

if (typeof context._renderGeneContextTab === 'function') {
  const ev = gata4Evidence();
  ev.clinvar_gene_landscape = {
    ok: true, gene: 'GATA4', total_classified: 3,
    tier_counts: { P: 1, LP: 0, VUS: 2, LB: 0, B: 0 },
    tier_counts_by_csq: { missense: { P: 1, LP: 0, VUS: 1, LB: 0, B: 0 },
                          synonymous: { P: 0, LP: 0, VUS: 1, LB: 0, B: 0 } },
    positions: [
      { gpos: 11750234, aa: 304, tier: 'VUS', csq: 'missense', stars: 2, count: 1 },
      { gpos: 11750236, aa: null, tier: 'VUS', csq: 'synonymous', stars: 1, count: 1 },
      { gpos: 11750213, aa: 297, tier: 'P', csq: 'missense', stars: 2, count: 4 },
    ],
    positions_truncated: false, g_min: 11750213, g_max: 11750236,
    plp_fraction_pct: 33.3, has_two_star_plp: true,
    condition_keywords: null, keyword_source: 'gencc',
    records_examined: 3, records_filtered_out: 0, plp_top_phenotypes: [],
  };
  let gcHtml;
  try {
    gcHtml = context._renderGeneContextTab(ev, {}, 'GATA4');
  } catch (e) {
    failures.push(`_renderGeneContextTab THREW ${e && e.name}: ${e && e.message}`);
  }
  if (typeof gcHtml === 'string' && gcHtml.length > 200) {

    if (!gcHtml.includes('p.Gly304Arg') || !gcHtml.includes('p.Gly304=')) {
      failures.push('gene context: the same-residue records are not rendered');
    }

    if (gcHtml.includes('ev-cls-ben">not PM5 evidence')) {
      failures.push('D5: "not PM5 evidence" still wears the benign class');
    }
    if (!gcHtml.includes('ev-dim">not PM5 evidence')) {
      failures.push('D5: "not PM5 evidence" chip missing its neutral class');
    }
    if (!gcHtml.includes("on ClinVar's MANE transcript")) {
      failures.push('the same-residue caption does not name both residue numbers');
    }

    if (!gcHtml.includes('ClinVar records at residue 303')) {
      failures.push('D1: the block heading is not the neutral list heading');
    }
    if (/PM5 \u2014 known pathogenic variant/.test(gcHtml)) {
      failures.push('D1: the asserting "PM5 — known pathogenic variant(s)" '
        + 'heading is still rendered');
    }

    const boxes = (gcHtml.match(/gc-pm5-callout-title/g) || []).length;
    if (boxes !== 1) {
      failures.push(`D1: expected 1 residue callout, found ${boxes}`);
    }

    if (!gcHtml.includes('gc-pm5-callout')) {
      failures.push('D1: the block did not render when PM5 has no candidates');
    }

    if (!gcHtml.includes('PM5 row of the Criteria tab')) {
      failures.push('D1: the block does not defer the verdict to the PM5 row');
    }
  } else if (!failures.length) {
    failures.push('_renderGeneContextTab returned no markup');
  }

  {
    const evFb = gata4Evidence();
    evFb.clinvar_gene_landscape = ev.clinvar_gene_landscape;
    delete evFb.same_site;
    evFb.clinvar_pm5_candidates = {
      ok: true, gene: 'GATA4', protein_position: 296,
      matched_protein_position: 297, numbering_differs: true,
      candidates: [
        { variation_id: 30098, name: 'NM_001308093.3(GATA4):c.889G>T (p.Gly297Cys)',
          clinical_significance: 'Pathogenic', tier: 'P',
          review_status: 'no assertion criteria provided', stars: 0 },
      ],
      count: 1, count_two_star: 0,
      ps1_candidates: [], ps1_count: 0, ps1_count_two_star: 0,
    };
    let h;
    try { h = context._renderGeneContextTab(evFb, {}, 'GATA4'); } catch (e) {
      failures.push(`D1 fallback THREW ${e && e.name}: ${e && e.message}`);
    }
    if (typeof h === 'string') {
      if (!h.includes('p.Gly297Cys')) {
        failures.push('D1: with same_site absent the block lists nothing, so '
          + 'the panel goes silent instead of degrading');
      }
      if (h.includes('PM5-eligible')) {
        failures.push('D1: the fallback path claims per-row PM5 eligibility it '
          + 'never received');
      }
      if (!h.includes('eligibility was not available')) {
        failures.push('D1: the fallback path does not say the verdicts are '
          + 'missing');
      }
    }
  }

  {
    const evNone = gata4Evidence();
    evNone.clinvar_gene_landscape = ev.clinvar_gene_landscape;
    evNone.same_site.clinvar_same_residue = {
      ok: true, gene: 'GATA4', display_position: 303, matched_position: 304,
      numbering_differs: true, records: [], count: 0, pm5_eligible_count: 0,
    };
    let h;
    try { h = context._renderGeneContextTab(evNone, {}, 'GATA4'); } catch (e) {
      failures.push(`D1 empty-residue THREW ${e && e.name}: ${e && e.message}`);
    }
    if (typeof h === 'string' && !h.includes('No other ClinVar record is recorded')) {
      failures.push('D1: an empty residue renders no block, so "nothing '
        + 'submitted" is indistinguishable from "block missing"');
    }
  }

  {
    const evSplice = gata4Evidence();
    evSplice.clinvar_gene_landscape = ev.clinvar_gene_landscape;
    evSplice.vep.most_severe_consequence = 'splice_donor_variant';
    evSplice.vep.hgvsp = null;
    delete evSplice.same_site;
    evSplice.clinvar_pm5_candidates = {
      ok: true, gene: 'GATA4', not_applicable: true, candidates: [], count: 0,
      ps1_candidates: [], ps1_count: 0, ps1_count_two_star: 0,
    };
    let h;
    try { h = context._renderGeneContextTab(evSplice, {}, 'GATA4'); } catch (e) {
      failures.push(`D1 no-residue THREW ${e && e.name}: ${e && e.message}`);
    }
    if (typeof h === 'string' && h.includes('ClinVar records at residue')) {
      failures.push('D1: a residue-scoped block rendered for a variant with '
        + 'no protein residue');
    }
  }
} else {
  failures.push('_renderGeneContextTab is not reachable on the script global');
}

if (failures.length) {
  console.error('FAIL\n' + failures.map(f => '  - ' + f).join('\n'));
  process.exit(1);
}
console.log(`OK — _renderEvidence executed and produced markup for ${cases.length} payloads`);
for (const [label] of cases) {
  console.log(`   ${String(results[label].length).padStart(7)} chars  ${label}`);
}
