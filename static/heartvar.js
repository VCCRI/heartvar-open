

    const API_ENDPOINT = '/api/curate/stream';

    let _TIER_POINTS = {};
    let _CRITERION_NAMES = {};

    let CRIT_POINTS = {};
    let HARD_CODED_CRITERIA_CODES = new Set();
    let CANONICAL_CRITERIA_ORDER = [];

    const SERVER_OWNED_CODES = new Set(['PM1', 'PP1', 'BS4']);

    const ACMG_READY = fetch('/static/acmg_constants.json')
      .then(r => {
        if (!r.ok) throw new Error(`HTTP ${r.status}`);
        return r.json();
      })
      .then(d => {
        _TIER_POINTS = d.tier_points;
        CRIT_POINTS = d.bare_code_points;
        _CRITERION_NAMES = d.criterion_names;
        HARD_CODED_CRITERIA_CODES = new Set(d.hard_coded_criteria_codes);
        CANONICAL_CRITERIA_ORDER = d.criteria_display_order;
        return d;
      })
      .catch(e => {
        console.error('Failed to load ACMG constants (static/acmg_constants.json):', e);
        throw e;
      });

    const COORD_INPUT_RE =
      /^(?:chr)?[\dXYM]+T?[:\-_]\d+[:\-_][ACGT]+[:\-_][ACGT]+$/i;

    function isCoordInput(s) {
      return !!s && COORD_INPUT_RE.test(s.trim());
    }

    function isGenomicHgvs(s) {
      s = String(s || '');
      return !isCoordInput(s) && /(^|[\s:])g\.\d/i.test(s);
    }

    let _geneValidateSeq = 0;
    function clearGeneSymbolMsg() {
      const el = document.getElementById('hvl-gene-msg');
      if (el) { el.textContent = ''; el.className = 'hvl-gene-msg'; }
    }
    async function validateGeneSymbol() {
      const input = document.getElementById('hvl-gene');
      const el = document.getElementById('hvl-gene-msg');
      if (!input || !el) return;
      const sym = input.value.trim();
      if (!sym) { clearGeneSymbolMsg(); return; }
      const seq = ++_geneValidateSeq;
      try {
        const r = await fetch('/api/gene/validate?symbol=' + encodeURIComponent(sym));
        if (!r.ok) { clearGeneSymbolMsg(); return; }
        const d = await r.json();
        if (seq !== _geneValidateSeq) return;
        if (d.status === 'alias' || d.status === 'unrecognized') {
          el.textContent = d.message || '';
          el.className = 'hvl-gene-msg hvl-gene-msg--warn show';
        } else {
          clearGeneSymbolMsg();
        }
      } catch (e) {
        clearGeneSymbolMsg();
      }
    }

    function linkifyPmids(s) {
      if (s == null) return '';
      return String(s).replace(
        /(PMIDs?:?\s*)(\d{4,9}(?:\s*,\s*\d{4,9})*)/gi,
        (m, label, ids) => {
          if (/pubmed\.ncbi/i.test(label)) return m;
          const linked = ids.replace(/\d{4,9}/g, id =>
            `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}/" target="_blank" rel="noopener">${id}</a>`);
          return label + linked;
        }
      );
    }

    function getSelectedGenomeBuild() {
      const row = document.getElementById('build-toggle-row');
      if (!row || !row.classList.contains('is-visible')) return null;
      return row.dataset.build || 'GRCh38';
    }

    function hvlSyncBuildToggle() {
      const row = document.getElementById('build-toggle-row');
      const field = document.getElementById('hvl-variant');
      if (!row || !field) return;
      if (isCoordInput(field.value)) row.classList.add('is-visible');
      else row.classList.remove('is-visible');
    }

    function hvlSetBuild(btn) {
      const row = document.getElementById('build-toggle-row');
      if (!row || !btn) return;
      row.dataset.build = btn.dataset.build || 'GRCh38';
      row.querySelectorAll('.hvl-build-opt').forEach(b =>
        b.classList.toggle('is-active', b === btn));
    }

    function hvlResetBuildToggle() {
      const row = document.getElementById('build-toggle-row');
      if (!row) return;
      row.dataset.build = 'GRCh38';
      row.classList.remove('is-visible');
      row.querySelectorAll('.hvl-build-opt').forEach(b =>
        b.classList.toggle('is-active', (b.dataset.build || 'GRCh38') === 'GRCh38'));
    }

    function getAiEnabled() {
      const cb = document.getElementById('hvl-ai-enable');
      return !!(cb && cb.checked);
    }

    let lastResult = null;
    let lastEvidence = null;
    let lastVariantId = null;

    let erepoVerdict = null;

    function renderErepoPanel() {
      const host = document.getElementById('summary-erepo-panel');
      if (!host) return;
      const v = erepoVerdict;
      if (!v || !v.found) { host.innerHTML = ''; return; }
      const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g,
        c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
      const crit = Array.isArray(v.criteria) ? v.criteria.map(esc).join(' · ') : '';
      const cls = esc(v.classification || '—');
      const vcep = esc(v.vcep || 'ClinGen VCEP');
      const link = v.url
        ? `<a class="erepo-ref-link" href="${esc(v.url)}" target="_blank" rel="noopener">View in eRepo ↗</a>`
        : '';
      host.innerHTML = `
        <div class="erepo-ref-box" role="note">
          <span class="erepo-ref-head">ClinGen VCEP classification (independent reference)</span>
          <div class="erepo-ref-class">${cls}<span class="erepo-ref-vcep"> · ${vcep}</span></div>
          ${crit ? `<div class="erepo-ref-crit">Criteria: ${crit}</div>` : ''}
          ${link}
          <div class="erepo-ref-note">For reference only — independent of HeartVar's analysis.</div>
        </div>`;
    }
    let lastVariant = {
      gene: '', hgvs_c: '', hpo: '', family: '', genome_build: '',

      zygosity: '', inheritance_input: '', proband_sex: '',
      trio_status: '', denovo_status: '', zygosity_inferred: false,

      seg_affected_carriers: 0, seg_affected_noncarriers: 0,
      denovo_confirmed_count: 0, denovo_unconfirmed_count: 0,
      seg_meioses: 0, in_trans_pathogenic: '',
      alt_cause_present: '', alt_cause_detail: '',

      family_history_summary: 'unknown',
    };

    window.heartvarProteinData = {
      hgvsp: null,
      hgvsc: null,
      gene: null,
      uniprotAccession: null,
      uniprotLength: null,
      uniprotFeatures: null,
      uniprotProteinName: null,
      uniprotUrl: null,
      uniprotOk: false,

      domainPlp: null,
    };

    function _populateProteinData(evidence, gene, hgvs_c) {
      const ev = evidence || {};

      const evVep = (ev && ev.vep) || {};
      gene = (gene || evVep.gene_symbol || evVep.derived_gene_symbol || '').trim();
      if (gene) window.heartvarProteinData.gene = gene;
      if (ev.vep) {
        const vep = ev.vep;

        window.heartvarProteinData.hgvsp = vep.hgvsp == null ? null : String(vep.hgvsp);
        window.heartvarProteinData.hgvsc = vep.hgvsc == null ? null : String(vep.hgvsc);
      }
      if (ev.uniprot) {
        const up = ev.uniprot;
        window.heartvarProteinData.uniprotAccession   = up.accession || null;
        window.heartvarProteinData.uniprotLength      = up.length    || null;
        window.heartvarProteinData.uniprotFeatures    = Array.isArray(up.features) ? up.features : null;
        window.heartvarProteinData.uniprotProteinName = up.protein_name || null;
        window.heartvarProteinData.uniprotUrl         = up.url || null;
        window.heartvarProteinData.uniprotOk          = !!up.ok;
      }
      if (ev.domain_plp) {
        window.heartvarProteinData.domainPlp = ev.domain_plp;
      }
    }

    function _mergeClinicalContextFromEvidence(evidence) {
      const cc = evidence && evidence.clinical_context;
      if (!cc || typeof cc !== 'object') return;

      if (cc.zygosity && !lastVariant.zygosity) {
        lastVariant.zygosity = cc.zygosity;
      }
      lastVariant.zygosity_inferred = !!cc.zygosity_inferred;

      if (cc.family_history_summary) {
        lastVariant.family_history_summary = cc.family_history_summary;
      }
    }

    function _resetProteinData(gene, hgvs_c) {
      window.heartvarProteinData.hgvsp = null;
      window.heartvarProteinData.hgvsc = null;
      window.heartvarProteinData.gene  = gene || null;
      window.heartvarProteinData.uniprotAccession   = null;
      window.heartvarProteinData.uniprotLength      = null;
      window.heartvarProteinData.uniprotFeatures    = null;
      window.heartvarProteinData.uniprotProteinName = null;
      window.heartvarProteinData.uniprotUrl         = null;
      window.heartvarProteinData.uniprotOk          = false;
      window.heartvarProteinData.domainPlp          = null;
    }

    const HV_PAGE_PATHS = { curate: '/', about: '/about', contact: '/contact' };

    function _hvPageForPath(pathname) {
      const p = (pathname || '/').replace(/\/+$/, '') || '/';
      if (p === '/about') return 'about';
      if (p === '/contact') return 'contact';
      return 'curate';
    }

    function _hvSyncUrl(page, section) {
      const url = (HV_PAGE_PATHS[page] || '/') + (section ? '#' + section : '');
      if (url === window.location.pathname + window.location.hash) return;

      try { history.pushState({ hvPage: page }, '', url); } catch (_e) {  }
    }

    window.hvSyncUrlToHome = function(){ _hvSyncUrl('curate'); };

    function hvRouteFromUrl() {
      const page = _hvPageForPath(window.location.pathname);
      const section = (window.location.hash || '').replace(/^#/, '');
      if (page === 'curate') {

        if (window.__hvRan) showPage('curate', { silent: true });
        else if (typeof returnToLanding === 'function') returnToLanding({ silent: true });
        return;
      }
      document.body.classList.remove('hv-landing-on');
      showPage(page, { silent: true });
      const el = section ? document.getElementById('about-' + section) : null;
      if (el) requestAnimationFrame(function(){ el.scrollIntoView({ block: 'start' }); });
      else window.scrollTo(0, 0);
    }

    window.addEventListener('popstate', function(){ hvRouteFromUrl(); });

    document.addEventListener('DOMContentLoaded', function(){
      if (_hvPageForPath(window.location.pathname) !== 'curate') hvRouteFromUrl();
    });

    function showPage(p, opts) {
      const silent = !!(opts && opts.silent);
      const split = document.querySelector('main.split');
      const about = document.getElementById('page-about');
      const contact = document.getElementById('page-contact');
      const admin = document.getElementById('page-admin');

      if (p === 'curate' && !window.__hvRan) {

        if (typeof returnToLanding === 'function') { returnToLanding({ silent: silent }); return; }
      }
      const onAbout = p === 'about';
      const onContact = p === 'contact';
      const onAdmin = p === 'admin';
      const onShell = !onAbout && !onContact && !onAdmin;
      if (split)   split.style.display   = onShell ? 'flex' : 'none';
      if (about)   about.style.display   = onAbout ? 'block' : 'none';
      if (contact) contact.style.display = onContact ? 'block' : 'none';
      if (admin)   admin.style.display   = onAdmin ? 'block' : 'none';
      document.body.classList.toggle('hv-page-about', onAbout);
      document.body.classList.toggle('hv-page-contact', onContact);
      document.body.classList.toggle('hv-page-admin', onAdmin);
      document.getElementById('nt-curate').classList.toggle('active', p === 'curate');
      document.getElementById('nt-about').classList.toggle('active', p === 'about');
      const adminTab = document.getElementById('nt-admin');
      if (adminTab) adminTab.classList.toggle('active', p === 'admin');

      if (onAdmin && typeof window.hvAdminLoadLogs === 'function') {
        window.hvAdminLoadLogs();
      }
      if (onAdmin && typeof window.hvAdminLoadDbStatus === 'function') {
        window.hvAdminLoadDbStatus();
      }

      if (!silent && !onAdmin) _hvSyncUrl(p, opts && opts.section);

      fadeInPage(onAbout ? about : onContact ? contact : onAdmin ? admin : split);
    }

    function fadeInPage(el) {
      if (!el) return;
      el.classList.remove('hv-page-fade');
      void el.offsetWidth;
      el.classList.add('hv-page-fade');
    }

    window.goAbout = function(section){
      document.body.classList.remove('hv-landing-on');
      showPage('about', { section: section || '' });
      if (section){
        var el = document.getElementById('about-' + section);
        if (el){ requestAnimationFrame(function(){ el.scrollIntoView({ behavior:'smooth', block:'start' }); }); return; }
      }
      window.scrollTo(0, 0);
    };

    window.goAboutSection = function(section){
      var el = document.getElementById('about-' + section);
      if (!el) return false;
      el.scrollIntoView({ behavior:'smooth', block:'start' });
      _hvSyncUrl('about', section);
      return false;
    };

    window.goContact = function(){
      document.body.classList.remove('hv-landing-on');
      showPage('contact');
      window.scrollTo(0, 0);
    };

    function showResultTab(name) {
      const panel = document.getElementById('result-panel');
      if (!panel) return;
      panel.querySelectorAll('.result-tab').forEach(t => {
        t.classList.toggle('active', t.dataset.tab === name);
      });
      panel.querySelectorAll('.result-tab-pane').forEach(p => {
        p.style.display = p.dataset.tab === name ? 'block' : 'none';
      });
    }

    function markCriteriaTabReady() {
      const panel = document.getElementById('result-panel');
      if (!panel) return;
      const dot = panel.querySelector('.result-tab-dot');
      if (dot) dot.setAttribute('data-state', 'ready');
    }

    function criterionPoints(code) {
      return Object.prototype.hasOwnProperty.call(CRIT_POINTS, code) ? CRIT_POINTS[code] : 0;
    }

    const CRITERIA_DEFINITIONS = {
      PVS1: { tier: 'Pathogenic Very Strong', direction: 'pathogenic',
        definition: 'Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a gene where loss-of-function is a known mechanism of disease.' },
      PS1:  { tier: 'Pathogenic Strong', direction: 'pathogenic',
        definition: 'Same amino acid change as a previously established pathogenic variant, regardless of nucleotide change.' },
      PS2:  { tier: 'Pathogenic Strong', direction: 'pathogenic',
        definition: 'De novo variant (confirmed with maternity and paternity) in a patient with the disease and no family history.' },
      PS3:  { tier: 'Pathogenic Strong', direction: 'pathogenic',
        definition: 'Well-established functional studies demonstrate a deleterious effect on the gene or gene product.' },
      PS4:  { tier: 'Pathogenic Strong', direction: 'pathogenic',
        definition: 'Prevalence of the variant in affected individuals is significantly increased compared to controls.' },
      PM1:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'Located in a mutational hotspot and/or well-established functional domain without benign variation.' },
      PM2:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'Absent from controls (or at extremely low frequency) in population databases such as gnomAD.' },
      PM3:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'For recessive disorders, detected in trans with a pathogenic variant.' },
      PM4:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'Protein length changes due to in-frame deletions/insertions in a non-repeat region, or stop-loss variants.' },
      PM5:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'Novel missense change at an amino acid residue where a different pathogenic missense change has been seen before.' },
      PM6:  { tier: 'Pathogenic Moderate', direction: 'pathogenic',
        definition: 'Assumed de novo, but without confirmation of paternity and maternity.' },
      PP1:  { tier: 'Pathogenic Supporting', direction: 'pathogenic',
        definition: 'Co-segregation with disease in multiple affected family members in a gene definitively known to cause the disease.' },
      PP2:  { tier: 'Pathogenic Supporting', direction: 'pathogenic',
        definition: 'Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease.' },
      PP3:  { tier: 'Pathogenic Supporting', direction: 'pathogenic',
        definition: 'Multiple lines of computational evidence support a deleterious effect (conservation, evolutionary, splicing impact, etc.).' },
      PP4:  { tier: 'Pathogenic Supporting', direction: 'pathogenic',
        definition: 'Patient\'s phenotype or family history is highly specific for a disease with a single genetic aetiology.' },
      PP5:  { tier: 'Pathogenic Supporting', direction: 'pathogenic',
        definition: 'Reputable source recently reports variant as pathogenic, but evidence not available to the laboratory to perform independent evaluation.' },
      BA1:  { tier: 'Benign Stand-alone', direction: 'benign',
        definition: 'Allele frequency >5% in a large population database such as gnomAD.' },
      BS1:  { tier: 'Benign Strong', direction: 'benign',
        definition: 'Allele frequency is greater than expected for the disorder.' },
      BS2:  { tier: 'Benign Strong', direction: 'benign',
        definition: 'Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder with full penetrance expected at an early age.' },
      BS3:  { tier: 'Benign Strong', direction: 'benign',
        definition: 'Well-established functional studies show no deleterious effect.' },
      BS4:  { tier: 'Benign Strong', direction: 'benign',
        definition: 'Lack of segregation in affected members of a family.' },
      BP1:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'Missense variant in a gene for which primarily truncating variants are known to cause disease.' },
      BP2:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder; or observed in cis with a pathogenic variant in any inheritance pattern.' },
      BP3:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'In-frame deletions/insertions in a repetitive region without a known function.' },
      BP4:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'Multiple lines of computational evidence suggest no impact on gene or gene product.' },
      BP5:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'Variant found in a case with an alternate molecular basis for disease.' },
      BP6:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'Reputable source recently reports variant as benign, but evidence not available to the laboratory to perform independent evaluation.' },
      BP7:  { tier: 'Benign Supporting', direction: 'benign',
        definition: 'A synonymous variant for which splicing prediction algorithms predict no impact on the splice consensus sequence nor the creation of a new splice site, and the nucleotide is not highly conserved.' },
    };

    function formatCritPoints(n) {
      if (n > 0) return `+${n}`;
      if (n < 0) return `−${Math.abs(n)}`;
      return '0';
    }

    function _scaleArrowPct(classification) {
      const m = { 'benign': 10, 'likely benign': 30, 'vus': 50, 'likely pathogenic': 70, 'pathogenic': 90 };
      const k = String(classification || '').trim().toLowerCase();
      return Object.prototype.hasOwnProperty.call(m, k) ? m[k] : null;
    }

    function renderSummaryCriteriaChips(criteria) {

      const heroPts = document.getElementById('summary-score-hero');
      if (heroPts && lastResult && Number.isFinite(lastResult.points_total)) {
        const p = lastResult.points_total;
        heroPts.textContent = `${p >= 0 ? '+' : '−'}${Math.abs(p)} pts`;
      }

      const scaleArrow = document.getElementById('summary-scale-arrow');
      if (scaleArrow && lastResult) {
        const pct = _scaleArrowPct(lastResult.classification);
        if (pct == null) { scaleArrow.style.display = 'none'; }
        else { scaleArrow.style.display = ''; scaleArrow.style.left = pct + '%'; }
      }

      const row = document.getElementById('summary-crit-row');
      if (!row) return;
      const met = (criteria || []).filter(c =>
        c && c.status === 'met' && (c.direction === 'pathogenic' || c.direction === 'benign')
      );
      if (!met.length) {
        row.innerHTML = `<span class="summary-crit-placeholder">No criteria met</span>`;
        return;
      }
      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

      const tagClassFor = (c) => {
        if (c.direction === 'benign') return 'crit-ben';
        let tier = strengthTier(c.criteria_strength) || strengthTier(c.code);
        if (!tier) {
          const code = String(c.code || '');
          if (/^PM\d/.test(code)) tier = 'Moderate';
          else if (/^PP\d/.test(code)) tier = 'Supporting';
          else tier = 'Strong';
        }
        if (tier === 'Moderate') return 'crit-mod';
        if (tier === 'Supporting') return 'crit-sup';
        return 'crit-path';
      };

      row.innerHTML = met.map(c => {
        const tagCls = tagClassFor(c);
        const pts = formatCritPoints(criterionPoints(c.code));
        const strength = (c.criteria_strength && c.criteria_strength !== c.code)
          ? c.criteria_strength : '';
        const titleAttr = strength ? ` title="${esc(strength).replace(/"/g, '&quot;')}"` : '';
        const ev = c.evidence ? esc(c.evidence) : '';
        const valHTML = ev
          ? `<span class="evrow__val"><span class="nb">${ev}</span></span>`
          : `<span class="evrow__val"><span class="nb">${esc(strength || 'met')}</span></span>`;
        return `<div class="evrow"${titleAttr}>
            <div class="evrow__main">
              <span class="evrow__src">${esc(c.code)}</span>
              ${valHTML}
            </div>
            <span class="evrow__crit ${tagCls}"><span class="cd" aria-hidden="true"></span>${esc(c.code)} ${pts}</span>
          </div>`;
      }).join('');
    }

    function showToast(msg) {
      const t = document.getElementById('toast');
      t.textContent = msg;
      t.classList.add('show');
      setTimeout(() => t.classList.remove('show'), 2200);
    }

    // ── Criteria-card toggle (event delegation) ─────────────────────────────
    // The 28 ACMG criteria cards are rendered into the DOM via innerHTML
    // every time the user generates a result. We CANNOT attach the click
    // handler inside renderResult() — a fresh listener would be added to
    // `document` on every render, and after N renders each click fires N
    // toggles on the same element (even N silently cancels itself out, which
    // looks like "clicks do nothing"). The listener lives here at module
    // scope so it's registered exactly once at page load and persists across
    // arbitrary numbers of renderResult() / innerHTML replacements.

    // Criteria accordion: each card toggles its own inline detail panel.
    window.toggleCritCard = function (hd) {
      const card = hd.closest('.crit-acc');
      if (!card) return;
      const detail = card.querySelector('.crit-acc__detail');
      const isOpen = card.classList.toggle('open');
      if (detail) detail.style.maxHeight = isOpen ? (detail.scrollHeight + 'px') : '0px';
    };
    document.addEventListener('click', function (e) {
      if (e.target.closest('a')) return;  // let source links work without toggling
      const hd = e.target.closest('.crit-acc__hd');
      if (hd) toggleCritCard(hd);
    });

    // Reformat the backend's gnomAD-style variant id (e.g. "3-12604200-C-T")
    // into a colon-separated coordinate like "chr3:12604200:C:T (GRCh38)" —
    // the same chr:pos:ref:alt shape accepted by the variant-input field.
    // Handles SNVs and indels: empty alleles arrive as a literal '-' from the
    // backend, so the regex accepts either nucleotides or '-' for ref/alt.
    // `assembly` is the genome build label from VEP (e.g. "GRCh38"); when
    // unknown the build suffix is omitted rather than guessed.
    // Lowercase the ACMG/AMP classification tier words when they appear
    // mid-sentence (preceded by whitespace, not at position 0). Avoids
    // Claude's "is Pathogenic driven by" mid-sentence capitalisation.
    // Order matters in the alternation: "Likely Pathogenic" first so it
    // wins over a bare "Pathogenic" inside the longer phrase.
    function lowercaseTiersMidSentence(s) {
      if (!s) return s;
      return String(s).replace(
        /(\s)(Likely Pathogenic|Likely pathogenic|Likely Benign|Likely benign|Pathogenic|Benign|VUS)\b/g,
        (_, ws, word) => ws + word.toLowerCase()
      );
    }

    // Map an `effort` enum value to a human label + leading emoji glyph.
    // Unknown values render with a neutral 📝 fallback so layout doesn't
    // break if the backend introduces a new category before the frontend
    // catches up.
    const EFFORT_BADGES = {
      literature_check: { icon: '📄', label: 'Literature' },
      clinical_data:    { icon: '🏥', label: 'Clinical data' },
      lab_work:         { icon: '🔬', label: 'Lab work' },
      family_study:     { icon: '👨‍👩‍👧', label: 'Family study' },
    };
    const PRIORITY_LABEL = { high: 'High priority', medium: 'Medium priority', low: 'Low priority' };

    // Extract the strength tier ("Strong", "Moderate", "Supporting",
    // "VeryStrong") from a criteria_strength string like "PS2_Moderate" or
    // "PVS1" (bare code → VeryStrong for PVS1, otherwise no tier).
    // Returns null when the input doesn't encode a tier worth chipping.
    function strengthTier(criteriaStrength) {
      if (!criteriaStrength || typeof criteriaStrength !== 'string') return null;
      if (criteriaStrength.endsWith('_Strong')) return 'Strong';
      if (criteriaStrength.endsWith('_Moderate')) return 'Moderate';
      if (criteriaStrength.endsWith('_Supporting')) return 'Supporting';
      // Bare codes: PVS1 is Very Strong by definition; BA1 is stand-alone.
      if (criteriaStrength === 'PVS1') return 'VeryStrong';
      return null;
    }

    // Format a gnomAD allele frequency for display.
    // AF<1e-4 → scientific notation (e.g. "6.16 × 10<sup>−6</sup>")
    // AF≥1e-4 → decimal with up to 4 sig figs, trailing zeros stripped
    // null / 0   → "absent"
    function formatAF(af) {
      if (af === null || af === undefined) return 'absent';
      const n = Number(af);
      if (!Number.isFinite(n) || n <= 0) return 'absent';
      if (n < 1e-4) {
        // toExponential(2) normalises the mantissa to [1, 10) AFTER rounding
        // so we never end up with "10.00 × 10⁻⁵" near power-of-ten boundaries.
        const [mant, expRaw] = n.toExponential(2).split('e');
        const exp = parseInt(expRaw, 10);
        const expStr = exp < 0 ? `−${Math.abs(exp)}` : String(exp);
        return `${mant} × 10<sup>${expStr}</sup>`;
      }
      return n.toPrecision(4).replace(/\.?0+$/, '');
    }

    function formatVariantCoord(vid, assembly) {
      if (!vid) return '';
      const m = vid.match(/^([^-]+)-(\d+)-([ACGTN]+|-)-([ACGTN]+|-)$/i);
      if (!m) return vid;  // unexpected shape — surface raw rather than blank
      const [, chrom, pos, ref, alt] = m;
      const build = assembly ? ` (${assembly})` : '';
      return `chr${chrom}:${pos}:${ref}:${alt}${build}`;
    }

    // Shared helpers for the export pipelines.
    // Empty values surface as a human-readable fallback rather than blanks.
    function _orNA(v, fallback) {
      const fb = fallback || 'Not available';
      if (v === null || v === undefined) return fb;
      if (typeof v === 'string' && v.trim() === '') return fb;
      if (Array.isArray(v) && v.length === 0) return fb;
      return v;
    }
    // Plain-text formatVariantCoord (HTML-free version of the on-page helper).
    function _coordPlain(vid, assembly) {
      if (!vid) return '';
      const m = vid.match(/^([^-]+)-(\d+)-([ACGTN]+|-)-([ACGTN]+|-)$/i);
      if (!m) return vid;
      const [, chrom, pos, ref, alt] = m;
      return `chr${chrom}:${pos}:${ref}:${alt}${assembly ? ' (' + assembly + ')' : ''}`;
    }
    function _afPlain(af) {
      if (af === null || af === undefined) return 'Absent';
      const n = Number(af);
      if (!Number.isFinite(n) || n <= 0) return 'Absent';
      if (n < 1e-4) return n.toExponential(2);
      return n.toPrecision(4).replace(/\.?0+$/, '');
    }
    function _fmtNum(x, dp) {
      if (x === null || x === undefined) return 'Not available';
      const n = Number(x);
      if (!Number.isFinite(n)) return 'Not available';
      return n.toFixed(dp);
    }
    // Padded key-value line: "  Label:          value".
    // For labels longer than the column width, falls back to a single
    // space after the colon so the value never abuts the label.
    function _kvLine(label, value, width) {
      const w = width || 18;
      const head = label + ':';
      const padded = head.length >= w ? head + ' ' : head.padEnd(w);
      return '  ' + padded + value;
    }
    const _HR = '──────────────────────────────────────────────';

    // Build the flat per-source list used by the evidence section of both
    // the TXT and CSV exports. Returns an array of section objects, each
    // with { category, source, label, kv: [[label, value], ...], extras }.
    // `extras` is a free-form trailing block (used for the per-record lists
    // in ClinVar / GenCC / PubMed where a single "Field/Value" row would
    // misrepresent the data).
    function _buildEvidenceSections(ev) {
      const sections = [];
      ev = ev || {};
      const sCase = s => {
        if (s === null || s === undefined || s === '') return '';
        return String(s).replace(/_/g, ' ').toLowerCase().replace(/^./, c => c.toUpperCase());
      };

      // ── FUNCTIONAL ANNOTATION ─────────────────────────────────────────
      const vep = ev.vep || {};
      const vepUserTxRow = vep.user_transcript
        ? [['User-specified transcript', vep.user_transcript]]
        : [];
      if (vep.ok) {
        const startStr = Number.isFinite(Number(vep.start)) ? String(Number(vep.start)) : (vep.start || 'Not available');
        const alleleColon = String(vep.allele_string || '').replace('/', ':');
        const coord = `chr${vep.seq_region_name || '?'}:${startStr}`
          + (alleleColon ? `:${alleleColon}` : '')
          + ` (${vep.assembly_name || 'GRCh38'})`;
        // Per-transcript consequence table → flat lines for the TXT/CSV/XLSX
        // export, so the curator's downloaded record carries every transcript's
        // consequence + HGVS, not just the scored one.
        const _txAll = vep.transcript_consequences_all || [];
        const _txExtras = _txAll.length ? [
          ...(vep.consequence_differs_significantly
            ? ['⚠ Consequence differs significantly between transcripts.'] : []),
          'Transcript consequences:',
          ..._txAll.map(t => {
            const flags = [
              t.is_mane_select ? 'MANE Select' : '',
              t.is_mane_plus_clinical ? 'MANE Plus Clinical' : '',
              t.is_picked ? 'scored' : '',
            ].filter(Boolean).join(', ');
            const cq = sCase(t.consequence) || '—';
            const hg = [
              t.hgvsc ? String(t.hgvsc).split(':').pop() : '',
              t.hgvsp ? String(t.hgvsp).split(':').pop() : '',
            ].filter(Boolean).join(' · ');
            return `  ${t.transcript_id}${flags ? ' [' + flags + ']' : ''}: ${cq}${hg ? ' — ' + hg : ''}`;
          }),
        ] : [];
        sections.push({
          category: 'Functional annotation', source: 'Ensembl VEP',
          kv: [
            ...vepUserTxRow,
            ['Coordinates',  coord],
            ['Transcript',   _orNA(vep.transcript_id)],
            ['Consequence',  _orNA(sCase(vep.most_severe_consequence))],
            ['Impact',       _orNA(vep.impact ? String(vep.impact).toUpperCase() : null)],
            ['REVEL',        (vep.revel_score !== null && vep.revel_score !== undefined) ? `${vep.revel_score} — drives PP3/BP4 (calibrated, capped Supporting)` : 'Not available (REVEL scores missense SNVs only)'],
            ['CADD (PHRED)', (vep.cadd_phred !== null && vep.cadd_phred !== undefined) ? `${vep.cadd_phred} — supporting in-silico context` : 'Not available'],
            ['PolyPhen',     vep.polyphen_prediction ? `${sCase(vep.polyphen_prediction)} (${vep.polyphen_score}) — supporting in-silico context` : 'Not available'],
            ['SIFT',         vep.sift_prediction ? `${sCase(vep.sift_prediction)} (${vep.sift_score}) — supporting in-silico context` : 'Not available'],
          ],
          extras: _txExtras,
        });
      } else {
        sections.push({ category: 'Functional annotation', source: 'Ensembl VEP',
          kv: [
            ...vepUserTxRow,
            ['Status', `Lookup failed: ${vep.error || 'unknown'}`],
          ],
        });
      }

      const sa = ev.spliceai || {};
      if (sa.skipped) {
        sections.push({ category: 'Functional annotation', source: 'SpliceAI',
          kv: [['Status', `Not applicable (indel format) — ${sa.reason || 'Broad API does not accept this allele format'}`]] });
      } else if (sa.ok) {
        const maxd = sa.max_delta;
        const t0 = (sa.scores_per_transcript || [])[0] || {};
        const interp = (maxd >= 0.8) ? 'High confidence splice impact'
                     : (maxd >= 0.5) ? 'Moderate confidence splice impact'
                     : (maxd >= 0.2) ? 'Low confidence splice impact'
                     : 'No predicted splice impact';
        sections.push({
          category: 'Functional annotation', source: 'SpliceAI',
          kv: [
            ['Max ΔScore',     _fmtNum(maxd, 2)],
            ['DS_AG',          _fmtNum(t0.DS_AG, 2)],
            ['DS_AL',          _fmtNum(t0.DS_AL, 2)],
            ['DS_DG',          _fmtNum(t0.DS_DG, 2)],
            ['DS_DL',          _fmtNum(t0.DS_DL, 2)],
            ['Interpretation', interp],
          ],
        });
      } else {
        const err = sa.error || 'not available';
        const isNoScores = /no .* scores|does not overlap|not annotated/i.test(err);
        sections.push({ category: 'Functional annotation', source: 'SpliceAI',
          kv: [['Status', isNoScores ? 'No scores for this variant position' : `Lookup failed: ${err}`]] });
      }

      // AlphaMissense — local pre-computed missense pathogenicity.
      const amEv = ev.alphamissense || {};
      if (amEv.available === false) {
        sections.push({ category: 'Functional annotation', source: 'AlphaMissense',
          kv: [['Status', `Not available — ${amEv.reason || 'ALPHAMISSENSE_PATH not configured'}`]] });
      } else if (amEv.not_applicable) {
        sections.push({ category: 'Functional annotation', source: 'AlphaMissense',
          kv: [['Status', 'Not applicable (non-missense)']] });
      } else {
        sections.push({
          category: 'Functional annotation', source: 'AlphaMissense',
          kv: [
            ['Score',           _fmtNum(amEv.score, 3)],
            ['Classification',  String(amEv.classification || '').replace('_', ' ')],
            ['Protein variant', amEv.protein_variant ? `p.${amEv.protein_variant}` : ''],
          ],
        });
      }

      // ── POPULATION FREQUENCY ──────────────────────────────────────────
      const g = ev.gnomad || {};
      if (g.ok) {
        const v = g.variant;
        const c = (g.gene && g.gene.gnomad_constraint) || {};
        const kv = [];
        if (v) {
          const ex = v.exome || {};
          const ge = v.genome || {};
          const exHas = ex.ac !== undefined && ex.ac !== null;
          const geHas = ge.ac !== undefined && ge.ac !== null;
          const af = ex.af ?? ge.af ?? null;
          const acTotal = (exHas ? ex.ac : 0) + (geHas ? ge.ac : 0);
          const anTotal = (exHas ? (ex.an || 0) : 0) + (geHas ? (ge.an || 0) : 0);
          const hom = (exHas ? (ex.ac_hom || 0) : 0) + (geHas ? (ge.ac_hom || 0) : 0);
          kv.push(['Variant ID',       _coordPlain(v.variantId, ev.vep?.assembly_name) || v.variantId]);
          kv.push(['Variant AF',       _afPlain(af)]);
          kv.push(['AC / AN (exome)',  exHas ? `${ex.ac} / ${ex.an}` : 'Not present']);
          kv.push(['AC / AN (genome)', geHas ? `${ge.ac} / ${ge.an}` : 'Not present']);
          kv.push(['AC / AN (total)',  (exHas || geHas)
            ? `${acTotal} / ${anTotal}`
              + ` (exome ${exHas ? `${ex.ac}/${ex.an}` : 'not in dataset'}`
              + ` + genome ${geHas ? `${ge.ac}/${ge.an}` : 'not in dataset'})`
            : 'Not available']);
          kv.push(['Homozygotes',      String(hom)]);
        } else {
          kv.push(['Variant AF',       'Absent (not found in gnomAD)']);
          kv.push(['AC / AN (total)',  'Not available']);
          kv.push(['Homozygotes',      'Not available']);
        }
        kv.push(['pLI',    _fmtNum(c.pLI, 2)]);
        kv.push(['LOEUF',  _fmtNum(c.oe_lof_upper, 2)]);
        kv.push(['mis_z',  _fmtNum(c.mis_z, 2)]);
        sections.push({ category: 'Population frequency', source: 'gnomAD v4', kv });
      } else {
        sections.push({ category: 'Population frequency', source: 'gnomAD v4',
          kv: [['Status', `Lookup failed: ${g.error || 'unknown'}`]] });
      }

      // ── CLINICAL DATABASES ────────────────────────────────────────────
      const cv = ev.clinvar || {};
      if (cv.ok && cv.found) {
        const recs = cv.records || [];
        const extras = recs.slice(0, 20).map(rec =>
          `  ${rec.accession}: ${rec.clinical_significance || '—'} [${rec.review_status || 'no review status'}]`
        );
        sections.push({ category: 'Clinical', source: 'ClinVar',
          kv: [['Submissions (this variant)', String(cv.total_submissions ?? recs.reduce((s, r) => s + (r.number_submitters || 1), 0))]],
          extras,
          csvRecords: recs.map(rec => [rec.accession, rec.clinical_significance || '—']),
        });
      } else if (cv.ok) {
        sections.push({ category: 'Clinical', source: 'ClinVar', kv: [['Status', 'No matching records']] });
      } else {
        sections.push({ category: 'Clinical', source: 'ClinVar',
          kv: [['Status', `Lookup failed: ${cv.error || 'unknown'}`]] });
      }

      // ── PROTEIN & EXPRESSION ──────────────────────────────────────────
      const up = ev.uniprot || {};
      if (up.ok && up.found) {
        const domains = (up.domains || []).map(d => `${d.name} (${d.start}-${d.end})`);
        sections.push({ category: 'Protein & expression', source: 'UniProt',
          kv: [
            ['Accession',        _orNA(up.accession)],
            ['Protein',          `${up.protein_name || 'Not available'} (${up.length} aa)`],
            ['Domains',          domains.length ? domains.join(', ') : 'None annotated'],
            ['Binding sites',    String((up.binding_sites || []).length || 'Not available')],
            ['Natural variants', String(up.natural_variant_count ?? 'Not available')],
          ],
        });
      } else if (up.ok) {
        sections.push({ category: 'Protein & expression', source: 'UniProt',
          kv: [['Status', `No reviewed human entry for ${up.gene}`]] });
      } else {
        sections.push({ category: 'Protein & expression', source: 'UniProt',
          kv: [['Status', `Lookup failed: ${up.error || 'unknown'}`]] });
      }

      // Domain pathogenicity context — PM1 evidence sourced from the
      // local ClinVar DB, scoped to the UniProt domain containing the
      // variant residue. Mirrors what the Protein-tab section renders so
      // the TXT/CSV export captures the PM1 verdict + counts. Top-5
      // variants land in `extras` so they appear in the per-record
      // detail block (same pattern as ClinVar's submission list).
      const dp = ev.domain_plp || {};
      if (!dp.ok) {
        sections.push({ category: 'Protein & expression', source: 'Domain P/LP (ClinVar)',
          kv: [['Status', `Lookup failed: ${dp.error || 'unknown'}`]] });
      } else if (dp.not_applicable) {
        sections.push({ category: 'Protein & expression', source: 'Domain P/LP (ClinVar)',
          kv: [
            ['Status',         'Not applicable'],
            ['Reason',         dp.reason || 'no domain context'],
            ['PM1 assessment', dp.pm1_assessment || 'not_applicable'],
          ] });
      } else {
        const top = dp.top_variants || [];
        sections.push({
          category: 'Protein & expression', source: 'Domain P/LP (ClinVar)',
          kv: [
            ['Domain',         `${dp.domain_name || '(unnamed)'} (aa ${dp.domain_start}-${dp.domain_end})`],
            ['Variant residue', String(dp.variant_position ?? '—')],
            ['P/LP in domain', `${dp.total_plp ?? 0} total (P=${dp.count_P ?? 0}, LP=${dp.count_LP ?? 0})`],
            ['≥2★ P/LP',       dp.has_two_star_plp ? 'Yes' : 'No'],
            ['PM1 status',     dp.pm1_status || 'not_applicable'],
            ['PM1 assessment', dp.pm1_assessment || ''],
          ],
          extras: top.slice(0, 5).map(v =>
            `  ${v.accession || ''}: ${v.tier} ${v.stars}★ — ${v.name || ''}`
          ),
          csvRecords: top.map(v => [v.accession || '', `${v.tier} ${v.stars}*`, v.name || '']),
        });
      }

      const gt = ev.gtex || {};
      // Short display labels for the cardiovascular tissues GTEx returns.
      const GTEX_SHORT = {
        'Heart Left Ventricle': 'Heart LV', 'Heart - Left Ventricle': 'Heart LV',
        'Heart Atrial Appendage': 'Heart AA', 'Heart - Atrial Appendage': 'Heart AA',
        'Artery Aorta': 'Aorta', 'Artery - Aorta': 'Aorta',
        'Artery Coronary': 'Coronary artery', 'Artery - Coronary': 'Coronary artery',
      };
      if (gt.ok && gt.found) {
        const kv = (gt.tissues || []).map(t => {
          const label = GTEX_SHORT[t.tissue_label] || t.tissue_label || 'Tissue';
          const tpm = t.median_tpm;
          return [label, tpm != null ? `${Number(tpm).toFixed(1)} TPM` : 'Not available'];
        });
        kv.push(['Gene ID', `${gt.gencode_id || 'Not available'} (${gt.dataset || 'gtex'})`]);
        sections.push({ category: 'Protein & expression', source: 'GTEx (cardiovascular)', kv });
      } else if (gt.ok) {
        sections.push({ category: 'Protein & expression', source: 'GTEx (cardiovascular)',
          kv: [['Status', `${gt.gene} not in GTEx`]] });
      } else {
        sections.push({ category: 'Protein & expression', source: 'GTEx (cardiovascular)',
          kv: [['Status', `Lookup failed: ${gt.error || 'unknown'}`]] });
      }

      const fhEv = ev.fetal_heart || {};
      if (fhEv.ok && fhEv.found) {
        const fhCells = fhEv.cell_types || [];
        // Surface the top-3 (cell_type × stage) rows by mean × pct so
        // the export captures the biologically dominant signal without
        // dumping the entire heatmap.
        const ranked = fhCells.map(c => Object.assign({}, c, {
          _score: (c.mean_expr || 0) * (c.pct_expressing || 0) / 100,
        })).sort((a, b) => b._score - a._score).slice(0, 3);
        const kvRows = ranked.map(c => [
          `${c.cell_type_label || c.cell_type} · ${c.stage}`,
          `mean log1p(CPM) ${Number(c.mean_expr).toFixed(2)} · ${Number(c.pct_expressing).toFixed(0)}% of ${c.n_cells} cells (${c.band})`,
        ]);
        kvRows.push(['Dataset', `${fhEv.dataset} · ${fhEv.url}`]);
        sections.push({ category: 'Protein & expression', source: 'Fetal heart (Farah 2024)',
          kv: kvRows,
        });
      } else if (fhEv.ok) {
        sections.push({ category: 'Protein & expression', source: 'Fetal heart (Farah 2024)',
          kv: [['Status', `${fhEv.gene || 'Gene'} not detected in fetal heart dataset`]] });
      } else {
        sections.push({ category: 'Protein & expression', source: 'Fetal heart (Farah 2024)',
          kv: [['Status', `Lookup failed: ${fhEv.error || 'unknown'}`]] });
      }

      const pv = ev.protvar || {};
      if (!pv.ok) {
        sections.push({ category: 'Protein & expression', source: 'ProtVar',
          kv: [['Status', `Lookup failed: ${pv.error || 'unknown'}`]] });
      } else if (!pv.applicable) {
        sections.push({ category: 'Protein & expression', source: 'ProtVar',
          kv: [['Status', 'Not applicable (non-missense variant)']] });
      } else if (!pv.found) {
        sections.push({ category: 'Protein & expression', source: 'ProtVar',
          kv: [['Status', 'No mapping returned']] });
      } else {
        const aa = `${pv.ref_aa || '?'}${pv.protein_position ?? '?'}${pv.alt_aa || '?'}`;
        // conservation_score arrives as either a number or {name, score}.
        let consText = 'Not available';
        const cs = pv.conservation_score;
        if (cs !== undefined && cs !== null) {
          if (typeof cs === 'object') {
            consText = cs.score !== undefined ? String(cs.score) + (cs.name ? ` (${cs.name})` : '') : 'Not available';
          } else {
            consText = String(cs);
          }
        }
        sections.push({ category: 'Protein & expression', source: 'ProtVar',
          kv: [
            ['UniProt',       _orNA(pv.uniprot)],
            ['Position',      aa],
            ['Consequence',   _orNA(sCase(pv.consequence))],
            ['Conservation',  consText],
            ['Features',      (pv.feature_types || []).length ? pv.feature_types.join(', ') : 'None'],
          ],
        });
      }

      // ── DISEASE DATABASES & INTERACTIONS ──────────────────────────────
      const cd = ev.chdgene || {};
      if (cd.ok && cd.listed) {
        sections.push({ category: 'Disease & interactions', source: 'CHDgene',
          kv: [
            ['Status',       'Listed (established CHD gene)'],
            ['CHD subtypes', (cd.chd_classification || []).join(', ') || 'Not available'],
            ['Inheritance',  (cd.inheritance || []).join(', ') || 'Not available'],
          ],
        });
      } else if (cd.ok) {
        sections.push({ category: 'Disease & interactions', source: 'CHDgene',
          kv: [['Status', 'Not listed in CHDgene curated list']] });
      } else {
        sections.push({ category: 'Disease & interactions', source: 'CHDgene',
          kv: [['Status', `Lookup failed: ${cd.error || 'unknown'}`]] });
      }

      const pa = ev.panelapp || {};
      if (pa.ok && pa.on_cardiovascular_panel) {
        const panels = pa.panels_found || [];
        const hpoTerms = pa.submitted_hpo || [];
        const kvRows = [
          ['Total panels',         String(pa.total_panels ?? panels.length)],
          ['Green / amber / red',  `${pa.green_panels || 0} / ${pa.amber_panels || 0} / ${pa.red_panels || 0}`],
          ['Submitted phenotype',  hpoTerms.length ? hpoTerms.join(', ') : 'None'],
          ['PP4 status',           pa.any_green_with_hpo_match ? 'Green panel + phenotype match — supports PP4'
                                  : pa.any_green_no_hpo_match    ? 'Green panel found, phenotype match required'
                                  : 'No green panel — PP4 not supported'],
        ];
        const extras = panels.map(m => {
          const pp4 = m.contributes_to_pp4 ? '✓ PP4'
                    : m.confidence !== 'green' ? '— ' + m.confidence
                    : '⚠ Phenotype required';
          return `  [${m.category}] ${m.panel_name}: ${m.confidence} (${m.confidence_level || '?'}/3) · MOI ${m.moi || 'unspecified'} · ${pp4}`;
        });
        // Per-panel CSV rows: Field = panel name, Value = headline.
        const csvRecords = panels.map(m => [
          m.panel_name || 'Panel',
          `${m.category} · ${m.confidence} · MOI ${m.moi || 'unspecified'} · ${m.contributes_to_pp4 ? 'contributes to PP4' : 'no PP4 contribution'}`,
        ]);
        sections.push({ category: 'Disease & interactions', source: 'PanelApp',
          kv: kvRows, extras, csvRecords });
      } else if (pa.ok) {
        sections.push({ category: 'Disease & interactions', source: 'PanelApp',
          kv: [['Status', `${pa.gene || 'Gene'} not on any cardiovascular panel`]] });
      } else {
        sections.push({ category: 'Disease & interactions', source: 'PanelApp',
          kv: [['Status', `Lookup failed: ${pa.error || 'unknown'}`]] });
      }

      const gc = ev.gencc || {};
      if (gc.ok && gc.found) {
        const subs = (gc.submissions || []).slice(0, 10);
        const extras = subs.map(s =>
          `  ${s.hpo_match ? '✓ ' : ''}${s.classification || '?'}: ${s.disease || '?'} (${s.moi || 'MoI ?'}) — ${s.submitter || '?'}`
        );
        const phenoCount = gc.phenotype_matched_count || 0;
        const phenoBest = gc.phenotype_matched_best_classification;
        const kvRows = [
          ['Phenotype-matched best',
            phenoCount
              ? `${phenoBest || '—'} (across ${phenoCount} matched submission(s))`
              : (gc.submitted_hpo
                  ? 'None — no curated disease matched the proband phenotype'
                  : 'Not assessed — no proband HPO submitted')],
          ['Overall best (all diseases)', _orNA(gc.best_classification)],
          ['Submissions',                  String(gc.submission_count ?? subs.length)],
          ['Disputed / refuted (matched)',
            phenoCount
              ? (gc.phenotype_matched_has_disputed_or_refuted ? 'Yes' : 'No')
              : '—'],
          ['Disputed / refuted (any disease)', gc.has_disputed_or_refuted ? 'Yes' : 'No'],
        ];
        sections.push({ category: 'Disease & interactions', source: 'GenCC',
          kv: kvRows,
          extras,
          csvRecords: subs.map(s => [s.classification || '?', `${s.disease || '?'} (${s.submitter || '?'})${s.hpo_match ? ' [phenotype match]' : ''}`]),
        });
      } else if (gc.ok) {
        sections.push({ category: 'Disease & interactions', source: 'GenCC',
          kv: [['Status', `${gc.gene || 'Gene'} not in GenCC`]] });
      } else {
        sections.push({ category: 'Disease & interactions', source: 'GenCC',
          kv: [['Status', `Lookup failed: ${gc.error || 'unknown'}`]] });
      }

      const mgi = ev.mgi || {};
      if (mgi.ok && mgi.found) {
        const cardiac = mgi.cardiac_phenotypes || [];
        sections.push({ category: 'Disease & interactions', source: 'MGI',
          kv: [
            ['Mouse orthologue',   _orNA(mgi.mouse_symbol)],
            ['Phenotypic alleles', String(mgi.phenotype_count ?? 'Not available')],
            ['Cardiac phenotypes', cardiac.length ? cardiac.join(', ') : 'None detected'],
          ],
        });
      } else if (mgi.ok) {
        sections.push({ category: 'Disease & interactions', source: 'MGI',
          kv: [['Status', mgi.reason || 'No mouse orthologue indexed']] });
      } else {
        sections.push({ category: 'Disease & interactions', source: 'MGI',
          kv: [['Status', `Lookup failed: ${mgi.error || 'unknown'}`]] });
      }

      const bg = ev.biogrid || {};
      if (bg.ok && bg.total_interactions) {
        const top = (bg.top_partners || []).map(p =>
          `${p.symbol} (n=${p.publication_count})`
        );
        const chd = (bg.chd_interactors_all || []).map(p => p.symbol || p);
        sections.push({ category: 'Disease & interactions', source: 'BioGRID',
          kv: [
            ['Total interactions', `${bg.total_interactions} curated, ${bg.unique_partners} unique partners`],
            ['Top partners',       top.length ? top.join(', ') : 'Not available'],
            ['CHD gene partners',  chd.length ? chd.join(', ') : 'None'],
          ],
        });
      } else if (bg.ok) {
        sections.push({ category: 'Disease & interactions', source: 'BioGRID',
          kv: [['Status', 'No curated interactions reported']] });
      } else {
        sections.push({ category: 'Disease & interactions', source: 'BioGRID',
          kv: [['Status', `Lookup failed: ${bg.error || 'unknown'}`]] });
      }

      // ── LITERATURE ────────────────────────────────────────────────────
      const pm = ev.pubmed || {};
      const fmtPaper = p => {
        const meta = `${p.first_author || '?'} ${p.year || ''}`.trim();
        return `  • "${p.title || '(no title)'}" — ${meta} (PMID ${p.pmid})`;
      };
      // Variant-specific PubMed (ev.pubmed.variant_papers).
      if (pm.ok) {
        const v = pm.variant_papers || [];
        const _pmTotal = (typeof pm.variant_total === 'number' && pm.variant_total >= v.length) ? pm.variant_total : v.length;
        sections.push({ category: 'Literature', source: 'PubMed — Variant-specific',
          kv: [['Total papers (PubMed)', String(_pmTotal)], ['Listed (top)', String(v.length)]],
          extras: v.length ? v.map(fmtPaper) : ['  No indexed variant-specific publications'],
          csvRecords: v.map(p => [`PMID ${p.pmid}`, `${p.first_author || '?'} ${p.year || ''} — ${p.title || ''}`]),
        });
      } else {
        sections.push({ category: 'Literature', source: 'PubMed — Variant-specific',
          kv: [['Status', `Lookup failed: ${pm.error || 'unknown'}`]] });
      }
      // Gene-level literature — the phenotype-matched gene papers shown in the
      // Gene tab's "Gene literature" section. Sourced from ev.gene_literature
      // (NOT ev.pubmed.gene_papers, which the backend never populates) so the
      // exported count matches what the UI displays.
      const glit = ev.gene_literature || {};
      if (glit.ok) {
        const gp = glit.papers || [];
        const diseaseLbl = (glit.phenotype_keywords || []).filter(Boolean).join(' / ') || 'Cardiac';
        sections.push({ category: 'Literature', source: 'PubMed — Gene-level',
          kv: [
            ['Phenotype filter', diseaseLbl],
            ['Indexed papers', String(gp.length)],
          ],
          extras: gp.length ? gp.map(fmtPaper) : ['  No phenotype-matched gene-level publications'],
          csvRecords: gp.map(p => [`PMID ${p.pmid}`, `${p.first_author || '?'} ${p.year || ''} — ${p.title || ''}`]),
        });
      } else {
        sections.push({ category: 'Literature', source: 'PubMed — Gene-level',
          kv: [['Status', glit.error ? `Lookup failed: ${glit.error}` : 'Not loaded']] });
      }

      return sections;
    }

    // Builds the plain-text report shared by Export-as-TXT and Copy-to-clipboard.
    function _buildTxtReport() {
      if (!lastResult) return '';
      const r = lastResult;
      const ev = lastEvidence || {};
      const lines = [];
      lines.push('Cardiac Variant Curation Report');
      lines.push(new Date().toLocaleString());
      lines.push('');
      lines.push(`Variant: ${lastVariant.gene} ${lastVariant.hgvs_c}`);
      // Query-context lines mirror the Summary tab's Query-details table
      // so screenshots and text reports stay in sync. Only emit a line
      // when the curator actually supplied the field.
      if (lastVariant.genome_build) lines.push(`Genome build: ${lastVariant.genome_build}`);
      if (lastVariant.hpo)          lines.push(`Phenotype / HPO: ${lastVariant.hpo}`);
      // Clinical-context lines — only emit fields the curator actually
      // supplied. Zygosity is annotated with "(inferred)" when the
      // backend filled it in from chrX + male sex.
      if (lastVariant.zygosity || lastVariant.zygosity_inferred) {
        const zVal = _clinicalLabel('zygosity', lastVariant.zygosity);
        const suffix = lastVariant.zygosity_inferred ? ' (inferred)' : '';
        lines.push(`Zygosity: ${zVal}${suffix}`);
      }
      if (lastVariant.inheritance_input) lines.push(`Inheritance: ${_clinicalLabel('inheritance', lastVariant.inheritance_input)}`);
      if (lastVariant.proband_sex)       lines.push(`Proband sex: ${_clinicalLabel('sex', lastVariant.proband_sex)}`);
      if (lastVariant.trio_status) {
        const tVal = _clinicalLabel('trio', lastVariant.trio_status);
        const dVal = lastVariant.denovo_status ? ` — ${_clinicalLabel('denovo', lastVariant.denovo_status)}` : '';
        lines.push(`Trio: ${tVal}${dVal}`);
      }
      if (lastVariant.family)       lines.push(`Family history: ${lastVariant.family}`);
      // Structured family & segregation evidence — only emit fields the
      // curator actually supplied (counts > 0 / non-empty selects).
      if (lastVariant.seg_affected_carriers)    lines.push(`Affected relatives carrying variant: ${lastVariant.seg_affected_carriers}`);
      if (lastVariant.seg_affected_noncarriers) lines.push(`Affected relatives NOT carrying: ${lastVariant.seg_affected_noncarriers}`);
      if (lastVariant.seg_meioses)              lines.push(`Informative meioses: ${lastVariant.seg_meioses}`);
      if (lastVariant.in_trans_pathogenic)      lines.push(`2nd pathogenic variant in trans: ${lastVariant.in_trans_pathogenic}`);
      if (lastVariant.alt_cause_present)        lines.push(`Alternate molecular cause: ${lastVariant.alt_cause_present}${lastVariant.alt_cause_detail ? ` (${lastVariant.alt_cause_detail})` : ''}`);
      lines.push(`Classification (HeartVar): ${r.classification}`);
      lines.push(`Confidence: ${r.confidence}`);
      if (Number.isFinite(r.points_total)) {
        lines.push(`Points (Tavtigian 2020): ${r.points_total >= 0 ? '+' : ''}${r.points_total}`);
      }
      // Curator overrides. Emitted whenever the Criteria tab was edited, and
      // NOT optional: a report stating a tier without saying a human forced it
      // would be a misleading clinical document. HeartVar's own call stays
      // above, labelled, so the two are never confusable.
      lines.push.apply(lines, _hvOverrideReportLines());
      // Patient-level carrier-status interpretation (informational; recessive /
      // X-linked / dual-MOI genes). Server-computed, not curator-supplied, so
      // emit unconditionally when present.
      if (ev.carrier_status && ev.carrier_status.detail) {
        lines.push(`Carrier status: ${ev.carrier_status.title} — ${ev.carrier_status.detail}`);
      }
      lines.push('');
      lines.push('Summary');
      if (Array.isArray(r.summary)) {
        r.summary.forEach(s => lines.push(`- ${s}`));
      } else if (r.summary) {
        lines.push(String(r.summary));
      }

      // ── OBSERVED EVIDENCE ──────────────────────────────────────────────
      lines.push('');
      lines.push('OBSERVED EVIDENCE');
      lines.push(_HR);

      const sections = _buildEvidenceSections(ev);
      const CATEGORIES = [
        ['Functional annotation',   'FUNCTIONAL ANNOTATION'],
        ['Population frequency',    'POPULATION FREQUENCY'],
        ['Clinical',                'CLINICAL DATABASES'],
        ['Protein & expression',    'PROTEIN & EXPRESSION'],
        ['Disease & interactions',  'DISEASE DATABASES & INTERACTIONS'],
        ['Literature',              'LITERATURE'],
      ];
      for (const [key, heading] of CATEGORIES) {
        const inCat = sections.filter(s => s.category === key);
        if (!inCat.length) continue;
        lines.push('');
        lines.push(heading);
        lines.push('');
        inCat.forEach((s, i) => {
          if (i > 0) lines.push('');
          lines.push(s.source);
          s.kv.forEach(([label, value]) => lines.push(_kvLine(label, value)));
          if (s.extras && s.extras.length) {
            lines.push('');
            s.extras.forEach(x => lines.push(x));
          }
        });
        lines.push('');
        lines.push(_HR);
      }

      // ── ACMG/AMP CRITERIA ──────────────────────────────────────────────
      lines.push('ACMG/AMP CRITERIA');
      lines.push('');
      const crits = r.criteria || [];
      const met   = crits.filter(c => c.status === 'met');
      const insuf = crits.filter(c => c.status === 'insufficient_data');
      const naMet = crits.filter(c => c.status === 'na' || c.status === 'not_met');
      met.forEach(c => {
        const tier = c.criteria_strength && c.criteria_strength !== c.code ? ` [${c.criteria_strength}]` : '';
        lines.push(`${c.code}${tier}: ${c.evidence}`);
      });
      if (insuf.length) {
        lines.push('');
        insuf.forEach(c => lines.push(`${c.code} [Insufficient data]: ${c.evidence}`));
      }
      if (naMet.length) {
        lines.push('');
        const codes = naMet.map(c => c.code).join(', ');
        lines.push(`Not met / not applicable: ${codes}`);
      }
      lines.push('');
      if (Number.isFinite(r.points_total)) {
        lines.push(`Total points: ${r.points_total >= 0 ? '+' : ''}${r.points_total}`);
      }
      lines.push('Applied guidelines: Richards et al. 2015; Tavtigian et al. 2020');

      return lines.join('\n');
    }

    function exportTXT() {
      if (!lastResult) return;
      try {
        const text = _buildTxtReport();
        const blob = new Blob([text], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${lastVariant.gene}_${lastVariant.hgvs_c}_curation.txt`
          .replace(/[^a-zA-Z0-9_.-]/g, '_');
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        showToast('Downloaded as TXT');
      } catch (err) {
        showToast('Download failed — could not build the report');
      }
    }

    async function copyReportToClipboard() {
      if (!lastResult) return;
      const text = _buildTxtReport();
      // Primary path: the async Clipboard API, which only exists in a secure
      // context (HTTPS or localhost). When the page is served over plain HTTP
      // navigator.clipboard is undefined, so fall back to a hidden-textarea
      // execCommand('copy') so the button still works everywhere.
      try {
        if (navigator.clipboard && navigator.clipboard.writeText) {
          await navigator.clipboard.writeText(text);
          showToast('Copied to clipboard');
          return;
        }
      } catch (err) { /* fall through to the legacy copy path */ }
      try {
        const ta = document.createElement('textarea');
        ta.value = text;
        ta.setAttribute('readonly', '');
        ta.style.position = 'fixed';
        ta.style.top = '-1000px';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        const ok = document.execCommand('copy');
        document.body.removeChild(ta);
        showToast(ok ? 'Copied to clipboard' : 'Copy failed — clipboard unavailable');
      } catch (err) {
        showToast('Copy failed — clipboard unavailable');
      }
    }

    // ── Minimal dependency-free .xlsx writer ──────────────────────────
    // Builds a valid OOXML SpreadsheetML workbook (multiple sheets, styled
    // headers, column widths, frozen header rows, colour-coded cells) and
    // packages it as an uncompressed ("stored") ZIP — no external library
    // or CDN required. Text cells use inline strings; pass {n:true} for a
    // numeric cell and {s:<index>} for a style from _XLSX_STYLES.

    // CRC-32 (IEEE) table + checksum over a Uint8Array — required for each
    // stored ZIP entry.
    const _XLSX_CRC_TABLE = (() => {
      const t = new Uint32Array(256);
      for (let n = 0; n < 256; n++) {
        let c = n;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
        t[n] = c >>> 0;
      }
      return t;
    })();
    function _crc32(bytes) {
      let c = 0xFFFFFFFF;
      for (let i = 0; i < bytes.length; i++) c = _XLSX_CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
      return (c ^ 0xFFFFFFFF) >>> 0;
    }

    // Package an array of {name, bytes:Uint8Array} into a stored (no
    // compression) ZIP and return the bytes. Little-endian throughout.
    function _zipStore(files) {
      const enc = new TextEncoder();
      const u16 = n => [n & 0xFF, (n >>> 8) & 0xFF];
      const u32 = n => [n & 0xFF, (n >>> 8) & 0xFF, (n >>> 16) & 0xFF, (n >>> 24) & 0xFF];
      const chunks = [];
      const central = [];
      let offset = 0;
      files.forEach(f => {
        const name = enc.encode(f.name);
        const data = f.bytes;
        const crc = _crc32(data);
        // Local file header (date fixed to 1980-01-01 → DOS date 33).
        const local = [].concat(
          u32(0x04034b50), u16(20), u16(0), u16(0), u16(0), u16(33),
          u32(crc), u32(data.length), u32(data.length),
          u16(name.length), u16(0));
        chunks.push(new Uint8Array(local), name, data);
        const cd = [].concat(
          u32(0x02014b50), u16(20), u16(20), u16(0), u16(0), u16(0), u16(33),
          u32(crc), u32(data.length), u32(data.length),
          u16(name.length), u16(0), u16(0), u16(0), u16(0), u32(0),
          u32(offset));
        central.push(new Uint8Array(cd), name);
        offset += local.length + name.length + data.length;
      });
      const cdStart = offset;
      let cdSize = 0;
      central.forEach(c => { chunks.push(c); cdSize += c.length; });
      const eocd = [].concat(
        u32(0x06054b50), u16(0), u16(0),
        u16(files.length), u16(files.length),
        u32(cdSize), u32(cdStart), u16(0));
      chunks.push(new Uint8Array(eocd));
      let total = 0;
      chunks.forEach(c => total += c.length);
      const out = new Uint8Array(total);
      let p = 0;
      chunks.forEach(c => { out.set(c, p); p += c.length; });
      return out;
    }

    function _xmlEsc(s) {
      return String(s == null ? '' : s)
        .replace(/[\u0000-\u0008\u000B\u000C\u000E-\u001F]/g, '')
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
    }

    // 0-based column index → spreadsheet letter (0→A, 26→AA …).
    function _colLetter(i) {
      let s = '';
      i += 1;
      while (i > 0) { const m = (i - 1) % 26; s = String.fromCharCode(65 + m) + s; i = Math.floor((i - 1) / 26); }
      return s;
    }

    // Shared style table. Indices referenced by cells via {s:n}:
    //   1 bold label · 2 title · 3 subtitle · 4 header (navy/white) ·
    //   5 wrapped text · 6/7/8 classification fill (patho/benign/VUS) ·
    //   9/10 criterion-met fill (pathogenic/benign) · 11 not-assessed.
    const _XLSX_STYLES =
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
      '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' +
        '<fonts count="6">' +
          '<font><sz val="11"/><name val="Calibri"/></font>' +
          '<font><b/><sz val="11"/><name val="Calibri"/></font>' +
          '<font><b/><sz val="15"/><color rgb="FF1F3A5F"/><name val="Calibri"/></font>' +
          '<font><i/><sz val="10"/><color rgb="FF6B6B6B"/><name val="Calibri"/></font>' +
          '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>' +
          '<font><i/><sz val="11"/><color rgb="FF8A8A8A"/><name val="Calibri"/></font>' +
        '</fonts>' +
        '<fills count="9">' +
          '<fill><patternFill patternType="none"/></fill>' +
          '<fill><patternFill patternType="gray125"/></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FF1F3A5F"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FFB23A2E"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FF2E7D4F"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FF2C5282"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FFFBE9E7"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FFE6F4EA"/><bgColor indexed="64"/></patternFill></fill>' +
          '<fill><patternFill patternType="solid"><fgColor rgb="FFF1F1F1"/><bgColor indexed="64"/></patternFill></fill>' +
        '</fills>' +
        '<borders count="2">' +
          '<border><left/><right/><top/><bottom/><diagonal/></border>' +
          '<border><left/><right/><top/><bottom style="thin"><color rgb="FFBFBFBF"/></bottom><diagonal/></border>' +
        '</borders>' +
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' +
        '<cellXfs count="12">' +
          '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>' +
          '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>' +
          '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/>' +
          '<xf numFmtId="0" fontId="3" fillId="0" borderId="0" xfId="0" applyFont="1"/>' +
          '<xf numFmtId="0" fontId="4" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1"><alignment vertical="center"/></xf>' +
          '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf>' +
          '<xf numFmtId="0" fontId="4" fillId="3" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment horizontal="left"/></xf>' +
          '<xf numFmtId="0" fontId="4" fillId="4" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment horizontal="left"/></xf>' +
          '<xf numFmtId="0" fontId="4" fillId="5" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment horizontal="left"/></xf>' +
          '<xf numFmtId="0" fontId="0" fillId="6" borderId="0" xfId="0" applyFill="1"/>' +
          '<xf numFmtId="0" fontId="0" fillId="7" borderId="0" xfId="0" applyFill="1"/>' +
          '<xf numFmtId="0" fontId="5" fillId="8" borderId="0" xfId="0" applyFont="1" applyFill="1"/>' +
        '</cellXfs>' +
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' +
      '</styleSheet>';

    // Render one sheet object → worksheet XML.
    //   sheet = { name, rows:[[cell,…],…], cols:[width,…],
    //             merges:["A1:B1",…], freezeHeader:bool }
    //   cell  = primitive | { v, s:<styleIdx>, n:<numeric?> }
    function _sheetXml(sheet) {
      const colsXml = (sheet.cols && sheet.cols.length)
        ? '<cols>' + sheet.cols.map((w, i) => `<col min="${i + 1}" max="${i + 1}" width="${w}" customWidth="1"/>`).join('') + '</cols>'
        : '';
      const rowsXml = (sheet.rows || []).map((row, ri) => {
        const r = ri + 1;
        const cells = (row || []).map((raw, ci) => {
          if (raw === null || raw === undefined) return '';
          const cell = (typeof raw === 'object') ? raw : { v: raw };
          const ref = _colLetter(ci) + r;
          const sAttr = cell.s ? ` s="${cell.s}"` : '';
          if (cell.n) return `<c r="${ref}"${sAttr}><v>${_xmlEsc(cell.v)}</v></c>`;
          if (cell.v === '' || cell.v === null || cell.v === undefined) return `<c r="${ref}"${sAttr}/>`;
          return `<c r="${ref}"${sAttr} t="inlineStr"><is><t xml:space="preserve">${_xmlEsc(cell.v)}</t></is></c>`;
        }).join('');
        return `<row r="${r}">${cells}</row>`;
      }).join('');
      const view = sheet.freezeHeader
        ? '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/><selection pane="bottomLeft"/></sheetView></sheetViews>'
        : '<sheetViews><sheetView workbookViewId="0"/></sheetViews>';
      const merges = (sheet.merges && sheet.merges.length)
        ? `<mergeCells count="${sheet.merges.length}">` + sheet.merges.map(m => `<mergeCell ref="${m}"/>`).join('') + '</mergeCells>'
        : '';
      return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' +
        view + '<sheetFormatPr defaultRowHeight="15"/>' + colsXml +
        '<sheetData>' + rowsXml + '</sheetData>' + merges + '</worksheet>';
    }

    // Assemble the full workbook package and return a Blob.
    function _xlsxBlob(sheets) {
      const enc = new TextEncoder();
      const files = [];
      const add = (name, str) => files.push({ name, bytes: enc.encode(str) });

      let ct = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' +
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' +
        '<Default Extension="xml" ContentType="application/xml"/>' +
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' +
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>';
      sheets.forEach((s, i) => {
        ct += `<Override PartName="/xl/worksheets/sheet${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>`;
      });
      ct += '</Types>';
      add('[Content_Types].xml', ct);

      add('_rels/.rels', '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' +
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>' +
        '</Relationships>');

      let wb = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>';
      sheets.forEach((s, i) => { wb += `<sheet name="${_xmlEsc(s.name)}" sheetId="${i + 1}" r:id="rId${i + 1}"/>`; });
      wb += '</sheets></workbook>';
      add('xl/workbook.xml', wb);

      let wr = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">';
      sheets.forEach((s, i) => {
        wr += `<Relationship Id="rId${i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet${i + 1}.xml"/>`;
      });
      wr += `<Relationship Id="rId${sheets.length + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>`;
      add('xl/_rels/workbook.xml.rels', wr);

      add('xl/styles.xml', _XLSX_STYLES);
      sheets.forEach((s, i) => add(`xl/worksheets/sheet${i + 1}.xml`, _sheetXml(s)));

      return new Blob([_zipStore(files)], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
    }

    // Export the current curation as a structured multi-sheet workbook:
    //   Summary       — variant identity, classification, clinical context,
    //                   interpretation notes (key/value layout).
    //   Evidence      — every data source as Source / Category / Field / Value.
    //   ACMG Criteria — one row per criterion, colour-coded by status.
    function exportXLSX() {
      if (!lastResult) return;
      try {
        const r = lastResult;
        const ev = lastEvidence || {};
        const v = lastVariant || {};
        const summaryArr = Array.isArray(r.summary) ? r.summary : (r.summary ? [String(r.summary)] : []);

        const pad = n => String(n).padStart(2, '0');
        const now = new Date();
        const dateStr = `${now.getFullYear()}-${pad(now.getMonth() + 1)}-${pad(now.getDate())} ${pad(now.getHours())}:${pad(now.getMinutes())}`;

        const protein = _stripHgvsPrefix(window.heartvarProteinData && window.heartvarProteinData.hgvsp)
          || ((ev.alphamissense && ev.alphamissense.protein_variant) ? 'p.' + ev.alphamissense.protein_variant : '');
        const cls = r.classification || '';
        const classFill = (() => {
          const c = cls.toLowerCase();
          if (c.includes('pathogenic')) return 6;
          if (c.includes('benign')) return 7;
          return 8;
        })();
        const stripStrength = s => (s ? String(s).split('_').pop() : '');
        const applied = (r.criteria || [])
          .filter(c => c.status === 'met')
          .map(c => c.code + (c.criteria_strength ? ` (${stripStrength(c.criteria_strength)})` : ''))
          .join(', ');

        // ── Sheet 1 — Summary ─────────────────────────────────────────
        const S = [];
        const merges1 = [];
        const hdr = label => { S.push([{ v: label, s: 4 }, { v: '', s: 4 }]); merges1.push(`A${S.length}:B${S.length}`); };
        const kvp = (label, value, vs) => {
          if (value === '' || value === null || value === undefined) return;
          S.push([{ v: label, s: 1 }, { v: value, s: vs || 0 }]);
        };
        const titleVariant = `${v.gene || ''} ${v.hgvs_c || ''}`.trim() + (protein ? `  (${protein})` : '');
        S.push([{ v: titleVariant, s: 2 }]); merges1.push('A1:B1');
        S.push([{ v: `HeartVar variant curation · generated ${dateStr}`, s: 3 }]); merges1.push('A2:B2');
        S.push([]);

        hdr('Classification');
        S.push([{ v: 'Result (HeartVar)', s: 1 }, { v: cls || '—', s: classFill }]);
        kvp('Confidence', r.confidence || '');
        if (Number.isFinite(r.points_total)) S.push([{ v: 'Points (Tavtigian 2020)', s: 1 }, { v: r.points_total, n: true }]);
        kvp('Criteria applied', applied, 5);
        // Curator overrides, when the Criteria tab was edited. Rows are added
        // rather than the HeartVar result being rewritten, so a reader of the
        // workbook can always tell the engine's call from the human's.
        {
          const _ovrList = (typeof __hvOverrideSummaryList === 'function') ? __hvOverrideSummaryList() : [];
          const _adj = (window.__hvOverrideState || {}).adjusted;
          if (_ovrList.length && _adj) {
            S.push([{ v: 'Result (curator-adjusted)', s: 1 }, { v: _adj.classification || '—', s: classFill }]);
            S.push([{ v: 'Points (curator-adjusted)', s: 1 }, { v: _adj.points_total, n: true }]);
            kvp('Curator changes', String(_ovrList.length), 5);
            _ovrList.forEach(o => kvp('  ' + o.code, o.text, 5));
          }
        }
        S.push([]);

        hdr('Variant');
        kvp('Gene', v.gene || '');
        kvp('HGVS (coding)', v.hgvs_c || '');
        kvp('Protein change', protein);
        kvp('Genome build', v.genome_build || '');
        const vep = ev.vep || {};
        if (vep.ok) {
          const startStr = Number.isFinite(Number(vep.start)) ? String(Number(vep.start)) : (vep.start || '');
          const alleleColon = String(vep.allele_string || '').replace('/', ':');
          kvp('Coordinates', `chr${vep.seq_region_name || '?'}:${startStr}${alleleColon ? ':' + alleleColon : ''} (${vep.assembly_name || 'GRCh38'})`);
          kvp('Transcript', vep.transcript_id || '');
          kvp('Consequence', vep.most_severe_consequence ? String(vep.most_severe_consequence).replace(/_/g, ' ').replace(/^./, c => c.toUpperCase()) : '');
          kvp('Impact', vep.impact ? String(vep.impact).toUpperCase() : '');
        }
        S.push([]);

        hdr('Clinical context');
        kvp('Phenotype / HPO', v.hpo || '');
        kvp('Zygosity', v.zygosity ? _clinicalLabel('zygosity', v.zygosity) + (v.zygosity_inferred ? ' (inferred)' : '') : '');
        kvp('Inheritance', v.inheritance_input ? _clinicalLabel('inheritance', v.inheritance_input) : '');
        kvp('Proband sex', v.proband_sex ? _clinicalLabel('sex', v.proband_sex) : '');
        kvp('Trio status', v.trio_status ? _clinicalLabel('trio', v.trio_status) : '');
        kvp('De novo status', v.denovo_status ? _clinicalLabel('denovo', v.denovo_status) : '');
        kvp('Family history', v.family || '', 5);
        kvp('Affected relatives carrying', v.seg_affected_carriers || '');
        kvp('Affected relatives not carrying', v.seg_affected_noncarriers || '');
        kvp('Informative meioses', v.seg_meioses || '');
        kvp('2nd pathogenic variant in trans', v.in_trans_pathogenic || '');
        kvp('Alternate molecular cause', v.alt_cause_present || '');
        kvp('Alternate cause detail', v.alt_cause_detail || '', 5);
        // Patient-level carrier-status interpretation (informational; recessive /
        // X-linked / dual-MOI genes). Server-computed, so emit when present.
        if (ev.carrier_status && ev.carrier_status.detail) {
          kvp('Carrier status', `${ev.carrier_status.title} — ${ev.carrier_status.detail}`, 5);
        }

        const bullets = summaryArr.filter(Boolean);
        if (bullets.length) {
          S.push([]);
          hdr('Interpretation notes');
          bullets.forEach(b => { S.push([{ v: String(b), s: 5 }]); merges1.push(`A${S.length}:B${S.length}`); });
        }

        // ── Sheet 2 — Evidence ────────────────────────────────────────
        const E = [[{ v: 'Source', s: 4 }, { v: 'Category', s: 4 }, { v: 'Field', s: 4 }, { v: 'Value', s: 4 }]];
        _buildEvidenceSections(ev).forEach(sec => {
          (sec.kv || []).forEach(([field, value]) => {
            E.push([{ v: sec.source || '' }, { v: sec.category || '' }, { v: field || '' }, { v: String(value == null ? '' : value), s: 5 }]);
          });
          (sec.csvRecords || []).forEach(([field, value]) => {
            E.push([{ v: sec.source || '' }, { v: sec.category || '' }, { v: String(field == null ? '' : field) }, { v: String(value == null ? '' : value), s: 5 }]);
          });
        });

        // ── Sheet 3 — ACMG criteria ───────────────────────────────────
        const C = [[{ v: 'Code', s: 4 }, { v: 'Status', s: 4 }, { v: 'Strength', s: 4 }, { v: 'Direction', s: 4 }, { v: 'Evidence', s: 4 }]];
        (r.criteria || []).forEach(c => {
          const st = c.status || '';
          let stStyle = 0;
          if (st === 'met') stStyle = (c.direction === 'benign') ? 10 : 9;
          else if (st === 'not_assessed' || st === 'not_evaluated') stStyle = 11;
          C.push([
            { v: c.code || '', s: 1 },
            { v: st.replace(/_/g, ' '), s: stStyle },
            { v: stripStrength(c.criteria_strength) },
            { v: c.direction || '' },
            { v: c.evidence || '', s: 5 },
          ]);
        });

        const blob = _xlsxBlob([
          { name: 'Summary', rows: S, cols: [30, 64], merges: merges1 },
          { name: 'Evidence', rows: E, cols: [22, 18, 28, 70], freezeHeader: true },
          { name: 'ACMG Criteria', rows: C, cols: [10, 14, 13, 12, 80], freezeHeader: true },
        ]);
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${v.gene || 'variant'}_${(v.hgvs_c || '').replace(/[^a-zA-Z0-9]/g, '_')}_curation.xlsx`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        showToast('Downloaded as Excel');
      } catch (err) {
        showToast('Download failed — could not build the workbook');
      }
    }

    // ── Streaming UI helpers ───────────────────────────────────────────────
    // Loading-screen variant. Switch to 'original' to restore the verbose
    // per-source checklist + simulator. The 'terminal' branch renders a
    // retro mono terminal scrolling through a pre-scripted pipeline log,
    // patched with real values as the SSE stream lands them.
    const LOADING_SCREEN = 'original'; // options: 'terminal' | 'original'

    // Ordered alphabetically by display name so the "Querying databases…"
    // loading panel reads top-to-bottom in predictable order.
    // Keys MUST match the backend db_pending/db_done `source` strings exactly
    // (see backend/evidence.py `sources`), or the row never resolves and the
    // completion gate (dbCompletion, below) hangs. Note the non-obvious keys:
    // Open Targets → 'opentargets_evidence', PMC Open Access → 'pmcoa',
    // Heart of Fetal Cells → 'fetal_heart'. AlphaFold reports whether a
    // bundled 3-D model exists (the model streams later on the protein tab).
    const DB_SOURCES = [
      { key: 'alphafold', label: 'AlphaFold' },
      { key: 'alphamissense', label: 'AlphaMissense' },
      { key: 'biogrid',  label: 'BioGRID' },
      { key: 'chdgene',  label: 'CHDgene' },
      { key: 'clinvar',  label: 'ClinVar' },
      { key: 'vep',      label: 'Ensembl VEP' },
      { key: 'gencc',    label: 'GenCC' },
      { key: 'gnomad',   label: 'gnomAD v4' },
      { key: 'gtex',     label: 'GTEx v10' },
      { key: 'fetal_heart', label: 'Heart of Fetal Cells' },
      { key: 'medgen',   label: 'MedGen' },
      { key: 'mgi',      label: 'MGI' },
      { key: 'opentargets_evidence', label: 'Open Targets' },
      { key: 'panelapp', label: 'PanelApp' },
      { key: 'pmcoa',    label: 'PMC Open Access' },
      { key: 'protvar',  label: 'ProtVar' },
      { key: 'pubmed',   label: 'PubMed' },
      { key: 'pubtator3', label: 'PubTator3' },
      { key: 'spliceai', label: 'SpliceAI' },
      { key: 'uniprot',  label: 'UniProt' },
    ];

    const PHASES = [
      { key: 'query',    label: 'Querying databases' },
      { key: 'annotate', label: 'Annotating variant' },
      { key: 'generate', label: 'Generating interpretation' },
    ];

    function renderProgressSkeleton(gene, hgvs_c) {
      const stepper = PHASES.map((p, i) => `
        <div class="phase-step" data-phase="${p.key}" data-state="${i === 0 ? 'active' : 'future'}">
          <div class="phase-num"><span class="phase-num-inner">${i + 1}</span></div>
          <div class="phase-label">${p.label}</div>
        </div>
        ${i < PHASES.length - 1 ? '<div class="phase-connector" data-phase-connector="' + p.key + '"></div>' : ''}
      `).join('');

      const rows = DB_SOURCES.map(s =>
        `<div class="progress-row" data-src="${s.key}" data-state="pending">
           <span class="progress-icon"><span class="spinner-sm"></span></span>
           <span class="progress-label">${s.label}</span>
           <span class="progress-status">pending</span>
         </div>`
      ).join('');
      // Body content depends on the active loading-screen variant. The
      // phase-stepper is shared (per spec it stays above the terminal in
      // terminal mode); the per-source checklist / sim-grid / Claude
      // progress are 'original'-only so their DOM never lands when the
      // terminal owns the loading surface. The setDbRowState / fakeCritSim
      // / claudeProgress helpers all early-return on missing elements, so
      // they remain safe to call from the SSE handler regardless.
      const terminalBody = `
        <div class="term-window" id="term-window">
          <div class="term-titlebar">
            <span class="term-dot red"></span>
            <span class="term-dot amber"></span>
            <span class="term-dot green"></span>
            <span class="term-title">heartvar — variant pipeline</span>
          </div>
          <div class="term-body" id="term-body" aria-live="polite"></div>
        </div>`;

      const originalBody = `
        <div class="progress-rows" id="progress-rows">${rows}</div>
        <!-- Step-3 fake-criteria grid — hidden until Step 3 becomes
             active, replaces the DB checklist visually. -->
        <div class="sim-grid-wrap" id="sim-grid-wrap" style="display:none">
          <div class="sim-grid-label" id="sim-grid-label">
            <span class="skel-spinner" aria-hidden="true"></span>
            <span id="sim-grid-label-text">Getting AI interpretation…</span>
          </div>
          <div class="sim-grid" id="sim-grid"></div>
        </div>
        <div class="claude-progress" id="claude-progress" style="display:none">
          <div class="claude-spinner-row" aria-live="polite">
            <span class="spinner-md"></span>
            <span class="claude-spinner-label">Generating interpretation…</span>
            <span class="claude-elapsed" id="claude-elapsed">0%</span>
          </div>
          <div class="claude-progress-bar">
            <div class="claude-progress-fill" id="claude-progress-fill"></div>
          </div>
        </div>`;

      const body = LOADING_SCREEN === 'terminal' ? terminalBody : originalBody;

      // The slim red top bar lives inside the progress card now (not the
      // result card). It runs through Step 3 only and is hidden during
      // Steps 1-2 via [data-state="done"] on first render.
      // Escape user-typed gene/HGVS before it is interpolated into innerHTML.
      const _escHero = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
      const heroGene = _escHero((gene || '').trim() || 'Variant');
      const heroHgvs = _escHero((hgvs_c || '').trim());
      return `
        <div class="card v11-loading" id="progress-card">
          <div class="gen-progress" id="gen-progress" data-state="done" aria-hidden="true">
            <div class="gen-progress-fill" id="gen-progress-fill"></div>
          </div>
          <div class="v11-loading__head">
            <span class="v11-loading__eyebrow"><span class="v11-loading__pulse" aria-hidden="true"></span>Analysing variant</span>
            <div class="v11-loading__variant">
              <span class="v11-loading__gene">${heroGene}</span>${heroHgvs ? '<span class="v11-loading__hgvs">' + heroHgvs + '</span>' : ''}
            </div>
            <svg class="v11-loading__ecg" viewBox="0 0 600 60" fill="none" preserveAspectRatio="none" aria-hidden="true">
              <path d="M0 30 H170 l7 -20 6 38 6 -32 7 14 H300 l7 -20 6 38 6 -32 7 14 H600" stroke="#8c1a1f" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>
            </svg>
          </div>
          <div class="phase-stepper">${stepper}</div>
          ${body}
        </div>`;
    }

    function setPhaseState(phaseKey, state) {
      const step = document.querySelector(`.phase-step[data-phase="${phaseKey}"]`);
      if (!step) return;
      step.setAttribute('data-state', state);
      if (state === 'done') {
        const conn = document.querySelector(`[data-phase-connector="${phaseKey}"]`);
        if (conn) conn.setAttribute('data-state', 'done');
      }
    }

    // Track and advance phases from SSE events.
    function makePhaseTracker() {
      let current = 'query';
      return {
        advance(next) {
          // Mark all phases up to and including `current` as done, set `next` active.
          const order = PHASES.map(p => p.key);
          for (const k of order) {
            if (k === next) { setPhaseState(k, 'active'); break; }
            setPhaseState(k, 'done');
          }
          current = next;
        },
        finish() {
          for (const p of PHASES) setPhaseState(p.key, 'done');
          const lastConn = document.querySelector('[data-phase-connector]:last-of-type');
          if (lastConn) lastConn.setAttribute('data-state', 'done');
          // Also mark the final phase's preceding connector as done if it isn't already
          const conns = document.querySelectorAll('[data-phase-connector]');
          conns.forEach(c => c.setAttribute('data-state', 'done'));
        },
      };
    }

    function setDbRowState(src, state) {
      const row = document.querySelector(`.progress-row[data-src="${src}"]`);
      if (!row) return;
      // pending_access and not_applicable both share the muted 'skipped'
      // visual treatment via CSS, but each carries its own status text:
      // pending_access = awaiting subscription / API access; not_applicable
      // = the source genuinely cannot score this variant (e.g. SpliceAI
      // on indels) and "not configured" would be misleading.
      const cssState = (state === 'pending_access' || state === 'not_applicable') ? 'skipped' : state;
      row.setAttribute('data-state', cssState);
      const iconCell = row.querySelector('.progress-icon');
      const statusCell = row.querySelector('.progress-status');
      if (state === 'done')    { iconCell.textContent = '✓'; statusCell.textContent = 'done'; }
      else if (state === 'error') { iconCell.textContent = '✗'; statusCell.textContent = 'failed'; }
      else if (state === 'skipped') { iconCell.textContent = '·'; statusCell.textContent = 'not configured'; }
      else if (state === 'pending_access') { iconCell.textContent = '·'; statusCell.textContent = 'pending API access'; }
      else if (state === 'not_applicable') { iconCell.textContent = '·'; statusCell.textContent = 'not applicable'; }
      else                     { iconCell.innerHTML  = '<span class="spinner-sm"></span>'; statusCell.textContent = 'pending'; }
    }

    const CLAUDE_PROGRESS_ESTIMATE_MS = 60000;
    const claudeProgress = {
      timer: null,
      startTime: 0,
      start() {
        const box = document.getElementById('claude-progress');
        if (!box) return;
        box.style.display = 'block';
        this.startTime = performance.now();
        const fill = document.getElementById('claude-progress-fill');
        const elapsed = document.getElementById('claude-elapsed');
        if (fill) fill.style.width = '0%';
        if (elapsed) elapsed.textContent = '0%';
        clearInterval(this.timer);
        this.timer = setInterval(() => {
          const dt = performance.now() - this.startTime;
          const pct = Math.min(90, (dt / CLAUDE_PROGRESS_ESTIMATE_MS) * 90);
          if (fill) fill.style.width = pct.toFixed(1) + '%';
          if (elapsed) elapsed.textContent = Math.round(pct) + '%';
        }, 200);
      },
      finish() {
        clearInterval(this.timer);
        this.timer = null;
        const box = document.getElementById('claude-progress');
        const fill = document.getElementById('claude-progress-fill');
        const elapsed = document.getElementById('claude-elapsed');
        if (fill) fill.style.width = '100%';
        if (elapsed) elapsed.textContent = '100%';
        setTimeout(() => { if (box) box.style.display = 'none'; }, 400);
      },
    };

    // ── Stage 2 generation progress bar (slim, top of result card) ──
    // Replaces the in-tab crit-progress bar. Fills to 90 % over a ~30 s
    // estimate, animates to 100 % on stage2_complete, then fades out so
    // the criteria reveal reads as a single motion.
    const GEN_PROGRESS_ESTIMATE_MS = 50000;
    const genProgress = {
      timer: null,
      startTime: 0,
      start() {
        const box = document.getElementById('gen-progress');
        const fill = document.getElementById('gen-progress-fill');
        if (!box || !fill) return;
        // Idempotent: if the bar is already animating, leave the existing
        // run alone so stage1→stage2 doesn't reset to 0 %.
        if (this.timer) return;
        box.removeAttribute('data-state');
        box.style.opacity = '1';
        fill.style.width = '0%';
        this.startTime = performance.now();
        this.timer = setInterval(() => {
          const dt = performance.now() - this.startTime;
          const pct = Math.min(90, (dt / GEN_PROGRESS_ESTIMATE_MS) * 90);
          fill.style.width = pct.toFixed(1) + '%';
        }, 200);
      },
      finish() {
        clearInterval(this.timer);
        this.timer = null;
        const fill = document.getElementById('gen-progress-fill');
        const box = document.getElementById('gen-progress');
        if (fill) fill.style.width = '100%';
        // Hold at 100 % for a beat, then fade out so the bar's completion
        // and the criteria reveal both register as a single transition.
        setTimeout(() => { if (box) box.setAttribute('data-state', 'done'); }, 300);
      },
    };
    // Back-compat alias — older code paths call critProgress.start/finish.
    const critProgress = genProgress;

    // ── Retro terminal loading screen (LOADING_SCREEN === 'terminal') ──
    // Self-paced reveal of a pre-scripted pipeline log. Real values from
    // SSE events are patched into {placeholders} as they land; lines
    // revealed before a value arrives display a muted "…" until the value
    // is available. On terminalAnim.finish() the script flushes any
    // remaining lines quickly and emits the green "done ✓" prompt.
    const TERMINAL_LINES = [
      // ── Stage 1: query (~20 s)
      {kind:'prompt', text:'heartvar$ query --variant {GENE}:{HGVS} --genome GRCh38'},
      {kind:'out',    text:'Resolving variant notation...'},
      {kind:'pad',    text:'Connecting to Ensembl VEP REST API...',  tail:'OK', tone:'ok'},
      {kind:'out',    text:'Consequence: {vep_consequence} | IMPACT: {impact}'},
      {kind:'out',    text:'Transcript: {transcript}'},
      {kind:'out',    text:'Protein change: {p_notation}'},
      {kind:'out',    text:'SpliceAI delta scores: {spliceai_result}'},
      {kind:'out',    text:'AlphaMissense score: {alphamissense_result}'},
      {kind:'pad',    text:'Connecting to gnomAD v4...',              tail:'OK', tone:'ok'},
      {kind:'out',    text:'Allele frequency: {af_result}'},
      {kind:'out',    text:'Popmax frequency: {popmax}'},
      {kind:'pad',    text:'Querying ClinVar local database...',      tail:'OK', tone:'ok'},
      {kind:'out',    text:'Submissions found: {clinvar_count}'},
      {kind:'out',    text:'Most recent classification: {clinvar_class}'},
      {kind:'pad',    text:'Querying UniProt {accession}...',         tail:'OK', tone:'ok'},
      {kind:'out',    text:'Protein: {protein_name} · {length} aa'},
      {kind:'out',    text:'Domain context: {domain_result}'},
      {kind:'pad',    text:'Fetching PubMed gene papers...',          tail:'OK', tone:'ok'},
      // ── Stage 2: annotate (~15 s)
      {kind:'spacer'},
      {kind:'prompt', text:'heartvar$ annotate --criteria ACMG-AMP-2015'},
      {kind:'out',    text:'Loading ACMG/AMP 2015 classification rules...'},
      {kind:'out',    text:'Evaluating loss-of-function criteria (PVS1)...'},
      {kind:'out',    text:'Checking population frequency thresholds (BA1, PM2)...'},
      {kind:'out',    text:'Evaluating functional evidence criteria (PS3, BS3)...'},
      {kind:'out',    text:'Checking de novo status (PS2, PM6)...'},
      {kind:'out',    text:'Cross-referencing computational predictors (PP3, BP4)...'},
      {kind:'out',    text:'Evaluating phenotype specificity (PP4)...'},
      {kind:'out',    text:'Checking reputable source evidence (PP5, BP6)...'},
      {kind:'out',    text:'Aggregating evidence weights...'},
      {kind:'pad',    text:'Preliminary score: {score} pts',          tail:'OK', tone:'warn'},
      // ── Stage 3: interpret (~25 s) — 5 progress-bar lines collapse
      //    into one row via replaceLast so the bar feels like it animates
      //    in place rather than scrolling 5 rows of [###].
      {kind:'spacer'},
      {kind:'prompt', text:'heartvar$ interpret --model claude-sonnet-4 --output structured'},
      {kind:'out',    text:'Preparing evidence summary for AI review...'},
      {kind:'pad',    text:'Sending {token_count} tokens to Claude...', tail:'OK', tone:'ok'},
      {kind:'out',    text:'Streaming response...'},
      {kind:'bar',    text:'[████░░░░░░░░░░░░░░░░] 20%'},
      {kind:'bar',    text:'[████████░░░░░░░░░░░░] 40%',   replaceLast:true},
      {kind:'bar',    text:'[████████████░░░░░░░░] 60%',   replaceLast:true},
      {kind:'bar',    text:'[████████████████░░░░] 80%',   replaceLast:true},
      {kind:'bar',    text:'[████████████████████] 100%',  replaceLast:true},
      {kind:'pad',    text:'Parsing structured output...',          tail:'OK', tone:'ok'},
      {kind:'pad',    text:'Validating ACMG criteria assignments...', tail:'OK', tone:'ok'},
      {kind:'out',    text:'Classification confirmed: {final_class}'},
      {kind:'pad',    text:'Generating clinical narrative...',      tail:'OK', tone:'ok'},
      {kind:'done',   text:'heartvar$ done ✓'},
    ];

    const terminalAnim = {
      values: {},
      timer: null,
      cursor: 0,
      lastLineEl: null,
      cursorEl: null,
      flushing: false,
      escHtml(s) {
        return String(s == null ? '' : s)
          .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      },
      // Resolve {placeholder} tokens against the running values map. Tokens
      // that have not yet been populated render as a muted "…" so the
      // script structure is preserved without misleading the curator.
      interpolate(text) {
        return text.replace(/\{([A-Z_]+)\}/g, (_, k) => {
          const v = this.values[k];
          if (v == null || v === '') return '<span class="term-pending">…</span>';
          return this.escHtml(v);
        });
      },
      mountCursor(target) {
        if (!this.cursorEl) {
          this.cursorEl = document.createElement('span');
          this.cursorEl.className = 'term-cursor';
        }
        if (target) target.appendChild(this.cursorEl);
      },
      buildLine(line) {
        const el = document.createElement('div');
        el.className = 'term-line';
        if (line.kind === 'spacer') {
          el.innerHTML = '&nbsp;';
          return el;
        }
        if (line.kind === 'prompt') {
          // Split into prompt + command at the first space after "$".
          const m = line.text.match(/^(\S+\$)\s+(.*)$/);
          const prompt = m ? m[1] : 'heartvar$';
          const cmd    = m ? m[2] : line.text;
          el.innerHTML =
            `<span class="term-prompt">${this.escHtml(prompt)}</span> ` +
            `<span class="term-cmd">${this.interpolate(cmd)}</span>`;
          return el;
        }
        if (line.kind === 'done') {
          el.innerHTML = `<span class="term-done">${this.escHtml(line.text)}</span>`;
          return el;
        }
        if (line.kind === 'bar') {
          el.innerHTML = `<span class="term-bar">${this.escHtml(line.text)}</span>`;
          return el;
        }
        if (line.kind === 'pad') {
          el.classList.add('pad');
          const toneClass = line.tone === 'ok'   ? 'term-ok'
                          : line.tone === 'warn' ? 'term-warn'
                          : 'term-out';
          el.innerHTML =
            `<span class="term-out">${this.interpolate(line.text)}</span>` +
            `<span class="term-tail ${toneClass}">${this.escHtml(line.tail || '')}</span>`;
          return el;
        }
        // Default 'out' line.
        const cls = line.tone === 'ok'   ? 'term-ok'
                  : line.tone === 'warn' ? 'term-warn'
                  : 'term-out';
        el.innerHTML = `<span class="${cls}">${this.interpolate(line.text)}</span>`;
        return el;
      },
      revealNext() {
        const body = document.getElementById('term-body');
        if (!body) return;
        if (this.cursor >= TERMINAL_LINES.length) return;
        const line = TERMINAL_LINES[this.cursor++];
        const el   = this.buildLine(line);
        if (line.replaceLast && this.lastLineEl && this.lastLineEl.parentNode === body) {
          body.replaceChild(el, this.lastLineEl);
        } else {
          body.appendChild(el);
        }
        this.lastLineEl = el;
        // Move the cursor block to sit at the end of the latest line; on
        // a 'done' line we drop the cursor entirely.
        if (this.cursorEl && this.cursorEl.parentNode) {
          this.cursorEl.parentNode.removeChild(this.cursorEl);
        }
        if (line.kind !== 'done' && line.kind !== 'spacer') {
          this.mountCursor(el);
        }
        // Auto-advance to the latest line. With overflow-y: hidden the
        // user can't scroll, so this just clamps to the bottom view.
        body.scrollTop = body.scrollHeight;
      },
      scheduleNext() {
        if (this.cursor >= TERMINAL_LINES.length) return;
        // Progress-bar rows feel right with ~2 s between them; everything
        // else uses a 600-1800 ms jittered delay for organic pacing.
        const next = TERMINAL_LINES[this.cursor];
        const isBar = next && next.kind === 'bar';
        const delay = isBar
          ? 1700 + Math.random() * 600
          : 600  + Math.random() * 1200;
        this.timer = setTimeout(() => {
          this.revealNext();
          this.scheduleNext();
        }, delay);
      },
      start({gene, hgvs_c}) {
        if (LOADING_SCREEN !== 'terminal') return;
        const body = document.getElementById('term-body');
        if (!body) return;
        this.values = {GENE: gene || '?', HGVS: hgvs_c || '?', token_count: '~2,400'};
        this.cursor = 0;
        this.lastLineEl = null;
        this.flushing = false;
        body.innerHTML = '';
        clearTimeout(this.timer);
        // First line lands immediately so the terminal isn't blank for a
        // beat at the start of the run.
        this.revealNext();
        this.scheduleNext();
      },
      // Patch one or more {PLACEHOLDER} values; any already-revealed line
      // that references them is re-rendered in place.
      setValues(patch) {
        if (LOADING_SCREEN !== 'terminal') return;
        Object.assign(this.values, patch);
        this.rerenderAll();
      },
      rerenderAll() {
        const body = document.getElementById('term-body');
        if (!body) return;
        // Walk every revealed line. We can't simply replace innerHTML
        // because the cursor element is appended to whichever line is
        // current; capture which line owns it, rebuild content, then
        // re-mount the cursor on the same line index.
        const cursorIdx = this.cursorEl && this.lastLineEl
          ? Array.prototype.indexOf.call(body.children, this.lastLineEl)
          : -1;
        const revealed = TERMINAL_LINES.slice(0, this.cursor);
        // Rebuild children one-for-one. If counts diverge (shouldn't, but
        // defend against it), fall back to a full rebuild.
        if (body.children.length !== revealed.length) {
          body.innerHTML = '';
          revealed.forEach(l => body.appendChild(this.buildLine(l)));
        } else {
          revealed.forEach((line, i) => {
            const fresh = this.buildLine(line);
            body.replaceChild(fresh, body.children[i]);
          });
        }
        this.lastLineEl = body.children[Math.max(0, cursorIdx)] || body.lastElementChild;
        if (this.cursorEl && this.cursorEl.parentNode) {
          this.cursorEl.parentNode.removeChild(this.cursorEl);
        }
        if (this.lastLineEl) this.mountCursor(this.lastLineEl);
      },
      // Pull values out of the per-source data payload arriving via the
      // db_done SSE event. Keeps the runCuration handler free of source-
      // specific destructuring; everything funnels through here.
      ingestSource(source, data) {
        if (LOADING_SCREEN !== 'terminal' || !data) return;
        const patch = {};
        if (source === 'vep' && data.ok) {
          if (data.most_severe_consequence) patch.vep_consequence = String(data.most_severe_consequence).replace(/_/g, ' ');
          if (data.impact)                  patch.impact = data.impact;
          if (data.transcript_id)           patch.transcript = data.transcript_id;
          if (data.hgvsp) {
            const i = String(data.hgvsp).indexOf(':');
            patch.p_notation = i >= 0 ? data.hgvsp.slice(i + 1) : data.hgvsp;
          } else {
            patch.p_notation = 'none (non-coding)';
          }
        } else if (source === 'spliceai') {
          if (data.ok && (data.delta_acceptor != null || data.delta_donor != null)) {
            const ds = ['delta_acceptor', 'delta_donor', 'delta_acceptor_gain', 'delta_donor_gain']
              .map(k => data[k]).filter(v => v != null);
            const maxDelta = ds.length ? Math.max(...ds.map(Number)) : null;
            patch.spliceai_result = maxDelta != null ? `Δ max = ${maxDelta.toFixed(2)}` : 'no impact';
          } else {
            patch.spliceai_result = 'not available';
          }
        } else if (source === 'alphamissense') {
          if (data.ok && data.am_pathogenicity != null) {
            patch.alphamissense_result = `${Number(data.am_pathogenicity).toFixed(3)} (${data.am_class || 'classified'})`;
          } else {
            patch.alphamissense_result = 'not available';
          }
        } else if (source === 'gnomad') {
          if (data.ok && data.variant_found) {
            patch.af_result = data.af != null ? Number(data.af).toExponential(2) : 'not reported';
            patch.popmax    = data.af_popmax != null ? Number(data.af_popmax).toExponential(2) : 'n/a';
          } else {
            patch.af_result = 'not found';
            patch.popmax    = 'n/a';
          }
        } else if (source === 'clinvar') {
          if (data.ok) {
            const recs = data.records || data.entries || [];
            patch.clinvar_count = String(recs.length);
            const first = recs[0] || {};
            patch.clinvar_class = first.classification || first.clinical_significance || 'no submissions';
          }
        } else if (source === 'uniprot') {
          if (data.ok && data.found) {
            patch.accession    = data.accession || '';
            patch.protein_name = data.protein_name || '';
            patch.length       = data.length != null ? String(data.length) : '?';
            const dn = (data.domains || []).length;
            patch.domain_result = dn ? `${dn} annotated domain${dn === 1 ? '' : 's'}` : 'no domains annotated';
          } else {
            patch.accession = 'n/a';
            patch.domain_result = 'no entry';
          }
        }
        if (Object.keys(patch).length) this.setValues(patch);
      },
      ingestStage1(stage1) {
        if (LOADING_SCREEN !== 'terminal' || !stage1) return;
        const patch = {};
        if (Number.isFinite(stage1.points_total)) {
          patch.score = (stage1.points_total >= 0 ? '+' : '−') + Math.abs(stage1.points_total);
        }
        if (stage1.classification) patch.final_class = stage1.classification;
        if (Object.keys(patch).length) this.setValues(patch);
      },
      // Called from stage2_complete / db_only_complete. Flushes any
      // remaining scripted lines at a fast pace so the terminal feels
      // synced to the real pipeline finish, then resolves on a 500 ms
      // delay so the caller can transition to the results panel.
      finish() {
        if (LOADING_SCREEN !== 'terminal') return Promise.resolve();
        if (this.flushing) return this.flushPromise || Promise.resolve();
        this.flushing = true;
        clearTimeout(this.timer);
        this.flushPromise = new Promise(resolve => {
          const tick = () => {
            if (this.cursor >= TERMINAL_LINES.length) {
              // Drop the cursor on the final "done ✓" line.
              if (this.cursorEl && this.cursorEl.parentNode) {
                this.cursorEl.parentNode.removeChild(this.cursorEl);
              }
              setTimeout(resolve, 500);
              return;
            }
            this.revealNext();
            setTimeout(tick, 110);
          };
          tick();
        });
        return this.flushPromise;
      },
    };

    // Step-3 fake-criteria simulator: 27 decorative boxes that complete
    // one at a time (~1.1 s each) so the user has something to watch
    // during Claude generation. Boxes hold no real meaning — the actual
    // ACMG results stream in via stage2_complete and are revealed by
    // revealResults() once Stage 2 finishes.
    // Pace the sim so it spans roughly the full Claude-generation wait
    // (typically ~45-55 s end-to-end). The last few boxes are held back
    // for rapidComplete() so something always animates on reveal, even
    // when the result lands ahead of the 27-tick schedule.
    //
    // The 28 entries below mirror the ACMG/AMP criteria the model
    // evaluates in the consolidated curation call. Boxes render code +
    // short label so the step-3 grid reads like a checklist of criteria
    // being ticked off rather than abstract placeholder bars.
    const SIM_CRITERIA = [
      { code: 'PVS1', name: 'Null variant' },
      { code: 'PS1',  name: 'Same AA change' },
      { code: 'PS2',  name: 'De novo (confirmed)' },
      { code: 'PS3',  name: 'Functional studies' },
      { code: 'PS4',  name: 'Affected prevalence' },
      { code: 'PM1',  name: 'Critical domain' },
      { code: 'PM2',  name: 'Absent in controls' },
      { code: 'PM3',  name: 'In trans w/ path' },
      { code: 'PM4',  name: 'Length change' },
      { code: 'PM5',  name: 'Novel missense' },
      { code: 'PM6',  name: 'De novo (assumed)' },
      { code: 'PP1',  name: 'Co-segregation' },
      { code: 'PP2',  name: 'Low-tolerance missense' },
      { code: 'PP3',  name: 'Computational' },
      { code: 'PP4',  name: 'Phenotype match' },
      { code: 'PP5',  name: 'Reputable source' },
      { code: 'BA1',  name: 'AF > 5%' },
      { code: 'BS1',  name: 'AF above expected' },
      { code: 'BS2',  name: 'Healthy adult' },
      { code: 'BS3',  name: 'Functional (benign)' },
      { code: 'BS4',  name: 'No segregation' },
      { code: 'BP1',  name: 'Missense in trunc gene' },
      { code: 'BP2',  name: 'In trans / cis' },
      { code: 'BP3',  name: 'In-frame repeat' },
      { code: 'BP4',  name: 'Computational (benign)' },
      { code: 'BP5',  name: 'Alt molecular cause' },
      { code: 'BP6',  name: 'Reputable (benign)' },
      { code: 'BP7',  name: 'Silent / no splice' },
    ];
    const SIM_BOX_COUNT = SIM_CRITERIA.length;
    const SIM_TICK_MS   = 1500;   // 27 × 1.5 s ≈ 40 s if every tick fires
    const SIM_RESERVE   = 4;      // last 4 boxes held for the reveal
    const SIM_FAST_MS   = 70;     // UPPER bound on the per-box finale pace
    // TOTAL time the rapid-complete finale may take, however many boxes are
    // left. Previously the finale was 70 ms x remaining, so it grew as the
    // backend got faster — a 3 s gather left ~25 boxes and cost ~1.75 s here.
    const SIM_FINALE_BUDGET_MS = 320;
    // Hold on "Analysis complete ✓" before revealing the card. Long enough to
    // register, short enough not to be the thing the curator is waiting for.
    const SIM_DONE_HOLD_MS = 180;
    const fakeCritSim = {
      timer: null,
      remaining: [],
      start() {
        const grid = document.getElementById('sim-grid');
        const wrap = document.getElementById('sim-grid-wrap');
        const rows = document.getElementById('progress-rows');
        if (!grid || !wrap) return;
        // Build boxes from the SIM_CRITERIA list so each cell carries the
        // ACMG code + short name. The "done" state is signalled by a
        // warm-sand background swap on the box itself (no tick / icon —
        // see the CSS comment by .sim-box[data-state="done"] for why
        // we deliberately avoid checkmarks here).
        grid.innerHTML = SIM_CRITERIA.map((c, i) =>
          `<div class="sim-box" data-state="pending" data-idx="${i}">
             <div class="sim-box-code">${c.code}</div>
             <div class="sim-box-name" title="${c.name}">${c.name}</div>
           </div>`
        ).join('');
        this.remaining = Array.from({ length: SIM_BOX_COUNT }, (_, i) => i);
        // Shuffle so boxes complete in random order rather than left-to-right.
        for (let i = this.remaining.length - 1; i > 0; i--) {
          const j = Math.floor(Math.random() * (i + 1));
          [this.remaining[i], this.remaining[j]] = [this.remaining[j], this.remaining[i]];
        }
        // Swap the DB checklist for the simulation grid.
        if (rows) rows.style.display = 'none';
        wrap.style.display = 'block';
        const labelText = document.getElementById('sim-grid-label-text');
        const label = document.getElementById('sim-grid-label');
        if (label) label.classList.remove('sim-done');
        if (labelText) labelText.textContent = 'Getting AI interpretation…';
        clearInterval(this.timer);
        this.timer = setInterval(() => this.tickOne(), SIM_TICK_MS);
        // Slim top progress bar runs alongside the simulation.
        genProgress.start();
      },
      tickOne() {
        // Stop short of completing every box — the reserve ensures the
        // reveal still shows a few boxes ticking even when generation
        // takes longer than the simulator's 27-tick schedule.
        if (this.remaining.length <= SIM_RESERVE) {
          clearInterval(this.timer);
          this.timer = null;
          return;
        }
        const idx = this.remaining.shift();
        const box = document.querySelector(`.sim-box[data-idx="${idx}"]`);
        if (box) box.setAttribute('data-state', 'done');
      },
      // Rapid-complete any remaining boxes, flip the label to "Analysis
      // complete ✓", then call `cb` after a short hold.
      //
      // ⚠ THE FINALE IS BUDGETED, NOT PER-BOX, and the reason is that the old
      // form got SLOWER as the backend got FASTER. It ticked every remaining
      // box at a fixed 70 ms and then held 500 ms. During loading the sim ticks
      // one box per SIM_TICK_MS (1500 ms) and stops with SIM_RESERVE left, so a
      // gather that finishes in 3 s leaves ~25 boxes untouched:
      //     25 x 70 ms + 500 ms = ~2.3 s of animation AFTER the data is ready.
      // Reported from production as "the loading screen still shows for a few
      // seconds after all databases are ticked green" — and it was, because
      // every second cut from the gather added ~0.5 s of leftover boxes here.
      //
      // A fixed TOTAL budget decouples the two: the finale always takes about
      // the same short time whether there are 4 boxes left or 25. The reserve
      // still exists so something visibly animates on reveal (see
      // SIM_RESERVE), it just no longer bills the curator for the backend
      // being quick.
      rapidComplete(cb) {
        clearInterval(this.timer);
        this.timer = null;
        const finishUp = () => {
          const label = document.getElementById('sim-grid-label');
          const labelText = document.getElementById('sim-grid-label-text');
          if (label) {
            label.classList.add('sim-done');
            // Replace the spinner with a green tick when complete.
            const spinner = label.querySelector('.skel-spinner');
            if (spinner) {
              const tick = document.createElement('span');
              tick.innerHTML = `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" class="sim-grid-tick"><polyline points="20 6 9 17 4 12"/></svg>`;
              spinner.replaceWith(tick.firstChild);
            }
          }
          if (labelText) labelText.textContent = 'Analysis complete';
          setTimeout(cb || (() => {}), SIM_DONE_HOLD_MS);
        };
        if (!this.remaining.length) {
          finishUp();
          return;
        }
        // Spread whatever is left across the budget rather than paying a
        // fixed cost per box. Floored at 8 ms so a large leftover still reads
        // as a sweep rather than an instant jump.
        const perBox = Math.max(
          8, Math.min(SIM_FAST_MS, SIM_FINALE_BUDGET_MS / this.remaining.length));
        const tickFast = () => {
          if (!this.remaining.length) {
            finishUp();
            return;
          }
          const idx = this.remaining.shift();
          const box = document.querySelector(`.sim-box[data-idx="${idx}"]`);
          if (box) box.setAttribute('data-state', 'done');
          setTimeout(tickFast, perBox);
        };
        tickFast();
      },
      reset() {
        clearInterval(this.timer);
        this.timer = null;
        this.remaining = [];
      },
    };

    // Common gate between "stage 2 has resolved" and "reveal the result
    // card". In terminal mode this waits for the scripted log to flush
    // its tail + emit the green "done ✓" line; in original mode it falls
    // back to fakeCritSim.rapidComplete so existing callers don't change.
    function chainLoadingFinish(cb) {
      if (LOADING_SCREEN === 'terminal') {
        terminalAnim.finish().then(cb || (() => {}));
      } else {
        fakeCritSim.rapidComplete(cb);
      }
    }

    // SSE parser: takes a ReadableStream<Uint8Array>, yields {event, data} objects.
    async function* parseSse(stream) {
      const reader = stream.getReader();
      const decoder = new TextDecoder();
      let buf = '';
      while (true) {
        const { value, done } = await reader.read();
        if (done) break;
        buf += decoder.decode(value, { stream: true });
        let idx;
        while ((idx = buf.indexOf('\n\n')) !== -1) {
          const block = buf.slice(0, idx);
          buf = buf.slice(idx + 2);
          let event = 'message', data = '';
          for (const line of block.split('\n')) {
            if      (line.startsWith('event: ')) event = line.slice(7).trim();
            else if (line.startsWith('data: '))  data += line.slice(6);
          }
          if (data) {
            try { yield { event, data: JSON.parse(data) }; }
            catch (e) { console.error('Bad SSE data:', e, data.slice(0, 200)); }
          }
        }
      }
    }

    // ── Client-side LLM helpers (browser-direct AI mode) ──
    // Strip leading code fences and extract the first complete JSON object
    // from a response string. Mirrors backend/claude.py parse_first_json_object.
    //
    // Unused: the browser-direct AI path this served was removed, and the
    // server parses its own JSON. Retained deliberately — see the release
    // harness notes before removing it.
    function parseFirstJsonObject(text) {
      const stripped = String(text).replace(/```(?:json)?/gi, '').trim();
      const start = stripped.indexOf('{');
      if (start < 0) return JSON.parse(stripped);
      // Walk braces accounting for string literals so we stop at the matching
      // close brace and ignore any trailing prose Claude/GPT may add.
      let depth = 0, inStr = false, esc = false, end = -1;
      for (let i = start; i < stripped.length; i++) {
        const ch = stripped[i];
        if (esc) { esc = false; continue; }
        if (ch === '\\') { esc = true; continue; }
        if (ch === '"') { inStr = !inStr; continue; }
        if (inStr) continue;
        if (ch === '{') depth++;
        else if (ch === '}') {
          depth--;
          if (depth === 0) { end = i; break; }
        }
      }
      if (end < 0) throw new SyntaxError('Unterminated JSON object in LLM response');
      return JSON.parse(stripped.slice(start, end + 1));
    }

    // _CRITERION_NAMES (code → short human-readable name) is populated from
    // static/acmg_constants.json by ACMG_READY. Used by normalizeCriteria when
    // the model emits the compact 4-tuple form (no `name`/`direction` fields).

    // Accept both legacy dict-per-criterion AND the compact 4-tuple form
    // [code, status, criteria_strength, evidence] that the prompt now
    // asks Claude to emit. Returns the full dict shape the rest of the
    // renderer expects (code/name/status/direction/criteria_strength/
    // evidence). Mirrors backend/app.py:_normalize_criteria so both the
    // server flow and the client-AI flow normalise identically.
    function normalizeCriteria(items) {
      const out = [];
      for (const it of items || []) {
        if (it && typeof it === 'object' && !Array.isArray(it)) {
          const code = (it.code || '').trim();
          const direction = it.direction || (code.startsWith('B') ? 'benign' : 'pathogenic');
          out.push({
            code,
            name: it.name || _CRITERION_NAMES[code] || code,
            status: it.status,
            direction,
            criteria_strength: it.criteria_strength || null,
            evidence: it.evidence || '',
            // Stored for downstream auditability — "hard_coded" vs "ai".
            // The renderer does not surface this field; it travels with
            // the criterion so server-side auditing / debugging tools
            // can tell which evaluator produced each entry.
            source: it.source || null,
          });
        } else if (Array.isArray(it)) {
          const code = (it[0] || '').toString().trim();
          out.push({
            code,
            name: _CRITERION_NAMES[code] || code,
            status: it[1] || null,
            direction: code.startsWith('B') ? 'benign' : 'pathogenic',
            criteria_strength: it[2] || null,
            evidence: it[3] || '',
            source: null,
          });
        }
      }
      return out;
    }

    // NO browser-side scorer/classifier lives here. The server is the single
    // source of the points total and the tier: backend/acmg/tiers.py owns
    // compute_points_total + classification_for_criteria (including the
    // ACMG-2015 benign combining floor), and the page renders the
    // `points_total` / `classification` it is handed. The old JS twins
    // (_pointsFor / computePointsTotal / classificationFor) were remnants of
    // the removed BYOK client-AI path, unreferenced since the server took the
    // key, and were deleted rather than kept in cross-language sync — a second
    // classifier that nothing calls can only drift. If a browser-direct AI path
    // ever returns, call the server for the tier instead of re-adding these.
      // Three ACMG helpers that mirrored app.py in the browser were removed on
      // 2026-09-08: mergeHardCodedAndAi, applyCrossCriterionExclusions and
      // gateCriteriaApplicability, with their _PM1_LOF_CONSEQUENCES table. They
      // existed for a bring-your-own-key path that ran the LLM call from the page
      // and so needed to combine criteria client-side. That path is gone, the
      // server owns the combine and the tier, and each function had exactly one
      // occurrence in the tree: its own definition. If a client-side flow ever
      // returns, call the server for the tier rather than re-adding these.

    // Did this response come from the PLATFORM rather than from HeartVar?
    //
    // Two independent signals, because either alone lets a case through:
    //   * the status codes App Service / a front door return when the container
    //     is not running (403 is Azure's "this web app is stopped"), and
    //   * an HTML document body on an endpoint that only ever returns JSON or
    //     an SSE stream, which catches any future platform page whatever its
    //     status.
    // A JSON body means the response came from US even on those statuses — the
    // app raises its own 403 with {detail:{error,message}} (backend/app.py's
    // admin guard), and that message must still reach the user rather than being
    // replaced by a maintenance notice.
    function _looksLikePlatformErrorPage(status, body) {
      const text = (body || '').trimStart();
      const head = text.slice(0, 200).toLowerCase();
      if (head.startsWith('<!doctype html') || head.startsWith('<html')) return true;
      if (![403, 502, 503, 504].includes(status)) return false;
      try {
        JSON.parse(text);
        return false;          // structured, so it is our own refusal
      } catch (e) {
        return true;           // no body, or a platform page we do not recognise
      }
    }

    async function runCuration() {
      // Guarantee the shared ACMG constants are loaded before any scoring or
      // rendering runs. In practice this fetch resolves long before the user
      // submits (and well before the API round-trip below), so this never
      // blocks — it just removes any theoretical race on a cold first submit.
      await ACMG_READY;
      const gene = document.getElementById('hvl-gene').value.trim();
      const hgvs_c = document.getElementById('hvl-variant').value.trim();
      // Coordinate input may omit the gene — the backend derives it
      // from VEP. HGVS input still requires both fields.
      const coordMode = isCoordInput(hgvs_c);
      if (!hgvs_c) { alert('Please enter a variant.'); return; }
      if (!gene && !coordMode) { alert('Please enter a gene symbol (required for HGVS input).'); return; }
      // Catch input formats the annotator can't resolve BEFORE running the
      // pipeline, so curators get a targeted hint instead of a downstream VEP
      // failure. (Resolving rsIDs / protein-only HGVS would need a backend
      // lookup; for now we guide the curator to a c. or genomic input.)
      if (/^\s*rs\d+\s*$/i.test(hgvs_c)) {
        alert('rsIDs aren’t supported yet. Please enter the variant as a coding HGVS (e.g. c.1988G>A) or genomic coordinates (e.g. 14-23424115-G-A).');
        return;
      }
      if (/^\s*p\.\(?[A-Za-z*]/i.test(hgvs_c)) {
        alert('Protein-only notation (e.g. p.Arg502Trp) isn’t supported yet — the variant can’t be annotated from a protein change alone. Please enter the coding HGVS (e.g. c.1504C>T) or genomic coordinates (e.g. 14-23424115-G-A).');
        return;
      }
      // Results present full-width.
      document.body.classList.add('hv-results-full');
      window.__hvRan = true;
      // Reconcile the build toggle with the current input before reading it,
      // so a build selected for an earlier coordinate entry can never ride
      // onto an HGVS submission (the toggle hides itself for non-coord input).
      hvlSyncBuildToggle();
      const genome_build = getSelectedGenomeBuild();

      const hpo = document.getElementById('hvl-hpo').value.trim();
      const family = document.getElementById('hvl-famhx').value.trim();
      // Structured clinical-context fields. All five are sent to the
      // backend on every submission — empty string for "unknown" so the
      // backend gets explicit unknowns rather than missing keys.
      const zygosity = document.getElementById('hvl-zygosity').value;
      const inheritance_input = document.getElementById('hvl-inheritance').value;
      const proband_sex = document.getElementById('hvl-sex').value;
      const trio_status = document.getElementById('hvl-trio').value;
      // 🔴 THE DE-NOVO VALUE IS NO LONGER DISCARDED WHEN TRIO IS UNKNOWN.
      // This read `trio_status ? dropdown.value : ''`, so selecting "Confirmed
      // de novo" and leaving Trio status alone threw the answer away IN THE
      // BROWSER — the server never saw it, PS2 could not fire, and the curator
      // was told "De novo not confirmed by parental testing" about an input
      // they had actually supplied. Reported from production 2026-09-07.
      //
      // The stale-value concern it was written for is now handled properly at
      // both ends: the dropdown has an explicit "Not assessed" default that
      // asserts nothing, and the backend decides what each combination earns
      // (app.py denovo_confirmed / hard_coded.py _eval_pm6) — "Confirmed" is
      // trusted on an unknown trio, "Unconfirmed" is not, because SVI awards no
      // points when the parents may never have been tested.
      const denovo_status = document.getElementById('hvl-denovo').value;

      // Structured family & segregation evidence. Counts are parsed to
      // non-negative integers (blank / NaN → 0 = "not provided"); the
      // selects pass through their raw "" / "yes" / "no" values. Every
      // key matches a CurationRequest field 1:1 and is sent on EVERY
      // submission so curator-entered family data always reaches the
      // deterministic engine + AI prompt (PP1/BS4/PM3/PS2/PM6/BS2/BP5).
      const _segInt = (id) => {
        const n = parseInt(document.getElementById(id).value, 10);
        return Number.isFinite(n) && n > 0 ? n : 0;
      };
      const seg_affected_carriers    = _segInt('hvl-seg-affected-carriers');
      const in_trans_pathogenic      = document.getElementById('hvl-in-trans').value;
      const denovo_confirmed_count   = _segInt('hvl-denovo-confirmed-count');
      const denovo_unconfirmed_count = _segInt('hvl-denovo-unconfirmed-count');
      // seg_affected_noncarriers is READ FROM THE FORM again as of 2026-09-07
      // (see the note on its input in index.html). It was BS4's primary input
      // and hardcoding it to 0 meant that path could never fire while the
      // carriers box above fed PP1 normally.
      //
      // As of 2026-09-08 these two counts are the ONLY inputs to BS4 and PP1.
      // Both criteria are derived in Python now and neither is sent to the
      // model, so the free-text route through the LLM is gone. A blank box is
      // silently worth zero evidence, which is why index.html says so under
      // each one.
      //
      // seg_meioses and the alternate-cause fields genuinely have no input:
      // their markup was removed and they stay at their literal defaults to
      // keep the request shape unchanged (the backend accepts every key).
      // seg_meioses gates PP1's strength ladder, so PP1 cannot reach Strong
      // from the form; BS2 / BP5 are still inferred from the free-text Family
      // history box rather than driven by dedicated inputs.
      const seg_affected_noncarriers = _segInt('hvl-seg-affected-noncarriers');
      const seg_meioses              = 0;
      const alt_cause_present        = '';
      const alt_cause_detail         = '';

      // Resolve the AI routing. AI runs on HeartVar's OWN key (server-side)
      // when the curator ticked "Include AI interpretation"; otherwise the
      // backend returns the deterministic evidence-only result. Reset
      // window.__hvAiUsed here — it is set true only when a real AI result
      // (stage2_complete) lands, so the results-only chatbot stays hidden when
      // AI is off OR when AI was requested but degraded to evidence-only.
      const aiMode = getAiEnabled() ? 'server' : 'none';
      window.__hvAiUsed = false;

      // Entry is the landing's Generate button, which already flipped the view
      // before runCuration ran. Guard so the spinner toggle no-ops if absent.
      const btn = document.getElementById('submitBtn');
      if (btn) {
        btn.disabled = true;
        btn.innerHTML = '<span class="spinner"></span>Interpreting variant...';
      }

      const rp = document.getElementById('result-panel');
      // Swap the right-column view from the empty placeholder to the live
      // progress skeleton. No scrollIntoView — the right column is already
      // visible in the split layout.
      const ph = document.getElementById('result-placeholder');
      if (ph) ph.style.display = 'none';
      rp.style.display = 'block';
      rp.innerHTML = renderProgressSkeleton(gene, hgvs_c);
      // Chatbot is results-only: clear any leftover .show/.open from a prior
      // run so "Ask about this variant" never appears over the loading screen.
      { const _hvcFab = document.getElementById('hvc-fab'); if (_hvcFab) _hvcFab.classList.remove('show'); if (typeof hvChatClose === 'function') hvChatClose(); }
      // Terminal mode owns the scrolling log in place of the per-source
      // checklist; start it immediately so the first prompt line lands
      // before the first SSE event arrives.
      terminalAnim.start({gene, hgvs_c});

      const phases = makePhaseTracker();
      const dbCompletion = Object.fromEntries(DB_SOURCES.map(s => [s.key, false]));

      let stage1Data = null;
      let stage2Done = false;
      let dbOnlyPayload = null;
      // Early-reveal state. `partialEvidence` accumulates from db_done so the
      // preliminary card renders a populated Evidence tab rather than an empty
      // one; `literaturePending` drives the "…" pill; `prelimRevealed` stops
      // db_only_complete replaying the reveal animation over a visible card.
      const partialEvidence = {};
      let literaturePending = false;
      let prelimRevealed = false;
      const stage2Chunks = [];
      // Reset cross-stage state so exports don't reuse the previous variant.
      lastResult = null;
      lastEvidence = null;
      lastVariantId = null;
      erepoVerdict = null;
      lastVariant = {
        gene, hgvs_c, hpo, family, genome_build,
        zygosity, inheritance_input, proband_sex, trio_status, denovo_status,
        // Structured family & segregation evidence — echoed into the
        // Summary tab + TXT/CSV exports so the curator can confirm the
        // family data that drove PP1/BS4/PM3/PS2/PM6/BS2/BP5.
        seg_affected_carriers, seg_affected_noncarriers,
        denovo_confirmed_count, denovo_unconfirmed_count,
        seg_meioses,
        in_trans_pathogenic, alt_cause_present, alt_cause_detail,
        zygosity_inferred: false,
        // Defaults to "unknown"; overwritten by the backend's classifier
        // result via _mergeClinicalContextFromEvidence when stage 2 or
        // db_only_complete lands.
        family_history_summary: 'unknown',
      };
      // Clear the Protein-tab store before the new run lands — partial
      // state from the previous variant must never leak into this one.
      _resetProteinData(gene, hgvs_c);
      // Same reasoning for the held-evidence token: a new run must never be
      // able to add AI to the PREVIOUS variant's gather.
      HV_INTERPRET.token = null;
      HV_INTERPRET.gene = '';
      HV_INTERPRET.hgvsC = '';

      try {
        const resp = await fetch(API_ENDPOINT, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json', 'Accept': 'text/event-stream' },
          body: JSON.stringify({
            gene, hgvs_c, hpo, family, ai_mode: aiMode, genome_build,
            zygosity, inheritance_input, proband_sex, trio_status, denovo_status,
            // Structured family & segregation evidence — field names match
            // CurationRequest exactly so the backend normaliser picks them
            // up directly. Always sent so curator family data is never
            // gated behind a flag.
            seg_affected_carriers, seg_affected_noncarriers,
        denovo_confirmed_count, denovo_unconfirmed_count,
            seg_meioses, in_trans_pathogenic,
            alt_cause_present, alt_cause_detail,
          })
        });
        if (!resp.ok) {
          // Auth/quota responses carry a structured body ({detail:{error,message}})
          // so they can be handled as states rather than shown as raw status codes.
          // Both only ever occur on an AI-bearing request — the evidence-only path
          // is public and unmetered.
          //   401 signin_required    — session lapsed between load and submit
          //   429 ai_quota_exhausted — the day's per-user AI allowance is spent
          let structured = null;
          if (resp.status === 401 || resp.status === 429) {
            try { structured = (await resp.json()).detail; } catch (e) {}
          }
          if (structured && structured.error === 'signin_required') {
            if (typeof hvAuthHandle401 === 'function') hvAuthHandle401();
            throw new Error(structured.message
              || 'Sign in to use the AI interpretation, or untick it to continue.');
          }
          if (structured && structured.error === 'ai_quota_exhausted') {
            throw new Error(structured.message);
          }
          const errText = await resp.text();
          // A stopped or restarting App Service does not serve OUR app at all —
          // the Azure platform front end returns its own HTML page, and
          // deploy/maintenance.html is served in its place where IT has wired
          // one up. Dumping that body into the error banner printed the whole
          // document, comments and CSS included, at a curator mid-curation
          // (reported 2026-08-31). The body is never useful to a user, so it
          // goes to the console for support and the banner gets the same
          // wording the maintenance page itself shows.
          //
          // ⚠ Like maintenance.html, this wording assumes a PLANNED refresh or
          // deploy. For an unplanned outage it under-states the problem — swap
          // for the commented line below, exactly as that file documents.
          if (_looksLikePlatformErrorPage(resp.status, errText)) {
            console.warn(
              `HeartVar: platform error page on ${resp.status}, body suppressed:`,
              errText.slice(0, 2000));
            throw new Error(
              // 'HeartVar is briefly offline for maintenance. Please try again in a few minutes.'
              'HeartVar is updating its reference databases. Please try again in 5 minutes.');
          }
          throw new Error(`Backend ${resp.status}: ${errText}`);
        }
        for await (const { event, data } of parseSse(resp.body)) {
          switch (event) {
            case 'stage1_start':
              // DB gather + stage 1 share the existing pipeline progress card.
              break;
            case 'db_pending':
              setDbRowState(data.source, 'pending');
              break;
            case 'db_done': {
              const payload = data.data || {};
              // Accumulated so the early reveal has real evidence to render
              // rather than an empty Evidence tab. The final db_only_complete
              // payload replaces it wholesale.
              if (data.source) partialEvidence[data.source] = payload;
              let rowState = 'done';
              // Subscription-gated sources return {available: false,
              // pending: true} — surface that as 'pending API access'
              // rather than a generic skip/failure.
              if (payload.pending === true) rowState = 'pending_access';
              else if (payload.skipped && payload.not_applicable) rowState = 'not_applicable';
              else if (payload.skipped) rowState = 'skipped';
              else if (payload.available === false) rowState = 'skipped';
              // Sources that distinguish transport failures from valid
              // "no result at this position" answers carry an explicit
              // lookup_failed flag (SpliceAI). When the flag
              // is False the API answered cleanly — that's a real
              // result and the row should still tick green even though
              // ok is False. Only flip to 'error' when the flag is
              // missing (legacy sources) or explicitly True.
              else if (payload.ok === false && payload.lookup_failed === false) rowState = 'done';
              else if (payload.ok === false) rowState = 'error';
              setDbRowState(data.source, rowState);
              // Patch any matching {placeholder} in the terminal script
              // with the just-arrived real value (no-op in original mode).
              terminalAnim.ingestSource(data.source, payload);
              // Eagerly land VEP + UniProt fields in the Protein-tab
              // store the moment each source completes, ahead of the
              // result. Helper merges fields per-source so arrivals
              // don't clobber each other.
              if (data.source === 'vep' || data.source === 'uniprot' || data.source === 'domain_plp') {
                _populateProteinData({ [data.source]: payload }, gene, hgvs_c);
              }
              dbCompletion[data.source] = true;
              if (Object.values(dbCompletion).every(Boolean)) {
                // In no-key mode there's no LLM work after the gather, so
                // we skip the "generating interpretation" phase — the
                // post-stream branch reveals the evidence-only card
                // straight from "DB lookups complete".
                if (aiMode !== 'none') {
                  phases.advance('generate');
                  // Hand off the DB checklist for the Step-3 simulation grid
                  // + slim top progress bar. The tabbed result card stays
                  // hidden until stage2_complete (see revealResults).
                  fakeCritSim.start();
                }
              }
              break;
            }
            case 'preliminary_classification': {
              // ── THE DETERMINISTIC RESULT, BEFORE LITERATURE LANDS ────────
              // The backend fires this the moment the last NON-literature
              // source is in. Literature is the critical path (pmcoa waits on
              // both pubmed and pubtator3, so it runs in series) and no ACMG
              // criterion reads any of it, so the tier here IS the tier that
              // arrives in db_only_complete — verified server-side by running
              // the identical scoring pipeline. Measured on a cold gene: 2.3 s
              // instead of 8.2 s (DSP), 5.6 s instead of 10.2 s (MYH7).
              //
              // ADDITIVE: if this event never arrives, or anything here
              // throws, the ordinary db_only_complete path still reveals the
              // full card exactly as before.
              try {
                literaturePending = true;
                window.__hvLiteraturePending = true;
                prelimRevealed = true;
                dbOnlyPayload = null;
                genProgress.finish();
                chainLoadingFinish(() => {
                  revealEvidenceOnlyResults(
                    Object.assign({}, data, { db_evidence: partialEvidence }),
                    gene, hgvs_c,
                  );
                });
              } catch (e) { prelimRevealed = false; literaturePending = false; }
              break;
            }
            case 'db_only_complete':
              // No-key mode: backend stopped after DB gather. Cache the
              // evidence; the post-stream branch builds the placeholder
              // card so Summary / Criteria / Recommendations show the
              // "add an API key" message while Evidence stays populated.
              dbOnlyPayload = data;
              lastEvidence = data.db_evidence;
              lastVariantId = data.variant_id;
              _populateProteinData(data.db_evidence, gene, hgvs_c);
              _mergeClinicalContextFromEvidence(data.db_evidence);
              // If the preliminary card is already on screen, replace its
              // contents in place — literature has landed, so the "…" pill
              // goes and the full evidence renders. Re-running the reveal
              // ANIMATION over a visible card would look like a second result
              // arriving, so the post-stream branch is disarmed by nulling
              // dbOnlyPayload. When there was no early reveal this is a no-op
              // and the original post-stream path runs untouched.
              if (prelimRevealed) {
                literaturePending = false;
                window.__hvLiteraturePending = false;
                dbOnlyPayload = null;
                try { revealEvidenceOnlyResults(data, gene, hgvs_c); }
                catch (e) { /* the card already shows the same tier */ }
              }
              break;
            case 'stage2_start':
              // Simulator already running; no-op.
              break;
            case 'stage2_chunk':
              if (typeof data.chunk === 'string') stage2Chunks.push(data.chunk);
              break;
            case 'stage2_complete': {
              // Single consolidated payload — classification, confidence,
              // summary, criteria and the server-checked Tavtigian
              // points_total all arrive together. Shaping is shared with the
              // post-hoc "Add AI interpretation" path (_shapeAiResult).
              const shaped = _shapeAiResult(data, gene, hgvs_c);
              const canonicalPoints = shaped.canonicalPoints;
              stage1Data = shaped.stage1Data;
              // Feed the terminal driver so any still-pending {score} /
              // {final_class} placeholders resolve when finish() flushes.
              terminalAnim.ingestStage1(stage1Data);
              genProgress.finish();
              chainLoadingFinish(() => {
                revealResults(stage1Data, gene, hgvs_c, canonicalPoints);
              });
              stage2Done = true;
              // A real AI result landed → enable the results-only chatbot
              // (which also runs on the server key). The evidence-only path
              // never reaches here, so the FAB stays hidden when AI is off or
              // was degraded to evidence-only.
              window.__hvAiUsed = true;
              break;
            }
            case 'erepo_verdict': {
              // ClinGen eRepo VCEP verdict — INDEPENDENT read-only reference.
              // Arrives after stage2_complete (the summary card may not be in
              // the DOM yet because revealResults animates in), so we cache it
              // and let renderErepoPanel() apply it both now and again after
              // the summary renders. Never influences HeartVar's own result.
              erepoVerdict = (data && data.found) ? data : null;
              renderErepoPanel();
              break;
            }
            case 'error': {
              throw new Error((data && data.message) || 'Unknown server error');
            }
          }
        }
        // Post-stream branching: AI-on lands stage1Data via the single
        // stage2_complete event (handled inline above); the evidence-only
        // path (AI off, or AI requested-but-unavailable) reveals its card
        // from the cached db_only_complete payload once the stream closes.
        if (dbOnlyPayload) {
          genProgress.finish();
          chainLoadingFinish(() => {
            revealEvidenceOnlyResults(dbOnlyPayload, gene, hgvs_c);
          });
        } else if (!stage1Data && !prelimRevealed) {
          // `prelimRevealed` matters: the early-reveal path deliberately nulls
          // dbOnlyPayload after re-rendering in place, so by the time the
          // stream closes there is a fully populated card on screen and
          // neither branch above should fire. Without this guard that state
          // is indistinguishable from "the stream died", and the card gets
          // replaced by an error.
          throw new Error('Stream ended before the model returned a result');
        }
      } catch (err) {
        clearInterval(critProgress.timer);
        critProgress.timer = null;
        fakeCritSim.reset();
        genProgress.finish();
        const _errMsg = String((err && err.message) || err || 'Something went wrong.')
          .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
        // Genomic-HGVS hint: a g. HGVS is assembly-specific and is NOT lifted
        // over (only chr:pos:ref:alt coordinates are), so a GRCh37 accession
        // is a common cause of the "reference base doesn't match" error here.
        // Surface the actionable hint on the error screen too — this is where
        // a wrong-build genomic HGVS most often lands (it never reaches the
        // Summary tab's flag).
        const _ghgvsErrHint = isGenomicHgvs(hgvs_c)
          ? `<div class="summary-ghgvs-note" role="note" style="margin-top:14px">
              <span class="summary-ghgvs-note__ic" aria-hidden="true">⚠</span>

              <span>You entered a genomic (<span class="ev-mono">g.</span>) HGVS, which is assembly-specific and is <strong>not lifted over</strong> (only genomic coordinates are). If this is a GRCh37 / hg19 accession (e.g. <span class="ev-mono">NC_000014.8</span>), re-enter the variant as genomic coordinates so it can be lifted over, or use the matching GRCh38 accession (<span class="ev-mono">NC_000014.9</span>).</span>

            </div>`

          : '';

        rp.innerHTML = `<div class="err">
          <strong>Error:</strong> ${_errMsg}
          ${_ghgvsErrHint}
          <div class="err-actions">
            <button type="button" class="err-back-btn" onclick="if(typeof returnToLanding==='function')returnToLanding();">← Back to search</button>
          </div>
        </div>`;
      } finally {
        clearInterval(critProgress.timer);
        critProgress.timer = null;
        clearInterval(claudeProgress.timer);
        claudeProgress.timer = null;
        if (btn) {
          btn.disabled = false;
          btn.innerHTML = 'Generate variant interpretation →';
        }
      }
    }

    function _shapeAiResult(data, gene, hgvs_c) {

      const criteria = normalizeCriteria(data && data.criteria);
      const canonicalPoints = Number.isFinite(data && data.points_total)
        ? data.points_total : null;
      const stage1Data = {
        classification: data && data.classification,
        confidence: data && data.confidence,
        summary: data && data.summary,
        points_total: canonicalPoints,

        borderline_reasoning: (data && data.borderline_reasoning) || null,
        vus_subclassification: (data && data.vus_subclassification) || null,

        gene_context: (data && data.gene_context) || {},
        db_evidence: data && data.db_evidence,
        variant_id: data && data.variant_id,
      };
      lastEvidence = stage1Data.db_evidence || lastEvidence;
      lastVariantId = stage1Data.variant_id || lastVariantId;
      if (stage1Data.db_evidence) {
        _populateProteinData(stage1Data.db_evidence, gene, hgvs_c);
        _mergeClinicalContextFromEvidence(stage1Data.db_evidence);
      }
      lastResult = {
        classification: stage1Data.classification,
        confidence: stage1Data.confidence,
        points_total: canonicalPoints,
        summary: stage1Data.summary,
        borderline_reasoning: stage1Data.borderline_reasoning,
        vus_subclassification: stage1Data.vus_subclassification,
        gene_context: stage1Data.gene_context,
        criteria,
      };
      return { stage1Data, criteria, canonicalPoints };
    }

    function revealResults(result, gene, hgvs_c, canonicalPoints, opts) {
      opts = opts || {};
      const panel = document.getElementById('result-panel');
      const progressCard = document.getElementById('progress-card');
      const swap = () => {

        renderStage1Card(result, gene, hgvs_c);

        renderErepoPanel();
        const criteria = (lastResult && lastResult.criteria) || [];
        if (criteria.length) {
          renderSummaryCriteriaChips(criteria);
          renderStage2Criteria(criteria, gene, hgvs_c, lastEvidence, lastVariantId);
        } else if (opts.criteriaErrorMessage) {

          const area = document.getElementById('criteria-area');
          if (area) {
            area.innerHTML = `<div class="crit-error">
              <strong>Could not load the per-criterion breakdown.</strong>
              ${opts.criteriaErrorMessage}
              The classification and summary remain valid.
            </div>`;
          }
          markCriteriaTabReady();
        }

        if (Number.isFinite(canonicalPoints)) {
          const scoreEl = document.getElementById('summary-score-value');
          if (scoreEl) {
            scoreEl.textContent = `${canonicalPoints >= 0 ? '+' : '−'}${Math.abs(canonicalPoints)} pts`;
            scoreEl.title = 'Tavtigian 2020 point sum — Pathogenic ≥10 · LP 6-9 · VUS 0-5 · LB -1 to -6 · B ≤-7';
          }
        }
        showResultTab('summary');
        showToast('Analysis complete');
      };
      if (progressCard) {

        progressCard.classList.add('load-fade-out');
        setTimeout(swap, 300);
      } else {
        swap();
      }
    }

    function _buildPreliminarySummary(classification, points, criteria) {
      const met = (criteria || []).filter(c => c.status === 'met');
      const metCodes = met.map(c => c.code).join(', ');
      const ptsStr = Number.isFinite(points)
        ? ` (${points >= 0 ? '+' : '−'}${Math.abs(points)} pts)` : '';
      const tier = classification || 'VUS';
      const lead = met.length
        ? `Preliminary classification: ${tier}${ptsStr}, derived from ${met.length} deterministic ACMG/AMP criteria (${metCodes}).`
        : `Preliminary classification: ${tier}${ptsStr}. No deterministic criteria reached their threshold for this variant.`;

      return `${lead} Judgement-based criteria were not assessed, so this tier may understate pathogenicity.`;
    }

    const HV_INTERPRET = { token: null, gene: '', hgvsC: '', busy: false };

    function _interpretStatus(text) {
      const el = document.getElementById('prelim-ai-status');
      if (el) el.textContent = text || '';
    }

    window.hvAddAiInterpretation = async function () {
      if (!HV_INTERPRET.token || HV_INTERPRET.busy) return;
      const btn = document.getElementById('prelim-ai-btn');
      HV_INTERPRET.busy = true;
      if (btn) { btn.disabled = true; btn.textContent = 'Interpreting…'; }
      _interpretStatus('Running the AI interpretation');
      try {
        const resp = await fetch('/api/curate/interpret', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json', 'Accept': 'text/event-stream' },
          body: JSON.stringify({ token: HV_INTERPRET.token }),
        });
        if (!resp.ok) {

          let detail = null;
          try { detail = (await resp.json()).detail; } catch (e) {}
          if (detail && detail.error === 'signin_required') {
            if (typeof hvAuthHandle401 === 'function') hvAuthHandle401();
          }
          if (resp.status === 410) HV_INTERPRET.token = null;
          throw new Error((detail && detail.message)
            || `The AI interpretation could not be started (${resp.status}).`);
        }
        let stage2 = null;
        for await (const { event, data } of parseSse(resp.body)) {
          switch (event) {
            case 'stage2_start':
            case 'stage2_chunk':

              break;
            case 'stage2_restart':
              _interpretStatus('The model stream dropped — retrying.');
              break;
            case 'stage2_complete':
              stage2 = data;
              break;
            case 'erepo_verdict':
              erepoVerdict = (data && data.found) ? data : null;
              break;
            case 'error':
              throw new Error((data && data.message) || 'Unknown server error');
          }
        }
        if (!stage2) throw new Error('The interpretation ended before a result arrived.');

        const shaped = _shapeAiResult(stage2, HV_INTERPRET.gene, HV_INTERPRET.hgvsC);

        window.__hvAiUsed = true;
        HV_INTERPRET.token = null;
        revealResults(
          shaped.stage1Data, HV_INTERPRET.gene, HV_INTERPRET.hgvsC,
          shaped.canonicalPoints,
        );
      } catch (err) {
        const msg = String((err && err.message) || err || 'Something went wrong.');
        _interpretStatus(msg);
        if (btn) {
          btn.textContent = HV_INTERPRET.token
            ? 'Retry AI interpretation' : 'Add AI interpretation';
          btn.disabled = !HV_INTERPRET.token;
        }
      } finally {
        HV_INTERPRET.busy = false;
      }
    };

    function revealEvidenceOnlyResults(payload, gene, hgvs_c) {
      const progressCard = document.getElementById('progress-card');
      const swap = () => {

        const criteria = normalizeCriteria(payload.criteria || payload.hard_coded_criteria || []);
        const points = Number.isFinite(payload.points_total) ? payload.points_total : null;
        const classification = payload.classification || null;
        const summaryText = _buildPreliminarySummary(classification, points, criteria);

        const fauxStage1 = {
          classification,
          confidence: null,
          summary: summaryText,
          points_total: points,

          gene_context: {},
          db_evidence: payload.db_evidence,
          variant_id: payload.variant_id,
          preliminary: true,
        };
        renderStage1Card(fauxStage1, gene, hgvs_c);

        lastResult = {
          classification, points_total: points, confidence: null,
          summary: summaryText,
          gene_context: {},
          criteria,
        };

        const evPill = document.getElementById('summary-class-pill');
        if (evPill) evPill.textContent = 'Preliminary · AI-free';

        const mainCol = document.querySelector('#result-panel .v11-summary__main');

        if (!window.__hvLiteraturePending) {
          const _stale = document.getElementById('lit-pending-pill');
          if (_stale) _stale.remove();
        }
        if (mainCol && !mainCol.querySelector('.prelim-banner')) {
          const banner = document.createElement('div');
          banner.className = 'prelim-banner';
          banner.setAttribute('style', 'margin:0 0 18px;padding:12px 14px;border:1px solid #e7d9c4;border-left:3px solid var(--maroon,#8a2310);background:#faf4ea;border-radius:10px;font-size:13.5px;line-height:1.55;color:var(--ink,#3a2f2a);');

          const _aiUnavail = payload.ai_unavailable;
          const _unavailNote = _aiUnavail === 'budget'
            ? '<strong>AI interpretation is temporarily unavailable</strong> — the shared daily AI limit has been reached, so this is the evidence-only result. Please try again tomorrow. '
            : (_aiUnavail === 'unconfigured'
                ? '<strong>AI interpretation is currently unavailable on this server</strong> — showing the evidence-only result. '
                : '');

          const _tok = payload.interpret_token || null;
          const _optInNote = (_aiUnavail || _tok)
            ? ''
            : 'Tick <em>“Include AI interpretation”</em> under the search box and re-run for a complete, more accurate interpretation. ';

          const _litPill = window.__hvLiteraturePending
            ? ' <span id="lit-pending-pill" class="lit-pending-pill">'
              + 'literature loading<span class="term-pending">…</span></span>'
            : '';
          const _bannerText = '<strong>Preliminary classification; generated without AI.</strong> ' + _unavailNote + _optInNote + _litPill;
          banner.innerHTML = _bannerText;
          if (_tok) {

            banner.innerHTML = '';
            const textCol = document.createElement('span');
            textCol.setAttribute('style', 'flex:1 1 auto');
            textCol.innerHTML = _bannerText;
            banner.appendChild(textCol);
            banner.setAttribute('style', banner.getAttribute('style')
              + ';display:flex;align-items:center;justify-content:space-between;gap:18px');
            HV_INTERPRET.token = _tok;
            HV_INTERPRET.gene = gene;
            HV_INTERPRET.hgvsC = hgvs_c;
            const btn = document.createElement('button');
            btn.id = 'prelim-ai-btn';
            btn.type = 'button';
            btn.textContent = 'Add AI interpretation';
            btn.setAttribute('style', 'appearance:none;border:1px solid var(--maroon,#8a2310);background:var(--maroon,#8a2310);color:#fdf8f1;font:inherit;font-size:13px;font-weight:600;letter-spacing:.01em;padding:7px 15px;border-radius:999px;cursor:pointer;white-space:nowrap');

            if (_aiUnavail) {
              btn.disabled = true;
              btn.setAttribute('style', btn.getAttribute('style') + ';opacity:.5;cursor:not-allowed');
              btn.title = _aiUnavail === 'budget'
                ? 'The shared daily AI limit has been reached \u2014 please try again tomorrow.'
                : 'AI interpretation is unavailable on this server.';
            }
            btn.addEventListener('click', function () { window.hvAddAiInterpretation(); });

            const actionCol = document.createElement('div');
            actionCol.setAttribute('style', 'flex:0 1 auto;display:flex;align-items:center;justify-content:flex-end;gap:11px');
            const status = document.createElement('span');
            status.id = 'prelim-ai-status';
            status.setAttribute('style', 'max-width:300px;font-size:12.5px;line-height:1.45;color:#6b5c52');
            actionCol.appendChild(status);
            actionCol.appendChild(btn);
            banner.appendChild(actionCol);
          }
          mainCol.insertBefore(banner, mainCol.firstChild);
        }

        if (criteria.length) {
          renderSummaryCriteriaChips(criteria);
          renderStage2Criteria(criteria, gene, hgvs_c, payload.db_evidence, payload.variant_id);
        } else {
          const critArea = document.getElementById('criteria-area');
          if (critArea) critArea.innerHTML = '<div class="crit-error">No criteria were returned for this variant.</div>';
          markCriteriaTabReady();
        }

        showResultTab('summary');
        showToast(payload.ai_unavailable
          ? 'Evidence-only result — AI was temporarily unavailable'
          : (payload.interpret_token
              ? 'Preliminary classification ready — press “Add AI interpretation” for the full AI result'
              : 'Preliminary classification ready — tick “Include AI interpretation” for the full AI result'));
      };
      if (progressCard) {
        progressCard.classList.add('load-fade-out');
        setTimeout(swap, 300);
      } else {
        swap();
      }
    }

    function _buildDbLinks(gene, hgvs_c, evidence, variantId) {

      const evVep = (evidence && evidence.vep) || {};
      gene = (gene || evVep.gene_symbol || evVep.derived_gene_symbol || '').trim();
      const gnomadUrl = variantId
        ? `https://gnomad.broadinstitute.org/variant/${variantId}?dataset=gnomad_r4`
        : `https://gnomad.broadinstitute.org/gene/${encodeURIComponent(gene)}?dataset=gnomad_r4`;
      const clinvarRecs = ((evidence && evidence.clinvar && evidence.clinvar.records) || []);
      const clinvarUrl = clinvarRecs.length
        ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${clinvarRecs[0].variation_id}/`
        : `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(gene + ' ' + hgvs_c)}`;
      const pubmedUrl = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(gene + ' ' + hgvs_c)}`;
      const uniprotUrl = `https://www.uniprot.org/uniprotkb?query=${encodeURIComponent('gene:' + gene + ' AND organism_id:9606')}`;
      const caddUrl = variantId
        ? `https://cadd.gs.washington.edu/snv?variants=${variantId.replace(/-/g, ' ')}`
        : 'https://cadd.gs.washington.edu/';
      const spliceaiUrl = variantId
        ? `https://spliceailookup.broadinstitute.org/#variant=${variantId}&hg=38`
        : 'https://spliceailookup.broadinstitute.org/';

      const vepUrl = gene
        ? `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(gene)}`
        : 'https://www.ensembl.org/Homo_sapiens/Tools/VEP';

      return {

        PVS1:  [{ label: 'Ensembl VEP', url: vepUrl }],
        PS1:   [{ label: 'ClinVar', url: clinvarUrl }],

        PS2:   [{ label: 'Clinical / family (curator-supplied)' }],
        PS3:   [{ label: 'PubMed', url: pubmedUrl }],

        PS4:   [{ label: 'PubMed', url: pubmedUrl }],
        PM1:   [{ label: 'UniProt', url: uniprotUrl }, { label: 'ClinVar (domain P/LP)', url: clinvarUrl }],
        PM2:   [{ label: 'gnomAD', url: gnomadUrl }],

        PM3:   [{ label: 'Clinical / family (curator-supplied)' }, { label: 'ClinVar (partner variant)', url: clinvarUrl }],

        PM4:   [{ label: 'Ensembl VEP', url: vepUrl }],
        PM5:   [{ label: 'ClinVar', url: clinvarUrl }],
        PM6:   [{ label: 'Clinical / family (curator-supplied)' }],

        PP1:   [{ label: 'Family segregation (curator-supplied)' }],
        PP2:   [{ label: 'gnomAD (gene constraint)', url: gnomadUrl }],

        PP3:   [{ label: 'in-silico (REVEL/CADD)', url: caddUrl }, { label: 'SpliceAI', url: spliceaiUrl }],
        PP4:   [{ label: 'Clinical phenotype (curator-supplied)' }],
        PP5:   [{ label: 'ClinVar', url: clinvarUrl }],
        BA1:   [{ label: 'gnomAD', url: gnomadUrl }],
        BS1:   [{ label: 'gnomAD', url: gnomadUrl }],
        BS2:   [{ label: 'gnomAD', url: gnomadUrl }],
        BS3:   [{ label: 'PubMed', url: pubmedUrl }],

        BS4:   [{ label: 'Family segregation (curator-supplied)' }],

        BP1:   [{ label: 'Gene mechanism' }, { label: 'UniProt', url: uniprotUrl }],

        BP2:   [{ label: 'Clinical / family (curator-supplied)' }, { label: 'ClinVar (partner variant)', url: clinvarUrl }],

        BP3:   [{ label: 'Ensembl VEP', url: vepUrl }],
        BP4:   [{ label: 'in-silico (REVEL/CADD)', url: caddUrl }, { label: 'SpliceAI', url: spliceaiUrl }],

        BP5:   [{ label: 'Clinical (curator-supplied)' }],
        BP6:   [{ label: 'ClinVar', url: clinvarUrl }],

        BP7:   [{ label: 'SpliceAI', url: spliceaiUrl }, { label: 'phyloP (UCSC)', url: 'https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way' }],
      };
    }

    function _ccTier(c) {
      let tw = strengthTier(c.criteria_strength);
      if (!tw) {
        const code = String(c.code || '');
        if (/^PVS/.test(code)) tw = 'VeryStrong';
        else if (/^BA/.test(code)) tw = 'Stand-alone';
        else if (/^(PS|BS)\d/.test(code)) tw = 'Strong';
        else if (/^[PB]M\d/.test(code)) tw = 'Moderate';
        else if (/^[PB]P\d/.test(code)) tw = 'Supporting';
      }
      return tw;
    }

    function _ccPoints(c) {
      const tw = _ccTier(c);
      let n = tw === 'Stand-alone' ? 8 : _TIER_POINTS[tw];
      if (n == null) n = Math.abs(criterionPoints(String(c.code || '')));
      return c.direction === 'benign' ? -n : n;
    }

    const HV_OVR = {
      base: [],
      effective: [],
      overrides: {},
      gene: '',
      inheritance: '',
      basePoints: null,
      baseTier: '',
      adjusted: null,
      pending: false,
      dbLinks: null,
    };
    window.__hvOverrideState = HV_OVR;

    const HV_TIER_LADDER = {
      pathogenic: ['VeryStrong', 'Strong', 'Moderate', 'Supporting'],
      benign: ['Strong', 'Moderate', 'Supporting'],
    };
    const HV_TIER_LABELS = {
      VeryStrong: 'Very Strong', Strong: 'Strong',
      Moderate: 'Moderate', Supporting: 'Supporting',
      'Stand-alone': 'Stand-alone',
    };

    function _hvDirectionOf(c) {
      if (c && (c.direction === 'benign' || c.direction === 'pathogenic')) return c.direction;
      return String((c && c.code) || '').charAt(0) === 'B' ? 'benign' : 'pathogenic';
    }

    function _hvStrengthOptions(c) {
      const code = String((c && c.code) || '').split('/')[0].trim();
      if (code === 'BA1') return [{ value: 'BA1', label: 'Stand-alone' }];
      const tiers = HV_TIER_LADDER[_hvDirectionOf(c)] || HV_TIER_LADDER.pathogenic;
      return tiers.map(t => ({ value: `${code}_${t}`, label: HV_TIER_LABELS[t] || t }));
    }

    function _hvSelectedStrength(c) {
      const code = String((c && c.code) || '').split('/')[0].trim();
      const ovr = HV_OVR.overrides[code];
      if (ovr && ovr.criteria_strength) return ovr.criteria_strength;
      if (c && c.criteria_strength) return c.criteria_strength;
      if (code === 'BA1') return 'BA1';
      const tier = _ccTier(c);
      return tier ? `${code}_${tier}` : '';
    }

    function _hvOverrideCount() { return Object.keys(HV_OVR.overrides).length; }

    function _hvInitOverrides(criteria, gene, hgvsC, inheritance, points, tier, dbLinks) {
      HV_OVR.base = (criteria || []).map(c => Object.assign({}, c));
      HV_OVR.effective = HV_OVR.base.map(c => Object.assign({}, c));
      HV_OVR.overrides = {};
      HV_OVR.gene = gene || '';
      HV_OVR.hgvsC = hgvsC || '';
      HV_OVR.inheritance = inheritance || '';
      HV_OVR.basePoints = (typeof points === 'number') ? points : null;
      HV_OVR.baseTier = tier || '';
      HV_OVR.adjusted = null;
      HV_OVR.dbLinks = dbLinks || null;
    }

    function _hvEditedCriteria() {
      return HV_OVR.base.map(c => {
        const code = String(c.code || '').split('/')[0].trim();
        const ovr = HV_OVR.overrides[code];
        return {
          code: code,
          status: ovr ? ovr.status : (c.status || 'not_met'),
          criteria_strength: ovr ? ovr.criteria_strength : (c.criteria_strength || null),
          direction: _hvDirectionOf(c),
          curator_override: !!ovr,
        };
      });
    }

    window.hvToggleCriterion = function(code, on) {
      const base = HV_OVR.base.find(c => String(c.code || '').split('/')[0].trim() === code);
      const wasMet = !!base && base.status === 'met';
      if (on === wasMet && !HV_OVR.overrides[code]) return;
      if (on === wasMet) {

        delete HV_OVR.overrides[code];
      } else {
        HV_OVR.overrides[code] = {
          status: on ? 'met' : 'not_met',
          criteria_strength: on ? _hvSelectedStrength(base || { code: code }) : null,
        };
      }
      _hvPushRescore();
    };

    window.hvSetCriterionStrength = function(code, value) {
      const base = HV_OVR.base.find(c => String(c.code || '').split('/')[0].trim() === code);
      const wasMet = !!base && base.status === 'met';

      const isRep = (code === 'PP5' || code === 'BP6');
      if (wasMet && !isRep && base.criteria_strength === value) {
        delete HV_OVR.overrides[code];
        _hvPushRescore();
        return;
      }

      HV_OVR.overrides[code] = { status: 'met', criteria_strength: value };
      _hvPushRescore();
    };

    window.hvResetOverrides = function() {
      HV_OVR.overrides = {};
      HV_OVR.adjusted = null;
      HV_OVR.effective = HV_OVR.base.map(c => Object.assign({}, c));
      _hvRenderCriteriaArea();
    };

    let _hvRescoreTimer = null;
    function _hvPushRescore() {
      if (_hvRescoreTimer) clearTimeout(_hvRescoreTimer);
      _hvRescoreTimer = setTimeout(_hvDoRescore, 120);

      _hvRenderCriteriaArea();
    }

    async function _hvDoRescore() {
      if (!_hvOverrideCount()) {
        HV_OVR.adjusted = null;
        HV_OVR.effective = HV_OVR.base.map(c => Object.assign({}, c));
        _hvRenderCriteriaArea();
        return;
      }
      HV_OVR.pending = true;
      try {
        const r = await fetch('/api/acmg/rescore', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            gene: HV_OVR.gene,
            inheritance_input: HV_OVR.inheritance,
            criteria: _hvEditedCriteria(),
          }),
        });
        if (!r.ok) throw new Error('rescore ' + r.status);
        const data = await r.json();

        const byCode = {};
        (data.criteria || []).forEach(c => { byCode[c.code] = c; });
        HV_OVR.effective = HV_OVR.base.map(c => {
          const code = String(c.code || '').split('/')[0].trim();
          const s = byCode[code];
          return s ? Object.assign({}, c, {
            status: s.status,
            criteria_strength: s.criteria_strength,
            curator_override: s.curator_override,
          }) : Object.assign({}, c);
        });
        HV_OVR.adjusted = {
          points_total: data.points_total,
          classification: data.classification,
        };
      } catch (_e) {

        HV_OVR.adjusted = HV_OVR.adjusted || null;
        showToast('Could not recompute the classification — check your connection.');
      } finally {
        HV_OVR.pending = false;
        _hvRenderCriteriaArea();
      }
    }

    function _hvRenderCriteriaArea() {
      const area = document.getElementById('criteria-area');
      if (!area || !HV_OVR.base.length) return;
      area.innerHTML = _hvCriteriaSectionsHTML();
    }

    function _hvOverrideSummaryList() {
      return Object.keys(HV_OVR.overrides).map(code => {
        const base = HV_OVR.base.find(c => String(c.code || '').split('/')[0].trim() === code);
        const ovr = HV_OVR.overrides[code];
        const wasMet = !!base && base.status === 'met';
        const wasStr = (base && base.criteria_strength) || '';
        const tier = strengthTier(ovr.criteria_strength) || ovr.criteria_strength;
        if (ovr.status !== 'met') return { code: code, text: `${code} switched OFF (HeartVar applied ${wasStr || 'it'})` };

        if (code === 'PP5' || code === 'BP6') {
          return { code: code, text: `${code} counted at ${tier} on your instruction — HeartVar scores reputable-source criteria 0 (ClinGen SVI 2018)` };
        }
        if (!wasMet) return { code: code, text: `${code} switched ON at ${tier}` };
        return { code: code, text: `${code} strength changed ${wasStr} → ${ovr.criteria_strength}` };
      });
    }
    window.__hvOverrideSummaryList = _hvOverrideSummaryList;

    function _hvOverrideReportLines() {
      const n = _hvOverrideCount();
      if (!n || !HV_OVR.adjusted) return [];
      const out = [
        '',
        `Classification (curator-adjusted): ${HV_OVR.adjusted.classification}`,
        `Points (curator-adjusted): ${formatCritPoints(HV_OVR.adjusted.points_total)}`,
        `Curator changes to the ACMG criteria (${n}):`,
      ];
      _hvOverrideSummaryList().forEach(o => out.push(`  - ${o.text}`));
      out.push(
        '  The adjusted tier was recomputed by the same ACMG engine, including',
        "  the cross-criterion exclusions and the ACMG-2015 benign combining floor.",
      );
      return out;
    }

    function _ccard(c, DB_LINKS, opts) {
      opts = opts || {};
      const esc = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const isMet = c.status === 'met';
      const cls = (isMet && c.direction === 'pathogenic') ? 'cp'
                : (isMet && c.direction === 'benign') ? 'cb' : 'cn';
      const codeRaw = String(c.code || '');
      const firstCode = (codeRaw.split('/')[0] || '').trim();
      const def = CRITERIA_DEFINITIONS[codeRaw] || CRITERIA_DEFINITIONS[firstCode];

      const isRepCode = (firstCode === 'PP5' || firstCode === 'BP6');
      const isRepAssertion = isRepCode && !HV_OVR.overrides[firstCode];

      const tw = isMet ? _ccTier(c) : null;
      const twDisplay = tw === 'VeryStrong' ? 'Very Strong' : (tw || '');
      let pillLabel, pillCls;
      if (isMet && isRepAssertion) {
        pillLabel = 'Reference'; pillCls = 'is-na';
      } else if (isMet) {
        pillLabel = twDisplay || 'Met';
        if (c.direction === 'benign') pillCls = 'is-ben';
        else if (tw === 'Strong' || tw === 'VeryStrong' || tw === 'Stand-alone') pillCls = 'is-strong';
        else if (tw === 'Moderate') pillCls = 'is-mod';
        else pillCls = 'is-supp';
      } else if (c.status === 'insufficient_data') {
        pillLabel = 'Not assessed'; pillCls = 'is-na';
      } else if (c.status === 'not_assessed') {

        pillLabel = 'Not assessed'; pillCls = 'is-na';
      } else {
        pillLabel = (c.status === 'na') ? 'N/A' : 'Not met'; pillCls = 'is-na';
      }

      const ptsLabel = (isMet && isRepAssertion) ? '0'
        : isMet ? formatCritPoints(_ccPoints(c)) : '—';
      const ptsCls = (isMet && isRepAssertion) ? 'is-none'
        : isMet ? (c.direction === 'benign' ? 'is-ben' : 'is-path') : 'is-none';

      const links = DB_LINKS[codeRaw] || DB_LINKS[firstCode] || [];
      const sourcesHTML = links.map(l =>
        l.url
          ? `<a class="crit-acc__src" href="${l.url}" target="_blank" rel="noopener">${esc(l.label)} ↗</a>`
          : `<span class="crit-acc__src is-nolink">${esc(l.label)}</span>`
      ).join('');
      let rationale;
      if (isMet && isRepAssertion) {
        const assertion = c.direction === 'benign'
          ? 'a benign / likely-benign assertion in ClinVar'
          : 'a pathogenic / likely-pathogenic assertion in ClinVar';
        rationale = `This variant carries ${assertion}. Shown for reference only — it contributes 0 points to the total. `
          + `Reusing ClinVar's own classification as scoring evidence would be circular (it would let ClinVar's verdict drive this one), `
          + `so per ClinGen SVI guidance (Biesecker & Harrison 2018) ${firstCode} is excluded from the point-based tier.`;
      }
      else if (isMet) rationale = `Applied at ${twDisplay || 'base'} strength (${ptsLabel} pts).`;
      else if (c.status === 'insufficient_data') rationale = 'The data required to evaluate this criterion was not available for this variant.';

      else if (c.status === 'not_assessed') rationale = SERVER_OWNED_CODES.has(firstCode)
        ? 'Derived in Python from curator-entered structured inputs, which were not supplied. Enabling AI interpretation will not produce this criterion.'
        : 'Requires AI / clinical-geneticist judgement — not evaluated in this preliminary (AI-free) classification. Enable AI interpretation for a complete interpretation.';
      else if (c.status === 'na') rationale = 'Not applicable to this variant under the ACMG/AMP framework.';
      else rationale = 'Evaluated, but the threshold for this criterion was not met.';

      const rows = [];
      if (def && def.definition) rows.push(`<span class="crit-acc__k">Definition</span><span class="crit-acc__v">${esc(def.definition)}</span>`);
      if (c.evidence) rows.push(`<span class="crit-acc__k">Evidence</span><span class="crit-acc__v">${linkifyPmids(esc(c.evidence))}</span>`);
      rows.push(`<span class="crit-acc__k">Rationale</span><span class="crit-acc__v">${esc(rationale)}</span>`);

      if (isMet && sourcesHTML) rows.push(`<span class="crit-acc__k">Sources</span><span class="crit-acc__v crit-acc__sources">${sourcesHTML}</span>`);

      const ovrHTML = _hvCritControlHTML(c, firstCode, isRepCode);
      const wasEdited = !!HV_OVR.overrides[firstCode];
      return `<div class="crit-acc ${cls}${wasEdited ? ' is-overridden' : ''}">
        <button type="button" class="crit-acc__hd">
          <span class="crit-acc__code">${esc(codeRaw)}</span>
          <span class="crit-acc__name">${esc(c.name || '')}</span>
          ${wasEdited ? '<span class="crit-acc__ovrbadge" title="You changed this criterion">edited</span>' : ''}
          <span class="crit-acc__pill ${pillCls}">${esc(pillLabel)}</span>
          <span class="crit-acc__pts ${ptsCls}">${ptsLabel}</span>
          <span class="crit-acc__chev" aria-hidden="true">›</span>
        </button>
        <div class="crit-acc__detail"><div class="crit-acc__inner">${rows.join('')}</div></div>
        ${ovrHTML}
      </div>`;
    }

    function _hvCritControlHTML(c, code, isRepCode) {
      if (!code || !HV_OVR.base.length) return '';
      const base = HV_OVR.base.find(b => String(b.code || '').split('/')[0].trim() === code);
      const ovr = HV_OVR.overrides[code];

      const on = ovr ? ovr.status === 'met' : (!!base && base.status === 'met');

      const withheld = !!ovr && ovr.status === 'met' && c.status !== 'met';
      const opts = _hvStrengthOptions(c).map(o =>
        `<option value="${_hvEsc(o.value)}"${o.value === _hvSelectedStrength(c) ? ' selected' : ''}>${_hvEsc(o.label)}</option>`
      ).join('');

      const svi = (isRepCode && ovr)
        ? `<span class="crit-ovr__warn is-on" title="ClinGen SVI (Biesecker &amp; Harrison 2018) recommends retiring the reputable-source criteria, and HeartVar scores them 0. You have chosen to count this one; the change is recorded in the exported report.">counted at your request — SVI 2018 retires this criterion</span>`
        : '';
      return `<div class="crit-ovr">
        <label class="crit-ovr__sw" title="Include or exclude this criterion from the classification">
          <input type="checkbox"${on ? ' checked' : ''} onchange="hvToggleCriterion('${_hvEsc(code)}', this.checked)" aria-label="Apply ${_hvEsc(code)}">
          <span class="crit-ovr__track" aria-hidden="true"></span>
          <span class="crit-ovr__lbl">${on ? 'Applied' : 'Excluded'}</span>
        </label>
        <label class="crit-ovr__strwrap">
          <span class="crit-ovr__strlbl">Strength</span>
          <select class="crit-ovr__str"${on ? '' : ' disabled'} onchange="hvSetCriterionStrength('${_hvEsc(code)}', this.value)" aria-label="${_hvEsc(code)} strength">${opts}</select>
        </label>
        ${svi}
        ${withheld ? '<span class="crit-ovr__withheld" title="A cross-criterion rule withheld this — see the Evidence line above">withheld by a conflict rule</span>' : ''}
      </div>`;
    }

    function _rsec(title, items, DB_LINKS, opts) {
      if (!items.length) return '';
      opts = opts || {};
      const kind = opts.kind || 'na';
      const dotCls = kind === 'met' ? 'is-met' : kind === 'insuf' ? 'is-insuf' : 'is-na';
      let count;
      if (kind === 'met') {

        const sum = items.reduce((a, c) => {
          const fc = (String(c.code || '').split('/')[0] || '').trim();
          if ((fc === 'PP5' || fc === 'BP6') && !HV_OVR.overrides[fc]) return a;
          return a + (_ccPoints(c) || 0);
        }, 0);
        count = `${items.length} · ${formatCritPoints(sum)} pts`;
      } else {
        count = String(items.length);
      }
      return `<div class="crit-group">
        <div class="crit-group__hd">
          <span class="crit-group__dot ${dotCls}" aria-hidden="true"></span>
          <span class="crit-group__title">${title}</span>
          <span class="crit-group__count">${count}</span>
        </div>
        ${items.map(c => _ccard(c, DB_LINKS, opts)).join('')}
      </div>`;
    }

    function toggleClingenMore(btn) {
      const rows = btn.closest('.gc-validity-rows');
      if (!rows) return;
      const expanded = rows.classList.toggle('expanded-clingen');
      const more = btn.dataset.moreCount || '0';
      btn.textContent = expanded ? 'Show fewer' : `+ ${more} more curation(s)`;
    }

    const _hvEsc = s => String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

    function _transcriptModelLabel(evidence, explicitId) {
      const vep = (evidence && evidence.vep) || {};
      const txId = explicitId || vep.selected_transcript_id || vep.transcript_id || '';
      if (!txId) return '';

      const mane = vep.is_mane_select
        ? 'MANE Select'
        : ((vep.is_mane_clinical || vep.is_mane_plus_clinical) ? 'MANE Plus Clinical' : '');
      return `<span class="hv-model" title="The gene model this view is drawn against — exon ranks and c. positions are relative to this transcript">Model <span class="hv-model__id">${_hvEsc(txId)}</span>${mane ? ` <span class="hv-model__badge">${mane}</span>` : ''}</span>`;
    }

    function _proteinModelLabel(opts) {
      opts = opts || {};
      const pd = window.heartvarProteinData || {};
      const acc = opts.accession || pd.uniprotAccession || '';
      if (!acc) return '';
      const len = opts.length || pd.uniprotLength || null;
      return `<span class="hv-model" title="The protein model this view is drawn against — residue numbers are relative to this UniProt sequence">Model <span class="hv-model__id">UniProt ${_hvEsc(acc)}</span>${len ? ` <span class="hv-model__badge">${_hvEsc(len)} aa</span>` : ''}</span>`;
    }

    const _CSQ_LABELS = {
      missense: 'Missense', truncating: 'Truncating', splice: 'Splice/intronic',
      inframe: 'In-frame', synonymous: 'Synonymous', utr: 'UTR', other: 'Other',
    };
    const _CSQ_ORDER = ['missense', 'truncating', 'splice', 'inframe', 'synonymous', 'utr', 'other'];
    function _csqLabel(k) { return _CSQ_LABELS[k] || k; }

    function _csqBucketFromVepTerm(term) {
      term = (term || '').toLowerCase();
      if (term.indexOf('missense') >= 0) return 'missense';
      if (term.indexOf('stop_gained') >= 0 || term.indexOf('frameshift') >= 0
        || term.indexOf('start_lost') >= 0 || term.indexOf('stop_lost') >= 0
        || term.indexOf('transcript_ablation') >= 0) return 'truncating';
      if (term.indexOf('splice') >= 0 || term.indexOf('intron') >= 0) return 'splice';
      if (term.indexOf('synonymous') >= 0) return 'synonymous';
      if (term.indexOf('inframe') >= 0) return 'inframe';
      if (term.indexOf('utr') >= 0) return 'utr';
      return 'other';
    }

    function _buildLandscapeChips(data, selected) {
      const allTotal = Object.values(data.tierCounts || {}).reduce((s, n) => s + (n || 0), 0);
      const chip = (key, label, n, on) =>
        `<button type="button" class="hv-ls-chip${on ? ' is-active' : ''}" data-csq="${key}" onclick="window.__hvSetLandscapeCsq && window.__hvSetLandscapeCsq('${key}')" style="font-size:11.5px;padding:3px 10px;border-radius:999px;border:1px solid var(--line,#e6ddd4);cursor:pointer;background:${on ? 'var(--maroon,#8c1a1f)' : '#fff'};color:${on ? '#fff' : 'var(--ink-soft,#6a5d54)'}">${label} ${n}</button>`;
      const chips = [chip('all', 'All types', allTotal, selected === 'all')];
      for (const k of _CSQ_ORDER) {
        const c = (data.tierByCsq || {})[k];
        if (!c) continue;
        const n = Object.values(c).reduce((s, x) => s + (x || 0), 0);
        if (!n) continue;
        chips.push(chip(k, _csqLabel(k), n, selected === k));
      }
      return `<div class="hv-ls-chips" style="display:flex;flex-wrap:wrap;gap:6px;margin-bottom:10px">${chips.join('')}</div>`;
    }

    function _buildLandscapeLollipop(data, csq, showVus) {
      let pts = data.positions || [];
      if (csq !== 'all') pts = pts.filter(p => p.csq === csq);
      const exons = (data.exons && data.exons.length)
        ? [...data.exons].sort((a, b) => a.rank - b.rank) : null;
      if (!pts.length && !exons) return '';

      const orient = (((exons && data.exonStrand != null) ? data.exonStrand : data.strand) < 0) ? -1 : 1;
      const W = 1000, padX = 16, midY = 80, reach = 54;
      const innerW = W - 2 * padX;

      let xOf;
      const exonX0 = [], exonPx = [];
      if (exons) {
        const nE = exons.length, nGaps = Math.max(0, nE - 1);
        const gapPx = Math.min(6, (innerW * 0.30) / Math.max(1, nGaps));
        const exonLens = exons.map(e => Math.max(1, Math.abs(e.end - e.start) + 1));
        const totalExonLen = exonLens.reduce((s, l) => s + l, 0) || 1;
        const exonPxTotal = Math.max(1, innerW - gapPx * nGaps);
        for (let i = 0; i < nE; i++) exonPx.push((exonLens[i] / totalExonLen) * exonPxTotal);
        let cx = padX;
        for (let i = 0; i < nE; i++) { exonX0.push(cx); cx += exonPx[i] + gapPx; }
        xOf = (g) => {
          for (let i = 0; i < nE; i++) {
            const e = exons[i], eLo = Math.min(e.start, e.end), eHi = Math.max(e.start, e.end);
            if (g >= eLo && g <= eHi) {
              const tx5 = orient < 0 ? e.end : e.start;
              const frac = Math.min(1, Math.max(0, Math.abs(g - tx5) / exonLens[i]));
              return exonX0[i] + frac * exonPx[i];
            }
          }
          for (let i = 0; i < nE - 1; i++) {
            const t3 = orient < 0 ? exons[i].start : exons[i].end;
            const t5 = orient < 0 ? exons[i + 1].end : exons[i + 1].start;
            if (g > Math.min(t3, t5) && g < Math.max(t3, t5)) {
              const frac = Math.min(1, Math.max(0, Math.abs(g - t3) / (Math.abs(t5 - t3) || 1)));
              return exonX0[i] + exonPx[i] + frac * gapPx;
            }
          }
          const first = exons[0], tx5f = orient < 0 ? first.end : first.start;
          return ((orient < 0) ? (g > tx5f) : (g < tx5f)) ? padX : (padX + innerW);
        };
      } else {
        let lo = data.gMin, hi = data.gMax;
        if (lo == null || hi == null) {
          const gs = pts.map(p => p.gpos);
          lo = gs.length ? Math.min.apply(null, gs) : 0;
          hi = gs.length ? Math.max.apply(null, gs) : 1;
        }
        if (data.probandGpos != null) { lo = Math.min(lo, data.probandGpos); hi = Math.max(hi, data.probandGpos); }
        const span = (hi - lo) || 1;
        xOf = (g) => { const f = (g - lo) / span; return padX + (orient < 0 ? (1 - f) : f) * innerW; };
      }

      const geneTop = 168, geneH = 11, geneMid = geneTop + geneH / 2;
      const H = exons ? 210 : 174;
      const endsY = exons ? 162 : 160;
      const pinBottom = exons ? (geneTop + geneH + 2) : 150;
      const probandLabelY = exons ? 202 : 166;
      const benignLabelY = exons ? 148 : 147;
      const TIER = {
        P: { c: '#c8341a', dir: -1 }, LP: { c: '#d96b3a', dir: -1 },
        B: { c: '#2a6e3a', dir: 1 }, LB: { c: '#8db888', dir: 1 },
        VUS: { c: '#f0ad4e', dir: 0 },
      };
      const ln = n => Math.log(1 + n);
      const maxCount = Math.max.apply(null, [1].concat(pts.map(p => p.count)));
      const hOf = c => 10 + (ln(c) / ln(maxCount + 1)) * (reach - 10);
      const rOf = c => Math.min(6, 2.2 + ln(c));
      let needles = '', vusMarks = '', vusCount = 0;

      const opOf = st => ((st || 0) >= 2 ? 0.95 : 0.5);

      const _clinvarTx = data.clinvarTxId || '';
      const _pResC = data.probandResidueClinvar;
      const _pResM = data.probandResidueModel;
      const _numDiffers = _pResC != null && _pResM != null && _pResC !== _pResM;
      const _codonSet = Array.isArray(data.codonPositions) ? data.codonPositions : null;
      const aaTipFor = (p) => {
        const inCodon = !!(_codonSet && _codonSet.indexOf(p.gpos) >= 0);
        const atProbandResidue = (p.aa != null && _pResC != null && p.aa === _pResC) || inCodon;
        if (p.aa == null && !inCodon) return '';
        if (atProbandResidue && _numDiffers) {

          return ` · residue ${_pResM} on this model`
            + ` (${_pResC} in ClinVar${_clinvarTx ? ' · ' + _clinvarTx : ''})`
            + (p.aa == null ? ' · residue from position; the ClinVar name has none' : '');
        }
        if (p.aa == null) {

          return ` · residue ${_pResM != null ? _pResM : _pResC}`
            + ' · from position; the ClinVar name has none';
        }

        return ` · residue ${p.aa}`
          + (_numDiffers ? ` (ClinVar numbering${_clinvarTx ? ' · ' + _clinvarTx : ''})` : '');
      };
      for (const p of pts) {
        const t = TIER[p.tier]; if (!t) continue;
        const x = xOf(p.gpos).toFixed(1);
        const aaTip = aaTipFor(p);
        if (t.dir === 0) {
          vusCount++;
          vusMarks += `<circle cx="${x}" cy="${midY}" r="2.1" fill="${t.c}" opacity="${opOf(p.stars)}"><title>VUS${aaTip} · ${p.count} record${p.count > 1 ? 's' : ''} · ${p.stars}★</title></circle>`;
          continue;
        }
        const op = opOf(p.stars);
        const y2 = (midY + t.dir * hOf(p.count)).toFixed(1);
        const r = rOf(p.count).toFixed(1);
        needles += `<line x1="${x}" y1="${midY}" x2="${x}" y2="${y2}" stroke="${t.c}" stroke-width="1.1" opacity="${op}"></line>`
          + `<circle cx="${x}" cy="${y2}" r="${r}" fill="${t.c}" opacity="${op}"><title>${p.tier}${aaTip} · ${p.count} record${p.count > 1 ? 's' : ''} · ${p.stars}★</title></circle>`;
      }
      const baseline = `<line x1="${padX}" y1="${midY}" x2="${W - padX}" y2="${midY}" stroke="#d8cfc6" stroke-width="1"></line>`;

      let geneModel = '';
      if (exons) {
        geneModel = `<line x1="${padX}" y1="${geneMid}" x2="${W - padX}" y2="${geneMid}" stroke="#cdbfb2" stroke-width="1.2"></line>`;
        for (const e of exons) {
          const xa = xOf(e.start), xb = xOf(e.end);
          const x1 = Math.min(xa, xb), w = Math.max(1.2, Math.abs(xb - xa));
          const isP = data.probandExonRank != null && e.rank === data.probandExonRank;
          const fill = isP ? '#cfe0f0' : '#e7dccb';
          const stroke = isP ? '#1a4f8a' : '#c2b29c';
          geneModel += `<rect x="${x1.toFixed(1)}" y="${geneTop}" width="${w.toFixed(1)}" height="${geneH}" rx="1.5" fill="${fill}" stroke="${stroke}" stroke-width="${isP ? 1.3 : 0.6}"><title>exon ${e.rank}</title></rect>`;
          if (w >= 15) {
            geneModel += `<text x="${(x1 + w / 2).toFixed(1)}" y="${geneTop + geneH + 8}" text-anchor="middle" font-size="7.5" fill="${isP ? '#1a4f8a' : '#8a7d72'}" font-weight="${isP ? 700 : 400}">${e.rank}</text>`;
          }
        }
      }
      const ends = `<text x="${padX}" y="${endsY}" font-size="10" fill="#8a7d72">5′</text>`
        + `<text x="${W - padX}" y="${endsY}" text-anchor="end" font-size="10" fill="#8a7d72">3′</text>`;
      let pin = '';
      if (data.probandGpos != null) {
        const px = xOf(data.probandGpos).toFixed(1);
        const exonTxt = data.probandExon ? ` · exon ${_hvEsc(data.probandExon)}` : '';
        pin = `<line x1="${px}" y1="6" x2="${px}" y2="${pinBottom}" stroke="#1a4f8a" stroke-width="1.4" stroke-dasharray="3 2"></line>`
          + `<text x="${px}" y="${probandLabelY}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#1a4f8a">▲ this variant${exonTxt}</text>`;
      }
      const vusGroup = vusCount ? `<g class="lolli-vus" style="display:${showVus ? '' : 'none'}">${vusMarks}</g>` : '';
      const svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" preserveAspectRatio="xMidYMid meet" role="img" aria-label="ClinVar positional variant landscape for ${_hvEsc(data.gene)}" style="display:block;max-width:100%">
        <text x="${padX}" y="14" font-size="10" fill="#a99a8d">↑ pathogenic</text>
        <text x="${padX}" y="${benignLabelY}" font-size="10" fill="#a99a8d">↓ benign</text>
        ${baseline}${geneModel}${ends}${vusGroup}${needles}${pin}
      </svg>`;
      const exonNote = exons ? ' Boxes are exons (5′→3′); the proband’s exon is highlighted in blue.' : '';
      const trunc = data.truncated
        ? `<div class="gc-bar-cap-note" style="margin-top:4px;color:var(--amber)">⚠ Dense gene — showing the highest-priority variants; some lower-confidence ones omitted.</div>` : '';
      const toggle = vusCount
        ? `<button type="button" class="gc-lolli-vus-toggle" style="margin-top:6px;font-size:12px;background:none;border:1px solid var(--line,#e6ddd4);border-radius:6px;padding:3px 9px;cursor:pointer;color:var(--ink-soft,#6a5d54)" onclick="window.__hvToggleLandscapeVus && window.__hvToggleLandscapeVus()">${showVus ? 'Hide' : 'Show'} ${vusCount} VUS</button>` : '';

      const modelNote = data.modelLabel
        ? `<div class="hv-model-note">${data.modelLabel}</div>` : '';

      const numberingNote = _numDiffers
        ? `<div class="hv-model-note">Residue numbers on the marks are ClinVar's`
          + `${_clinvarTx ? ` (${_hvEsc(_clinvarTx)})` : ''}. This variant is`
          + ` residue ${_hvEsc(_pResC)} there and residue ${_hvEsc(_pResM)} on the`
          + ` model above — the same residue.</div>`
        : '';
      return `<div class="gc-lolli-wrap" style="margin-top:12px">
        <div class="gc-bar-cap-note" style="margin-bottom:2px">Positional landscape — each mark is a ClinVar variant at its position along the gene (5′→3′; needle size ∝ records; faded = &lt;2★). The blue line marks this variant.${exonNote}</div>
        ${modelNote}
        ${numberingNote}
        ${svg}
        ${toggle}
        ${trunc}
      </div>`;
    }

    function _buildLandscapeBody(data, csq, showVus) {
      const counts = csq === 'all'
        ? (data.tierCounts || {})
        : ((data.tierByCsq || {})[csq] || { P: 0, LP: 0, VUS: 0, LB: 0, B: 0 });
      const segs = [
        { count: counts.P || 0, label: 'P', cls: 'gc-seg-p', color: '#c8341a' },
        { count: counts.LP || 0, label: 'LP', cls: 'gc-seg-lp', color: '#d96b3a' },
        { count: counts.VUS || 0, label: 'VUS', cls: 'gc-seg-vus', color: '#f0ad4e' },
        { count: counts.LB || 0, label: 'LB', cls: 'gc-seg-lb', color: '#8db888' },
        { count: counts.B || 0, label: 'B', cls: 'gc-seg-b', color: '#2a6e3a' },
      ];
      const total = segs.reduce((s, x) => s + x.count, 0);
      if (!total) {
        return `<div class="gc-empty">No ${csq === 'all' ? '' : _csqLabel(csq).toLowerCase() + ' '}ClinVar records for ${_hvEsc(data.gene)}.</div>`;
      }
      const bar = segs.map(s => {
        if (!s.count) return '';
        const pct = (s.count / total) * 100;
        const txt = pct >= 8 ? `${s.label} ${s.count}` : '';
        return `<div class="gc-bar-seg ${s.cls}" style="flex:${s.count}" title="${s.label}: ${s.count} (${pct.toFixed(1)}%)">${txt}</div>`;
      }).join('');
      const legend = segs.map(s =>
        `<span><span class="gc-bar-legend-swatch" style="background:${s.color}"></span>${s.label} (${s.count})</span>`
      ).join('');
      return `<div class="gc-bar">${bar}</div>
        <div class="gc-bar-legend">${legend}</div>
        ${_buildLandscapeLollipop(data, csq, showVus)}`;
    }

    window.__hvSetLandscapeCsq = function (csq) {
      const s = window.__hvLandscapeState; if (!s) return;
      s.csq = csq; s.showVus = false;
      const el = document.getElementById('hv-ls-chart');
      if (el) el.innerHTML = _buildLandscapeBody(s.data, s.csq, s.showVus);
      document.querySelectorAll('.hv-ls-chip').forEach(c => {
        const on = c.dataset.csq === csq;
        c.classList.toggle('is-active', on);
        c.style.background = on ? 'var(--maroon,#8c1a1f)' : '#fff';
        c.style.color = on ? '#fff' : 'var(--ink-soft,#6a5d54)';
      });
    };
    window.__hvToggleLandscapeVus = function () {
      const s = window.__hvLandscapeState; if (!s) return;
      s.showVus = !s.showVus;
      const el = document.getElementById('hv-ls-chart');
      if (el) el.innerHTML = _buildLandscapeBody(s.data, s.csq, s.showVus);
    };

    function _renderGeneContextTab(ev, geneContext, gene) {
      ev = ev || {};
      geneContext = geneContext || {};

      gene = (gene || (ev.vep && (ev.vep.gene_symbol || ev.vep.derived_gene_symbol)) || '').trim();

      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
      const tierPillClass = (raw) => {
        const t = String(raw || '').toLowerCase();
        if (t.includes('definitive') || t.startsWith('well')) return 'gc-pill gc-pill-definitive';
        if (t === 'strong' || t.includes('strong evidence')) return 'gc-pill gc-pill-strong';
        if (t === 'moderate' || t.includes('moderate evidence')) return 'gc-pill gc-pill-moderate';
        if (t.includes('supportive')) return 'gc-pill gc-pill-supportive';
        if (t.includes('disputed') || t.includes('refuted')) return 'gc-pill gc-pill-disputed';
        return 'gc-pill gc-pill-limited';
      };

      const claudeNote = (label, sentence) => {
        if (sentence && String(sentence).trim()) {
          return `<div class="gc-claude-note">
            <span class="gc-claude-note-label">${label}</span>${esc(sentence)}
          </div>`;
        }
        return `<div class="gc-claude-note gc-claude-placeholder">
          <span class="gc-claude-note-label">${label}</span>Enable AI interpretation for AI-derived gene-level synthesis.
        </div>`;
      };

      const sCase = s => {
        if (s === null || s === undefined || s === '') return '';
        const t = String(s).replace(/_/g, ' ').toLowerCase();
        return t.charAt(0).toUpperCase() + t.slice(1);
      };
      const kv = (label, value) =>
        `<div class="ev-kv"><div class="ev-kv-label">${label}</div><div class="ev-kv-value">${value}</div></div>`;
      const kvGroup = items => `<div class="ev-kv-group">${items.filter(Boolean).join('')}</div>`;
      const clsClass = raw => {
        const s = String(raw || '').toLowerCase();
        if (!s) return '';
        if (/conflict/.test(s)) return 'ev-cls ev-cls-conflict';
        if (/likely pathogenic/.test(s)) return 'ev-cls ev-cls-lp';
        if (/pathogenic/.test(s))        return 'ev-cls ev-cls-path';
        if (/likely benign/.test(s))     return 'ev-cls ev-cls-lb';
        if (/benign/.test(s))            return 'ev-cls ev-cls-ben';
        if (/uncertain|vus/.test(s))     return 'ev-cls ev-cls-vus';
        if (/definitive/.test(s))        return 'ev-cls ev-cls-definitive';
        if (/strong/.test(s))            return 'ev-cls ev-cls-strong';
        if (/moderate/.test(s))          return 'ev-cls ev-cls-moderate';
        if (/limited|animal/.test(s))    return 'ev-cls ev-cls-limited';
        if (/disputed|refuted/.test(s))  return 'ev-cls ev-cls-disputed';
        return 'ev-cls';
      };

      const gc = ev.gencc || {};
      const clingenSubs = (gc && gc.clingen_submissions) || [];

      let clingenContent;
      if (clingenSubs.length) {

        const hpoSubmitted = Boolean(gc.submitted_hpo);
        const phenoMatched = (gc.clingen_phenotype_matched_count || 0) > 0;
        const displayBest = phenoMatched
          ? gc.clingen_phenotype_matched_best_classification
          : (hpoSubmitted ? null : gc.best_clingen_classification);
        const displayDisputed = phenoMatched
          ? gc.clingen_phenotype_matched_has_disputed_or_refuted
          : false;
        const flag = displayDisputed
          ? ' <span title="Disputed/Refuted record present" style="color:#c8341a">⚠</span>'
          : '';
        const metaText = phenoMatched
          ? `Phenotype match · ${gc.clingen_phenotype_matched_count} of ${clingenSubs.length}`
          : (hpoSubmitted
              ? 'No phenotype-matched curation'
              : 'All curations (no proband HPO)');
        const headline = displayBest
          ? `<div class="gc-validity-headline"><span class="${tierPillClass(displayBest)}">${esc(displayBest)}</span>${flag}<span class="gc-validity-meta">${esc(metaText)}</span></div>`
          : `<div class="gc-validity-empty">No match for proband phenotype</div>`;
        const renderRow = (s, extra) => {
          const mark = s.hpo_match
            ? ' <span class="gc-pheno-tick" title="Matches proband phenotype">✓</span>'
            : '';
          const cls = `gc-validity-list-row${extra ? ' gc-validity-extra' : ''}${s.hpo_match ? ' is-match' : ''}`;
          return `<div class="${cls}"><span class="${tierPillClass(s.classification)}">${esc(s.classification)}</span><span>${esc(s.disease)}${s.moi ? ` (${esc(s.moi)})` : ''}${mark}</span></div>`;
        };

        const sortedClingen = [...clingenSubs].sort(
          (a, b) => Number(Boolean(b.hpo_match)) - Number(Boolean(a.hpo_match))
        );
        const firstRows = sortedClingen.slice(0, 3).map(s => renderRow(s, false));
        const extraRows = sortedClingen.slice(3).map(s => renderRow(s, true));
        const moreBtn = extraRows.length
          ? `<button type="button" class="gc-validity-more" data-more-count="${extraRows.length}" onclick="toggleClingenMore(this)">+ ${extraRows.length} more curation(s)</button>`
          : '';
        clingenContent = `${headline}
          <div class="gc-validity-list">${firstRows.join('')}${extraRows.join('')}</div>
          ${moreBtn}`;
      } else {
        clingenContent = `<div class="gc-validity-empty">No GCEP curation — gene has not been formally curated by a ClinGen Gene Curation Expert Panel.</div>`;
      }

      let genccContent;
      if (gc && gc.found && gc.best_classification) {
        const hpoSubmitted = Boolean(gc.submitted_hpo);
        const phenoMatched = (gc.phenotype_matched_count || 0) > 0;
        const displayBest = phenoMatched
          ? gc.phenotype_matched_best_classification
          : (hpoSubmitted ? null : gc.best_classification);
        const displayDisputed = phenoMatched
          ? gc.phenotype_matched_has_disputed_or_refuted
          : false;
        const flag = displayDisputed
          ? ' <span title="Disputed/Refuted submitter(s) present in matched subset" style="color:#c8341a">⚠</span>'
          : '';
        const metaText = phenoMatched
          ? `Phenotype match · ${gc.phenotype_matched_count} of ${gc.submission_count}`
          : (hpoSubmitted
              ? 'No phenotype-matched submission'
              : 'All submitters (no proband HPO)');

        const scopeSubs = phenoMatched
          ? (gc.submissions || []).filter(s => s.hpo_match)
          : (gc.submissions || []);
        const tierRank = [
          'Definitive', 'Strong', 'Moderate', 'Limited',
          'Animal Model Only', 'Disputed Evidence', 'Refuted Evidence',
        ];
        const tierOrder = new Map(tierRank.map((t, i) => [t, i]));
        const groups = new Map();
        scopeSubs.forEach(s => {
          const tier = (s && s.classification || '').trim();
          const submitter = (s && s.submitter || '').trim();
          const disease = (s && s.disease || '').trim();
          if (!tier || !submitter) return;
          const key = `${tier} ${disease}`;
          if (!groups.has(key)) groups.set(key, { tier, disease, submitters: new Set() });
          groups.get(key).submitters.add(submitter);
        });
        let breakdownHtml = '';
        if (groups.size) {
          const entries = Array.from(groups.values()).sort((a, b) => {
            const ar = tierOrder.has(a.tier) ? tierOrder.get(a.tier) : 99;
            const br = tierOrder.has(b.tier) ? tierOrder.get(b.tier) : 99;
            if (ar !== br) return ar - br;
            return a.disease.localeCompare(b.disease);
          });
          const rows = entries.map(({ tier, disease, submitters }) => {
            const names = Array.from(submitters);
            const shown = names.slice(0, 3).map(esc).join(', ');
            const more = names.length - 3;
            const suffix = more > 0 ? `, + ${more} more` : '';
            const diseaseHtml = disease ? `${esc(disease)} · ` : '';
            return `<div class="gc-validity-list-row"><span class="${tierPillClass(tier)}">${esc(tier)} (${names.length})</span><span>${diseaseHtml}${shown}${suffix}</span></div>`;
          });
          breakdownHtml = `<div class="gc-validity-list">${rows.join('')}</div>`;
        } else {
          breakdownHtml = `<div class="gc-validity-list"><div class="gc-validity-list-row"><span>${gc.submission_count || 0} total submissions across all curators.</span></div></div>`;
        }
        const headline = displayBest
          ? `<div class="gc-validity-headline"><span class="${tierPillClass(displayBest)}">${esc(displayBest)}</span>${flag}<span class="gc-validity-meta">${esc(metaText)}</span></div>`
          : `<div class="gc-validity-empty">No match for proband phenotype</div>`;
        genccContent = `${headline}${breakdownHtml}`;
      } else {
        genccContent = `<div class="gc-validity-empty">No GenCC submissions for this gene.</div>`;
      }

      const clingenSourceUrl = gene
        ? `https://search.clinicalgenome.org/kb/genes/${encodeURIComponent(gene)}`
        : '';
      const genccSourceUrl = (gc && gc.url) || (gene
        ? `https://search.thegencc.org/genes/HGNC?q=${encodeURIComponent(gene)}`
        : '');
      const clingenLabel = clingenSourceUrl
        ? `<a href="${clingenSourceUrl}" target="_blank" rel="noopener">ClinGen ↗</a>`
        : 'ClinGen';
      const genccLabel = genccSourceUrl
        ? `<a href="${genccSourceUrl}" target="_blank" rel="noopener">GenCC ↗</a>`
        : 'GenCC';

      const pa = ev.panelapp || {};
      const panelConfPill = c => {
        const s = String(c || '').toLowerCase();
        if (s === 'green') return 'gc-pill gc-pill-green';
        if (s === 'amber') return 'gc-pill gc-pill-amber';
        if (s === 'red')   return 'gc-pill gc-pill-red';
        return 'gc-pill gc-pill-limited';
      };
      let panelappContent;
      if (!pa.ok) {
        panelappContent = `<div class="gc-validity-empty"><em>Lookup failed: ${esc(pa.error || 'unknown')}</em></div>`;
      } else if (!pa.on_cardiovascular_panel) {
        panelappContent = `<div class="gc-validity-empty">${esc(pa.gene || 'Gene')} not on any cardiovascular panel</div>`;
      } else {
        const panels = pa.panels_found || [];
        const hpoSubmitted = (pa.submitted_hpo || []).length > 0;

        const confRank = { green: 3, amber: 2, red: 1 };
        let best = null, bestRank = 0;
        for (const p of panels) {
          const r = confRank[p.confidence] || 0;
          if (r > bestRank) { bestRank = r; best = p.confidence; }
        }
        const anyPp4 = panels.some(p => p.contributes_to_pp4);
        const meta = `${panels.length} cardiovascular panel${panels.length === 1 ? '' : 's'}`
          + (best === 'green' ? (anyPp4 ? ' · counts toward PP4' : ' · PP4 needs phenotype match') : '');
        const headline = best
          ? `<div class="gc-validity-headline"><span class="${panelConfPill(best)}">${sCase(best)}</span><span class="gc-validity-meta">${esc(meta)}</span></div>`
          : '';
        const renderPanel = m => {
          const ppLabel = (m.confidence === 'green')
            ? (m.contributes_to_pp4
                ? ` <span class="ev-panel-flag ev-cls ev-cls-ben">✓ counts toward PP4</span>`
                : ` <span class="ev-panel-flag ev-cls ev-cls-moderate">⚠ ${hpoSubmitted ? 'phenotype does not match this panel' : 'enter a phenotype/HPO term for PP4'}</span>`)
            : '';
          const moi = m.moi ? ` <span class="ev-dim">· MOI ${esc(m.moi)}</span>` : '';
          return `<div class="gc-validity-list-row">`
            + `<span class="${panelConfPill(m.confidence)}">${sCase(m.confidence)} (${esc(m.confidence_level)}/3)</span>`
            + `<span><a href="${m.panel_url}" target="_blank" rel="noopener">${esc(m.panel_name)}</a>${moi}${ppLabel}</span>`
            + `</div>`;
        };
        const list = `<div class="gc-validity-list">${panels.map(renderPanel).join('')}</div>`;
        const unrec = pa.unrecognised_hpo_tokens || [];
        const unrecLine = unrec.length
          ? `<div class="ev-note" style="margin-top:6px;background:var(--amber-light);border:1px solid var(--amber-border);border-left:3px solid var(--amber);color:var(--amber);padding:6px 10px;border-radius:4px">⚠ ${unrec.length} phenotype term${unrec.length === 1 ? '' : 's'} not recognised: ${unrec.map(t => `<code>${esc(t)}</code>`).join(', ')} — did not contribute to panel matching.</div>`
          : '';
        panelappContent = `${headline}${list}${unrecLine}`;
      }
      const panelappSourceUrl = gene
        ? `https://panelapp-aus.org/panels/entities/${encodeURIComponent(gene)}`
        : '';
      const panelappLabel = panelappSourceUrl
        ? `<a href="${panelappSourceUrl}" target="_blank" rel="noopener">PanelApp ↗</a>`
        : 'PanelApp';

      const cd = ev.chdgene || {};
      let chdgeneContent;
      if (!cd.ok) {
        chdgeneContent = `<div class="gc-validity-empty"><em>Lookup failed: ${esc(cd.error || 'unknown')}</em></div>`;
      } else if (!cd.listed) {
        chdgeneContent = `<div class="gc-validity-empty">${esc(cd.gene || 'Gene')} not in CHDgene curated list</div>`;
      } else {
        const chd = (cd.chd_classification || []).join(', ') || '—';
        const inh = (cd.inheritance || []).join(', ') || '—';
        const extra = cd.extra_cardiac_phenotype
          ? ' <span class="ev-dim">· extra-cardiac features reported</span>' : '';
        const headline = `<div class="gc-validity-headline"><span class="gc-pill gc-pill-definitive">Listed</span><span class="gc-validity-meta">Established CHD gene</span></div>`;
        chdgeneContent = `${headline}
          <div class="gc-validity-list">
            <div class="gc-validity-list-row"><span>CHD subtypes: ${esc(chd)}</span></div>
            <div class="gc-validity-list-row"><span>Inheritance: ${esc(inh)}${extra}</span></div>
          </div>`;
      }
      const chdgeneSourceUrl = (cd && cd.url) ? cd.url : '';
      const chdgeneLabel = chdgeneSourceUrl
        ? `<a href="${chdgeneSourceUrl}" target="_blank" rel="noopener">CHDgene ↗</a>`
        : 'CHDgene';

      const audit = ev.phenotype_audit || {};

      const interp = audit.interpretations || [];
      const interpLine = interp.length
        ? `<div class="gc-validity-readas">Read <b>${interp.length === 1 ? 'term' : 'terms'}</b> as: `
          + interp.map(t => `<code>${esc(t.token)}</code> → ${esc((t.parts || []).filter(Boolean).join(' + ') || t.label || t.token)}`).join('; ')
          + `</div>`
        : '';
      const unmatchedTerms = audit.unmatched || [];
      const nomatchLine = unmatchedTerms.length
        ? `<div class="gc-validity-nomatch">⚠ ${unmatchedTerms.length} submitted phenotype term${unmatchedTerms.length === 1 ? '' : 's'}`
          + ` match${unmatchedTerms.length === 1 ? 'es' : ''} none of ${gene ? esc(gene) + '’s' : 'this gene’s'} curated gene-disease records: `
          + unmatchedTerms.map(t => `<b>${esc(t.label || t.token)}</b>${(t.label && t.label !== t.token) ? ` <span class="ev-dim">(${esc(t.token)})</span>` : ''}`).join(', ')
          + ` — ${unmatchedTerms.length === 1 ? 'it did' : 'they did'} not contribute to phenotype matching or PP4.</div>`
        : '';
      const section1 = `<div class="gc-collapsible" id="gc-coll-validity" style="--gc-accent: var(--blue);">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Gene-disease validity</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            <div class="gc-validity-rows">
              <div class="gc-validity-row">
                <div class="gc-validity-source">${clingenLabel}</div>
                <div class="gc-validity-content">${clingenContent}</div>
              </div>
              <div class="gc-validity-row">
                <div class="gc-validity-source">${genccLabel}</div>
                <div class="gc-validity-content">${genccContent}</div>
              </div>
              <div class="gc-validity-row">
                <div class="gc-validity-source">${panelappLabel}</div>
                <div class="gc-validity-content">${panelappContent}</div>
              </div>
              <div class="gc-validity-row">
                <div class="gc-validity-source">${chdgeneLabel}</div>
                <div class="gc-validity-content">${chdgeneContent}</div>
              </div>
            </div>
            ${interpLine}${nomatchLine}
          </div>
        </div>
      </div>`;

      const _pm5Num = ev.clinvar_pm5_candidates || {};
      const _ssNum = ev.same_site || {};
      const _probandResClinvar = (_pm5Num.matched_protein_position != null)
        ? _pm5Num.matched_protein_position
        : (_ssNum.matched_residue != null ? _ssNum.matched_residue : null);
      const _probandResModel = (_pm5Num.protein_position != null)
        ? _pm5Num.protein_position
        : (_ssNum.display_residue != null ? _ssNum.display_residue : null);
      const _maneTxRow = ((ev.vep || {}).transcript_consequences_all || [])
        .find(r => r && r.is_mane_select) || null;
      const _clinvarNumberingTx = _maneTxRow ? (_maneTxRow.transcript_id || '') : '';
      const _codonPositions = (ev.same_site && ev.same_site.codon_positions) || null;
      const landscape = ev.clinvar_gene_landscape || {};
      const total = Number(landscape.total_classified) || 0;
      let section2Body;
      if (!landscape.ok || total === 0) {
        section2Body = `<div class="gc-empty">${landscape.ok
          ? 'No classified ClinVar records indexed for this gene.'
          : 'ClinVar gene-landscape lookup unavailable.'}</div>`;
      } else {

        const keywordNote =
          landscape.condition_keywords && landscape.condition_keywords.length
            && landscape.keyword_source !== 'gencc'
              ? `<div class="gc-bar-cap-note">Filtered to cardiac submissions.</div>`
              : '';

        const plpPhens = landscape.plp_top_phenotypes || [];
        const plpPhenBlock = plpPhens.length
          ? `<div class="gc-plp-phen-wrap">
              <div class="gc-plp-phen-label">P/LP linked to</div>
              <div class="gc-plp-phen-pills">${
                plpPhens.map(p =>
                  `<span class="gc-plp-phen-pill">${esc(p.phenotype)} (${p.count})</span>`
                ).join('')
              }</div>
            </div>`
          : '';

        const scopeCap = `<div class="gc-bar-cap-note" style="margin-bottom:8px">Gene-wide — the classification mix across <strong>all ${total} classified ClinVar record${total === 1 ? '' : 's'} for ${esc(gene)}</strong>, not this variant. Use the chips to filter by variant type.</div>`;

        const _vepBucket = _csqBucketFromVepTerm((ev.vep && ev.vep.most_severe_consequence) || '');
        const _tierByCsq = landscape.tier_counts_by_csq || {};
        const _exonData = ev.transcript_exons || null;

        const _probandExonRank = (function () {
          const m = /^(\d+)\s*\//.exec((ev.vep && ev.vep.exon) || '');
          return m ? parseInt(m[1], 10) : null;
        })();
        const lsData = {
          positions: Array.isArray(landscape.positions) ? landscape.positions : [],
          tierByCsq: _tierByCsq,
          tierCounts: landscape.tier_counts || {},
          gMin: landscape.g_min, gMax: landscape.g_max,
          strand: (ev.vep && ev.vep.strand) || 1,
          probandGpos: (ev.vep && ev.vep.start != null) ? ev.vep.start : null,
          probandExon: (ev.vep && ev.vep.exon) || '',
          probandExonRank: _probandExonRank,
          exons: (_exonData && _exonData.ok && Array.isArray(_exonData.exons)) ? _exonData.exons : null,
          exonStrand: _exonData && _exonData.strand,
          gene: gene,
          truncated: !!landscape.positions_truncated,

          modelLabel: _transcriptModelLabel(ev, _exonData && _exonData.transcript_id),

          clinvarTxId: _clinvarNumberingTx,
          probandResidueClinvar: _probandResClinvar,
          probandResidueModel: _probandResModel,
          codonPositions: _codonPositions,
        };

        const _codonPos = _codonPositions;
        const _atCodon = _codonPos
          ? (lsData.positions || []).filter(pt => _codonPos.indexOf(pt.gpos) >= 0)
          : [];
        const _vusAtCodon = _atCodon.some(pt => pt.tier === 'VUS');
        let defaultCsq = _tierByCsq[_vepBucket] ? _vepBucket : 'all';

        if (defaultCsq !== 'all' && _atCodon.some(pt => pt.csq !== defaultCsq)) {
          defaultCsq = 'all';
        }
        window.__hvLandscapeState = { data: lsData, csq: defaultCsq, showVus: _vusAtCodon };
        section2Body = `${scopeCap}
          ${_buildLandscapeChips(lsData, defaultCsq)}
          <div id="hv-ls-chart">${_buildLandscapeBody(lsData, defaultCsq, _vusAtCodon)}</div>
          ${keywordNote}
          ${plpPhenBlock}`;
      }

      const pm5 = ev.clinvar_pm5_candidates || {};
      let sameResidueBlock = '';
      {
        const sr = (ev.same_site && ev.same_site.clinvar_same_residue) || {};
        const recs = (sr.ok && Array.isArray(sr.records)) ? sr.records : [];
        const stars = n => '★'.repeat(n || 0) + '☆'.repeat(Math.max(0, 4 - (n || 0)));

        const dispPos = (sr.display_position != null) ? sr.display_position
                      : (pm5.protein_position != null) ? pm5.protein_position : null;
        const matchPos = (sr.matched_position != null) ? sr.matched_position
                       : (pm5.matched_protein_position != null) ? pm5.matched_protein_position : null;
        const differs = (sr.numbering_differs != null) ? !!sr.numbering_differs
                      : !!pm5.numbering_differs;
        const cap = (dispPos == null)
          ? 'this amino acid'
          : (differs && matchPos != null
              ? `residue ${esc(dispPos)} `
                + `<span class="ev-dim">(= residue ${esc(matchPos)} on ClinVar's MANE transcript)</span>`
              : `residue ${esc(dispPos)}`);

        if (dispPos != null) {
          let rows, foot;
          if (recs.length) {
            rows = recs.map(r => {
              const flag = r.pm5_eligible
                ? '<span class="ev-cls ev-cls-path">PM5-eligible</span>'
                : '<span class="ev-cls ev-dim">not PM5 evidence</span>';

              const reason = r.pm5_eligible ? '' : (r.pm5_ineligible_reason || '');
              const cls = (r.clinical_significance
                           && reason.indexOf(r.clinical_significance) < 0)
                ? ` <span class="ev-dim">${esc(r.clinical_significance)}</span>` : '';

              const why = reason
                ? ` <span class="ev-dim">: ${esc(reason)}</span>`
                : '';
              return `<li><code>${esc(r.name)}</code>${cls} `
                + `<span title="${esc(r.review_status || '')}">${stars(r.stars)}</span> `
                + `${flag}${why}</li>`;
            }).join('');
            foot = `PM5 requires a DIFFERENT MISSENSE change at this residue that is
              itself established Pathogenic/Likely pathogenic. Rows marked
              "not PM5 evidence" must not be cited for PM5 or PS1. Whether PM5
              actually applies is shown on the PM5 row of the Criteria tab.`;
          } else {

            const fb = []
              .concat(Array.isArray(pm5.candidates) ? pm5.candidates : [])
              .concat(Array.isArray(pm5.ps1_candidates) ? pm5.ps1_candidates : []);
            if (fb.length) {
              rows = fb.map(c => `<li><code>${esc(c.name)}</code>`
                + ` <span class="ev-dim">${esc(c.clinical_significance || c.tier || '')}</span> `
                + `<span title="${esc(c.review_status || '')}">${stars(c.stars)}</span></li>`).join('');
              foot = `Per-record PM5 eligibility was not available for this
                curation, so no row here is marked as qualifying. Whether PM5
                applies is shown on the PM5 row of the Criteria tab.`;
            } else {
              rows = '';
              foot = `No other ClinVar record is recorded at this residue. This is
                an absence of records, not an absence of evidence for the
                variant itself.`;
            }
          }
          sameResidueBlock = `<div class="gc-pm5-callout">
            <div class="gc-pm5-callout-title">ClinVar records at ${cap}</div>
            <div class="hv-model-note">${foot}</div>
            ${rows ? `<ul class="gc-pm5-callout-list">${rows}</ul>` : ''}
          </div>`;
        }
      }

      const clinvarBrowseLink = gene
        ? `<div class="gc-source-link"><a href="https://www.ncbi.nlm.nih.gov/clinvar/?gene=${encodeURIComponent(gene)}&term=${encodeURIComponent(`"${gene}"[GENE]`)}" target="_blank" rel="noopener">View all ${esc(gene)} variants on ClinVar ↗</a></div>`
        : '';
      const section2 = `<div class="gc-collapsible" id="gc-coll-clinvar" style="--gc-accent: #c26514;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">ClinVar variant landscape</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            ${section2Body}
            ${sameResidueBlock}
            ${clinvarBrowseLink}
          </div>
        </div>
      </div>`;

      const moiNormalize = (raw) => {
        const t = String(raw || '').trim();
        if (!t) return null;
        const lc = t.toLowerCase();
        const map = {
          ad: 'Autosomal dominant',
          ar: 'Autosomal recessive',
          xl: 'X-linked',
          xld: 'X-linked dominant',
          xlr: 'X-linked recessive',
          yl: 'Y-linked',
          mt: 'Mitochondrial',
          'semi-d': 'Semidominant',
          'sem-d': 'Semidominant',
          mu: 'Multifactorial',
          dig: 'Digenic',
        };
        if (map[lc]) return map[lc];

        return t.charAt(0).toUpperCase() + t.slice(1).toLowerCase();
      };
      const moiSet = new Set();

      (((ev.chdgene || {}).inheritance) || []).forEach(m => {
        const n = moiNormalize(m);
        if (n) moiSet.add(n);
      });
      ((gc && gc.submissions) || []).forEach(s => {
        const n = moiNormalize(s.moi);
        if (n) moiSet.add(n);
      });

      const moiList = Array.from(moiSet).sort();
      const moiText = moiList.length
        ? moiList.join(', ')
        : 'Not specified by curated sources.';

      const constraint = ((ev.gnomad && ev.gnomad.gene && ev.gnomad.gene.gnomad_constraint) || {});
      const pLI = constraint.pLI;
      const loeuf = constraint.oe_lof_upper;
      const misZ = constraint.mis_z;

      let constraintType = null;
      let constraintSub = '';
      let constraintNote = '';
      if (pLI != null && pLI >= 0.9) {
        constraintType = 'LoF intolerant';
        const bits = [`pLI ${Number(pLI).toFixed(2)}`];
        if (loeuf != null) bits.push(`LOEUF ${Number(loeuf).toFixed(2)}`);
        constraintSub = bits.join(' · ');
      } else if (misZ != null && misZ >= 3 && (pLI == null || pLI < 0.5)) {
        constraintType = 'Missense intolerant';
        const bits = [
          `mis_z ${Number(misZ).toFixed(2)}`,
          `pLI ${pLI == null ? 'n/a' : Number(pLI).toFixed(2)}`,
        ];
        constraintSub = bits.join(' · ');
      } else if (pLI != null || loeuf != null || misZ != null) {
        const bits = [];
        if (pLI != null)  bits.push(`pLI ${Number(pLI).toFixed(2)}`);
        if (loeuf != null) bits.push(`LOEUF ${Number(loeuf).toFixed(2)}`);
        if (misZ != null) bits.push(`mis_z ${Number(misZ).toFixed(2)}`);
        constraintNote = `Constraint values do not strongly favour a single mechanism (${bits.join(', ')}).`;
      } else {
        constraintNote = 'gnomAD constraint not available.';
      }

      const inhTile = `<div class="gc-mech-tile${constraintType ? '' : ' gc-mech-tile-full'}">
        <div class="gc-mech-tile-label">Inheritance</div>
        <div class="gc-mech-tile-value">${esc(moiText)}</div>
      </div>`;
      const conTile = constraintType
        ? `<div class="gc-mech-tile">
            <div class="gc-mech-tile-label">Constraint</div>
            <div class="gc-mech-tile-value">${esc(constraintType)}</div>
            <div class="gc-mech-tile-sub">${esc(constraintSub)}</div>
          </div>`
        : `<div class="gc-mech-tile gc-mech-tile-full">
            <div class="gc-mech-tile-label">Constraint</div>
            <div class="gc-mech-tile-note">${esc(constraintNote)}</div>
          </div>`;

      const gnomadGeneUrl = gene
        ? `https://gnomad.broadinstitute.org/gene/${encodeURIComponent(gene)}?dataset=gnomad_r4`
        : '';
      const constraintSourceLink = gnomadGeneUrl
        ? `<div class="gc-source-link">Constraint values: <a href="${gnomadGeneUrl}" target="_blank" rel="noopener">gnomAD gene page ↗</a></div>`
        : '';

      const mech = ev.gene_mechanism || null;

      const MECH_SHORT = {
        haploinsufficiency: 'Loss-of-function',
        recessive_lof: 'Loss-of-function (recessive)',
        dominant_negative: 'Dominant-negative',
        gof_or_dn: 'Gain-of-function / DN',
        mixed: 'Mixed mechanism',
        undetermined: 'Not established',
      };
      const mechChipStyle = (cls) => {
        if (cls === 'haploinsufficiency' || cls === 'recessive_lof')
          return 'background:var(--blue-light);color:var(--blue);border:1px solid var(--blue-border)';
        if (cls === 'gof_or_dn' || cls === 'dominant_negative')
          return 'background:#f6e9ea;color:var(--maroon);border:1px solid #e3c4c6';
        if (cls === 'mixed')
          return 'background:var(--amber-light);color:var(--amber);border:1px solid var(--amber-border)';
        return 'background:#f1f1ee;color:var(--muted);border:1px solid var(--line)';
      };
      const CONF_WORD = { established: 'curated evidence', emerging: 'emerging evidence', curated: 'HeartVar-curated', undetermined: '' };
      let mechBlockHtml = '';
      if (mech && mech.mechanism) {
        const chip = `<span class="gc-mech-chip" style="display:inline-block;padding:3px 10px;border-radius:999px;font-size:12.5px;font-weight:600;${mechChipStyle(mech.mechanism)}">${esc(mech.label || MECH_SHORT[mech.mechanism] || mech.mechanism)}</span>`;
        const conf = CONF_WORD[mech.confidence] ? `<span style="font-size:12px;color:var(--muted)">${esc(CONF_WORD[mech.confidence])}</span>` : '';

        const sentence = (mech.mechanism === 'mixed' || mech.mechanism === 'undetermined') && mech.summary
          ? `<div style="font-size:13px;line-height:1.5;color:var(--ink);margin-top:7px">${esc(mech.summary)}</div>`
          : '';
        const srcRows = (mech.sources || []).map(s => {
          const nm = s.url
            ? `<a href="${esc(s.url)}" target="_blank" rel="noopener">${esc(s.name)} ↗</a>`
            : `<strong>${esc(s.name)}</strong>`;
          return `<div style="font-size:12px;line-height:1.5;color:var(--muted)">${nm}${s.detail ? ' — ' + esc(s.detail) : ''}</div>`;
        }).join('');
        const srcBlock = srcRows
          ? `<div style="margin-top:8px;border-top:1px dashed var(--line);padding-top:7px">${srcRows}</div>`
          : '';
        mechBlockHtml = `<div class="gc-mech-established" style="background:var(--card);border:1px solid var(--line);border-radius:8px;padding:11px 13px;margin-bottom:12px">
          <div style="text-transform:uppercase;font-size:11px;letter-spacing:.08em;color:var(--muted);font-weight:600;margin-bottom:7px">Established disease mechanism</div>
          <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap">${chip}${conf}</div>
          ${sentence}
          ${srcBlock}
        </div>`;
      }

      const mechFlag = ev.mechanism_flag || null;
      const isWarn = !!(mechFlag && mechFlag.severity === 'warning');
      const flagStyle = isWarn
        ? 'background:var(--amber-light);border:1px solid var(--amber-border);border-left:3px solid var(--amber);color:var(--amber)'
        : 'background:var(--blue-light);border:1px solid var(--blue-border);border-left:3px solid var(--blue);color:var(--blue)';
      const mechFlagHtml = (mechFlag && mechFlag.detail)
        ? `<div class="gc-mech-flag" role="note" style="margin-bottom:12px;${flagStyle};padding:9px 12px;border-radius:6px;font-size:13px;line-height:1.5"><strong>${isWarn ? '⚠' : 'ℹ'} ${esc(mechFlag.title || 'Mechanism note')}.</strong> ${esc(mechFlag.detail)}</div>`
        : '';

      const mechSummary = isWarn
        ? '<span class="gc-coll-summary" style="color:var(--amber);font-weight:600">⚠ re-assess variant type</span>'
        : '<span class="gc-coll-summary"></span>';
      const section3 = `<div class="gc-collapsible${isWarn ? ' open' : ''}" id="gc-coll-mech" style="--gc-accent: var(--amber);">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Disease mechanism</span>
          ${mechSummary}
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            ${mechFlagHtml}
            ${mechBlockHtml}
            <div class="gc-mech-grid">
              ${inhTile}
              ${conTile}
            </div>
            ${claudeNote('Variant fit', geneContext.mechanism_consistency)}
            ${constraintSourceLink}
          </div>
        </div>
      </div>`;

      const gtex = ev.gtex || {};
      const tissues = (gtex && gtex.tissues) || [];
      const wantedTissues = ['Heart_Left_Ventricle', 'Heart_Atrial_Appendage', 'Artery_Aorta', 'Artery_Coronary'];

      const exprData = wantedTissues.map(key => {
        const t = tissues.find(x =>
          x && (x.tissue === key || x.tissueSiteDetailId === key || (x.tissue_label && x.tissue_label.replace(/\s+/g, '_') === key))
        );
        const tpm = t && (t.median_tpm != null ? t.median_tpm : (t.median != null ? t.median : null));
        return { key, label: key.replace(/_/g, ' '), tpm };
      });
      const exprCards = exprData.map(({ key, label, tpm }) => {
        if (tpm == null) {
          return `<div class="gc-expr-card">
            <div class="gc-expr-tissue">${esc(label)}</div>
            <div class="gc-expr-tpm gc-card-dim">—</div>
            <div class="gc-expr-band gc-band-low">no data</div>
          </div>`;
        }
        let band, bandText;
        if (tpm >= 10)      { band = 'gc-band-high';     bandText = 'High expression'; }
        else if (tpm >= 1)  { band = 'gc-band-moderate'; bandText = 'Moderate expression'; }
        else                { band = 'gc-band-low';      bandText = 'Low expression'; }
        return `<div class="gc-expr-card">
          <div class="gc-expr-tissue">${esc(label)}</div>
          <div class="gc-expr-tpm">${Number(tpm).toFixed(2)}<span class="gc-expr-tpm-unit">TPM</span></div>
          <div class="gc-expr-band ${band}">${bandText}</div>
        </div>`;
      }).join('');

      const gtexSourceLink = gtex.url
        ? `<div class="gc-source-link"><a href="${gtex.url}" target="_blank" rel="noopener">View gene on GTEx Portal ↗</a></div>`
        : '';
      const section4 = `<div class="gc-collapsible" id="gc-coll-expr" style="--gc-accent: #0F6E56;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Cardiac &amp; vascular expression (GTEx)</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            <div class="gc-expr-grid">${exprCards}</div>
            ${gtexSourceLink}
          </div>
        </div>
      </div>`;

      const fh = ev.fetal_heart || {};
      let section4bBody;
      if (!fh.ok) {
        section4bBody = `<div class="gc-empty">${
          fh.error
            ? esc('Fetal heart unavailable: ' + fh.error)
            : 'Fetal heart data not loaded.'}</div>`;
      } else if (!fh.found) {
        section4bBody = `<div class="gc-empty">${esc(fh.gene || 'Gene')} not detected in the Farah 2024 fetal heart dataset.</div>`;
      } else {
        const stages = fh.stages || [];
        const cells = fh.cell_types || [];

        const byType = new Map();
        const typeOrder = [];
        for (const c of cells) {
          if (!byType.has(c.cell_type)) {
            byType.set(c.cell_type, { label: c.cell_type_label || c.cell_type, n: c.n_cells, by_stage: {} });
            typeOrder.push(c.cell_type);
          }
          byType.get(c.cell_type).by_stage[c.stage] = c;
        }
        const headerCells = stages.map(s => `<th>${esc(s)}</th>`).join('');
        const rows = typeOrder.map(t => {
          const r = byType.get(t);
          const cellsHtml = stages.map(s => {
            const c = r.by_stage[s];
            if (!c) return `<td class="gc-fetal-cell fb-absent">—</td>`;
            const meanFmt = Number(c.mean_expr).toFixed(2);
            const pctFmt = Number(c.pct_expressing).toFixed(0);
            const bandClass = 'fb-' + (c.band || 'low');
            return `<td class="gc-fetal-cell ${bandClass}" title="${esc(t)} · ${esc(s)} · mean log1p(CPM) ${meanFmt} · ${pctFmt}% of ${c.n_cells} cells">
              <span class="gc-fetal-mean">${meanFmt}</span>
              <span class="gc-fetal-pct">${pctFmt}%</span>
            </td>`;
          }).join('');

          const rowNRaw = Math.max(...Object.values(r.by_stage).map(x => x ? x.n_cells : 0));
          const rowNFmt = isFinite(rowNRaw) && rowNRaw > 0 ? rowNRaw.toLocaleString() : '—';
          return `<tr>
            <td class="gc-fetal-rowlabel">${esc(r.label)} <span class="gc-fetal-rowlabel-sub">${esc(t)} · ${rowNFmt} cells max/stage</span></td>
            ${cellsHtml}
          </tr>`;
        }).join('');
        const legend = ['absent','low','moderate','high','very_high'].map(b => {
          const lbl = b.replace('_',' ');
          return `<span><span class="fl-swatch fb-${b}"></span>${esc(lbl)}</span>`;
        }).join('');
        const sourceUrl = fh.url || 'https://cells.ucsc.edu/?ds=hoc';
        section4bBody = `<div class="gc-fetal-wrap">
          <table class="gc-fetal-table">
            <thead><tr><th class="gc-fetal-rowlabel">Cell type</th>${headerCells}</tr></thead>
            <tbody>${rows}</tbody>
          </table>
        </div>
        <div class="gc-fetal-legend"><strong>Band:</strong> ${legend}</div>
        <div class="gc-fetal-footnote">
          Each cell shows mean log1p(CPM) and the % of cells expressing the gene at that stage.
          Bands are dataset-relative quantiles of non-zero expression.
          PCW = post-conceptional weeks (measured from conception, ≈ clinical gestational age − 2);
          9–15 PCW spans first to early second trimester, when chamber septation completes and the
          ventricular wall compacts.
          Source: <a href="${esc(sourceUrl)}" target="_blank" rel="noopener">Farah 2024 · "Heart of Cells" (UCSC Cell Browser) ↗</a>.
        </div>`;
      }
      const section4b = `<div class="gc-collapsible" id="gc-coll-fetal" style="--gc-accent: #B7335A;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Fetal cardiac expression (Farah 2024)</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">${section4bBody}</div>
        </div>
      </div>`;

      const ot = ev.opentargets_evidence || {};
      let section5Body;
      if (!ot.ok) {
        section5Body = `<div class="gc-empty">${
          ot.error
            ? esc('Open Targets unavailable: ' + ot.error)
            : 'Open Targets data not loaded.'}</div>`;
      } else {
        const matched = ot.matched_disease || null;
        const topDiseases = ot.top_diseases || [];

        const rows = matched
          ? [{ name: matched.name || '?', score: ot.overall_association_score, matched: true }]
          : topDiseases.map(d => ({ name: d.name || '?', score: d.score, matched: false }));
        if (!rows.length) {
          section5Body = `<div class="gc-empty">No association data.</div>`;
        } else {
          const items = rows.map(r => {
            const num = r.score == null ? null : Number(r.score);
            const pct = num == null ? 0 : Math.max(0, Math.min(1, num)) * 100;
            const sc = num == null ? 'n/a' : num.toFixed(2);
            const pill = r.matched ? `<span class="gc-ot-hpo-pill">HPO match</span>` : '';
            return `<li class="gc-ot-list-row">
              <span class="gc-ot-list-name">${esc(r.name)}${pill}</span>
              <span class="gc-ot-list-bar"><span class="gc-ot-list-fill" style="width:${pct.toFixed(1)}%"></span></span>
              <span class="gc-ot-list-score">${sc}</span>
            </li>`;
          }).join('');
          section5Body = `<ol class="gc-ot-list">${items}</ol>`;
        }
      }

      const otEnsembl = ot.ensembl_id || (ev.vep && ev.vep.gene_id) || '';
      const otSourceLink = otEnsembl
        ? `<div class="gc-source-link"><a href="https://platform.opentargets.org/target/${encodeURIComponent(otEnsembl)}" target="_blank" rel="noopener">View gene on Open Targets ↗</a></div>`
        : '';
      const section5 = `<div class="gc-collapsible" id="gc-coll-ot" style="--gc-accent: #6b6a14;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Open Targets</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            ${section5Body}
            ${claudeNote('Open Targets read', geneContext.opentargets_interpretation)}
            ${otSourceLink}
          </div>
        </div>
      </div>`;

      const glit = ev.gene_literature || {};
      const litPapers = (glit.ok ? (glit.papers || []) : []).slice(0, 5);
      let section6Body;
      if (!glit.ok) {
        section6Body = `<div class="gc-empty">${
          glit.error
            ? esc('Gene literature unavailable: ' + glit.error)
            : 'Gene literature not loaded.'}</div>`;
      } else if (!litPapers.length) {
        section6Body = `<div class="gc-empty">No phenotype-matched gene-level papers found.</div>`;
      } else {

        const items = litPapers.map(p => {
          const pmid = String(p.pmid || '');
          const url = pmid
            ? `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(pmid)}/`
            : '';
          const titleText = p.title || '(untitled)';
          const titleHtml = `<span class="gc-lit-title">${esc(titleText)}</span>`;
          const metaParts = [p.first_author, p.journal, p.year].filter(Boolean).map(esc);
          const metaHtml = metaParts.length
            ? `<span class="gc-lit-meta"> · ${metaParts.join(' · ')}</span>`
            : '';
          const pmidLink = pmid
            ? `<a class="gc-lit-pmid" href="${url}" target="_blank" rel="noopener">PMID ${esc(pmid)} ↗</a>`
            : '';
          const titleAttr = p.title ? ` title="${esc(p.title)}"` : '';
          return `<li class="gc-lit-item"${titleAttr}>
            <span class="gc-lit-line">${titleHtml}${metaHtml}</span>
            ${pmidLink}
          </li>`;
        }).join('');
        section6Body = `<ul class="gc-lit-list">${items}</ul>`;
      }

      const geneTokenLit = (ev.vep && ev.vep.gene_symbol) || glit.gene || gene || '';
      const phenoPhrases = (glit.phenotype_keywords || []).filter(Boolean);
      const diseaseLabel = phenoPhrases.length ? phenoPhrases.join(' / ') : 'Cardiac';
      const geneCHDSearchURLLit = glit.search_url || 'https://pubmed.ncbi.nlm.nih.gov/';
      const geneSearchLink = geneTokenLit
        ? `<div class="ev-search-links">
             <a class="ev-search-link" href="${geneCHDSearchURLLit}" target="_blank" rel="noopener">Search PubMed for ${esc(geneTokenLit)} in ${esc(diseaseLabel)} ↗</a>
           </div>`
        : '';
      const section6 = `<div class="gc-collapsible" id="gc-coll-lit" style="--gc-accent: #534AB7;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Gene literature</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            ${claudeNote('Synthesis', geneContext.gene_literature_summary)}
            ${section6Body}
            ${geneSearchLink}
          </div>
        </div>
      </div>`;

      const mgi = ev.mgi || {};
      let section9Body;
      if (!mgi.ok) {
        section9Body = `<div class="gc-empty"><em>Lookup failed: ${esc(mgi.error || 'unknown')}</em></div>`;
      } else if (!mgi.found) {
        section9Body = `<div class="gc-empty"><em>${esc(mgi.reason || 'No mouse orthologue indexed')}</em></div>`;
      } else {
        const cardiac = mgi.cardiac_phenotypes || [];
        const nC = cardiac.length;
        const orthoValue = `<a href="${mgi.url}" target="_blank" rel="noopener" class="ev-mono">${esc(mgi.mouse_symbol)}</a> <span class="ev-dim">(${esc(mgi.mgi_id)} · ortholog confidence: ${esc(mgi.ortholog_confidence || 'n/a')})</span>`;
        const topGroup = kvGroup([
          kv('Mouse orthologue', orthoValue),
          kv('Phenotypic alleles', `<strong>${mgi.phenotype_count || 0}</strong>`),
        ]);
        let cardiacBlock;
        if (nC) {
          const tags = cardiac.slice(0, 12).map(t => `<span class="ev-tag">${esc(t)}</span>`).join('');
          cardiacBlock = `<div class="ev-kv-group"><div class="ev-kv-label" style="flex:none;margin-bottom:4px">Cardiac/cardiovascular phenotypes (${nC})</div><div class="ev-tag-row">${tags}</div></div>`;
        } else {
          cardiacBlock = `<div class="ev-note">No cardiac MP terms detected</div>`;
        }
        section9Body = topGroup + cardiacBlock;
      }
      const section9 = `<div class="gc-collapsible" id="gc-coll-mgi" style="--gc-accent: #7c3d8a;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">MGI (mouse models)</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">${section9Body}</div>
        </div>
      </div>`;

      const bg = ev.biogrid || {};
      let section10Body;
      if (!bg.ok) {
        section10Body = `<div class="gc-empty"><em>Lookup failed: ${esc(bg.error || 'unknown')}</em></div>`;
      } else if (!bg.total_interactions) {
        section10Body = `<div class="gc-empty"><em>No curated interactions reported</em></div>`;
      } else {
        const totalValue = `<strong>${bg.total_interactions}</strong> curated <span class="ev-dim">·</span> <strong>${bg.unique_partners}</strong> unique partners`;
        const topHtml = (bg.top_partners || []).map(p => {
          const cls = p.is_chdgene ? 'ev-partner-chd' : 'ev-partner';
          return `<span class="${cls}">${esc(p.symbol)}</span> <span class="ev-partner-count">(n=${p.publication_count})</span>`;
        }).join(', ') || '<span class="ev-dim">—</span>';
        section10Body = kvGroup([
          kv('Total interactions', totalValue),
          kv('Top partners', topHtml),
          kv('Link', `<a href="${bg.url}" target="_blank" rel="noopener" class="ev-dim">↗ open in BioGRID</a>`),
        ]);
      }
      const section10 = `<div class="gc-collapsible" id="gc-coll-biogrid" style="--gc-accent: #2a6e3a;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">BioGRID</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">${section10Body}</div>
        </div>
      </div>`;

      const md = ev.medgen || {};
      let section11Body;
      if (!md.ok) {
        section11Body = `<div class="gc-empty">MedGen lookup unavailable — ${md.url
          ? `<a href="${md.url}" target="_blank" rel="noopener">search MedGen ↗</a>`
          : '<em>no fallback link available</em>'}</div>`;
      } else if (!Array.isArray(md.conditions) || md.conditions.length === 0) {
        section11Body = `<div class="gc-empty">No MedGen conditions indexed for ${esc(md.gene || gene || 'this gene')}.</div>`;
      } else {
        const rows = md.conditions.map(c => {
          const heart = !!c.heart_related;
          const nameLink = c.medgen_url
            ? `<a href="${c.medgen_url}" target="_blank" rel="noopener">${esc(c.name)}</a>`
            : esc(c.name);
          const cardiacPill = heart
            ? `<span class="gc-medgen-cardiac-pill" title="Heart-related condition">Cardiac</span>`
            : '';
          const nameCell = `${nameLink}${cardiacPill}`;
          const mimCell = c.mim
            ? `<span class="gc-medgen-mim"><a href="${c.omim_url || `https://www.omim.org/entry/${esc(c.mim)}`}" target="_blank" rel="noopener" title="Open OMIM entry ${esc(c.mim)}">${esc(c.mim)} ↗</a></span>`
            : '<span class="ev-dim">—</span>';
          const moiCell = c.moi ? esc(c.moi) : '<span class="ev-dim">—</span>';
          return `<tr class="${heart ? 'is-cardiac' : ''}"><td>${nameCell}</td><td>${mimCell}</td><td>${moiCell}</td></tr>`;
        }).join('');
        const totalShown = md.conditions.length;
        const totalFound = Number(md.total_found) || totalShown;
        const overflowNote = totalFound > totalShown && md.url
          ? `<div class="gc-medgen-footnote">Showing ${totalShown} of ${totalFound} associated conditions. <a href="${md.url}" target="_blank" rel="noopener">View all on MedGen ↗</a></div>`
          : '';
        section11Body = `<table class="gc-medgen-table">
            <thead><tr><th>Condition</th><th>MIM</th><th>Inheritance</th></tr></thead>
            <tbody>${rows}</tbody>
          </table>${overflowNote}`;
      }

      const medgenSourceLink = md.url
        ? `<div class="gc-source-link"><a href="${md.url}" target="_blank" rel="noopener">View gene on MedGen ↗</a></div>`
        : '';
      const section11 = `<div class="gc-collapsible" id="gc-coll-medgen" style="--gc-accent: #4a5664;">
        <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
          <span class="gc-coll-title">Gene-disease associations (MedGen)</span>
          <span class="gc-coll-summary"></span>
          <span class="gc-coll-chevron" aria-hidden="true">▶</span>
        </button>
        <div class="gc-coll-body">
          <div class="gc-coll-inner">
            ${section11Body}
            ${medgenSourceLink}
          </div>
        </div>
      </div>`;

      const _gcxBest = (gc.clingen_phenotype_matched_count ? gc.clingen_phenotype_matched_best_classification : gc.best_clingen_classification)
        || (gc.phenotype_matched_count ? gc.phenotype_matched_best_classification : gc.best_classification) || '';
      const _gcxPheno = (gc.clingen_phenotype_matched_count || gc.phenotype_matched_count) > 0;
      const _gcxPillCls = raw => { const t = String(raw || '').toLowerCase();
        if (t.includes('definitive') || t.startsWith('well')) return 'is-def';
        if (t.includes('strong')) return 'is-strong'; if (t.includes('moderate')) return 'is-mod';
        if (t.includes('supportive')) return 'is-sup'; if (t.includes('disputed') || t.includes('refuted')) return 'is-na';
        return 'is-na'; };
      const headValidityPill = _gcxBest
        ? `<span class="gcx-pill ${_gcxPillCls(_gcxBest)}"><span class="gcx-dot"></span>${esc(_gcxBest)}</span>` + (_gcxPheno ? `<span class="gcx-stat__sub-text">pheno match</span>` : '')
        : `<span class="gcx-pill is-na"><span class="gcx-dot"></span>No GCEP curation</span>`;
      const headInh = (moiList && moiList.length) ? moiList.join(', ') : 'Not specified';
      const headMech = constraintType || 'Not informative';
      const headMechSub = constraintType ? constraintSub : '';
      const headLinkBits = [
        clingenSourceUrl && `<a href="${clingenSourceUrl}" target="_blank" rel="noopener">ClinGen ↗</a>`,
        genccSourceUrl && `<a href="${genccSourceUrl}" target="_blank" rel="noopener">GenCC ↗</a>`,
        gnomadGeneUrl && `<a href="${gnomadGeneUrl}" target="_blank" rel="noopener">gnomAD ↗</a>`,
        gene && `<a href="https://www.ncbi.nlm.nih.gov/clinvar/?gene=${encodeURIComponent(gene)}" target="_blank" rel="noopener">ClinVar ↗</a>`,
        otEnsembl && `<a href="https://platform.opentargets.org/target/${encodeURIComponent(otEnsembl)}" target="_blank" rel="noopener">Open Targets ↗</a>`,
      ].filter(Boolean);

      const gcxTiles = [
        `<div class="ev-snap__tile" data-tone="${_gcxBest ? 'info' : 'muted'}"><div class="ev-snap__k">Gene–disease validity</div><div class="ev-snap__v is-text">${esc(_gcxBest || 'No GCEP curation')}</div><div class="ev-snap__sub">${_gcxPheno ? 'Phenotype match' : 'ClinGen GCEP'}</div></div>`,
        `<div class="ev-snap__tile" data-tone="info"><div class="ev-snap__k">Inheritance</div><div class="ev-snap__v is-text">${esc(headInh)}</div><div class="ev-snap__sub">Mode of inheritance</div></div>`,
        `<div class="ev-snap__tile" data-tone="info"><div class="ev-snap__k">Constraint mechanism</div><div class="ev-snap__v is-text">${esc(headMech)}</div>${headMechSub ? `<div class="ev-snap__sub">${esc(headMechSub)}</div>` : ''}</div>`,
      ];
      const gcxSnapStrip = `<div class="ev-snap">${gcxTiles.join('')}</div>`;

      const geneExtHTML = _buildExternalToolsHTML(ev, 'gene');
      const geneExtSection = geneExtHTML
        ? `<section class="ev-exttools v11-section"><div class="v11-eyebrow"><h3>Open external tools</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>${geneExtHTML}</section>`
        : '';

      return `
        <div class="tab-overview">
          <div class="tab-overview__head">
            <div class="tab-overview__title"><span class="tab-overview__gene ${gene ? '' : 'is-empty'}">${esc(gene || 'Unknown gene')}</span></div>
            <div class="tab-overview__sub">Gene-level evidence &amp; disease validity</div>
          </div>
          ${gcxSnapStrip}
        </div>
        <div class="gcx-acc">
          <div class="gcx-group">
            <div class="gcx-group__hd"><h3>Disease association</h3><span class="gcx-group__ln"></span><button type="button" class="ev-expand-toggle" onclick="toggleAllGeneSections(this)" data-state="collapsed">Expand all</button></div>
            ${section1}${section2}${section5}${section11}
          </div>
          <div class="gcx-group">
            <div class="gcx-group__hd"><h3>Biological context</h3><span class="gcx-group__ln"></span><button type="button" class="ev-expand-toggle" onclick="toggleAllGeneSections(this)" data-state="collapsed">Expand all</button></div>
            ${section3}${section4}${section4b}${section9}${section10}${section6}
          </div>
        </div>
        ${geneExtSection}`;
    }

    function _renderTranscriptTable(vep) {
      const rows = (vep && vep.transcript_consequences_all) || [];
      if (!rows.length) return '';
      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
      const cCase = s => esc(String(s || '').replace(/_/g, ' ').replace(/^./, c => c.toUpperCase()));
      const impactCls = imp => {
        const i = String(imp || '').toUpperCase();
        return i === 'HIGH' ? 'ev-cls ev-cls-path'
             : i === 'MODERATE' ? 'ev-cls ev-cls-lp'
             : i === 'LOW' ? 'ev-cls ev-cls-vus' : 'ev-dim';
      };
      const stripV = id => String(id || '').split('.')[0];
      const txHref = id => {
        const s = stripV(id);
        if (/^ENST\d/i.test(s)) return `https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=${encodeURIComponent(s)}`;
        if (/^(?:NM|NR|XM|XR)_\d/i.test(s)) return `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(s)}`;
        return '';
      };
      const lastTok = s => esc(String(s).split(':').pop());

      const TX_VISIBLE = 6;
      const txPrio = r => (r.is_mane_select ? 4 : 0) + (r.is_mane_plus_clinical ? 2 : 0) + (r.is_picked ? 1 : 0);
      const ranked = rows
        .map((r, i) => [r, i])
        .sort((a, b) => (txPrio(b[0]) - txPrio(a[0])) || (a[1] - b[1]))
        .map(pair => pair[0]);
      const bodyRows = ranked.map((r, i) => {
        const href = txHref(r.transcript_id);
        const txCell = href
          ? `<a href="${href}" target="_blank" rel="noopener" class="ev-mono">${esc(r.transcript_id)}</a>`
          : `<span class="ev-mono">${esc(r.transcript_id)}</span>`;
        const badges = [];
        if (r.is_mane_select) badges.push('<span class="ev-tx-badge ev-tx-badge--select" title="MANE Select transcript">MANE Select</span>');
        if (r.is_mane_plus_clinical) badges.push('<span class="ev-tx-badge ev-tx-badge--clinical" title="MANE Plus Clinical transcript">MANE Plus Clinical</span>');
        if (r.is_picked) badges.push('<span class="ev-tx-badge ev-tx-badge--scored" title="Transcript used for the ACMG classification">Scored</span>');
        const srcLine = r.source ? `<div class="ev-tx-src ev-dim">${esc(r.source)}</div>` : '';
        const conseq = `<span class="${impactCls(r.impact)}">${cCase(r.consequence)}</span>`;
        const hgvsc = r.hgvsc ? `<span class="ev-mono">${lastTok(r.hgvsc)}</span>` : '<span class="ev-dim">—</span>';
        const hgvsp = r.hgvsp ? `<span class="ev-mono">${lastTok(r.hgvsp)}</span>` : '<span class="ev-dim">—</span>';
        const cls = [r.is_picked ? 'ev-tx-row--scored' : '', i >= TX_VISIBLE ? 'ev-tx-row--extra' : ''].filter(Boolean).join(' ');
        const rowCls = cls ? ` class="${cls}"` : '';
        return `<tr${rowCls}><td>${txCell}${badges.length ? ' ' + badges.join(' ') : ''}${srcLine}</td>`
          + `<td>${conseq}</td><td>${hgvsc}</td><td>${hgvsp}</td></tr>`;
      }).join('');
      const moreBtn = rows.length > TX_VISIBLE
        ? `<button type="button" class="ev-tx-more" onclick="var w=this.closest('.ev-tx-wrap');var s=w.classList.toggle('show-all');this.textContent=s?'Show fewer transcripts':'Show all ${rows.length} transcripts'">Show all ${rows.length} transcripts</button>`
        : '';
      const setNote = vep.transcript_set
        ? `<div class="ev-tx-setnote ev-dim">Transcripts from the ${esc(vep.transcript_set === 'GENCODE' ? 'Ensembl/GENCODE' : vep.transcript_set)} set returned by this VEP query.</div>`
        : '';
      const alert = vep.consequence_differs_significantly
        ? `<div class="ev-tx-alert" role="note">⚠ <strong>Consequence differs across transcripts.</strong> The predicted molecular consequence is not the same on all transcripts (e.g. loss-of-function on one isoform but a milder change on another). Confirm the transcript used for classification is the clinically relevant one.</div>`
        : '';
      const table = `<table class="ev-tx-table"><thead><tr>`
        + `<th>Transcript</th><th>Consequence</th><th>HGVS c.</th><th>HGVS p.</th>`
        + `</tr></thead><tbody>${bodyRows}</tbody></table>`;
      return `<div class="ev-tx-wrap"><div class="ev-tx-title">Transcript consequences `
        + `<span class="ev-dim">(${rows.length})</span></div>${alert}${table}${moreBtn}${setNote}</div>`;
    }

    function _renderEvidence(ev) {
      if (!ev) return '';

      const sCase = s => {
        if (s === null || s === undefined || s === '') return '';
        const t = String(s).replace(/_/g, ' ').toLowerCase();
        return t.charAt(0).toUpperCase() + t.slice(1);
      };

      const kv = (label, value, opts) => {
        const cls = (opts && opts.bold) ? 'ev-kv-value ev-kv-bold' : 'ev-kv-value';
        return `<div class="ev-kv"><div class="ev-kv-label">${label}</div><div class="${cls}">${value}</div></div>`;
      };
      const kvGroup = items => `<div class="ev-kv-group">${items.filter(Boolean).join('')}</div>`;

      const evTile = (k, v, sub, tone, isText) =>
        `<div class="ev-snap__tile" data-tone="${tone || 'muted'}">`
        + `<div class="ev-snap__k">${k}</div>`
        + `<div class="ev-snap__v${isText ? ' is-text' : ''}">${v}</div>`
        + (sub ? `<div class="ev-snap__sub">${sub}</div>` : '')
        + `</div>`;
      const evTiles = tiles => `<div class="ev-tilerow">${tiles.filter(Boolean).join('')}</div>`;

      const evScoreBar = (val, opts = {}) => {
        const v = Number(val);
        const pct = Number.isFinite(v) ? Math.max(0, Math.min(1, v)) * 100 : 0;
        const zoneCls = opts.neutral ? 'ev-bar--neutral'
                      : opts.dir === 'hi-good' ? 'ev-bar--rev' : 'ev-bar--zones';
        const ticks = (opts.ticks || []).map(t =>
          `<span class="ev-bar-tick" style="left:${(t.at * 100).toFixed(1)}%" title="${t.label || ''}"></span>`).join('');
        const dot = Number.isFinite(v)
          ? `<span class="ev-bar-dot" style="left:${pct.toFixed(1)}%"></span>` : '';
        const scale = (opts.scale || []).length
          ? `<div class="ev-bar-scale">${opts.scale.map((s, i) =>
              `<span class="${(i === 1 && opts.scale.length === 3) ? 'mid' : ''}">${s}</span>`).join('')}</div>` : '';
        return `<div class="ev-barwrap">`
          + (opts.label ? `<div class="ev-bar-lab">${opts.label}</div>` : '')
          + `<div class="ev-bar ${zoneCls}">${ticks}${dot}</div>${scale}</div>`;
      };

      const afPos = af => {
        if (!(af > 0)) return 0;
        const lo = -6, hi = -2, l = Math.log10(af);
        return Math.max(0, Math.min(1, (l - lo) / (hi - lo)));
      };

      const clsClass = raw => {
        const s = String(raw || '').toLowerCase();
        if (!s) return '';
        if (/conflict/.test(s)) return 'ev-cls ev-cls-conflict';
        if (/likely pathogenic/.test(s)) return 'ev-cls ev-cls-lp';
        if (/pathogenic/.test(s))        return 'ev-cls ev-cls-path';
        if (/likely benign/.test(s))     return 'ev-cls ev-cls-lb';
        if (/benign/.test(s))            return 'ev-cls ev-cls-ben';
        if (/uncertain|vus/.test(s))     return 'ev-cls ev-cls-vus';
        if (/definitive/.test(s))        return 'ev-cls ev-cls-definitive';
        if (/strong/.test(s))            return 'ev-cls ev-cls-strong';
        if (/moderate/.test(s))          return 'ev-cls ev-cls-moderate';
        if (/limited|animal/.test(s))    return 'ev-cls ev-cls-limited';
        if (/disputed|refuted/.test(s))  return 'ev-cls ev-cls-disputed';
        return 'ev-cls';
      };
      const rows = [];

      const vep = ev.vep || {};
      const userTxRow = vep.user_transcript
        ? kv('User-specified transcript', `<span class="ev-mono">${vep.user_transcript}</span>`)
        : '';

      const _maneRow = (vep.transcript_consequences_all || [])
        .find(r => r && r.is_mane_select) || null;
      const _maneAcc = vep.mane_select_accession
        || (_maneRow && (_maneRow.mane_select_accession || _maneRow.transcript_id))
        || '';
      const nonManeRow = (vep.ok && !vep.is_mane_select && !vep.is_mane_clinical)
        ? kv('Numbering caveat',
            `<span class="ev-cls ev-dim">Not MANE Select</span> `
            + (_maneAcc
                ? `ClinVar names its records on <span class="ev-mono">${_hvEsc(_maneAcc)}</span>, `
                  + `not on <span class="ev-mono">${_hvEsc(vep.selected_transcript_id || vep.transcript_id || '')}</span>. `
                : `ClinVar names its records on this gene's MANE Select transcript, not on `
                  + `<span class="ev-mono">${_hvEsc(vep.selected_transcript_id || vep.transcript_id || '')}</span>. `)
            + `Residue and c. numbers shown here are on your transcript; ClinVar's `
            + `may differ. HeartVar matches ClinVar in MANE numbering internally.`)
        : '';
      if (vep.ok) {
        const imp = String(vep.impact || '').toUpperCase();
        const pill =
          imp === 'HIGH'     ? { tone: 'warn',  text: 'High impact' }
          : imp === 'MODERATE' ? { tone: 'amber', text: 'Moderate impact' }
          : imp === 'LOW'    ? { tone: 'info',  text: 'Low impact' }
          : { tone: 'muted', text: 'Modifier' };
        const consText = sCase(vep.most_severe_consequence) || 'Consequence unknown';
        const summary = consText;

        const txId = vep.transcript_id || '';
        const txStripped = txId.split('.')[0];
        let txHref = '';
        if (/^ENST\d/i.test(txStripped)) {
          txHref = `https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=${encodeURIComponent(txStripped)}`;
        } else if (/^(?:NM|NR|XM|XR)_\d/i.test(txStripped)) {
          txHref = `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(txStripped)}`;
        }
        const txLink = txId
          ? (txHref
              ? `<a href="${txHref}" target="_blank" rel="noopener" class="ev-mono">${txId}</a>`
              : `<span class="ev-mono">${txId}</span>`)
          : '';

        const rsid = (ev.gnomad && ev.gnomad.variant && ev.gnomad.variant.rsid) || '';
        const rsidPill = rsid
          ? `<a class="ev-rsid-pill" href="https://www.ensembl.org/id/${encodeURIComponent(rsid)}" target="_blank" rel="noopener" title="Open ${rsid} on Ensembl Variation">${rsid} ↗</a>`
          : '';
        const startStr = Number.isFinite(Number(vep.start)) ? String(Number(vep.start)) : (vep.start || '');
        const alleleColon = String(vep.allele_string || '').replace('/', ':');
        const coordValue = `chr${vep.seq_region_name || '?'}:${startStr}`
          + (alleleColon ? `:${alleleColon}` : '')
          + ` (${vep.assembly_name || 'GRCh38'})`;
        const hgvsLine = (vep.hgvsc || vep.hgvsp)
          ? `<div class="ev-mono ev-dim" style="margin-top:2px">${vep.hgvsc || ''}${vep.hgvsp ? ' · ' + vep.hgvsp : ''}</div>`
          : '';
        const txValue = (txLink || '—') + rsidPill + hgvsLine;
        const impactCls = imp === 'HIGH' ? 'ev-cls ev-cls-path'
                       : imp === 'MODERATE' ? 'ev-cls ev-cls-lp'
                       : imp === 'LOW' ? 'ev-cls-vus' : 'ev-dim';
        const impactValue = imp ? `<span class="${impactCls}">${imp}</span>` : '<span class="ev-dim">—</span>';
        const naDim = '<span class="ev-dim">not available</span>';
        rows.push({
          source: 'vep',
          group: 'functional',
          label: 'Ensembl VEP',
          summary,
          pill,
          body: kvGroup([
            userTxRow,
            nonManeRow,
            kv('Coordinates', coordValue),
            kv('Transcript', txValue),
            kv('Consequence', consText),
            kv('Impact', impactValue),

            kv('REVEL', (vep.revel_score !== null && vep.revel_score !== undefined)
              ? `${_fmtNum(vep.revel_score, 3)} <span class="ev-dim">→ drives PP3/BP4 (calibrated, capped Supporting)</span>`
              : `${naDim} <span class="ev-dim">· REVEL scores missense SNVs only</span>`),
            kv('CADD (PHRED)', (vep.cadd_phred !== null && vep.cadd_phred !== undefined)
              ? `${_fmtNum(vep.cadd_phred, 1)} <span class="ev-dim">· supporting in-silico context</span>` : naDim),
            kv('PolyPhen', vep.polyphen_prediction
              ? `${sCase(vep.polyphen_prediction)} (${vep.polyphen_score}) <span class="ev-dim">· supporting in-silico context</span>` : naDim),
            kv('SIFT', vep.sift_prediction
              ? `${sCase(vep.sift_prediction)} (${vep.sift_score}) <span class="ev-dim">· supporting in-silico context</span>` : naDim),
          ]) + _renderTranscriptTable(vep),
        });
      } else {
        rows.push({
          source: 'vep', group: 'functional', label: 'Ensembl VEP',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Error' },
          body: kvGroup([
            userTxRow,
            nonManeRow,
            kv('Status', `<em>Lookup failed: ${vep.error || 'unknown'}</em>`),
          ]),
        });
      }

      const _sameSiteRows = () => {
        const ss = ev.same_site || {};
        if (!ss.ok) return [];
        if (!ss.gnomad_available) {

          const why = ss.gnomad_note || ss.gnomad_reason;
          return why
            ? [kv('Other alleles at this site',
                  `<span class="ev-dim">Not available — ${_hvEsc(why)}</span>`)]
            : [];
        }
        const fmtAllele = (a) => {
          const href = `https://gnomad.broadinstitute.org/variant/${encodeURIComponent(a.variant_id)}?dataset=gnomad_r4`;

          const obs = !!a.ac;
          return `<span style="display:inline-block;margin-right:14px;white-space:nowrap">`
            + `<a href="${href}" target="_blank" rel="noopener" class="ev-mono">${_hvEsc(a.ref || '')}&gt;${_hvEsc(a.alt || '')}</a>`
            + (obs
                ? ` ${formatAF(a.af)} <span class="ev-dim">(AC ${a.ac})</span>`
                : ` <span class="ev-dim">not observed (AC 0)</span>`)
            + (a.rsid ? ` <span class="ev-dim">· ${_hvEsc(a.rsid)}</span>` : '')
            + `</span>`;
        };
        const out = [];
        const nt = ss.same_nucleotide || [];

        const _base = ss.chrom
          ? `chr${String(ss.chrom).replace(/^chr/i, '')}:${ss.variant_position || ''}`
          : '';
        out.push(kv(
          `Other alleles at this base${_base ? ` <span class="ev-dim">(${_hvEsc(_base)})</span>` : ''}`,
          nt.length
            ? nt.map(fmtAllele).join('')
            : '<span class="ev-dim">None reported in gnomAD v4</span>'
        ));
        if (ss.codon_span_available) {
          const cd = ss.same_codon || [];

          const confirmed = ss.codon_span_confirmed_on_picked_transcript;
          const note = ss.codon_span_note || '';
          const label = confirmed
            ? 'Elsewhere in this codon'
            : `<span title="${_hvEsc(note)}">Elsewhere in this codon</span>`;
          const prov = (!confirmed && ss.codon_span_source)
            ? ` <span class="ev-dim">· codon span from ${_hvEsc(ss.codon_span_source)}</span>`
            : '';
          out.push(kv(
            label,
            (cd.length
              ? cd.map(a => `<span class="ev-mono ev-dim">${a.position}</span> ${fmtAllele(a)}`).join('')
              : '<span class="ev-dim">None reported in gnomAD v4</span>') + prov
          ));
        }

        const _probandRsid = ((ev.gnomad || {}).variant || {}).rsid || null;
        const _rsidKnown = nt.map(a => !!a.rsid).concat([!!_probandRsid]);
        if (nt.length && _rsidKnown.some(x => !x) && _rsidKnown.some(x => x)) {
          out.push(kv('rsID note',
            '<span class="ev-dim">Alleles at one base can differ in whether they '
            + 'carry an rsID — these are matched on position and allele, not rsID.</span>'));
        }
        return out;
      };
      const g = ev.gnomad || {};
      if (g.ok) {
        const v = g.variant;
        let pill;
        let freqGroup;
        if (v) {
          const ex = v.exome || {};
          const ge = v.genome || {};
          const popmaxVals = [ex.faf95?.popmax, ge.faf95?.popmax].filter(x => x !== null && x !== undefined);
          const popmax = popmaxVals.length ? Math.max(...popmaxVals) : null;
          const af = ex.af ?? ge.af ?? null;
          if (popmax !== null && popmax >= 0.05) {
            pill = { tone: 'warn', text: 'Common' };
          } else if (popmax !== null && popmax >= 0.01) {
            pill = { tone: 'amber', text: `AF ${formatAF(popmax)}` };
          } else {
            pill = { tone: 'info', text: af !== null ? `AF ${formatAF(af)}` : 'Rare' };
          }
          const gnomadHref = `https://gnomad.broadinstitute.org/variant/${encodeURIComponent(v.variantId)}?dataset=gnomad_r4`;
          const variantValue = `<a href="${gnomadHref}" target="_blank" rel="noopener" class="ev-mono">${formatVariantCoord(v.variantId, ev.vep?.assembly_name)}</a>`
            + (v.rsid ? ` <span class="ev-dim">· ${v.rsid}</span>` : '');

          const exHasCounts = ex.ac !== undefined && ex.ac !== null;
          const geHasCounts = ge.ac !== undefined && ge.ac !== null;
          const anyCounts = exHasCounts || geHasCounts;
          const ac = (exHasCounts ? ex.ac : 0) + (geHasCounts ? ge.ac : 0);
          const anTotal = (exHasCounts ? (ex.an || 0) : 0) + (geHasCounts ? (ge.an || 0) : 0);
          const hom = (exHasCounts ? (ex.ac_hom || 0) : 0) + (geHasCounts ? (ge.ac_hom || 0) : 0);
          const het = Math.max(0, ac - 2 * hom);
          const exHemi = ex.ac_hemi;
          const geHemi = ge.ac_hemi;
          const hemiAvailable = (exHemi !== null && exHemi !== undefined) || (geHemi !== null && geHemi !== undefined);
          const hemi = (exHemi || 0) + (geHemi || 0);

          const onSexChrom = /^(X|Y|chrX|chrY)-/i.test(v.variantId || '');

          const popLabels = {
            afr: 'African / African-American',
            amr: 'Latino / Admixed American',
            asj: 'Ashkenazi Jewish',
            eas: 'East Asian',
            fin: 'Finnish European',
            mid: 'Middle Eastern',
            nfe: 'Non-Finnish European',
            sas: 'South Asian',
            remaining: 'Remaining',
          };

          const popMap = new Map();
          for (const src of [ex, ge]) {
            for (const p of (src && src.populations) || []) {
              if (!p.id || !(p.id in popLabels)) continue;
              const cur = popMap.get(p.id) || { ac: 0, hom: 0, hemi: 0 };
              cur.ac   += p.ac      || 0;
              cur.hom  += p.ac_hom  || 0;
              cur.hemi += p.ac_hemi || 0;
              popMap.set(p.id, cur);
            }
          }
          const carriersByPop = [...popMap.entries()]
            .filter(([, v]) => v.ac > 0)
            .sort((a, b) => b[1].ac - a[1].ac);
          let carriersPopHtml = '';
          if (carriersByPop.length) {
            const parts = carriersByPop.map(([id, v]) => {
              const label = popLabels[id] || id;
              const homNote = v.hom > 0 ? ` <span class="ev-dim">(${v.hom} hom)</span>` : '';
              const hemiNote = (onSexChrom && v.hemi > 0) ? ` <span class="ev-dim">(${v.hemi} hemi)</span>` : '';
              return `<span style="display:inline-block;margin-right:14px;white-space:nowrap"><strong>${label}:</strong> ${v.ac}${homNote}${hemiNote}</span>`;
            });
            carriersPopHtml = parts.join('');
          }

          const topPop = carriersByPop.length ? carriersByPop[0] : null;
          freqGroup =
            evTiles([
              evTile('Genome AF', geHasCounts ? formatAF(ge.af) : '—',
                     geHasCounts ? `${ge.ac} / ${ge.an}` : 'not in genomes',
                     geHasCounts ? 'info' : 'muted', !geHasCounts),
              evTile('Exome AF', exHasCounts ? formatAF(ex.af) : '—',
                     exHasCounts ? `${ex.ac} / ${ex.an}` : 'not in exomes',
                     exHasCounts ? 'info' : 'muted', !exHasCounts),
              evTile('Homozygotes', anyCounts ? String(hom) : '—',
                     hom > 0 ? 'observed' : (anyCounts ? 'none observed' : ''),
                     hom > 0 ? 'warn' : 'good'),
              topPop ? evTile('Top population', String(topPop[0]).toUpperCase(),
                     `${popLabels[topPop[0]] || topPop[0]} · ${topPop[1].ac}`, 'muted', true) : '',
            ])
            + (af > 0 ? evScoreBar(afPos(af), {
                label: 'Allele frequency · gnomAD (log scale)',
                scale: ['10⁻⁶', '10⁻⁴', '10⁻²'],
                neutral: true,
              }) : '')
            + kvGroup([
              kv('Variant ID', variantValue),

              kv('Total AC / AN', anyCounts
                ? `${ac} / ${anTotal} <span class="ev-dim">(exome ${
                    exHasCounts ? `${ex.ac}/${ex.an}` : 'not in dataset'
                  } · genome ${
                    geHasCounts ? `${ge.ac}/${ge.an}` : 'not in dataset'
                  })</span>` : '—'),
              kv('Carriers', anyCounts
                ? `${het} heterozygous <span class="ev-dim">·</span> ${hom} homozygous`
                  + ((hemiAvailable && onSexChrom) ? ` <span class="ev-dim">·</span> ${hemi} hemizygous` : '')
                : '—'),
              (carriersByPop.length > 1 && carriersPopHtml) ? kv('Carriers by population', carriersPopHtml) : '',
              ..._sameSiteRows(),
            ]);
        } else {
          pill = { tone: 'good', text: 'Absent' };
          freqGroup = evTiles([ evTile('gnomAD v4', 'Absent', 'Not observed in gnomAD v4 (exomes + genomes)', 'good', true) ])
            + kvGroup(_sameSiteRows());
        }
        const summary = '';
        rows.push({ source: 'gnomad', group: 'population', label: 'gnomAD v4', summary, pill, body: freqGroup });
      } else {
        rows.push({
          source: 'gnomad', group: 'population', label: 'gnomAD v4',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Error' },
          body: `<em>Lookup failed: ${g.error || 'unknown'}</em>`,
        });
      }

      const sa = ev.spliceai || {};
      if (sa.skipped) {
        rows.push({
          source: 'spliceai', group: 'functional', label: 'SpliceAI',
          summary: '',
          pill: { tone: 'muted', text: 'Not applicable' },
          body: `<em>Not applicable (indel)${sa.reason ? ' — ' + sa.reason : ''}</em>`,
        });
      } else if (sa.ok) {
        const fmt = (x) => (x === null || x === undefined ? '—' : Number(x).toFixed(2));
        const maxd = sa.max_delta;

        let pill;
        if (maxd === null || maxd === undefined) {
          pill = { tone: 'muted', text: 'No score' };
        } else if (maxd >= 0.8) {
          pill = { tone: 'warn',  text: 'High impact' };
        } else if (maxd >= 0.5) {
          pill = { tone: 'amber', text: 'Moderate impact' };
        } else if (maxd >= 0.2) {
          pill = { tone: 'muted', text: 'Low impact' };
        } else {
          pill = { tone: 'good',  text: 'No impact' };
        }
        const summary = '';
        const t0 = (sa.scores_per_transcript || [])[0];
        const maxColour = (maxd >= 0.8) ? 'ev-cls ev-cls-path'
                       : (maxd >= 0.5) ? 'ev-cls ev-cls-lp'
                       : (maxd >= 0.2) ? 'ev-cls ev-cls-vus'
                       : 'ev-cls ev-cls-ben';
        const interp = (maxd >= 0.8) ? 'High confidence splice impact'
                     : (maxd >= 0.5) ? 'Moderate confidence splice impact'
                     : (maxd >= 0.2) ? 'Low confidence splice impact'
                     : 'No predicted splice impact';
        const txLine = (t0 && (t0.gene || t0.transcript_id))
          ? kv('Transcript', `<span class="ev-mono ev-dim">${t0.gene || ''}${t0.gene && t0.transcript_id ? ' · ' : ''}${t0.transcript_id || ''}</span>`)
          : '';

        const saLink = sa.variant_id
          ? kv('Link', `<a href="https://spliceailookup.broadinstitute.org/#variant=${encodeURIComponent(sa.variant_id)}&hg=38&bc=basic&distance=500&mask=1" target="_blank" rel="noopener" class="ev-dim">SpliceAI Lookup ↗</a>`)
          : '';
        const dsLabels = { DS_AG: 'Acceptor gain', DS_AL: 'Acceptor loss', DS_DG: 'Donor gain', DS_DL: 'Donor loss' };
        const dsTone = x => { const n = Number(x); return !Number.isFinite(n) ? 'muted' : n >= 0.8 ? 'warn' : n >= 0.5 ? 'amber' : n >= 0.2 ? 'muted' : 'good'; };
        const dsTiles = t0 ? evTiles(['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'].map(k =>
          evTile(k.replace('DS_', ''), fmt(t0[k]), dsLabels[k], dsTone(t0[k]), false))) : '';
        rows.push({ source: 'spliceai', group: 'functional', label: 'SpliceAI', summary, pill, body:
          ((maxd !== null && maxd !== undefined) ? evScoreBar(maxd, {
            label: `Max Δscore · ${interp}`,
            ticks: [{ at: 0.2, label: 'low ≥ 0.2' }, { at: 0.5, label: 'moderate ≥ 0.5' }, { at: 0.8, label: 'high ≥ 0.8' }],
            scale: ['0', '0.5', '1'],
          }) : '')
          + dsTiles
          + kvGroup([
            kv('Max ΔScore', `<span class="${maxColour}">${fmt(maxd)}</span>`
              + (sa.model_message ? ` <span class="ev-dim">— ${sa.model_message}</span>` : '')),
            kv('Interpretation', interp),
            txLine,
            saLink,
          ])
        });
      } else if (sa.lookup_failed) {

        rows.push({
          source: 'spliceai', group: 'functional', label: 'SpliceAI',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Lookup failed' },
          body: `<em>SpliceAI lookup did not complete — ${sa.error || 'unknown error'}. The score for this variant is unknown; this is <strong>not</strong> evidence of "no impact".</em>`,
        });
      } else {

        const reason = sa.reason || sa.error || '';
        rows.push({
          source: 'spliceai', group: 'functional', label: 'SpliceAI',
          summary: '',
          pill: { tone: 'muted', text: 'No score at position' },
          body: `<em>SpliceAI has no score at this position — the variant does not overlap an annotated splice site${reason ? ` (${reason})` : ''}.</em>`,
        });
      }

      const am = ev.alphamissense || {};
      if (am.available === false) {
        rows.push({
          source: 'alphamissense', group: 'functional', label: 'AlphaMissense',
          summary: 'Not available',
          pill: { tone: 'muted', text: 'Not configured' },
          body: `<em>Not available — ${am.reason || 'ALPHAMISSENSE_PATH not configured'}</em>`,
        });
      } else if (am.not_applicable) {
        rows.push({
          source: 'alphamissense', group: 'functional', label: 'AlphaMissense',
          summary: '',
          pill: { tone: 'muted', text: 'Not applicable' },
          body: `<em>Not applicable (non-missense) — AlphaMissense covers only single missense substitutions. Expected outcome for frameshift, stop-gain, splice, and in-frame indel variants.</em>`,
        });
      } else {
        const score = am.score;
        const cls = am.classification || 'ambiguous';

        const scoreClass = cls === 'likely_pathogenic' ? 'ev-cls ev-cls-path'
                         : cls === 'likely_benign'     ? 'ev-cls ev-cls-ben'
                         : 'ev-cls ev-cls-vus';
        const pill = cls === 'likely_pathogenic' ? { tone: 'warn', text: 'Likely pathogenic' }
                   : cls === 'likely_benign'     ? { tone: 'good', text: 'Likely benign' }
                   : { tone: 'muted', text: 'Ambiguous' };
        const readableCls = cls === 'likely_pathogenic' ? 'Likely pathogenic'
                          : cls === 'likely_benign'     ? 'Likely benign'
                          : 'Ambiguous';
        const interp = cls === 'likely_pathogenic'
            ? 'Supports PP3 (concordant deleterious in-silico evidence)'
            : cls === 'likely_benign'
              ? 'Supports BP4 (concordant benign in-silico evidence)'
              : 'Neither PP3 nor BP4 — score within ambiguous band (0.340–0.564)';
        rows.push({
          source: 'alphamissense', group: 'functional', label: 'AlphaMissense',
          summary: '',
          pill,
          body:
            evScoreBar(score, {
              label: `AlphaMissense pathogenicity · ${readableCls}`,
              ticks: [{ at: 0.340, label: 'likely-benign ≤ 0.340' }, { at: 0.564, label: 'likely-pathogenic ≥ 0.564' }],
              scale: ['0 · benign', '0.34–0.56', '1 · pathogenic'],
            })
            + kvGroup([
              kv('Score', `<span class="${scoreClass}">${Number(score).toFixed(3)}</span>`
                         + ` <span class="ev-dim">— ${readableCls}</span>`),
              am.protein_variant
                ? kv('Protein variant', `<span class="ev-mono">p.${am.protein_variant}</span>`)
                : '',
              kv('Interpretation', interp),
            ]),
        });
      }

      const cv = ev.clinvar || {};
      if (cv.ok) {
        if (cv.found) {
          const renderStars = n => {
            const k = Math.max(0, Math.min(4, n || 0));
            return `<span class="ev-stars" title="ClinVar review status: ${k}/4 stars">${'★'.repeat(k)}${'☆'.repeat(4 - k)}</span>`;
          };
          const recs = cv.records || [];

          const totalSubs = (cv.total_submissions != null)
            ? cv.total_submissions
            : recs.reduce((s, r) => s + (r.number_submitters || 1), 0);
          const matchedSubs = (cv.phenotype_matched_submissions != null)
            ? cv.phenotype_matched_submissions
            : totalSubs;
          const subsGap = Math.max(0, totalSubs - matchedSubs);

          const isPath = r => /pathogenic|likely pathogenic/i.test(r.clinical_significance || '');
          const isBenign = r => /benign/i.test(r.clinical_significance || '') && !isPath(r);
          const isConflict = r => /conflict/i.test(r.clinical_significance || '') || /conflict/i.test(r.review_status || '');
          const sorted = [...recs].sort((a, b) => (b.stars || 0) - (a.stars || 0));
          const best = sorted[0] || {};
          let pill;
          if (sorted.some(isPath)) {
            const r = sorted.find(isPath);
            pill = { tone: 'warn', text: r.clinical_significance.split(/[,/]/)[0].trim() || 'Pathogenic' };
          } else if (sorted.some(isBenign)) {
            const r = sorted.find(isBenign);
            pill = { tone: 'good', text: r.clinical_significance.split(/[,/]/)[0].trim() || 'Benign' };
          } else if (sorted.some(isConflict)) {
            pill = { tone: 'amber', text: 'Conflicting' };
          } else if (/uncertain|vus/i.test(best.clinical_significance || '')) {
            pill = { tone: 'amber', text: 'Uncertain significance' };
          } else {
            pill = { tone: 'info', text: `${totalSubs} submission${totalSubs === 1 ? '' : 's'}` };
          }
          const submissionsLabel = subsGap > 0
            ? `${totalSubs} total <span class="ev-dim">(${matchedSubs} matching submitted phenotype)</span>`
            : `${totalSubs} submission${totalSubs === 1 ? '' : 's'}`;
          const summary = subsGap > 0
            ? `${totalSubs} total (${matchedSubs} matching submitted phenotype)`
            : `${totalSubs} submission${totalSubs === 1 ? '' : 's'}`;
          const trs = recs.slice(0, 5).map(rec => {
            const cvHref = `https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(rec.variation_id)}/`;
            return `<tr>`
              + `<td><a href="${cvHref}" target="_blank" rel="noopener" class="ev-mono">${rec.accession}</a></td>`
              + `<td><span class="${clsClass(rec.clinical_significance)}">${rec.clinical_significance || '—'}</span></td>`
              + `<td><span class="ev-dim">${rec.review_status || ''}</span> ${renderStars(rec.stars)}</td>`
              + `</tr>`;
          }).join('');
          const tableBody = `<table class="ev-detail-table">`
            + `<thead><tr><th>Accession</th><th>Classification</th><th>Review status</th></tr></thead>`
            + `<tbody>${trs}</tbody></table>`;

          const allConditions = (cv.all_conditions && cv.all_conditions.length)
            ? cv.all_conditions
            : (() => {
                const seen = new Set();
                const out = [];
                for (const r of recs) {
                  for (const c of (r.conditions || [])) {
                    if (!seen.has(c)) { seen.add(c); out.push(c); }
                  }
                }
                return out;
              })();
          const conditionsHtml = allConditions.length
            ? allConditions.map(c => `<span class="ev-dim" style="margin-right:6px">${c}</span>`).join('<span class="ev-dim">·</span> ')
            : '<em class="ev-dim">No conditions listed</em>';
          const gapNoteHtml = subsGap > 0
            ? `<div class="ev-dim" style="margin-top:6px">Note: ${subsGap} additional submission${subsGap === 1 ? '' : 's'} classified under related conditions (e.g. ${allConditions[0] || 'a related diagnosis'}) — these may represent the same clinical entity as the submitted phenotype.</div>`
            : '';

          const mix = { path: 0, lp: 0, vus: 0, lb: 0, ben: 0, other: 0 };
          for (const r of recs) {
            const s = String(r.clinical_significance || '').toLowerCase();
            const n = r.number_submitters || 1;
            if (/conflict/.test(s)) mix.other += n;
            else if (/likely pathogenic/.test(s)) mix.lp += n;
            else if (/pathogenic/.test(s)) mix.path += n;
            else if (/likely benign/.test(s)) mix.lb += n;
            else if (/benign/.test(s)) mix.ben += n;
            else if (/uncertain|vus/.test(s)) mix.vus += n;
            else mix.other += n;
          }
          const mixSegs = [
            ['path', 'P', mix.path], ['lp', 'LP', mix.lp], ['vus', 'VUS', mix.vus],
            ['lb', 'LB', mix.lb], ['ben', 'B', mix.ben], ['other', 'Other', mix.other],
          ].filter(s => s[2] > 0);
          const mixTotal = mixSegs.reduce((s, x) => s + x[2], 0);
          const mixBar = (mixTotal > 0 && mixSegs.length > 0)
            ? `<div class="ev-clbar">` + mixSegs.map(s => `<span class="ev-clbar-seg ev-clbar-${s[0]}" style="flex:${s[2]}" title="${s[1]}: ${s[2]}"></span>`).join('') + `</div>`
              + `<div class="ev-clbar-legend">` + mixSegs.map(s => `<span><i class="ev-clbar-sw ev-clbar-${s[0]}"></i>${s[1]} ${s[2]}</span>`).join('') + `</div>`
            : '';
          const classMain = (best.clinical_significance || pill.text || '').split(/[,/]/)[0].trim() || pill.text;
          const toneColor = { warn: '#c8341a', good: '#2a6e3a', amber: '#926010', info: 'var(--ink)', muted: 'var(--muted)' }[pill.tone] || 'var(--ink)';
          const classLine = `<div class="ev-clline">`
            + `<span class="ev-clline-main" style="color:${toneColor}">${classMain}</span>`
            + renderStars(best.stars)
            + `<span class="ev-dim">${submissionsLabel}</span>`
            + `</div>`;
          const condGroup = kvGroup([ kv('Conditions in ClinVar', conditionsHtml) ]);
          rows.push({ source: 'clinvar', group: 'clinical', label: 'ClinVar', summary, pill, body: classLine + mixBar + tableBody + condGroup + gapNoteHtml });
        } else {
          rows.push({
            source: 'clinvar', group: 'clinical', label: 'ClinVar',
            summary: '',
            pill: { tone: 'muted', text: 'No records' },
            body: '<em>No matching records</em>',
          });
        }
      } else {
        rows.push({
          source: 'clinvar', group: 'clinical', label: 'ClinVar',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Error' },
          body: `<em>Lookup failed: ${cv.error || 'unknown'}</em>`,
        });
      }

      const pm = ev.pubmed || {};
      if (!pm.ok) {
        rows.push({
          source: 'pubmed', group: 'literature', label: 'PubMed',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Error' },
          body: `<em>Lookup failed: ${pm.error || 'unknown'}</em>`,
        });
      } else {
        const renderPaper = p => {
          const url = p.url || `https://pubmed.ncbi.nlm.nih.gov/${p.pmid}/`;
          const title = p.title || '(no title)';
          const meta = `${p.first_author || '?'} ${p.year || ''}`.trim();
          const titleAttr = title.replace(/"/g, '&quot;');
          return `<div class="ev-paper">
            <a class="ev-paper-title" href="${url}" target="_blank" rel="noopener" title="${titleAttr}">${title}</a>
            <span class="ev-paper-meta">— ${meta}</span>
          </div>`;
        };
        const vPapers = pm.variant_papers || [];
        const variantSearchURL = pm.variant_search_url || 'https://pubmed.ncbi.nlm.nih.gov/';
        const nShown = vPapers.length;

        const nTotal = (typeof pm.variant_total === 'number' && pm.variant_total >= nShown) ? pm.variant_total : nShown;
        const pill = nTotal > 0
          ? { tone: 'good', text: `${nTotal} variant paper${nTotal === 1 ? '' : 's'}` }
          : { tone: 'muted', text: 'No variant publications' };
        const segs = [];
        segs.push(nShown
          ? vPapers.map(renderPaper).join('') + (nTotal > nShown
              ? `<div class="ev-dim" style="margin-top:6px">Showing the top ${nShown} of ${nTotal} — open the PubMed search for the full list.</div>`
              : '')
          : '<em>No indexed variant-specific publications found</em>');
        segs.push(`<div class="ev-search-links">`
          + `<a class="ev-search-link" href="${variantSearchURL}" target="_blank" rel="noopener">Search PubMed for this variant ↗</a>`
          + `</div>`);
        rows.push({
          source: 'pubmed', group: 'literature', label: 'PubMed',
          summary: '', pill, body: segs.join(''),
        });
      }

      const pv = ev.protvar || {};
      if (!pv.ok) {
        rows.push({
          source: 'protvar', group: 'protein_expression', label: 'ProtVar',
          summary: 'Lookup failed',
          pill: { tone: 'warn', text: 'Error' },
          body: `<em>Lookup failed: ${pv.error || 'unknown'}</em>`,
        });
      } else if (!pv.applicable) {
        rows.push({
          source: 'protvar', group: 'protein_expression', label: 'ProtVar',
          summary: 'Non-missense variant',
          pill: { tone: 'muted', text: 'Not applicable' },
          body: `<em>Not applicable (non-missense variant)</em>`,
        });
      } else if (!pv.found) {
        rows.push({
          source: 'protvar', group: 'protein_expression', label: 'ProtVar',
          summary: 'No mapping returned',
          pill: { tone: 'muted', text: 'No mapping' },
          body: `<em>No mapping returned</em>`,
        });
      } else {
        const aa = `${pv.ref_aa || '?'}${pv.protein_position ?? '?'}${pv.alt_aa || '?'}`;
        const summary = aa;
        rows.push({
          source: 'protvar', group: 'protein_expression', label: 'ProtVar',
          summary,
          pill: { tone: 'info', text: pv.uniprot || 'Mapped' },
          body: kvGroup([
            kv('UniProt', `<a href="${pv.url || '#'}" target="_blank" rel="noopener" class="ev-mono">${pv.uniprot}</a>`),
            kv('Residue change', `<span class="ev-mono">${aa}</span> <span class="ev-dim">(${sCase(pv.consequence) || '—'})</span>`),
            (pv.conservation_score !== undefined && pv.conservation_score !== null)
              ? kv('Conservation', String(pv.conservation_score)) : '',
            (pv.feature_types || []).length ? kv('Features', pv.feature_types.join(', ')) : '',
            pv.colocated_total ? kv('Co-located variants', `${pv.colocated_total} reported`) : '',
          ]),
        });
      }

      if (vep.seq_region_name && vep.start) {
        const ucscChrRaw = String(vep.seq_region_name);
        const ucscChr    = /^chr/i.test(ucscChrRaw) ? ucscChrRaw : `chr${ucscChrRaw}`;
        const ucscPos    = Number(vep.start);
        const ucscStart  = Math.max(1, ucscPos - 200);
        const ucscEnd    = ucscPos + 200;
        const ucscUrl    = `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${encodeURIComponent(`${ucscChr}:${ucscStart}-${ucscEnd}`)}`;
        const ucscRegion = `${ucscChr}:${ucscStart.toLocaleString()}–${ucscEnd.toLocaleString()}`;
        rows.push({
          source: 'ucsc', group: 'visualisation', label: 'UCSC Genome Browser',
          summary: 'Genome region viewer · ±200 bp context',
          pill: { tone: 'info', text: 'View region' },
          body: kvGroup([
            kv('Variant locus',    `<span class="ev-mono">${ucscChr}:${ucscPos.toLocaleString()}</span>`),
            kv('Region (±200 bp)', `<a class="ev-mono" href="${ucscUrl}" target="_blank" rel="noopener">${ucscRegion} ↗</a>`),
          ]),
        });
      }

      const GROUPS = [
        { key: 'functional', label: 'Functional annotation',          members: ['functional', 'protein_expression'] },
        { key: 'population', label: 'Population & clinical databases', members: ['population', 'clinical'] },
        { key: 'literature', label: 'Literature & visualisation',      members: ['literature', 'visualisation'] },
      ];

      const fmtSnap = x => (x === null || x === undefined ? '—' : Number(x).toFixed(2));
      let snapTiles = [];
      if (vep.ok) {
        const imp = String(vep.impact || '').toUpperCase();
        const tTone = imp === 'HIGH' ? 'warn' : imp === 'MODERATE' ? 'amber' : imp === 'LOW' ? 'info' : 'muted';
        const impLabel = imp ? sCase(imp) + ' impact' : 'Modifier';
        snapTiles.push(`<div class="ev-snap__tile" data-tone="${tTone}"><div class="ev-snap__k">Consequence</div><div class="ev-snap__v is-text">${sCase(vep.most_severe_consequence) || 'Unknown'}</div><div class="ev-snap__sub">${impLabel}</div></div>`);
      }
      {
        const g2 = ev.gnomad || {};
        if (g2.ok && g2.variant) {
          const ex = g2.variant.exome || {}, ge = g2.variant.genome || {};
          const pmv = [ex.faf95 && ex.faf95.popmax, ge.faf95 && ge.faf95.popmax].filter(x => x != null);
          const pmax = pmv.length ? Math.max(...pmv) : null;
          const af2 = (ex.af != null ? ex.af : (ge.af != null ? ge.af : null));
          const exC = ex.ac != null, geC = ge.ac != null;
          const ac2 = (exC ? ex.ac : 0) + (geC ? ge.ac : 0);
          const hom2 = (exC ? (ex.ac_hom || 0) : 0) + (geC ? (ge.ac_hom || 0) : 0);
          const het2 = Math.max(0, ac2 - 2 * hom2);
          const tone2 = pmax != null && pmax >= 0.05 ? 'warn' : pmax != null && pmax >= 0.01 ? 'amber' : 'info';
          const shown = (pmax != null ? pmax : af2);
          snapTiles.push(`<div class="ev-snap__tile" data-tone="${tone2}"><div class="ev-snap__k">gnomAD AF</div><div class="ev-snap__v">${shown != null ? formatAF(shown) : '—'}</div><div class="ev-snap__sub">${het2} het <span class="nb">·</span> ${hom2} hom</div></div>`);
        } else if (g2.ok) {
          snapTiles.push(`<div class="ev-snap__tile" data-tone="good"><div class="ev-snap__k">gnomAD AF</div><div class="ev-snap__v is-text">Absent</div><div class="ev-snap__sub">Not in gnomAD v4</div></div>`);
        }
      }
      {
        const sa2 = ev.spliceai || {};
        if (sa2.ok && sa2.max_delta != null) {
          const md = sa2.max_delta;
          const tone3 = md >= 0.8 ? 'warn' : md >= 0.5 ? 'amber' : md >= 0.2 ? 'muted' : 'good';
          const lab3 = md >= 0.8 ? 'High impact' : md >= 0.5 ? 'Moderate impact' : md >= 0.2 ? 'Low impact' : 'No impact';
          snapTiles.push(`<div class="ev-snap__tile" data-tone="${tone3}"><div class="ev-snap__k">SpliceAI Δ</div><div class="ev-snap__v">${fmtSnap(md)}</div><div class="ev-snap__sub">${lab3}</div></div>`);
        } else if (!sa2.skipped && !sa2.ok && (sa2.lookup_failed || sa2.error || sa2.reason)) {
          snapTiles.push(`<div class="ev-snap__tile" data-tone="muted"><div class="ev-snap__k">SpliceAI Δ</div><div class="ev-snap__v is-text">—</div><div class="ev-snap__sub">No score</div></div>`);
        }
      }
      {
        const am2 = ev.alphamissense || {};
        if (am2 && am2.available !== false && !am2.not_applicable && am2.score != null) {
          const c2 = am2.classification || 'ambiguous';
          const tone4 = c2 === 'likely_pathogenic' ? 'warn' : c2 === 'likely_benign' ? 'good' : 'muted';
          const lab4 = c2 === 'likely_pathogenic' ? 'Likely pathogenic' : c2 === 'likely_benign' ? 'Likely benign' : 'Ambiguous';
          snapTiles.push(`<div class="ev-snap__tile" data-tone="${tone4}"><div class="ev-snap__k">AlphaMissense</div><div class="ev-snap__v">${Number(am2.score).toFixed(2)}</div><div class="ev-snap__sub">${lab4}</div></div>`);
        }
      }
      {
        const cv2 = ev.clinvar || {};
        if (cv2.ok && cv2.found) {
          const recs2 = cv2.records || [];
          const sorted2 = [...recs2].sort((a, b) => (b.stars || 0) - (a.stars || 0));
          const best2 = sorted2[0] || {};
          const totalSubs2 = (cv2.total_submissions != null) ? cv2.total_submissions : recs2.reduce((s, r) => s + (r.number_submitters || 1), 0);
          const isP = r => /pathogenic|likely pathogenic/i.test(r.clinical_significance || '');
          const isB = r => /benign/i.test(r.clinical_significance || '') && !isP(r);
          const isC = r => /conflict/i.test(r.clinical_significance || '') || /conflict/i.test(r.review_status || '');
          let tone5, lab5;
          if (sorted2.some(isP)) { tone5 = 'warn'; lab5 = sorted2.find(isP).clinical_significance.split(/[,/]/)[0].trim(); }
          else if (sorted2.some(isB)) { tone5 = 'good'; lab5 = sorted2.find(isB).clinical_significance.split(/[,/]/)[0].trim(); }
          else if (sorted2.some(isC)) { tone5 = 'amber'; lab5 = 'Conflicting'; }
          else if (/uncertain|vus/i.test(best2.clinical_significance || '')) { tone5 = 'amber'; lab5 = 'Uncertain'; }
          else { tone5 = 'info'; lab5 = `${totalSubs2} submission${totalSubs2 === 1 ? '' : 's'}`; }
          const k5 = Math.max(0, Math.min(4, best2.stars || 0));
          const stars5 = `<span class="ev-stars">${'★'.repeat(k5)}${'☆'.repeat(4 - k5)}</span>`;
          snapTiles.push(`<div class="ev-snap__tile" data-tone="${tone5}"><div class="ev-snap__k">ClinVar</div><div class="ev-snap__v is-text">${lab5}</div><div class="ev-snap__sub">${stars5} <span class="nb">·</span> ${totalSubs2} sub</div></div>`);
        } else if (cv2.ok) {
          snapTiles.push(`<div class="ev-snap__tile" data-tone="muted"><div class="ev-snap__k">ClinVar</div><div class="ev-snap__v is-text">No records</div><div class="ev-snap__sub">Not submitted</div></div>`);
        }
      }

      const snapStrip = snapTiles.length ? `<div class="ev-snap">${snapTiles.join('')}</div>` : '';

      const _ovEsc = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const _ovSym = (vep.gene_symbol || vep.derived_gene_symbol || '').trim();
      const _ovC   = _stripHgvsPrefix(vep.hgvsc || '');
      const _ovP   = _stripHgvsPrefix(vep.hgvsp || '');
      const _ovHead = (_ovSym || _ovC)
        ? `<div class="tab-overview__head">
             <div class="tab-overview__title"><span class="tab-overview__gene">${_ovEsc(_ovSym)}</span>${_ovC ? ' ' + _ovEsc(_ovC) : ''}</div>
             <div class="tab-overview__sub">${_ovP ? _ovEsc(_ovP) : 'Variant annotation &amp; population frequency'}</div>
           </div>`
        : '';

      const escapeAttr = s => String(s || '').replace(/"/g, '&quot;');
      const renderRow = r => {
        const pillHTML = r.pill
          ? `<span class="ev-pill" data-tone="${r.pill.tone}">${r.pill.text}</span>`
          : '<span aria-hidden="true"></span>';
        return `<details class="ev-collapsible ev-acc">
          <summary>
            <span class="ev-row-name">${r.label}</span>
            ${pillHTML}
            <span class="ev-row-finding" title="${escapeAttr(r.summary)}">${r.summary || ''}</span>
            <span class="ev-chevron" aria-hidden="true">▶</span>
          </summary>
          <div class="ev-row-detail ev-acc__detail">${r.body || ''}</div>
        </details>`;
      };
      const groupsHTML = GROUPS.map(grp => {
        const inGroup = rows.filter(r => grp.members.includes(r.group));
        if (!inGroup.length) return '';
        return `<section class="ev-accgroup ev-group" data-group="${grp.key}">
          <div class="v11-eyebrow"><span class="v11-eyebrow__mark" aria-hidden="true">§</span><h3>${grp.label}</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span><button type="button" class="ev-expand-toggle" onclick="toggleAllEvidence(this)" data-state="collapsed">Expand all</button></div>
          ${inGroup.map(renderRow).join('')}
        </section>`;
      }).join('');

      const vepFailed = ev.vep && ev.vep.ok === false;
      const userTxNote = (ev.vep && ev.vep.user_transcript)
        ? ` The user-specified transcript was <span class="ev-mono">${ev.vep.user_transcript}</span>.`
        : '';
      const warningBanner = vepFailed
        ? `<div class="ev-warning-banner" role="alert">
             <svg class="ev-warning-icon" viewBox="0 0 24 24" fill="none" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
               <path d="M12 9v4M12 17h.01M10.29 3.86 1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/>
             </svg>
             <div>
               <span class="ev-warning-title">Variant annotation failed</span>
               — functional consequence could not be determined. Results below may be incomplete.${userTxNote}
             </div>
           </div>`
        : '';

      const extToolsHTML = _buildExternalToolsHTML(ev, 'variant');
      const extToolsSection = extToolsHTML
        ? `<section class="ev-exttools v11-section">
             <div class="v11-eyebrow"><h3>Open external tools</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
             ${extToolsHTML}
           </section>`
        : '';

      return `<div class="rsec">
        ${warningBanner}
        <div class="tab-overview">
          ${_ovHead}
          ${snapStrip}
        </div>
        <div class="evidence-groups">
          ${groupsHTML}
        </div>
        ${extToolsSection}
      </div>`;
    }

    const PROTEIN_FEATURE_COLOURS = {
      'Domain':              '#3b6ea8',
      'Region (disordered)': '#5b8c5a',
      'Region (other)':      '#888888',
      'Transmembrane':       '#7d4ea0',
      'Signal peptide':      '#c9893a',
      'Coiled coil':         '#2e8b7a',
      'Motif':               '#b07a2a',
      'Repeat':              '#a86fc4',
      'Zinc finger':         '#5fa49a',
      'DNA binding':         '#264d80',
      'Active site':         '#b03a2e',
      'Binding site':        '#a13a3a',
    };

    const _DISORDER_KEYWORDS = ['disorder', 'low complexity', 'low-complexity', 'compositional'];

    function _featureCategory(f) {
      if (!f || !f.type) return 'Region (other)';
      if (f.type === 'Region') {
        const d = (f.description || '').toLowerCase();
        return _DISORDER_KEYWORDS.some(k => d.includes(k))
          ? 'Region (disordered)'
          : 'Region (other)';
      }

      if (f.type === 'Signal') return 'Signal peptide';
      return f.type;
    }

    function parseAminoAcidPosition(pNotation) {
      if (!pNotation) return null;
      let s = String(pNotation).trim();
      const colon = s.indexOf(':');
      if (colon >= 0) s = s.slice(colon + 1);

      s = s.replace(/^p\.\(/, 'p.').replace(/\)$/, '');

      let m = s.match(/^p\.([A-Z][a-z]{2})(\d+)/);
      if (m) return parseInt(m[2], 10);

      m = s.match(/^p\.([A-Z*])(\d+)/);
      if (m) return parseInt(m[2], 10);

      console.warn('[heartvar] parseAminoAcidPosition returned null:', pNotation);
      return null;
    }

    const AA3TO1 = {
      ALA:'A', ARG:'R', ASN:'N', ASP:'D', CYS:'C', GLN:'Q', GLU:'E', GLY:'G',
      HIS:'H', ILE:'I', LEU:'L', LYS:'K', MET:'M', PHE:'F', PRO:'P', SER:'S',
      THR:'T', TRP:'W', TYR:'Y', VAL:'V', SEC:'U', PYL:'O',
    };

    function parseAminoAcidRef(pNotation) {
      if (!pNotation) return null;
      let s = String(pNotation).trim();
      const colon = s.indexOf(':');
      if (colon >= 0) s = s.slice(colon + 1);
      s = s.replace(/^p\.\(/, 'p.').replace(/\)$/, '');
      let m = s.match(/^p\.([A-Z][a-z]{2})\d/);
      if (m) return AA3TO1[m[1].toUpperCase()] || null;
      m = s.match(/^p\.([A-Z*])\d/);
      if (m) return m[1] === '*' ? '*' : m[1];
      return null;
    }

    function _hvStructSafeMax(acc) {
      if (!window.__hvStructManifest) {
        window.__hvStructManifest = fetch('/structures/manifest.json')
          .then(r => (r.ok ? r.json() : {}))
          .catch(() => ({}));
      }
      return window.__hvStructManifest.then(m => {
        for (const k in m) {
          const e = m[k];
          if (e && e.accession === acc && typeof e.safe_max_residue === 'number') {
            return e.safe_max_residue;
          }
        }
        return null;
      });
    }

    function _pdbModelAcc(pdb) {
      let title = '';
      const lines = String(pdb).slice(0, 4000).split('\n');
      for (const ln of lines) {
        if (ln.startsWith('TITLE')) title += ln.slice(10).replace(/\s+$/, '');
        else if (ln.startsWith('ATOM') || ln.startsWith('HETATM')) break;
      }
      const m = title.match(/\(([A-Za-z0-9]+(?:-\d+)?)\)\s*$/);
      return m ? m[1] : null;
    }

    function _stripHgvsPrefix(hgvsp) {
      if (!hgvsp) return '';
      const idx = String(hgvsp).indexOf(':');
      return idx >= 0 ? String(hgvsp).slice(idx + 1) : String(hgvsp);
    }

    function _ensureProteinTooltip() {
      if (window.__proteinTooltipReady) return;
      const tip = document.createElement('div');
      tip.className = 'protein-tooltip';
      tip.setAttribute('role', 'tooltip');
      document.body.appendChild(tip);
      document.addEventListener('mousemove', (e) => {
        const target = e.target && e.target.closest
          ? e.target.closest('.protein-track rect.feature')
          : null;
        if (!target || !target.getAttribute('data-tip')) {
          tip.style.display = 'none';
          return;
        }
        tip.textContent = target.getAttribute('data-tip');
        tip.style.display = 'block';

        const OFFSET = 12;
        const r = tip.getBoundingClientRect();
        let left = e.clientX + OFFSET;
        let top  = e.clientY + OFFSET;
        if (left + r.width  > window.innerWidth  - 8) left = e.clientX - r.width  - OFFSET;
        if (top  + r.height > window.innerHeight - 8) top  = e.clientY - r.height - OFFSET;
        tip.style.left = `${Math.max(4, left)}px`;
        tip.style.top  = `${Math.max(4, top)}px`;
      });
      window.__proteinTooltipReady = true;
    }

    function _ensure3Dmol() {
      if (window.$3Dmol) return Promise.resolve();
      if (window.__hv3DmolPromise) return window.__hv3DmolPromise;
      window.__hv3DmolPromise = new Promise((resolve, reject) => {
        const s = document.createElement('script');
        s.src = '/static/vendor/3Dmol-min.js';
        s.async = true;
        s.onload = () => resolve();
        s.onerror = () => reject(new Error('Failed to load 3Dmol viewer'));
        document.head.appendChild(s);
      });
      return window.__hv3DmolPromise;
    }

    const HV_PROBAND_COLOR = '#ffffff';

    const HV_PROBAND_HALO_COLOR = '#26221f';
    const HV_PROBAND_HALO_OPACITY = 0.65;

    const HV_PROBAND_HALO_PAD = 0.60;

    const HV_PROBAND_SPHERE_SCALE = 0.55;

    const HV_PROBAND_LABEL_COLOR = '#26221f';

    const HV_VDW_RADII = { C: 1.7, N: 1.55, O: 1.52, S: 1.8, H: 1.2, P: 1.8 };

    function _hvElem(a) {
      let el = String(a && a.elem || '').toUpperCase();
      if (!el) el = String(a && a.atom || '').replace(/[^A-Za-z]/g, '').charAt(0).toUpperCase();
      return el;
    }

    function _hvDrawProbandHalo(viewer, pos) {
      for (const a of (viewer.selectedAtoms({ resi: pos }) || [])) {
        viewer.addSphere({
          center: { x: a.x, y: a.y, z: a.z },
          radius: (HV_VDW_RADII[_hvElem(a)] || HV_VDW_RADII.C) * HV_PROBAND_SPHERE_SCALE
            + HV_PROBAND_HALO_PAD,
          color: HV_PROBAND_HALO_COLOR,
          opacity: HV_PROBAND_HALO_OPACITY,
        });
      }
    }

    function _hvPolar(a) {
      if (!a) return false;
      let el = String(a.elem || '').toUpperCase();
      if (!el) el = String(a.atom || '').replace(/[^A-Za-z]/g, '').charAt(0).toUpperCase();
      return el === 'N' || el === 'O';
    }

    function _hvRenderNeighbourhood(viewer, pos, labelColor) {
      labelColor = labelColor || HV_PROBAND_LABEL_COLOR;
      const NEIGHBOUR_A = 5.0, HBOND_MAX = 3.5, HBOND_MIN = 2.2;

      const shell = viewer.selectedAtoms({ within: { distance: NEIGHBOUR_A, sel: { resi: pos } } }) || [];
      const neighbourResis = Array.from(new Set(shell.map(a => a.resi))).filter(r => r !== pos);
      if (neighbourResis.length) {

        viewer.addStyle({ resi: neighbourResis }, { stick: { radius: 0.12, colorscheme: 'whiteCarbon' } });
      }

      const BACKBONE = { N: 1, O: 1, OXT: 1 };
      const selfPolar = (viewer.selectedAtoms({ resi: pos }) || [])
        .filter(_hvPolar)
        .filter(a => !BACKBONE[String(a.atom || '').toUpperCase()]);
      const shellPolar = shell.filter(a => a.resi !== pos).filter(_hvPolar);
      let contacts = 0;
      for (const a of selfPolar) {
        for (const b of shellPolar) {
          const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
          const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
          if (d >= HBOND_MIN && d <= HBOND_MAX) {
            viewer.addLine({ dashed: true, color: '#1769aa', start: { x: a.x, y: a.y, z: a.z }, end: { x: b.x, y: b.y, z: b.z } });
            contacts++;
          }
        }
      }

      const self = viewer.selectedAtoms({ resi: pos }) || [];
      if (self.length) {
        const ca = self.find(a => String(a.atom || '').toUpperCase() === 'CA') || self[0];
        const one = AA3TO1[String(self[0].resn || '').toUpperCase()] || '';
        viewer.addLabel((one || String(self[0].resn || '')) + pos, {
          position: { x: ca.x, y: ca.y, z: ca.z },
          inFront: true, fontSize: 12, fontColor: labelColor,
          backgroundColor: 'white', backgroundOpacity: 0.78,
          borderColor: labelColor, borderThickness: 0.5,
        });
      }
      return { neighbours: neighbourResis.length, contacts };
    }

    function _hvOverlayClinVar(viewer, probandPos, isCanonical, safeMax) {
      const list = window.__hv3DClinVarPLP || [];
      const TIERCOL = { P: '#c8341a', LP: '#d96b3a' };

      const tierByResi = {};
      for (const p of list) {
        if (p == null || p.aa == null) continue;
        if (p.tier !== 'P' && p.tier !== 'LP') continue;
        const prev = tierByResi[p.aa];
        tierByResi[p.aa] = (prev === 'P' || p.tier === 'P') ? 'P' : 'LP';
      }
      let plotted = 0, droppedNumbering = 0, droppedModel = 0;
      const resis = Object.keys(tierByResi);
      for (const key of resis) {
        const resi = parseInt(key, 10);

        if (probandPos != null && resi === probandPos) continue;
        const numberingTrusted = isCanonical || (safeMax != null && resi <= safeMax);
        if (!numberingTrusted) { droppedNumbering++; continue; }
        const atoms = viewer.selectedAtoms({ resi: resi });
        if (!atoms || !atoms.length) { droppedModel++; continue; }
        viewer.addStyle({ resi: resi }, { sphere: { color: TIERCOL[tierByResi[key]], scale: 0.42, opacity: 0.65 } });
        plotted++;
      }
      return { plotted: plotted, droppedNumbering: droppedNumbering, droppedModel: droppedModel, total: resis.length };
    }

    function _hv3DClinVarResiCount() {
      const uniq = {};
      for (const p of (window.__hv3DClinVarPLP || [])) { if (p && p.aa != null) uniq[p.aa] = 1; }
      return Object.keys(uniq).length;
    }

    window.__hvToggle3DClinVar = function () {
      try {
        const st = window.__hv3DState;
        if (!st || typeof st.drawScene !== 'function') return;
        st.show = !st.show;
        const res = st.drawScene(st.show);
        const btn = document.getElementById('hv3d-cv-toggle');
        if (btn) {
          btn.textContent = st.show ? 'Hide ClinVar P/LP' : ('Show ClinVar P/LP (' + _hv3DClinVarResiCount() + ')');
          btn.style.background = st.show ? 'var(--maroon,#8c1a1f)' : '#fff';
          btn.style.color = st.show ? '#fff' : 'var(--ink-soft,#6a5d54)';
        }
        const note = document.getElementById('hv3d-cv-note');
        if (note) {
          const o = res && res.overlay;
          if (st.show && o) {
            const bits = [o.plotted + ' P/LP residue' + (o.plotted === 1 ? '' : 's') + ' shown'];
            if (st.probandTier) bits.push('your variant is itself ClinVar ' + (st.probandTier === 'P' ? 'pathogenic — shown red' : 'likely-pathogenic — shown orange'));
            const dropped = (o.droppedNumbering || 0) + (o.droppedModel || 0);
            if (dropped) bits.push(dropped + ' omitted (outside the model’s reliably-numbered region)');
            note.innerHTML = bits.join(' · ');
            note.style.display = '';
          } else {
            note.innerHTML = '';
            note.style.display = 'none';
          }
        }
      } catch (e) {  }
    };

    window.__hvLoad3D = function (collEl) {
      try {
        if (!collEl || !collEl.classList.contains('open')) return;
        const mount = collEl.querySelector('.hv3d-mount');
        if (!mount || mount.dataset.init) return;
        mount.dataset.init = '1';
        const acc = mount.dataset.acc;
        const pos = mount.dataset.pos ? parseInt(mount.dataset.pos, 10) : null;
        const refAA = mount.dataset.refaa || null;
        const status = mount.querySelector('.hv3d-status');
        const setStatus = (html, show) => {
          if (!status) return;
          status.style.display = show === false ? 'none' : '';
          status.innerHTML = html;
        };
        const linksHTML = () => {
          const afUrl = 'https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(acc);
          return '<a href="' + afUrl + '" target="_blank" rel="noopener" style="color:var(--blue)">Open in AlphaFold ↗</a>';
        };
        const fallback = (msg) => {
          setStatus((msg || 'No bundled 3-D model for this protein.') + '<br>' + linksHTML());
        };

        const flagNote = (html) => {
          const note = document.createElement('div');
          note.style.cssText = 'position:absolute;left:8px;right:8px;bottom:8px;background:rgba(255,255,255,.95);padding:6px 10px;border-radius:6px;font-size:11.5px;line-height:1.5;color:var(--ink-2);border:1px solid var(--line,#e6ddd4)';
          note.innerHTML = html;
          mount.appendChild(note);
        };
        setStatus('Loading 3-D structure…');

        Promise.all([
          _ensure3Dmol()
            .then(() => fetch('/structures/' + encodeURIComponent(acc) + '.pdb'))
            .then(r => (r.ok ? r.text() : null)),
          _hvStructSafeMax(acc),
        ])
          .then(([pdb, safeMax]) => {
            if (pdb == null) { fallback(); return; }
            setStatus('', false);

            const modelAcc = _pdbModelAcc(pdb);
            const isCanonical = !!modelAcc && modelAcc === acc;
            const isIsoform = !!modelAcc && modelAcc !== acc;
            const viewer = window.$3Dmol.createViewer(mount, { backgroundColor: 'white' });
            viewer.addModel(pdb, 'pdb');

            if (isIsoform) {
              const cap = mount.parentElement && mount.parentElement.querySelector('.gc-bar-cap-note');
              if (cap) {
                const reg = (safeMax != null) ? ('canonical residues 1–' + safeMax) : 'an unverified region';
                cap.innerHTML = '<span style="color:var(--amber)">⚠ AlphaFold has no model of this protein’s '
                  + 'canonical sequence; the structure shown is isoform <b>' + modelAcc + '</b>, reliable only for '
                  + reg + '.</span><br>' + cap.innerHTML;
              }
            }

            let inModel = false, modelAA = null, aaOk = null, numberingTrusted = false, pinOk = false;
            if (pos != null) {
              const atoms = viewer.selectedAtoms({ resi: pos });
              inModel = !!(atoms && atoms.length);
              modelAA = inModel ? (AA3TO1[String(atoms[0].resn || '').toUpperCase()] || null) : null;

              aaOk = (refAA && modelAA) ? (modelAA === refAA) : null;

              numberingTrusted = isCanonical || (safeMax != null && pos <= safeMax);
              pinOk = inModel && numberingTrusted && aaOk !== false;
            }

            let probandTier = null;
            if (pos != null) {
              for (const p of (window.__hv3DClinVarPLP || [])) {
                if (p && p.aa === pos && (p.tier === 'P' || p.tier === 'LP')) {
                  probandTier = (probandTier === 'P' || p.tier === 'P') ? 'P' : 'LP';
                }
              }
            }

            const drawScene = (showClinVar) => {
              viewer.removeAllShapes();
              viewer.removeAllLabels();

              viewer.setStyle({}, { cartoon: { colorscheme: { prop: 'b', gradient: 'roygb', min: 50, max: 90 } } });
              let env = { neighbours: 0, contacts: 0 };
              if (pinOk) {

                const probColor = (showClinVar && probandTier)
                  ? (probandTier === 'P' ? '#c8341a' : '#d96b3a')
                  : HV_PROBAND_COLOR;

                try { _hvDrawProbandHalo(viewer, pos); } catch (e) {  }
                viewer.addStyle({ resi: pos }, { stick: { color: probColor, radius: 0.35 } });
                viewer.addStyle({ resi: pos }, { sphere: { color: probColor, scale: HV_PROBAND_SPHERE_SCALE } });

                const probLabelColor = (showClinVar && probandTier)
                  ? probColor
                  : HV_PROBAND_LABEL_COLOR;
                try { env = _hvRenderNeighbourhood(viewer, pos, probLabelColor) || env; } catch (e) {  }
              }
              let overlay = null;
              if (showClinVar) {
                try { overlay = _hvOverlayClinVar(viewer, pos, isCanonical, safeMax); } catch (e) {  }
              }
              viewer.render();
              viewer.resize();
              return { env: env, overlay: overlay };
            };

            const _init = drawScene(false);

            if (pinOk) { viewer.zoomTo({ resi: pos }); viewer.zoom(0.6); }
            else { viewer.zoomTo(); }
            viewer.render();

            window.__hv3DState = { drawScene: drawScene, show: false, probandTier: probandTier };

            if (pinOk) {
              const _cap = mount.parentElement && mount.parentElement.querySelector('.gc-bar-cap-note');
              if (_cap && !_cap.dataset.envNoted && (_init.env.neighbours || _init.env.contacts)) {
                _cap.dataset.envNoted = '1';
                const _bits = [];
                if (_init.env.neighbours) _bits.push('neighbouring residues within ~5 Å are shown as thin sticks');
                if (_init.env.contacts) _bits.push(_init.env.contacts + ' candidate side-chain polar contact' + (_init.env.contacts === 1 ? '' : 's') + ' shown as dashed blue line' + (_init.env.contacts === 1 ? '' : 's'));
                _cap.innerHTML += '<br><span style="color:var(--ink-2)">' + _bits.join('; ').replace(/^./, c => c.toUpperCase())
                  + '. Distance-based, from this single static wild-type model — not the mutant structure.</span>';
              }
            }

            if (pos != null && !pinOk) {
              let why;
              if (aaOk === false) {
                why = '⚠ The model’s residue ' + pos + ' (' + (modelAA || '?')
                  + ') doesn’t match the variant’s reference amino acid (' + refAA
                  + '), so the numbering can’t be confirmed — not highlighting, to avoid showing the wrong residue.';
              } else if (!numberingTrusted) {
                const reg = (safeMax != null) ? ('residues 1–' + safeMax) : 'an unverified region';
                why = '⚠ This model’s numbering is reliable only for ' + reg
                  + ', so residue ' + pos + ' (your variant) can’t be located reliably here.';
              } else {
                why = '⚠ Residue ' + pos + ' (your variant) isn’t present in this AlphaFold model, '
                  + 'so its position can’t be shown here.';
              }
              flagNote(why + '<br>' + linksHTML());
            }
          })
          .catch(() => fallback('3-D viewer could not load.'));
      } catch (e) {  }
    };

    function _renderProteinTab(evidence, gene, hgvs_c) {

      const store = window.heartvarProteinData || {};
      const evVep = (evidence && evidence.vep) || {};
      const evUp  = (evidence && evidence.uniprot) || {};

      gene = (gene || evVep.gene_symbol || evVep.derived_gene_symbol || '').trim();

      if (!store.hgvsp && !evVep.hgvsp) {
        console.warn('Protein tab: hgvsp missing from heartvarProteinData', {
          gene: store.gene || gene,
          hgvs_c,
          store_hgvsp: store.hgvsp,
          evidence_hgvsp: evVep.hgvsp,
        });
      }
      if (!store.uniprotFeatures && !Array.isArray(evUp.features)) {
        console.warn('Protein tab: uniprotFeatures missing from heartvarProteinData', {
          gene: store.gene || gene,
          store_features: store.uniprotFeatures,
          evidence_features: evUp.features,
        });
      }

      const rawHgvsp    = store.hgvsp    || evVep.hgvsp    || '';
      const hgvsP       = _stripHgvsPrefix(rawHgvsp);
      const aaPos       = parseAminoAcidPosition(rawHgvsp);
      const seqLen      = store.uniprotLength      || evUp.length || 0;
      const features    = (store.uniprotFeatures && store.uniprotFeatures.slice())
                       || (Array.isArray(evUp.features) ? evUp.features.slice() : []);
      const proteinName = store.uniprotProteinName || evUp.protein_name || '';
      const accession   = store.uniprotAccession   || evUp.accession    || '';
      const uniprotOk   = store.uniprotOk || !!evUp.ok;
      const uniprotUrl  = store.uniprotUrl
        || evUp.url
        || (accession ? `https://www.uniprot.org/uniprotkb/${encodeURIComponent(accession)}` : '');
      const up = { ok: uniprotOk, length: seqLen, features };

      const _bgEsc = s => String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const _bg = (evidence && evidence.biogrid) || {};
      let _bgBody = '';
      if (_bg.ok && _bg.total_interactions) {
        const _bgParts = (_bg.top_partners || []).map(p =>
          `<span class="${p.is_chdgene ? 'ev-partner-chd' : 'ev-partner'}">${_bgEsc(p.symbol)}</span> <span class="ev-partner-count">(n=${p.publication_count})</span>`
        ).join(', ') || '<span class="ev-dim">—</span>';
        _bgBody = `<div style="font-size:14px;color:var(--ink);margin-bottom:8px"><b>${_bg.total_interactions}</b> curated interactions <span class="ev-dim">·</span> <b>${_bg.unique_partners}</b> unique partners</div>
          <div style="font-size:13.5px;line-height:1.7;color:var(--ink-2)"><span style="font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:var(--muted);font-weight:600;margin-right:8px">Top partners</span>${_bgParts}</div>
          ${_bg.url ? `<a class="pv2-id__link" style="margin-top:12px;display:inline-block" href="${_bg.url}" target="_blank" rel="noopener">Open in BioGRID ↗</a>` : ''}`;
      } else if (_bg.ok) {
        _bgBody = '<div class="protein-empty">No curated interactions reported.</div>';
      } else if (_bg.error) {
        _bgBody = `<div class="protein-empty">Lookup failed: ${_bgEsc(_bg.error)}</div>`;
      }
      const _bgSection = _bgBody
        ? `<div class="pv2-section"><div class="v11-eyebrow"><h3>Protein interactions</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div><div class="pv2-trackcard" style="padding:18px 20px">${_bgBody}</div></div>`
        : '';

      const _protColl = (title, accent, body, opts) => {
        if (!body) return '';
        const open = (opts && opts.open) ? ' open' : '';
        const summary = (opts && opts.summary) || '';
        return `<div class="gc-collapsible${open}" style="--gc-accent:${accent};">
          <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open')">
            <span class="gc-coll-title">${title}</span>
            <span class="gc-coll-summary">${summary}</span>
            <span class="gc-coll-chevron" aria-hidden="true">▶</span>
          </button>
          <div class="gc-coll-body"><div class="gc-coll-inner">${body}</div></div>
        </div>`;
      };

      const _protSectionsHead = `<div class="ev-snap-head">
            <div class="v11-eyebrow"><h3>Annotations</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
            <button type="button" class="ev-expand-toggle" onclick="toggleAllGeneSections(this)" data-state="collapsed">Expand all</button>
          </div>`;

      const _protExtHTML = _buildExternalToolsHTML(evidence, 'protein');
      const _protExtSection = _protExtHTML
        ? `<section class="ev-exttools v11-section"><div class="v11-eyebrow"><h3>Open external tools</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>${_protExtHTML}</section>`
        : '';

      if (window.__hvDebug) {
        console.log('[protein-tab]', {
          gene, hgvs_c,
          vep_hgvsp_raw: rawHgvsp || null,
          parsed_aa_pos: aaPos,
          seq_len: seqLen,
          uniprot_ok: uniprotOk,
          feature_count: features.length,
          store_source: store.hgvsp ? 'window.heartvarProteinData' : 'evidence param',
        });
      }

      const escAttr = s => String(s == null ? '' : s).replace(/"/g, '&quot;');

      const escText = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');

      if (!up.ok || !seqLen) {
        const _lc = seqLen ? `${escText(seqLen)} aa` : '<span class="pv2-stat__v is-dim">—</span>';
        const _vc = hgvsP ? escText(hgvsP) : '<span class="pv2-stat__v is-dim">No protein consequence</span>';
        const _lk = uniprotUrl ? `<a class="pv2-id__link" href="${uniprotUrl}" target="_blank" rel="noopener">UniProt ${escText(accession || '')} ↗</a>` : '';
        return `<div class="protein-tab">
          <div class="pv2-id">
            <div class="pv2-id__head">
              <div class="pv2-id__title"><span class="pv2-id__gene">${escText(gene || 'Protein')}</span>${proteinName ? ' · ' + escText(proteinName) : ''}</div>
              <div class="pv2-id__sub">UniProt protein identity &amp; domain architecture</div>
            </div>
            ${_lk}
            <div class="pv2-stats">
              <div class="pv2-stat"><div class="pv2-stat__k">Length of protein</div><div class="pv2-stat__v">${_lc}</div></div>
              <div class="pv2-stat pv2-stat--variant"><div class="pv2-stat__k">Variant</div><div class="pv2-stat__v">${_vc}</div></div>
              <div class="pv2-stat"><div class="pv2-stat__k">Domain context</div><div class="pv2-stat__v is-dim">—</div></div>
            </div>
          </div>
          ${_protSectionsHead}
          <div class="gcx-acc pv2-acc">
            ${_protColl('Domain architecture', 'var(--blue)', `<div class="pv2-trackcard"><div class="protein-empty">No domain annotations available for this protein.</div></div>`)}
            ${_protColl('Domain pathogenicity context', 'var(--maroon)', _renderDomainPlpSection())}
            ${_protColl('Protein interactions', 'var(--green)', _bgBody)}
          </div>
          ${_protExtSection}
        </div>`;
      }

      const W = 900, H = 110;
      const PAD_L = 30, PAD_R = 20;
      const TRACK_Y = 44, TRACK_H = 24;
      const AXIS_Y = 84;
      const innerW = W - PAD_L - PAD_R;
      const xFor = (pos) => PAD_L + ((pos - 1) / Math.max(1, seqLen - 1)) * innerW;

      features.sort((a, b) => (b.end - b.start) - (a.end - a.start));

      const containingFeatures = (aaPos == null) ? [] :
        features.filter(f => f.start <= aaPos && aaPos <= f.end);

      const LABEL_CH_PX = 5.6;
      const LABEL_PAD_PX = 6;
      const LABEL_MIN_BLOCK_PX = 50;

      const rects = features.map(f => {
        const x1 = xFor(f.start);
        const x2raw = xFor(Math.max(f.start, f.end));

        const w = Math.max(3, x2raw - x1 - 1);
        const category = _featureCategory(f);
        const colour = PROTEIN_FEATURE_COLOURS[category] || '#888';
        const lenAa = (f.end - f.start + 1);
        const containsVariant = containingFeatures.includes(f);
        const tip = `${category}${f.description ? ': ' + f.description : ''} (${f.start}–${f.end}, ${lenAa} aa)`;

        let labelTxt = '';
        if (f.description && w >= LABEL_MIN_BLOCK_PX) {
          const maxChars = Math.max(1, Math.floor((w - LABEL_PAD_PX) / LABEL_CH_PX));
          if (f.description.length <= maxChars) {
            labelTxt = f.description;
          }
        }
        return `
          <g>
            <rect class="feature${containsVariant ? ' contains-variant' : ''}"
                  x="${x1.toFixed(2)}" y="${TRACK_Y}"
                  width="${w.toFixed(2)}" height="${TRACK_H}"
                  rx="2" ry="2"
                  fill="${colour}"
                  data-tip="${escAttr(tip)}"
                  aria-label="${escAttr(tip)}"></rect>
            ${labelTxt
              ? `<text class="feat-label" x="${(x1 + w / 2).toFixed(2)}" y="${TRACK_Y + TRACK_H / 2 + 3}" text-anchor="middle" pointer-events="none">${escText(labelTxt)}</text>`
              : ''}
          </g>`;
      }).join('');

      _ensureProteinTooltip();

      const pinSvg = (aaPos != null && aaPos >= 1 && aaPos <= seqLen)
        ? (() => {
            const x = xFor(aaPos);
            const labelText = hgvsP || '';

            const labelAnchor = (x > W - PAD_R - 80) ? 'end' : 'middle';

            const labelY = 11;
            const circleCY = 28;
            const circleR = 4;
            const connectorTop = labelY + 3;
            const connectorBot = circleCY - circleR;
            const pinBodyTop = circleCY + circleR;
            const pinBodyBot = TRACK_Y;
            return `
              <text class="pin-label" x="${x.toFixed(2)}" y="${labelY}" text-anchor="${labelAnchor}">${escText(labelText)}</text>
              <line class="pin-connector" x1="${x.toFixed(2)}" y1="${connectorTop}" x2="${x.toFixed(2)}" y2="${connectorBot.toFixed(2)}"/>
              <line class="pin-line" x1="${x.toFixed(2)}" y1="${pinBodyTop.toFixed(2)}" x2="${x.toFixed(2)}" y2="${pinBodyBot}"/>
              <circle class="pin-head" cx="${x.toFixed(2)}" cy="${circleCY}" r="${circleR}"/>`;
          })()
        : '';

      const axis = `
        <text class="axis-tick" x="${PAD_L}" y="${AXIS_Y}" text-anchor="start">1</text>
        <text class="axis-tick" x="${W - PAD_R}" y="${AXIS_Y}" text-anchor="end">${seqLen} aa</text>`;

      const trackBar = `<rect class="track-bar" x="${PAD_L}" y="${TRACK_Y}" width="${innerW}" height="${TRACK_H}" rx="3" ry="3"/>`;

      const svg = `
        <svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Protein domain track for ${escAttr(gene)} with variant pin">
          ${trackBar}
          ${rects}
          ${pinSvg}
          ${axis}
        </svg>`;

      const callout = containingFeatures.length
        ? `<div class="protein-callout">Variant falls within ${escText(containingFeatures.map(f => {
              const _c = _featureCategory(f);
              return _c === 'Region (disordered)' ? 'a disordered region' : (f.description ? `${f.description} (${_c.toLowerCase()})` : _c);
            }).join(', '))}.</div>`
        : (aaPos != null && features.length
            ? (() => {

                const before = features
                  .filter(f => f.end < aaPos)
                  .sort((a, b) => b.end - a.end)[0];
                const after = features
                  .filter(f => f.start > aaPos)
                  .sort((a, b) => a.start - b.start)[0];
                const label = f => f.description ? f.description : _featureCategory(f);
                if (before && after) return `<div class="protein-callout" style="color:var(--text-tertiary);font-weight:400">Variant sits between ${escText(label(before))} and ${escText(label(after))}.</div>`;
                if (after) return `<div class="protein-callout" style="color:var(--text-tertiary);font-weight:400">Variant lies in the N-terminal region (before ${escText(label(after))}).</div>`;
                if (before) return `<div class="protein-callout" style="color:var(--text-tertiary);font-weight:400">Variant lies in the C-terminal region (after ${escText(label(before))}).</div>`;
                return '';
              })()
            : '');

      const PROTEIN_LEGEND_ORDER = [
        'Domain', 'Region (disordered)', 'Region (other)',
        'Zinc finger', 'Active site', 'Binding site', 'Motif',
      ];
      const presentCategories = new Set(features.map(_featureCategory));
      const legend = `<div class="protein-legend">${
        PROTEIN_LEGEND_ORDER.map(c => {
          const cls = presentCategories.has(c)
            ? 'protein-legend-item'
            : 'protein-legend-item is-absent';
          const colour = PROTEIN_FEATURE_COLOURS[c] || '#888';
          return `<span class="${cls}"><span class="protein-legend-swatch" style="background:${colour}"></span>${c}</span>`;
        }).join('')
      }</div>`;

      const domainContext = containingFeatures.length
        ? containingFeatures.map(f => f.description || _featureCategory(f)).join(', ')
        : (aaPos != null
            ? (() => {
                if (!features.length) return 'Unannotated region';
                const before = features.filter(f => f.end < aaPos).sort((a, b) => b.end - a.end)[0];
                const after  = features.filter(f => f.start > aaPos).sort((a, b) => a.start - b.start)[0];
                const label = f => f.description || _featureCategory(f);
                if (before && after) return `Unannotated region (between ${label(before)} and ${label(after)})`;
                if (after)  return `Unannotated region (N-terminal, before ${label(after)})`;
                if (before) return `Unannotated region (C-terminal, after ${label(before)})`;
                return 'Unannotated region';
              })()
            : '—');

      const _lc = seqLen ? `${escText(seqLen)} aa` : '<span class="pv2-stat__v is-dim">—</span>';
      const _vc = hgvsP ? escText(hgvsP) : '<span class="pv2-stat__v is-dim">No protein consequence</span>';
      const _lk = uniprotUrl ? `<a class="pv2-id__link" href="${uniprotUrl}" target="_blank" rel="noopener">UniProt ${escText(accession || '')} ↗</a>` : '';
      const _idHeader = `
        <div class="pv2-id">
          <div class="pv2-id__head">
            <div class="pv2-id__title"><span class="pv2-id__gene">${escText(gene || 'Protein')}</span>${proteinName ? ' · ' + escText(proteinName) : ''}</div>
            <div class="pv2-id__sub">UniProt protein identity &amp; domain architecture</div>
          </div>
          ${_lk}
          <div class="pv2-stats">
            <div class="pv2-stat"><div class="pv2-stat__k">Length of protein</div><div class="pv2-stat__v">${_lc}</div></div>
            <div class="pv2-stat pv2-stat--variant"><div class="pv2-stat__k">Variant</div><div class="pv2-stat__v">${_vc}</div></div>
            <div class="pv2-stat"><div class="pv2-stat__k">Domain context</div><div class="pv2-stat__v">${escText(domainContext || '—')}</div></div>
          </div>
        </div>`;
      const _inFeat = containingFeatures.length > 0;
      const _calloutBody = _inFeat
        ? `Variant falls within <b>${escText(containingFeatures.map(f => { const _c = _featureCategory(f); return _c === 'Region (disordered)' ? 'a disordered region' : (f.description ? `${f.description} (${_c.toLowerCase()})` : _c); }).join(', '))}</b>.`
        : (aaPos != null && features.length
            ? (() => {
                const before = features.filter(f => f.end < aaPos).sort((a, b) => b.end - a.end)[0];
                const after  = features.filter(f => f.start > aaPos).sort((a, b) => a.start - b.start)[0];
                const lab = f => f.description ? f.description : _featureCategory(f);
                if (before && after) return `Variant sits between <b>${escText(lab(before))}</b> and <b>${escText(lab(after))}</b> (unannotated region).`;
                if (after)  return `Variant lies in the N-terminal region, before <b>${escText(lab(after))}</b>.`;
                if (before) return `Variant lies in the C-terminal region, after <b>${escText(lab(before))}</b>.`;
                return '';
              })()
            : (hgvsP ? '' : 'No protein-level position — splice / intronic / synonymous variant.'));
      const _calloutBand = _calloutBody
        ? `<div class="pv2-callout${_inFeat ? ' is-hit' : ''}"><span class="pv2-callout__icon">${_inFeat ? '◆' : '◇'}</span><span>${_calloutBody}</span></div>`
        : '';

      const _modelNote = _proteinModelLabel({ accession: accession, length: seqLen });
      const _archBody = `<div class="pv2-trackcard">
            <div class="protein-track">${svg}</div>
            ${_modelNote ? `<div class="hv-model-note">${_modelNote}</div>` : ''}
            ${_calloutBand}
            ${legend}
          </div>`;

      const _lsPositions = (evidence && evidence.clinvar_gene_landscape
        && Array.isArray(evidence.clinvar_gene_landscape.positions))
        ? evidence.clinvar_gene_landscape.positions : [];
      const _plpPos = _lsPositions.filter(p => p && (p.tier === 'P' || p.tier === 'LP') && p.aa != null);
      window.__hv3DClinVarPLP = _plpPos;
      const _plpResiSet = {};
      for (const _p of _plpPos) _plpResiSet[_p.aa] = 1;
      const _plpResiCount = Object.keys(_plpResiSet).length;
      const _plpToggle = (accession && _plpResiCount) ? `
            <div class="hv3d-cv-controls" style="display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-top:8px">
              <button type="button" id="hv3d-cv-toggle" onclick="window.__hvToggle3DClinVar && window.__hvToggle3DClinVar()" style="font-size:11.5px;padding:3px 10px;border-radius:999px;border:1px solid var(--line,#e6ddd4);cursor:pointer;background:#fff;color:var(--ink-soft,#6a5d54)">Show ClinVar P/LP (${_plpResiCount})</button>
              <span id="hv3d-cv-note" style="display:none;font-size:11px;color:var(--ink-2)"></span>
            </div>
            <div class="gc-bar-cap-note" style="margin-top:4px">Overlays every residue carrying a ClinVar <b style="color:#c8341a">pathogenic</b> / <b style="color:#d96b3a">likely-pathogenic</b> missense variant (canonical-transcript numbering) as a translucent sphere, so you can see whether this variant sits within a 3-D cluster of known pathogenic residues. Your queried variant is highlighted and labelled: it is shown as a <b>pale sphere ringed in charcoal</b> by default, turning <b style="color:#c8341a">red</b> / <b style="color:#d96b3a">orange</b> when the overlay is on if it is itself a ClinVar pathogenic / likely-pathogenic residue — the ring stays either way, so the queried residue is always the ringed one. Visual aid only — it does not feed the ACMG classification.</div>` : '';
      const _3dSection = accession ? `<div class="gc-collapsible" style="--gc-accent: var(--green);">
          <button type="button" class="gc-coll-header" onclick="this.parentElement.classList.toggle('open'); window.__hvLoad3D && window.__hvLoad3D(this.parentElement);">
            <span class="gc-coll-title">3-D structure</span>
            <span class="gc-coll-summary"></span>
            <span class="gc-coll-chevron" aria-hidden="true">▶</span>
          </button>
          <div class="gc-coll-body"><div class="gc-coll-inner">
            <div class="hv3d-mount" data-acc="${escText(accession)}" data-pos="${aaPos != null ? aaPos : ''}" data-refaa="${escAttr(parseAminoAcidRef(rawHgvsp) || '')}" style="width:100%;height:360px;position:relative;background:#faf7f3;border:1px solid var(--line,#e6ddd4);border-radius:8px;overflow:hidden">
              <div class="hv3d-status" style="position:absolute;inset:0;display:flex;align-items:center;justify-content:center;text-align:center;padding:16px;color:var(--muted);font-size:13px;line-height:1.6">Expand to load the interactive 3-D structure.</div>
            </div>
            ${_modelNote ? `<div class="hv-model-note">${_modelNote}</div>` : ''}
            <div class="gc-bar-cap-note" style="margin-top:6px">AlphaFold predicted model (DeepMind/EBI · CC-BY 4.0), coloured by pLDDT confidence (blue = high, red = low).${aaPos != null ? ` The altered residue <b>p.${aaPos}</b> is highlighted and labelled when it lies within the model’s reliably-numbered region (otherwise a note explains why).` : ''} Drag to rotate · scroll to zoom.</div>
            ${_plpToggle}
          </div></div>
        </div>` : '';

      return `<div class="protein-tab">
        ${_idHeader}
        ${_protSectionsHead}
        <div class="gcx-acc pv2-acc">
          ${_protColl('Domain architecture', 'var(--blue)', _archBody)}
          ${_3dSection}
          ${_protColl('Domain pathogenicity context', 'var(--maroon)', _renderDomainPlpSection())}
          ${_protColl('Protein interactions', 'var(--green)', _bgBody)}
        </div>
        ${_protExtSection}
      </div>`;
    }

    function _renderDomainPlpSection() {
      const dp = (window.heartvarProteinData && window.heartvarProteinData.domainPlp) || null;
      if (!dp) return '';
      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const escAttr = s => String(s == null ? '' : s).replace(/"/g, '&quot;');

      if (dp.ok === false) {
        return `<div class="protein-domain-plp">
          <div class="v11-eyebrow"><h3>Domain pathogenicity context</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
          <div class="protein-domain-plp-empty">Lookup failed: ${esc(dp.error || 'unknown error')}</div>
        </div>`;
      }
      if (dp.not_applicable) {
        return `<div class="protein-domain-plp">
          <div class="v11-eyebrow"><h3>Domain pathogenicity context</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
          <div class="protein-domain-plp-empty">${esc(dp.pm1_assessment || 'PM1 not applicable.')}</div>
        </div>`;
      }

      const cP  = Number(dp.count_P)  || 0;
      const cLP = Number(dp.count_LP) || 0;
      const total = cP + cLP;
      const dom = dp.domain_name || '(unnamed domain)';
      const dStart = dp.domain_start;
      const dEnd = dp.domain_end;
      const headerLine = `ClinVar P/LP variants in ${esc(dom)} (aa ${esc(dStart)}–${esc(dEnd)})`;

      if (total === 0) {
        return `<div class="protein-domain-plp">
          <div class="v11-eyebrow"><h3>Domain pathogenicity context</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
          <div class="protein-domain-plp-sub">${headerLine}</div>
          <div class="protein-domain-plp-empty">0 P/LP variants found in ${esc(dom)} — domain is not an established mutational hotspot.</div>
          ${_renderPdpPm1(dp)}
        </div>`;
      }

      const segs = [
        { count: cP,  label: 'P',  cls: 'gc-seg-p',  swatch: '#c8341a' },
        { count: cLP, label: 'LP', cls: 'gc-seg-lp', swatch: '#d96b3a' },
      ];
      const bar = segs.map(s => {
        if (!s.count) return '';
        const pct = (s.count / total) * 100;
        const txt = pct >= 12 ? `${s.label} ${s.count}` : '';
        return `<div class="gc-bar-seg ${s.cls}" style="flex:${s.count}" title="${s.label}: ${s.count} (${pct.toFixed(1)}%)">${txt}</div>`;
      }).join('');
      const legend = segs.map(s => `<span>
        <span class="gc-bar-legend-swatch" style="background:${s.swatch}"></span>${s.label} (${s.count})
      </span>`).join('');

      const variantRows = (dp.top_variants || []).map(v => {
        const stars = Number(v.stars) || 0;
        const starStr = '★'.repeat(stars) + '☆'.repeat(Math.max(0, 4 - stars));
        const tierCls = v.tier === 'P' ? 'is-p' : 'is-lp';
        const tier = esc(v.tier || '');
        const name = esc(v.name || '(unnamed)');
        const conds = (v.conditions && v.conditions.length)
          ? v.conditions.slice(0, 2).join('; ')
          : '';

        const vid = v.variation_id;
        const url = vid
          ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(vid)}/`
          : null;
        const inner = `
          <span class="pdp-stars" title="${stars}-star review"><span>${'★'.repeat(stars)}</span><span class="pdp-stars-empty">${'☆'.repeat(Math.max(0, 4 - stars))}</span></span>
          <span class="pdp-tier ${tierCls}">${tier}</span>
          <span class="pdp-name" title="${escAttr(v.name || '')}">${name}</span>
          ${conds ? `<span class="pdp-conditions" title="${escAttr(conds)}">${esc(conds)}</span>` : ''}`;
        return url
          ? `<a class="pdp-variant-row" href="${url}" target="_blank" rel="noopener">${inner}</a>`
          : `<div class="pdp-variant-row">${inner}</div>`;
      }).join('');
      const variantsBlock = variantRows
        ? `<div class="pdp-variants">${variantRows}</div>` : '';

      return `<div class="protein-domain-plp">
        <div class="v11-eyebrow"><h3>Domain pathogenicity context</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
        <div class="protein-domain-plp-sub">${headerLine}</div>
        <div class="gc-bar">${bar}</div>
        <div class="gc-bar-legend">${legend}</div>
        ${variantsBlock}
        ${_renderPdpPm1(dp)}
      </div>`;
    }

    function _renderPdpPm1(dp) {

      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
      const status = dp.pm1_status || 'not_applicable';
      const sentence = dp.pm1_assessment || '';
      const iconMap = { met: '✓', supporting: '⚡', not_met: '✗', not_applicable: '·' };
      const icon = iconMap[status] || '·';
      return `<div class="pdp-pm1 is-${esc(status)}">
        <span class="pdp-pm1-icon">${icon}</span>
        <span>${esc(sentence)}</span>
      </div>`;
    }

    function _buildExternalToolsHTML(ev, scope) {
      ev = ev || {};
      scope = scope || 'all';
      const enc2 = encodeURIComponent;
      const vepEv = ev.vep || {};

      const extGene = vepEv.gene_symbol
        || ((vepEv.queried_as || '').split(/[\s:]/)[0])
        || '';
      const extHgvsC = vepEv.hgvsc || '';

      let extChr = vepEv.seq_region_name || '';
      let extPos = vepEv.start || '';
      let extRef = '';
      let extAlt = '';
      const fwdId = vepEv.forward_variant_id || '';
      if (fwdId) {
        const parts = fwdId.split('-');
        if (parts.length === 4) {
          extChr = parts[0];
          extPos = parts[1];
          extRef = parts[2];
          extAlt = parts[3];
        }
      } else {
        const _COMP = { A: 'T', T: 'A', C: 'G', G: 'C', a: 't', t: 'a', c: 'g', g: 'c' };
        const complementAllele = a => (a || '').split('').reverse().map(c => _COMP[c] || c).join('');
        const alleleParts = (vepEv.allele_string || '').split('/');
        extRef = alleleParts[0] || '';
        extAlt = alleleParts[1] || '';
        if (vepEv.strand === -1) {
          extRef = complementAllele(extRef);
          extAlt = complementAllele(extAlt);
        }
      }
      const haveCoords = extChr && extPos && extRef && extAlt;
      const haveGene = !!extGene;

      const uniprotAcc = (ev.uniprot && ev.uniprot.accession) || '';

      const variant = [];
      const gene = [];
      const protein = [];

      if (haveGene) {

        const cvChange = (extHgvsC || '').replace(/^[^:]*:/, '');
        const cvTerm = cvChange
          ? enc2(`${extGene}[gene] AND ${cvChange}`)
          : enc2(`${extGene}[gene]`);
        variant.push({ label: 'ClinVar', url: `https://www.ncbi.nlm.nih.gov/clinvar/?term=${cvTerm}` });
      }
      if (haveCoords) {
        variant.push({
          label: 'gnomAD',
          url: `https://gnomad.broadinstitute.org/variant/${extChr}-${extPos}-${extRef}-${extAlt}?dataset=gnomad_r4`,
        });
      }
      if (haveCoords) {
        variant.push({
          label: 'SpliceAI Lookup',
          title: 'SpliceAI Lookup (Broad Institute)',
          url: `https://spliceailookup.broadinstitute.org/#variant=${extChr}-${extPos}-${extRef}-${extAlt}&hg=38&distance=500&mask=1&precomputed=0`,
        });
      }

      if (extChr && extPos) {
        const ucscChr = /^chr/i.test(extChr) ? extChr : `chr${extChr}`;
        const ucscStart = Math.max(1, Number(extPos) - 200);
        const ucscEnd = Number(extPos) + 200;
        variant.push({
          label: 'UCSC',
          title: 'UCSC Genome Browser (hg38) — ±200 bp around the variant',
          url: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${enc2(`${ucscChr}:${ucscStart}-${ucscEnd}`)}`,
        });
      }
      if (haveCoords) {
        variant.push({
          label: 'Franklin',
          url: `https://franklin.genoox.com/clinical-db/variant/snp/${extChr}-${extPos}-${extRef}-${extAlt}-hg38`,
        });
      }

      if (haveCoords) {
        const pvCoords = `${extChr} ${extPos} ${extRef} ${extAlt}`;
        variant.push({
          label: 'ProtVar',
          title: 'ProtVar (EBI) — functional, conservation, and co-located variant annotations',
          url: `https://www.ebi.ac.uk/ProtVar/query?search=${enc2(pvCoords)}`,
        });
      }

      const pmEv = ev.pubmed || {};
      if (pmEv.variant_search_url) {
        variant.push({
          label: 'PubMed',
          title: 'PubMed — variant-specific search',
          url: pmEv.variant_search_url,
        });
      }

      if (haveGene) {

        gene.push({
          label: 'ClinGen',
          url: `https://search.clinicalgenome.org/kb/genes/${enc2(extGene)}`,
        });

        gene.push({
          label: 'PanelApp',
          title: 'PanelApp Australia — gene panels',
          url: `https://panelapp-aus.org/panels/entities/${enc2(extGene)}`,
        });

        gene.push({
          label: 'OMIM',
          title: 'OMIM — gene-symbol search',
          url: `https://www.omim.org/search?search=${enc2(extGene)}`,
        });
      }

      const mdEv = ev.medgen || {};
      if (mdEv.url) {
        gene.push({
          label: 'MedGen',
          title: 'MedGen (NCBI) — gene-disease associations',
          url: mdEv.url,
        });
      } else if (haveGene) {
        gene.push({
          label: 'MedGen',
          title: 'MedGen (NCBI) — gene-disease associations',
          url: `https://www.ncbi.nlm.nih.gov/medgen/?term=${enc2(extGene)}%5Bgene%5D`,
        });
      }

      const gcEv = ev.gencc || {};
      if (gcEv.url) {
        gene.push({ label: 'GenCC', url: gcEv.url });
      } else if (haveGene) {
        gene.push({
          label: 'GenCC',
          url: `https://search.thegencc.org/genes/HGNC?q=${enc2(extGene)}`,
        });
      }

      const cdEv = ev.chdgene || {};
      if (cdEv.url) {
        gene.push({
          label: 'CHDgene',
          title: 'CHDgene (Victor Chang) — curated CHD gene list',
          url: cdEv.url,
        });
      }

      const otEv = ev.opentargets_evidence || {};
      const otEnsembl = otEv.ensembl_id || (ev.vep && ev.vep.gene_id) || '';
      if (otEnsembl) {
        gene.push({
          label: 'Open Targets',
          title: 'Open Targets Platform — gene-disease association profile',
          url: `https://platform.opentargets.org/target/${enc2(otEnsembl)}`,
        });
      }

      const mgiEv = ev.mgi || {};
      if (mgiEv.url) {
        gene.push({
          label: 'MGI',
          title: 'MGI — mouse ortholog phenotypes',
          url: mgiEv.url,
        });
      }

      const gtxEv = ev.gtex || {};
      if (gtxEv.url) {
        gene.push({
          label: 'GTEx',
          title: 'GTEx Portal — bulk-tissue median expression',
          url: gtxEv.url,
        });
      }

      if (haveGene) {
        gene.push({
          label: 'Franklin Gene',
          title: 'Franklin (Genoox) — gene-disease associations and mechanism',
          url: `https://franklin.genoox.com/clinical-db/gene/${enc2(extGene)}`,
        });
      }

      const glitEv = ev.gene_literature || {};
      if (glitEv.search_url) {
        gene.push({
          label: 'PubMed (gene)',
          title: 'PubMed — gene-level cardiac literature search',
          url: glitEv.search_url,
        });
      }

      if (haveGene) {

        protein.push({
          label: 'UniProt',
          url: uniprotAcc
            ? `https://www.uniprot.org/uniprotkb/${enc2(uniprotAcc)}/entry`
            : `https://www.uniprot.org/uniprotkb?query=${enc2(extGene)}+AND+organism_id:9606`,
        });
      }

      if (uniprotAcc) {
        protein.push({
          label: 'AlphaFold',
          title: 'AlphaFold DB — predicted 3-D structure (DeepMind/EBI)',
          url: `https://alphafold.ebi.ac.uk/entry/${enc2(uniprotAcc)}`,
        });
      }

      const bgEv = ev.biogrid || {};
      if (bgEv.url) {
        protein.push({
          label: 'BioGRID',
          title: 'BioGRID — curated protein–protein interactions',
          url: bgEv.url,
        });
      }

      const renderPill = l => {
        const t = l.title ? ` title="${l.title.replace(/"/g, '&quot;')}"` : '';
        return `<a class="db-link" href="${l.url}" target="_blank" rel="noopener"${t}>${l.label} ↗</a>`;
      };
      const renderGroup = (label, links, category) => links.length
        ? `<div class="summary-ext-group" data-category="${category}">
             <span class="summary-ext-group-label">${label}</span>
             <div class="summary-ext-pills">${links.map(renderPill).join('')}</div>
           </div>`
        : '';
      // Unlabelled flat pill list for a single scope (Variant / Gene / Protein tabs).
      const renderFlat = (links, category) => links.length
        ? `<div class="summary-ext-section"><div class="summary-ext-group" data-category="${category}"><div class="summary-ext-pills">${links.map(renderPill).join('')}</div></div></div>`
        : '';

      if (scope === 'variant') return renderFlat(variant, 'variant');
      if (scope === 'gene')    return renderFlat(gene, 'gene');
      if (scope === 'protein') return renderFlat(protein, 'protein');

      // scope === 'all' (Export tab): every group, labelled.
      if (!variant.length && !gene.length && !protein.length) return '';
      return `<div class="summary-ext-section">
        ${renderGroup('Variant', variant, 'variant')}
        ${renderGroup('Gene', gene, 'gene')}
        ${renderGroup('Protein', protein, 'protein')}
      </div>`;
    }

    // Expand-all / collapse-all toggle for the observed-evidence card.
    // Looks up every <details class="ev-collapsible"> inside the same
    // .rsec as the button, flips them in unison.
    function toggleAllEvidence(btn) {
      // A per-group button (inside a .ev-accgroup heading) flips only that
      // group's rows; fall back to the whole panel for any legacy caller.
      const sec = btn.closest('.ev-accgroup') || btn.closest('.rsec');
      if (!sec) return;
      const items = sec.querySelectorAll('details.ev-collapsible');
      const state = btn.getAttribute('data-state') || 'collapsed';
      const expand = state === 'collapsed';
      items.forEach(d => { d.open = expand; });
      btn.setAttribute('data-state', expand ? 'expanded' : 'collapsed');
      btn.textContent = expand ? 'Collapse all' : 'Expand all';
    }

    // Gene + Protein tabs use .gc-collapsible <div>s toggled via an `open`
    // CLASS (not native <details>), so they share this expand-all that flips
    // the class. On expand, any 3-D box opened this way kicks off its own
    // lazy viewer load (the box's manual onclick does the same).
    function toggleAllGeneSections(btn) {
      // A per-group button (inside a .gcx-group heading) flips only that
      // group's sections; the Protein tab's single top button has no enclosing
      // group, so it falls back to flipping every section in the pane.
      const scope = btn.closest('.gcx-group') || btn.closest('.result-tab-pane') || document;
      const items = scope.querySelectorAll('.gc-collapsible');
      const expand = (btn.getAttribute('data-state') || 'collapsed') === 'collapsed';
      items.forEach(d => {
        d.classList.toggle('open', expand);
        if (expand && window.__hvLoad3D && d.querySelector('.hv3d-mount')) window.__hvLoad3D(d);
      });
      btn.setAttribute('data-state', expand ? 'expanded' : 'collapsed');
      btn.textContent = expand ? 'Collapse all' : 'Expand all';
    }

    // Human-readable labels for the structured clinical-context dropdown
    // values. Kept as a single map so the Summary tab, TXT report, and
    // CSV export share one source of truth for the display strings —
    // never duplicate the option labels from the HTML directly.
    const _CLINICAL_LABELS = {
      zygosity:    { het: 'Heterozygous', hom: 'Homozygous', hemi: 'Hemizygous' },
      inheritance: {
        AD: 'Autosomal dominant', AR: 'Autosomal recessive',
        XLD: 'X-linked dominant', XLR: 'X-linked recessive',
        MT: 'Mitochondrial', DN: 'De novo',
      },
      sex:         { female: 'Female', male: 'Male' },
      trio:        { duo: 'Duo — one parent tested', trio: 'Full trio — both parents tested' },
      denovo:      {
        unconfirmed: 'Unconfirmed / assumed',
        confirmed: 'Confirmed de novo',
        inherited_affected: 'Inherited from affected parent',
        inherited_unaffected: 'Inherited from unaffected parent',
      },
      family_history: {
        no_family_history:       'No family history',
        positive_family_history: 'Positive family history',
        segregation_data:        'Segregation data available',
        unknown:                 'Unknown',
      },
    };
    function _clinicalLabel(group, value) {
      if (!value) return '';
      return (_CLINICAL_LABELS[group] && _CLINICAL_LABELS[group][value]) || value;
    }

    // Format the stored genome_build value as "GRCh38 / hg38" /
    // "GRCh37 / hg19". Falls back to the raw value when it's not one
    // of the two known builds (e.g. an upstream label change) so we
    // never blank out a real value just because the formatter missed.
    function _formatGenomeBuild(b) {
      const v = String(b || '').trim();
      if (!v) return '';
      if (/^grch38$/i.test(v) || /^hg38$/i.test(v)) return 'GRCh38 / hg38';
      if (/^grch37$/i.test(v) || /^hg19$/i.test(v)) return 'GRCh37 / hg19';
      return v;
    }

    // Builds a "chrN:pos:ref:alt" position string from the VEP block for
    // the Query Details "Variant (genomic)" cell. The colon form matches
    // the input-field regex (round-trippable) and uses one shape for
    // SNVs, deletions (alt='-'), insertions (ref='-') and MNVs alike.
    // Uses `allele_string` exactly as VEP returns it — VEP keeps
    // HGVS-input alleles in transcript orientation, which is what the
    // curator typed in the c. field, so a c.770C>T input renders as
    // ...:C:T here. Coordinate inputs land on the forward strand from
    // the start so the same code path is correct for both flows.
    function _formatHgvsG(vep) {
      if (!vep || !vep.ok) return '';
      const chrom = vep.seq_region_name;
      const pos = vep.start;
      const alleles = String(vep.allele_string || '').split('/');
      if (!chrom || !pos || alleles.length !== 2) return '';
      const [ref, alt] = alleles;
      if (!ref || !alt) return '';
      return `chr${chrom}:${pos}:${ref}:${alt}`;
    }

    // ── Summary right-rail "Case" recap (R2 style) ────────────────────
    // Compact, editorial recap of the CLINICAL / case context only — the
    // variant identity already lives in the verdict hero, so it is NOT
    // repeated here. Derives its fields from lastVariant + the
    // clinical-label maps; every
    // field is individually guarded and any value that is empty /
    // "Not provided" / "Unknown" / "N/A" is OMITTED so the box stays
    // tight. Returns '' when nothing case-related was supplied, so the
    // aside collapses cleanly. `evidence` is accepted for call-site symmetry
    // (HPO/phenotype actually lives on lastVariant); the param is unused.
    // Minimal HP code → readable term map (covers the example cases + common
    // CHD / cardiomyopathy terms). Unmapped codes fall back to the raw code.
    const HPO_NAMES = {
      'HP:0001638':'Cardiomyopathy','HP:0001639':'Hypertrophic cardiomyopathy',
      'HP:0001644':'Dilated cardiomyopathy','HP:0011664':'Left ventricular noncompaction',
      'HP:0001712':'Left ventricular hypertrophy','HP:0001645':'Sudden cardiac death',
      'HP:0011675':'Arrhythmia','HP:0001962':'Palpitations',
      'HP:0001629':'Ventricular septal defect','HP:0001631':'Atrial septal defect',
      'HP:0001674':'Atrioventricular septal defect','HP:0006695':'Atrioventricular canal defect',
      'HP:0001636':'Tetralogy of Fallot','HP:0001643':'Patent ductus arteriosus',
      'HP:0001680':'Coarctation of the aorta','HP:0001650':'Aortic valve stenosis',
      'HP:0010886':'Pulmonary artery stenosis','HP:0000589':'Coloboma',
      'HP:0000496':'Abnormal eye movements','HP:0001702':'Abnormal heart valve morphology'
    };
    function _phenotypeToText(s) {
      if (!s) return '';
      const seen = {}, out = [];
      String(s).split(',').forEach(function(tok){
        let t = tok.trim(); if (!t) return;
        if (/^HP:\d+$/i.test(t)) t = HPO_NAMES[t.toUpperCase()] || t;
        const k = t.toLowerCase();
        if (!seen[k]) { seen[k] = 1; out.push(t); }
      });
      return out.join(', ');
    }

    function _buildCaseRecapHTML(evidence) {
      const esc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
      const lv = lastVariant || {};

      // Resolved display strings.
      const geneVal  = (lv.gene || '').trim();
      const variantVal = (lv.hgvs_c || '').trim();
      const buildVal = (lv.genome_build || '').trim();
      const hpoVal   = _phenotypeToText((lv.hpo || '').trim());
      const inhVal   = lv.inheritance_input ? _clinicalLabel('inheritance', lv.inheritance_input) : '';
      const zygVal   = lv.zygosity ? _clinicalLabel('zygosity', lv.zygosity) : '';
      const zygExtra = (zygVal && lv.zygosity_inferred) ? ' (inferred)' : '';
      const sexVal   = lv.proband_sex
        ? (lv.proband_sex.charAt(0).toUpperCase() + lv.proband_sex.slice(1))
        : '';
      const trioVal   = lv.trio_status ? _clinicalLabel('trio', lv.trio_status) : '';
      const denovoVal = (lv.trio_status && lv.denovo_status)
        ? _clinicalLabel('denovo', lv.denovo_status)
        : '';
      const famSummary = lv.family_history_summary || '';
      const famLabel   = (famSummary && famSummary !== 'unknown')
        ? _clinicalLabel('family_history', famSummary)
        : (lv.family || '');
      const famTitle   = lv.family ? ` title="${esc(lv.family)}"` : '';

      // Structured segregation counts — only render when actually supplied.
      const _segNum = (v) => (v && Number(v) > 0) ? String(v) : '';
      const _yesNo  = (v) => v ? (v.charAt(0).toUpperCase() + v.slice(1)) : '';
      const segCarriersVal    = _segNum(lv.seg_affected_carriers);
      const segNoncarriersVal = _segNum(lv.seg_affected_noncarriers);
      const segMeiosesVal     = _segNum(lv.seg_meioses);
      const inTransVal        = _yesNo(lv.in_trans_pathogenic);
      const altCauseVal       = lv.alt_cause_present === 'yes'
        ? (lv.alt_cause_detail ? `Yes — ${lv.alt_cause_detail}` : 'Yes')
        : _yesNo(lv.alt_cause_present);

      const OMIT = new Set(['', 'not provided', 'unknown', 'n/a']);
      const has = (v) => { const s = String(v == null ? '' : v).trim(); return s && !OMIT.has(s.toLowerCase()); };
      const NP = '<span class="case-recap__np">Not provided</span>';
      const rows = [];
      // Core fields ALWAYS show (with "Not provided" when blank); advanced
      // segregation fields show only when the curator supplied them.
      const add = (label, value, opts) => {
        opts = opts || {};
        const ok = has(value);
        const valCls = (opts.mono && ok) ? 'case-recap__v case-recap__v--mono' : 'case-recap__v';
        const valHtml = ok ? `${esc(value)}${opts.suffix || ''}` : NP;
        rows.push(`<div class="case-recap__row"${opts.title || ''}>
            <span class="case-recap__k">${esc(label)}</span>
            <span class="${valCls}">${valHtml}</span>
          </div>`);
      };
      const addIf = (label, value, opts) => { if (has(value)) add(label, value, opts); };

      // Everything the curator entered for the search.
      add('Gene', geneVal, { mono: true });
      add('Variant', variantVal, { mono: true });
      addIf('Genome build', buildVal);
      add('Phenotype', hpoVal);
      add('Inheritance', inhVal);
      add('Zygosity', zygVal, { suffix: zygExtra ? `<span class="case-recap__note">${esc(zygExtra)}</span>` : '' });
      add('Proband sex', sexVal);
      add('Trio status', trioVal);
      add('De novo', denovoVal);
      add('Family history', famLabel, { title: famTitle });
      // Advanced segregation counts — only when present.
      addIf('Affected carriers', segCarriersVal);
      addIf('Affected non-carriers', segNoncarriersVal);
      addIf('Informative meioses', segMeiosesVal);
      addIf('2nd pathogenic in trans', inTransVal);
      addIf('Alternate molecular cause', altCauseVal);

      return `<div class="case-recap">
          <div class="case-recap__body">${rows.join('')}</div>
        </div>`;
    }

    // ════════ "Ask about this variant" chatbot ════════
    // Runs server-side on HeartVar's own key (POST /api/chat) and grounds
    // answers in the loaded classification + criteria + clinical context.
    // The pill only appears when an AI result is loaded (window.__hvAiUsed) —
    // i.e. the curator ticked "Include AI interpretation" and it succeeded.
    // Suggestion chips. `mode` routes to the server-side clinical-report
    // preamble (POST /api/chat mode="report"); omitted means a normal short
    // chat answer.
    const HVC_PROMPTS = [
      { q: 'Generate a clinical report', mode: 'report' },
      { q: 'Why was this classification reached?' },
      { q: 'What evidence would change the classification?' },
      { q: 'Summarise the population & ClinVar evidence' },
      { q: 'Does the clinical context affect any criteria?' },
    ];
    // A curator who types the request rather than clicking the chip should get
    // the same artefact. Narrow on purpose — it matches the phrase, not the
    // word "report" on its own, so "does ClinVar report this?" stays a chat
    // question.
    const HVC_REPORT_RE = /\bclinical report\b/i;
    let hvcHistory = [];
    function hvcEsc(s){ return String(s == null ? '' : s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }
    // Render a SHORT assistant reply as safe HTML. The model is told to keep
    // replies terse and heading-free, but defends both ways: HTML is escaped
    // first (model output is untrusted), then a minimal markdown subset is
    // applied — bold, italics, inline code, bullet/numbered lists, links, and
    // clickable PMIDs. Markdown headings ('#') are demoted to a bold line so a
    // stray '##' never leaks as a literal hashtag into the chat bubble.
    function hvcRenderMd(src){
      const inline = s => s
        .replace(/`([^`]+)`/g, (_, c) => '<code>' + c + '</code>')
        .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
        .replace(/(^|[^*])\*(?!\s)([^*\n]+?)\*(?!\*)/g, '$1<em>$2</em>')
        .replace(/\[([^\]]+)\]\((https?:\/\/[^)\s]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>')
        .replace(/\bPMID:?\s*(\d{4,9})\b/gi, '<a href="https://pubmed.ncbi.nlm.nih.gov/$1/" target="_blank" rel="noopener">PMID $1</a>');
      const lines = hvcEsc(src).replace(/\r\n?/g, '\n').split('\n');
      const out = []; let list = null; let para = [];
      const flushPara = () => { if (para.length){ out.push('<p>' + para.join('<br>') + '</p>'); para = []; } };
      const closeList = () => { if (list){ out.push('</' + list + '>'); list = null; } };
      for (const raw of lines){
        const line = raw.replace(/\s+$/, '');
        if (!line.trim()){ flushPara(); closeList(); continue; }
        const h = line.match(/^\s*#{1,6}\s+(.*)$/);
        if (h){ flushPara(); closeList(); out.push('<p class="hvc-h"><strong>' + inline(h[1]) + '</strong></p>'); continue; }
        const ul = line.match(/^\s*[-*•]\s+(.*)$/);
        const ol = line.match(/^\s*\d+[.)]\s+(.*)$/);
        if (ul){ flushPara(); if (list !== 'ul'){ closeList(); out.push('<ul>'); list = 'ul'; } out.push('<li>' + inline(ul[1]) + '</li>'); continue; }
        if (ol){ flushPara(); if (list !== 'ol'){ closeList(); out.push('<ol>'); list = 'ol'; } out.push('<li>' + inline(ol[1]) + '</li>'); continue; }
        closeList(); para.push(inline(line));
      }
      flushPara(); closeList();
      return out.join('');
    }
    function hvChatSyncFab(){
      const fab = document.getElementById('hvc-fab'); if (!fab) return;
      const ready = !!(window.__hvRan && window.__hvAiUsed && !document.body.classList.contains('hv-landing-on'));
      fab.classList.toggle('show', ready);
      if (!ready) hvChatClose();
    }
    function hvChatReset(){
      hvcHistory = [];
      const m = document.getElementById('hvc-msgs');
      if (m){ m.innerHTML = ''; delete m.dataset.init; }
      hvChatClose();
    }
    // Facts the "Generate a clinical report" mode needs and a chat answer never
    // did: the accessions for the identity line, the gene's mechanism and
    // conditions with MIM numbers, comparable ClinVar variants, and — carefully
    // — whether the absence of segregation or functional evidence is a CHECKED
    // absence or simply unassessed. The report is instructed to state a checked
    // absence and omit an unchecked one, so this is where that distinction is
    // established; getting it wrong here would put "no functional evidence has
    // been identified" into a report where nothing was ever looked at.
    // Everything is optional — each line appears only when its source resolved.
    function buildReportFacts(){
      const lv = (typeof lastVariant !== 'undefined' && lastVariant) || {};
      const ev = (typeof lastEvidence !== 'undefined' && lastEvidence) || {};
      const r  = (typeof lastResult  !== 'undefined' && lastResult)  || {};
      const vep = ev.vep || {};
      const out = [];
      const push = (k, v) => { if (v !== null && v !== undefined && String(v).trim() !== '') out.push(`${k}: ${v}`); };

      // ── Identity line inputs. The report wants the RefSeq NM_ accession the
      // example uses; mane_select_accession is exactly that, with the Ensembl
      // ENST as the fallback when a gene has no MANE Select.
      push('RefSeq transcript (MANE Select)', vep.mane_select_accession);
      push('Ensembl transcript', vep.selected_transcript_id || vep.transcript_id);
      push('HGVS coding', _stripHgvsPrefix(vep.hgvsc || lv.hgvs_c || ''));
      push('HGVS protein', _stripHgvsPrefix(vep.hgvsp || ''));
      push('Molecular consequence', (vep.most_severe_consequence || '').replace(/_/g, ' '));
      push('Exon', vep.exon);
      if (vep.phylop100way != null) push('phyloP100way conservation', vep.phylop100way);

      // A report-ready molecular-observation sentence, composed here rather than
      // left to the model. The report preamble asks for this bullet first, and
      // the model omitted it on two successive runs because its three
      // components (consequence, repeat status, conservation) sat on separate
      // context lines and had to be assembled. Pre-composing turns it into a
      // fact to echo instead of prose to construct. The conservation band is a
      // deterministic read of phyloP — the same >2.0 threshold BP7 uses — and
      // the number is always quoted beside it.
      {
        const csq = (vep.most_severe_consequence || '').replace(/_/g, ' ');
        const pp = vep.phylop100way;
        const bits = [];
        if (csq) bits.push(csq);
        // Repeat status only when a criterion actually assessed it (BP3/PM4
        // evaluate repeat regions); never asserted from silence.
        const rep = (Array.isArray(r.criteria) ? r.criteria : [])
          .find(c => /^(BP3|PM4)$/.test(String(c.code || '').split('/')[0].trim())
                     && /repeat|repetitive/i.test(c.evidence || ''));
        if (rep) {
          bits.push(/not\s+(?:in|located)|non-repetitive|NOT a repeat/i.test(rep.evidence)
            ? 'in a non-repetitive region' : 'in a repetitive region');
        }
        if (typeof pp === 'number') {
          bits.push(`${pp >= 4 ? 'highly conserved' : pp >= 2 ? 'conserved' : 'poorly conserved'} (phyloP ${pp})`);
        }
        if (bits.length > 1) push('Molecular observation (use as the first pathogenic-evidence bullet)', bits.join(', '));
      }

      // ── Population frequency. "Absent" is a strong report statement, so only
      // say it when gnomAD actually answered and returned no variant record.
      const g = ev.gnomad || {};
      if (g.ok) {
        const gv = g.variant;
        if (!gv) {
          push('gnomAD v4', 'variant ABSENT from gnomAD (exomes and genomes) — lookup succeeded and returned no record');
        } else {
          const ex = gv.exome || {}, ge = gv.genome || {};
          const parts = [];
          if (ex.af != null) parts.push(`exome AF ${ex.af}`);
          if (ge.af != null) parts.push(`genome AF ${ge.af}`);
          const fafs = [ex.faf95 && ex.faf95.popmax, ge.faf95 && ge.faf95.popmax].filter(x => x != null);
          if (fafs.length) parts.push(`FAF95 popmax ${Math.max.apply(null, fafs)}`);
          if (parts.length) push('gnomAD v4', parts.join(', '));
        }
      }

      // ── Protein domain containing the residue (the "well-established
      // functional domain" bullet). Read off the same UniProt features the
      // domain track draws.
      const pd = window.heartvarProteinData || {};
      const aa = (typeof parseAminoAcidPosition === 'function')
        ? parseAminoAcidPosition(vep.hgvsp || '') : null;
      if (Array.isArray(pd.uniprotFeatures) && aa != null) {
        const hits = pd.uniprotFeatures
          .filter(f => f && f.start <= aa && aa <= f.end)
          .map(f => `${f.description || f.type}${f.type ? ' (' + f.type + ')' : ''}`);
        if (hits.length) push(`UniProt features containing residue ${aa}`, hits.join('; '));
      }
      if (pd.uniprotAccession) push('UniProt accession', pd.uniprotAccession);

      // ── Gene-level: mechanism, conditions + MIM numbers, inheritance.
      const mech = ev.gene_mechanism || null;
      if (mech && mech.mechanism) {
        push('Gene disease mechanism', `${mech.label || mech.mechanism}${mech.summary ? ' — ' + mech.summary : ''}`);
      }
      const md = ev.medgen || {};
      if (md.ok && Array.isArray(md.conditions) && md.conditions.length) {
        const conds = md.conditions.slice(0, 8).map(c =>
          `${c.name}${c.mim ? ` (MIM#${c.mim})` : ''}${c.moi ? ` [${c.moi}]` : ''}`
        );
        push('Gene-associated conditions (MedGen/OMIM)', conds.join('; '));
      }

      // ── ClinVar: this variant, then comparable variants at the same residue
      // (the report's "a missense variant p.(ArgNNNHis) has been reported…"
      // bullet comes from here).
      const cv = ev.clinvar || {};
      if (cv.ok) {
        const recs = Array.isArray(cv.records) ? cv.records : [];
        if (recs.length) {
          push('ClinVar (this variant)', recs.slice(0, 4).map(rec => [
            rec.clinical_significance || 'no classification',
            rec.review_status || '',
            rec.accession || '',
            (rec.number_submitters != null ? `${rec.number_submitters} submitter(s)` : ''),
          ].filter(Boolean).join(' · ')).join(' | '));
        } else {
          push('ClinVar (this variant)', 'lookup succeeded — NO record found for this variant');
        }
      }
      const pm5 = ev.clinvar_pm5_candidates || {};
      if (pm5.ok && Array.isArray(pm5.candidates) && pm5.candidates.length) {
        push(`ClinVar variants at residue ${pm5.protein_position}`, pm5.candidates.slice(0, 6)
          .map(c => `${c.name} — ${c.tier}${c.stars != null ? ` (${c.stars}★)` : ''}`).join('; '));
      }

      // ── Checked-vs-unchecked absences. A criterion HeartVar EVALUATED and
      // did not meet is a checked absence the report may state; one marked
      // not_assessed (the AI-free flow) or insufficient_data was never looked
      // at, and claiming absence there would be a false negative in a clinical
      // document. Only the evaluated case is emitted.
      const crit = Array.isArray(r.criteria) ? r.criteria : [];
      const byCode = {};
      crit.forEach(c => { byCode[String(c.code || '').split('/')[0].trim()] = c; });
      const absence = (code, label) => {
        const c = byCode[code];
        if (!c) return;
        if (c.status === 'not_met' || c.status === 'na') {
          push(`${label} (EVALUATED, none found)`, c.evidence || 'threshold not met');
        }
      };
      absence('PS3', 'Published functional evidence');
      absence('PP1', 'Published segregation evidence');
      absence('PS4', 'Case/proband enrichment evidence');
      absence('PM5', 'Comparable variants at this residue');

      // ── How the variant was inherited in this family (last report bullet).
      if (lv.trio_status && lv.denovo_status) {
        push('Inheritance in this family', `${_clinicalLabel('denovo', lv.denovo_status)} (${_clinicalLabel('trio', lv.trio_status)})`);
      }
      if (lv.seg_affected_carriers) push('Affected relatives carrying the variant', lv.seg_affected_carriers);
      if (lv.in_trans_pathogenic) push('Second P/LP allele in trans', lv.in_trans_pathogenic);
      return out;
    }

    // `mode` shapes ONE thing: how met criteria are written.
    //
    // Chat mode labels them "PM1 [PM1_Moderate]: <evidence>" because a curator
    // asking why a call was reached wants the codes. Report mode drops the codes
    // and emits the evidence statements alone — the requested report format is
    // plain clinical prose with no code prefixes, and instructing that in the
    // preamble was NOT enough on its own: the model mirrored the context's own
    // "CODE [STRENGTH]:" shape and prefixed every bullet with it. The in-context
    // pattern beat the instruction, so the pattern is what changed. Verified in
    // the browser against a real report call.
    function buildChatContext(mode){
      const forReport = mode === 'report';
      const lv = (typeof lastVariant !== 'undefined' && lastVariant) || {};
      const r  = (typeof lastResult  !== 'undefined' && lastResult)  || {};
      const lines = [];
      const id = `${lv.gene || ''} ${lv.hgvs_c || ''}`.trim();
      if (id) lines.push('Variant: ' + id);
      if (lv.gene) lines.push('Gene: ' + lv.gene);
      if (r.classification) lines.push('HeartVar classification: ' + r.classification + (Number.isFinite(r.points_total) ? ` (ACMG/AMP ${r.points_total >= 0 ? '+' : ''}${r.points_total} pts)` : ''));
      const crit = Array.isArray(r.criteria) ? r.criteria : [];
      const fmt = c => `${c.code}${c.criteria_strength && c.criteria_strength !== c.code ? ' [' + c.criteria_strength + ']' : ''}: ${c.evidence || ''}`;
      const met = crit.filter(c => c.status === 'met');
      if (met.length) {
        const pathMet = met.filter(c => c.direction === 'pathogenic');
        const benMet = met.filter(c => c.direction === 'benign');
        if (forReport) {
          // Evidence statements only, grouped by direction so the report can
          // fill its two evidence sections without ever seeing a code to copy.
          if (pathMet.length) lines.push('Evidence supporting pathogenicity (statements, no codes):\n'
            + pathMet.map(c => '- ' + (c.evidence || '')).join('\n'));
          if (benMet.length) lines.push('Evidence supporting a benign reading (statements, no codes):\n'
            + benMet.map(c => '- ' + (c.evidence || '')).join('\n'));
        } else {
          lines.push('Criteria met:\n' + met.map(c => '- ' + fmt(c)).join('\n'));
        }
      }
      if (!forReport) {
        const insuf = crit.filter(c => c.status === 'insufficient_data');
        if (insuf.length) lines.push('Insufficient data to assess: ' + insuf.map(c => c.code).join(', '));
        const na = crit.filter(c => c.status === 'na' || c.status === 'not_met');
        if (na.length) lines.push('Not applicable / not met: ' + na.map(c => c.code).join(', '));
      }
      const cc = [];
      if (lv.hpo) cc.push('phenotype: ' + lv.hpo);
      if (lv.inheritance_input) cc.push('inheritance: ' + lv.inheritance_input);
      if (lv.zygosity) cc.push('zygosity: ' + lv.zygosity);
      if (lv.proband_sex) cc.push('proband sex: ' + lv.proband_sex);
      if (lv.trio_status) cc.push('trio: ' + lv.trio_status);
      if (lv.family) cc.push('family history: ' + lv.family);
      if (cc.length) lines.push('Clinical context — ' + cc.join('; '));
      // Annotation-level facts. Included for every question, not just the
      // report: the same accessions and gene-level detail make ordinary answers
      // better too, and one context builder cannot drift from itself.
      const facts = buildReportFacts();
      if (facts.length) lines.push('Annotation & gene-level facts:\n' + facts.map(f => '- ' + f).join('\n'));
      return lines.filter(Boolean).join('\n');
    }
    function hvcAddMsg(html, role){
      const m = document.getElementById('hvc-msgs');
      const d = document.createElement('div'); d.className = 'hvc-msg ' + role;
      const b = document.createElement('div'); b.className = 'hvc-bub'; b.innerHTML = html;
      d.appendChild(b); m.appendChild(d); m.scrollTop = m.scrollHeight; return d;
    }
    window.hvChatOpen = function(){
      if (!window.__hvAiUsed) return;   // only when an AI result is loaded
      const panel = document.getElementById('hvc-panel');
      // Cancel any in-flight collapse-to-bubble animation before re-opening.
      panel.classList.remove('closing');
      panel.classList.add('open');
      document.getElementById('hvc-fab').classList.add('hidden-by-panel');
      const lv = (typeof lastVariant !== 'undefined' && lastVariant) || {};
      const msgs = document.getElementById('hvc-msgs');
      if (msgs && !msgs.dataset.init){
        msgs.dataset.init = '1';
        hvcAddMsg(`I have the full interpretation for <strong>${hvcEsc(`${lv.gene || ''} ${lv.hgvs_c || ''}`.trim())}</strong> loaded. Ask me about the evidence, ACMG/AMP criteria, or classification.`, 'a');
        const p = document.createElement('div'); p.className = 'hvc-prompts'; p.id = 'hvc-prompts';
        HVC_PROMPTS.forEach(pr => {
          const b = document.createElement('button');
          b.type = 'button';
          b.className = 'hvc-prompt' + (pr.mode === 'report' ? ' hvc-prompt--report' : '');
          b.textContent = pr.q;
          b.onclick = () => hvChatSend(pr.q, pr.mode);
          p.appendChild(b);
        });
        const wrap = document.createElement('div'); wrap.className = 'hvc-msg a'; wrap.appendChild(p); msgs.appendChild(wrap);
      }
      setTimeout(() => { const i = document.getElementById('hvc-input'); if (i) i.focus(); }, 320);
    };
    window.hvChatClose = function(){
      const p = document.getElementById('hvc-panel');
      const f = document.getElementById('hvc-fab');
      // Reveal the bubble first so the panel visibly collapses *into* it.
      if (f) f.classList.remove('hidden-by-panel');
      if (!p) return;
      if (!p.classList.contains('open')){
        // Not open (e.g. reset on a new run) — nothing to animate.
        p.classList.remove('closing');
        return;
      }
      // Collapse down into the bubble (genie effect) rather than sliding off
      // to the right. `.closing` runs a scale/translate keyframe toward the
      // bottom-right bubble; we drop it (and `.open`) once the animation ends
      // so the panel returns to its parked off-screen state, ready to re-open.
      p.classList.add('closing');
      p.classList.remove('open');
      let cleaned = false;
      const done = () => {
        if (cleaned) return;
        cleaned = true;
        p.removeEventListener('animationend', done);
        // Snap back to the parked (off-screen) state with the base transition
        // suppressed — otherwise removing `.closing` would transition the
        // (now full-opacity) panel from the bubble back out to the right,
        // undoing the collapse we just played.
        p.style.transition = 'none';
        p.classList.remove('closing');
        void p.offsetWidth; // commit the parked transform before re-enabling
        p.style.transition = '';
      };
      p.addEventListener('animationend', done);
      // Fallback in case animationend never fires (reduced-motion / display
      // swap mid-animation): matches the .3s keyframe duration.
      setTimeout(done, 360);
    };
    window.hvChatSend = async function(preset, mode){
      const inp = document.getElementById('hvc-input');
      const text = String(preset != null ? preset : (inp ? inp.value : '')).trim();
      if (!text) return;
      // Chip-supplied mode wins; otherwise recognise a typed request for the
      // same artefact so the phrasing route and the button route agree.
      const sendMode = mode || (HVC_REPORT_RE.test(text) ? 'report' : 'chat');
      if (inp && preset == null){ inp.value = ''; inp.style.height = 'auto'; }
      const prompts = document.getElementById('hvc-prompts'); if (prompts){ const w = prompts.closest('.hvc-msg'); (w || prompts).remove(); }
      hvcAddMsg(hvcEsc(text), 'u');
      const typing = hvcAddMsg('<div class="hvc-typing"><span></span><span></span><span></span></div>', 'a');
      try {
        // Server-side chat: HeartVar runs the (scoped) call on its own key. We
        // send only the evidence CONTEXT we already have loaded + the question +
        // a little history; the server OWNS the "only this variant / refuse
        // off-topic" system prompt, so it can't be stripped client-side.
        const resp = await fetch('/api/chat', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            context: buildChatContext(sendMode),
            question: text,
            history: hvcHistory.slice(-8),
            mode: sendMode,
          }),
        });
        if (!resp.ok) {
          let msg = 'The assistant is unavailable right now (' + resp.status + ').';
          try {
            const j = await resp.json();
            // `detail` is a plain string for most failures but an object
            // ({error, message}) for the sign-in and per-user-quota cases —
            // reading it blindly would render "[object Object]".
            if (j && j.detail) {
              msg = (typeof j.detail === 'object') ? (j.detail.message || msg) : j.detail;
              if (j.detail.error === 'signin_required'
                  && typeof hvAuthHandle401 === 'function') {
                hvAuthHandle401();
              }
            }
          } catch (e) {}
          throw new Error(msg);
        }
        const body = await resp.json();
        const answer = (body && body.answer) || '';
        typing.remove();
        if (sendMode === 'report') {
          // A report exists to be pasted into a lab system, so it renders as
          // preformatted text (its line structure IS the format — collapsing it
          // through the markdown renderer would lose the identity line and the
          // section breaks) with a copy control that yields the plain text.
          const wrap = hvcAddMsg('', 'a');
          const bub = wrap.querySelector('.hvc-bub');
          bub.classList.add('hvc-bub--report');
          const pre = document.createElement('pre');
          pre.className = 'hvc-report';
          pre.textContent = answer;
          const bar = document.createElement('div');
          bar.className = 'hvc-report__bar';
          const btn = document.createElement('button');
          btn.type = 'button';
          btn.className = 'hvc-report__copy';
          btn.textContent = 'Copy report';
          btn.onclick = () => {
            navigator.clipboard.writeText(answer).then(
              () => { btn.textContent = 'Copied'; setTimeout(() => { btn.textContent = 'Copy report'; }, 1600); },
              () => { btn.textContent = 'Copy failed'; },
            );
          };
          bar.appendChild(btn);
          bub.appendChild(pre);
          bub.appendChild(bar);
        } else {
          hvcAddMsg(hvcRenderMd(answer), 'a');
        }
        hvcHistory.push({ role: 'u', text }); hvcHistory.push({ role: 'a', text: answer });
      } catch (e){
        typing.remove();
        hvcAddMsg("Sorry — " + hvcEsc(String((e && e.message) || e)).slice(0, 200), 'a');
      }
    };
    window.hvChatKey = function(e){ if (e.key === 'Enter' && !e.shiftKey){ e.preventDefault(); hvChatSend(); } };

    function renderStage1Card(stage1, gene, hgvs_c) {
      const evidence = stage1.db_evidence || {};
      const variantId = stage1.variant_id;

      // Stage-1 summary: array of three short sentences in the current
      // schema, single string in the legacy fallback. Either way we join
      // into a flowing paragraph for the Summary tab.
      let summaryProse = '';
      if (Array.isArray(stage1.summary)) {
        summaryProse = stage1.summary
          .map(s => lowercaseTiersMidSentence(String(s).trim()))
          .filter(Boolean)
          .join(' ');
      } else if (stage1.summary) {
        summaryProse = lowercaseTiersMidSentence(String(stage1.summary));
      }

      // Summary right-rail "Case" recap — compact clinical/case context
      // only (variant identity already lives in the verdict hero; the
      // external-tool deeplinks live at the bottom of the Variant tab,
      // built inside _renderEvidence). Returns '' when no case context was
      // supplied, in which case the rail collapses (handled by the grid below).
      const caseRecapHTML = _buildCaseRecapHTML(evidence);
      const classText = stage1.classification || '—';
      // Score-card content. Stage 1 sometimes omits points_total; the
      // canonical value lands at stage-2 reveal (see revealResults). Until
      // then we render an em dash so the card chrome stays consistent.
      const scoreText = Number.isFinite(stage1.points_total)
        ? `${stage1.points_total >= 0 ? '+' : '−'}${Math.abs(stage1.points_total)} pts`
        : '—';

      // ── Verdict-hero variant identity line ───────────────────────────
      // Derives the same identity fields (gene /
      // HGVS c. / protein / transcript + MANE / genomic g. / build) so the
      // editorial hero shows the same provenance the Query Details grid
      // does, without re-fetching anything. Every part is individually
      // guarded — a missing field is simply dropped from the line rather
      // than rendering an empty separator.
      const _heroEsc = s => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      const _hvep = (evidence && evidence.vep) || {};
      const heroGene = (lastVariant && lastVariant.gene) || _hvep.gene_symbol || gene || '';
      const _heroUserHgvs = (lastVariant && lastVariant.hgvs_c) || hgvs_c || '';
      const _heroResolvedC = _stripHgvsPrefix(_hvep.hgvsc || '');
      const heroHgvsC = isCoordInput(_heroUserHgvs)
        ? (_heroResolvedC || _heroUserHgvs)
        : (_heroUserHgvs || _heroResolvedC);
      const heroProtein = _stripHgvsPrefix(_hvep.hgvsp || '');
      const heroTranscript = _hvep.selected_transcript_id || _hvep.transcript_id || '';
      const heroManeHTML = heroTranscript
        ? (_hvep.is_mane_select
            ? '<span class="v11-mane">MANE Select</span>'
            : (_hvep.is_mane_clinical ? '<span class="v11-mane">MANE Plus Clinical</span>' : ''))
        : '';
      const heroGenomic = _formatHgvsG(_hvep);
      const heroBuild = _formatGenomeBuild(
        _hvep.input_build
        || (lastVariant && lastVariant.genome_build)
        || _hvep.assembly_name
        || ''
      );
      // Assemble the identity line from the parts that resolved, joining
      // present parts with a muted dot separator.
      const _heroSep = '<span class="sep" aria-hidden="true">·</span>';
      const _heroParts = [];
      if (heroHgvsC)     _heroParts.push(`<span class="v11-mono">${_heroEsc(heroHgvsC)}</span>`);

      if (heroProtein)   _heroParts.push(`<span class="v11-mono">${_heroEsc(heroProtein)}</span>`);
      if (heroTranscript)_heroParts.push(`<span class="v11-mono">${_heroEsc(heroTranscript)}</span>${heroManeHTML}`);
      if (heroGenomic)   _heroParts.push(`<span class="v11-mono">${_heroEsc(heroGenomic)}</span>`);
      if (heroBuild)     _heroParts.push(`<span class="v11-mono">${_heroEsc(heroBuild)}</span>`);
      const verdictVariantHTML = (heroGene || _heroParts.length)
        ? `<div class="v11-verdict__variant">${heroGene ? `<span class="gene">${_heroEsc(heroGene)}</span>` : ''}${_heroParts.join(_heroSep)}</div>`
        : '';

      const _fmtLiftCoord = (c) => {
        const p = String(c || '').split('-');
        return p.length === 4 ? `chr${p[0]}:${p[1]}:${p[2]}:${p[3]}` : _heroEsc(c || '');
      };
      const liftoverNoteHTML = (_hvep && _hvep.input_build === 'GRCh37')
        ? `<div class="summary-liftover-note" role="note">
            <span class="summary-liftover-note__ic" aria-hidden="true">⤴</span>
            <span>Genome build converted — you entered <strong>GRCh37 / hg19</strong> coordinates${_hvep.input_coords ? ` (<span class="v11-mono">${_fmtLiftCoord(_hvep.input_coords)}</span>)` : ''}, which were lifted over to <strong>GRCh38 / hg38</strong>${_hvep.lifted_coords_grch38 ? ` (<span class="v11-mono">${_fmtLiftCoord(_hvep.lifted_coords_grch38)}</span>)` : ''}. All annotations and the classification below use the GRCh38 coordinates.</span>
          </div>`
        : '';

      const genomicHgvsNoteHTML = isGenomicHgvs(_heroUserHgvs)
        ? `<div class="summary-ghgvs-note" role="note">
            <span class="summary-ghgvs-note__ic" aria-hidden="true">⚠</span>
            <span>Genomic HGVS entered — genomic (<span class="v11-mono">g.</span>) notation is assembly-specific (the build is set by the reference accession, e.g. <span class="v11-mono">NC_000014.9</span> = GRCh38, <span class="v11-mono">NC_000014.8</span> = GRCh37). HeartVar annotated this on GRCh38 <strong>without liftover</strong>. Confirm your accession is a GRCh38 build; for a GRCh37 / hg19 variant, enter it as genomic coordinates instead so it is lifted over.</span>
          </div>`
        : '';

      const carrierFlag = (evidence && evidence.carrier_status) || null;
      const carrierBandHTML = (carrierFlag && carrierFlag.detail) ? (() => {
        const warn = carrierFlag.severity === 'warning';
        const style = warn
          ? 'background:var(--amber-light);border:1px solid var(--amber-border);border-left:3px solid var(--amber);color:var(--amber)'
          : 'background:var(--blue-light);border:1px solid var(--blue-border);border-left:3px solid var(--blue);color:var(--blue)';
        const _b = carrierFlag.basis || {};

        const _rec = Array.isArray(_b.recessive_conditions) ? _b.recessive_conditions : [];
        const _dom = Array.isArray(_b.dominant_conditions) ? _b.dominant_conditions : [];
        const _condLine = (label, arr) => arr.length
          ? `<div style="margin-top:4px"><span style="font-weight:600">${label}:</span> ${_heroEsc(arr.join('; '))}</div>`
          : '';
        const condHTML = (_rec.length || _dom.length)
          ? `<div class="carrier-band__cond" style="margin-top:8px;font-size:12px;line-height:1.5;opacity:.95">
               ${_condLine('Recessive condition(s)', _rec)}${_condLine('Dominant condition(s)', _dom)}
             </div>`
          : '';

        const scopeNote = _b.phenotype_scoped
          ? 'Scoped to the entered phenotype.'
          : (_b.phenotype_provided
              ? 'Gene-wide (entered phenotype did not match a curated disease).'
              : 'Gene-wide (no phenotype entered).');
        const srcHTML = (Array.isArray(_b.sources) && _b.sources.length)
          ? `<div class="carrier-band__src" style="margin-top:7px;font-size:11.5px;opacity:.8">${_heroEsc(scopeNote)} Mode-of-inheritance source(s): ${_heroEsc(_b.sources.join(', '))}${_b.zygosity ? ` · zygosity: ${_heroEsc(_b.zygosity)}${_b.zygosity_inferred ? ' (inferred)' : ''}` : ''}</div>`
          : '';
        return `<section class="v11-section">
                <div class="v11-eyebrow"><h3>Carrier status</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
                <div class="carrier-band" role="note" style="${style};padding:11px 13px;border-radius:6px;font-size:13px;line-height:1.55">
                  <strong>${warn ? '⚠' : 'ℹ'} ${_heroEsc(carrierFlag.title || 'Carrier status')}.</strong> ${_heroEsc(carrierFlag.detail)}
                  ${condHTML}
                  ${srcHTML}
                </div>
              </section>`;
      })() : '';

      const _kf = [];
      const _kfTile = (k, v, sub, tone) => `<div class="ev-snap__tile" data-tone="${tone || 'info'}"><div class="ev-snap__k">${k}</div><div class="ev-snap__v is-text">${v}</div>${sub ? `<div class="ev-snap__sub">${sub}</div>` : ''}</div>`;
      try {
        const _gn = evidence.gnomad || {};
        if (_gn.ok) {
          const _gv = _gn.variant || null;
          if (_gv) {
            const _ex = _gv.exome || {}, _ge = _gv.genome || {};
            const _pmv = [_ex.faf95 && _ex.faf95.popmax, _ge.faf95 && _ge.faf95.popmax].filter(x => x != null);
            const _pmax = _pmv.length ? Math.max(..._pmv) : null;
            const _af = (_ex.af != null ? _ex.af : (_ge.af != null ? _ge.af : null));
            const _shown = (_pmax != null ? _pmax : _af);
            const _tone = _pmax != null && _pmax >= 0.05 ? 'warn' : _pmax != null && _pmax >= 0.01 ? 'amber' : 'info';
            _kf.push(_kfTile('gnomAD AF', _shown != null ? (typeof formatAF === 'function' ? formatAF(_shown) : Number(_shown).toExponential(1)) : '—', 'Filtering AF (popmax)', _tone));
          } else {
            _kf.push(_kfTile('gnomAD AF', 'Absent', 'Not in gnomAD v4', 'good'));
          }
        }
        const _cv = evidence.clinvar || {};
        if (_cv.ok && _cv.found) {
          const _recs = _cv.records || [];
          const _best = [..._recs].sort((a, b) => (b.stars || 0) - (a.stars || 0))[0] || {};
          const _sig = (_best.clinical_significance || '').split(/[,/]/)[0].trim() || 'Reported';
          const _tone = /pathogenic/i.test(_sig) ? 'warn' : (/benign/i.test(_sig) ? 'good' : 'amber');
          const _k = Math.max(0, Math.min(4, _best.stars || 0));
          _kf.push(_kfTile('ClinVar', _heroEsc(_sig), `${'★'.repeat(_k)}${'☆'.repeat(4 - _k)}`, _tone));
        } else if (_cv.ok) {
          _kf.push(_kfTile('ClinVar', 'No records', 'Not submitted', 'muted'));
        }
        const _sa = evidence.spliceai || {};
        if (_sa.ok && _sa.max_delta != null) {
          const _md = _sa.max_delta;
          const _tone = _md >= 0.8 ? 'warn' : _md >= 0.5 ? 'amber' : _md >= 0.2 ? 'muted' : 'good';
          _kf.push(_kfTile('SpliceAI Δ', Number(_md).toFixed(2), _md >= 0.5 ? 'Splice-altering' : 'Low impact', _tone));
        }
        const _am = evidence.alphamissense || {};
        if (_am && _am.available !== false && !_am.not_applicable && _am.score != null) {
          const _c = _am.classification || 'ambiguous';
          const _tone = _c === 'likely_pathogenic' ? 'warn' : (_c === 'likely_benign' ? 'good' : 'muted');
          const _lab = _c === 'likely_pathogenic' ? 'Likely pathogenic' : (_c === 'likely_benign' ? 'Likely benign' : 'Ambiguous');
          _kf.push(_kfTile('AlphaMissense', Number(_am.score).toFixed(2), _lab, _tone));
        }

        if (_hvep.most_severe_consequence) {
          const _csq = String(_hvep.most_severe_consequence).replace(/_/g, ' ');
          const _csqLabel = _csq.charAt(0).toUpperCase() + _csq.slice(1);
          const _imp = String(_hvep.impact || '').toUpperCase();
          const _impTone = _imp === 'HIGH' ? 'warn' : (_imp === 'MODERATE' ? 'amber' : (_imp === 'LOW' ? 'info' : 'muted'));
          _kf.push(_kfTile('Consequence', _heroEsc(_csqLabel), _imp ? (_imp.charAt(0) + _imp.slice(1).toLowerCase()) + ' impact' : '', _impTone));
        }
      } catch (e) {  }

      const _kfShown = _kf.slice(0, 4);
      const keyFactsHTML = _kfShown.length
        ? `<section class="v11-section"><div class="v11-eyebrow"><h3>At a glance</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div><div class="ev-snap ev-snap--summary" style="grid-template-columns:repeat(${_kfShown.length},minmax(0,1fr))">${_kfShown.join('')}</div></section>`
        : '';

      document.getElementById('result-panel').innerHTML = `
        <div class="card" id="result-header">
          <!-- Tab bar is the first thing in the result card. Criteria carries
               a loading dot until stage 2 completes (see renderStage2Criteria
               → markCriteriaTabReady). -->
          <div class="result-tabs" role="tablist">
            <button class="result-tab active" data-tab="summary"           onclick="showResultTab('summary')">Summary</button>
            <button class="result-tab"        data-tab="variant"           onclick="showResultTab('variant')">Variant</button>
            <button class="result-tab"        data-tab="gene-context"      onclick="showResultTab('gene-context')">Gene</button>
            <button class="result-tab"        data-tab="protein"           onclick="showResultTab('protein')">Protein</button>
            <button class="result-tab"        data-tab="criteria"          onclick="showResultTab('criteria')">ACMG Criteria<span class="result-tab-dot" data-state="loading" title="Criteria still loading"></span></button>
            <button class="result-tab"        data-tab="export"            onclick="showResultTab('export')">Export</button>
            <!-- Result actions, right-aligned in the tab ribbon. -->
            <div class="result-tab-actions">
              <button type="button" class="result-tab-action" onclick="if(typeof returnToLanding==='function')returnToLanding();" title="Edit the query and re-run">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M11 4H4a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14a2 2 0 0 0 2-2v-7"/><path d="M18.5 2.5a2.12 2.12 0 0 1 3 3L12 15l-4 1 1-4 9.5-9.5z"/></svg>
                <span class="lbl">Edit query</span>
              </button>
              <button type="button" class="result-tab-action result-tab-action--primary" onclick="if(typeof newVariant==='function')newVariant();" title="Start a new variant">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M12 5v14M5 12h14"/></svg>
                <span class="lbl">New variant</span>
              </button>
            </div>
          </div>

          <div class="result-tab-pane" data-tab="summary">
            <!-- ════════ R2 EDITORIAL SUMMARY ════════
                 TWO-COLUMN report (R2): a main editorial column (verdict
                 hero → Clinical summary → Key evidence) and a compact,
                 sticky right-rail "Case" recap. The full query-details
                 grid + external-tool links no longer live here — case
                 context is condensed into the rail and the external-tool
                 deeplinks moved to the bottom of the Variant tab. All tier
                 colour / score / criteria IDs preserved for the post-render
                 JS that patches them at stage-2 reveal. -->
            <div class="v11-summary" data-class="${classText}">

              <!-- ───────── MAIN COLUMN ───────── -->
              <div class="v11-summary__main">

              <!-- (a) VERDICT HERO — ported from R2 .verdict (R2 markup
                   ~650-704). #summary-class-card keeps the data-class
                   tier-colour mechanism; a faint maroon ECG trace draws
                   along the top hairline on load (.v11-verdict__trace),
                   and a small ECG sits to the RIGHT of the tier. -->
              <div class="v11-verdict summary-class-card" id="summary-class-card" data-class="${classText}">
                <span class="v11-verdict__accent" aria-hidden="true"></span>
                <div class="v11-verdict__body">
                  <div class="v11-verdict__topline">
                    <span class="v11-vchip v11-vchip--tier"><span class="v11-pulse" aria-hidden="true"></span>HeartVar classification</span>
                  </div>

                  <div class="v11-verdict__line">
                    <div>
                      <!-- Large tier — coloured by tier via the data-class
                           selectors on #summary-class-card. -->
                      <div class="summary-class-tier v11-verdict__tier" id="summary-class-tier">${classText}</div>

                      <!-- Confidence pill + dot + ACMG/AMP score sub-line. -->
                      <div class="v11-verdict__sub">
                        <span class="summary-class-pill" id="summary-class-pill">${stage1.confidence || '—'} confidence</span>
                        <span class="dotsep" aria-hidden="true">·</span>
                        <span class="acmg-lbl">ACMG/AMP</span>
                        <span class="v11-verdict__pts"><span class="summary-score-hero" id="summary-score-hero">${scoreText}</span></span>
                      </div>
                    </div>
                  </div>

                  <!-- Variant identity line — gene + HGVS c./p. +
                       transcript (MANE badge) + genomic g. + build. -->
                  ${verdictVariantHTML}
                </div>
              </div>

              <!-- hg19 → hg38 liftover note (only when GRCh37 input was
                   lifted over; empty string otherwise). -->
              ${liftoverNoteHTML}

              <!-- Genomic-HGVS caution (only when a g. HGVS was entered;
                   empty string otherwise). Mutually exclusive with the
                   liftover note above (HGVS vs coordinate input). -->
              ${genomicHgvsNoteHTML}

              ${stage1.borderline_reasoning ? `
              <!-- Borderline VUS-vs-LP reasoning surface — collapsed by
                   default, labelled as an internal assessment. -->
              <details class="summary-borderline v11-borderline">
                <summary class="summary-borderline-head">
                  <span class="summary-borderline-icon" aria-hidden="true">⚖</span>
                  <span>Borderline assessment</span>
                </summary>
                <div class="summary-borderline-body">${String(stage1.borderline_reasoning).replace(/&/g, '&amp;').replace(/</g, '&lt;')}</div>
              </details>
              ` : ''}

              <!-- (b) CLINICAL SUMMARY — narrative prose. (Key-evidence
                   section removed — the met criteria live in the Criteria tab.) -->
              <section class="v11-section">
                <div class="v11-eyebrow"><h3>Clinical summary</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
                <p class="summary-prose">${linkifyPmids(summaryProse)}</p>
              </section>

              <!-- (b2) CARRIER STATUS — patient-level genotype interpretation
                   for recessive / X-linked / dual-MOI genes (empty string when
                   not applicable). Informational; never changes the tier. -->
              ${carrierBandHTML}

              ${keyFactsHTML}

              <!-- ClinGen eRepo VCEP verdict placeholder — INDEPENDENT
                   read-only reference, populated by the erepo_verdict SSE
                   event. Kept in the DOM but hidden (CSS) so the handler can
                   still find it. -->
              <div id="summary-erepo-panel"></div>

              </div><!-- /.v11-summary__main -->

              <!-- ───────── RIGHT RAIL (R2 .rail, sticky) ─────────
                   (1) ACMG/AMP point-score gauge — carries the canonical
                       #summary-score-value the stage-2 reveal writes into.
                   (2) "Query details" railcard = the full case recap built by
                       _buildCaseRecapHTML (always present). -->
              <aside class="v11-summary__aside" aria-label="At a glance">
                <div class="v11-railcard v11-gaugecard" id="summary-gauge-card">
                  <div class="v11-railcard__hd"><span class="rk">ACMG/AMP point score</span></div>
                  <div class="v11-gauge">
                    <div class="v11-gauge__num"><span class="summary-score-value" id="summary-score-value">${scoreText}</span></div>
                    <div class="v11-scaleband">
                      ${(function(){ const p = _scaleArrowPct(classText); return p == null ? '' : `<div class="v11-scaleband__arrow" id="summary-scale-arrow" style="left:${p}%"></div>`; })()}
                      <div class="v11-scaleband__track" aria-hidden="true">
                        <div class="v11-scaleband__seg v11-seg-b"></div>
                        <div class="v11-scaleband__seg v11-seg-lb"></div>
                        <div class="v11-scaleband__seg v11-seg-vus"></div>
                        <div class="v11-scaleband__seg v11-seg-lp"></div>
                        <div class="v11-scaleband__seg v11-seg-p"></div>
                      </div>
                      <div class="v11-scaleband__ticks"><span>Benign</span><span>VUS</span><span>Pathogenic</span></div>
                    </div>
                  </div>
                </div>
                ${caseRecapHTML ? `<div class="v11-railcard"><div class="v11-railcard__hd"><span class="rk">Query details</span></div>${caseRecapHTML}</div>` : ''}
              </aside>

            </div>
          </div>

          <div class="result-tab-pane" data-tab="variant" style="display:none">
            ${_renderEvidence(evidence)}
          </div>

          <div class="result-tab-pane" data-tab="gene-context" style="display:none">
            ${_renderGeneContextTab(evidence, stage1.gene_context || {}, gene)}
          </div>

          <div class="result-tab-pane" data-tab="protein" style="display:none">
            ${_renderProteinTab(evidence, gene, hgvs_c)}
          </div>

          <div class="result-tab-pane" data-tab="criteria" style="display:none">
            <!-- Criteria sections injected by renderStage2Criteria at reveal. -->
            <div id="criteria-area"></div>
          </div>

          <div class="result-tab-pane" data-tab="export" style="display:none">
            <div class="export-tab">
              <p class="export-intro">Download the full curation report or copy a formatted summary for your case notes.</p>
              <div class="export-cards">
                <button class="export-card" onclick="exportTXT()">
                  <span class="export-card__ic"><svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4"><path d="M9 1.5H4a1 1 0 0 0-1 1v11a1 1 0 0 0 1 1h8a1 1 0 0 0 1-1V5.5L9 1.5Z"/><path d="M9 1.5V5.5H13"/></svg></span>
                  <span class="export-card__title">Download TXT</span>
                  <span class="export-card__desc">Plain-text curation report</span>
                </button>
                <button class="export-card" onclick="exportXLSX()">
                  <span class="export-card__ic"><svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4"><rect x="2" y="2.5" width="12" height="11" rx="1.2"/><path d="M2 6.5h12M6 6.5v7M2 10h12"/></svg></span>
                  <span class="export-card__title">Download Excel</span>
                  <span class="export-card__desc">Workbook — Summary · Evidence · ACMG tabs</span>
                </button>
                <button class="export-card" onclick="copyReportToClipboard()">
                  <span class="export-card__ic"><svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4"><rect x="5.5" y="5.5" width="8.5" height="8.5" rx="1.2"/><path d="M10.5 5.5V3a1 1 0 0 0-1-1H3a1 1 0 0 0-1 1v6.5a1 1 0 0 0 1 1h2.5"/></svg></span>
                  <span class="export-card__title">Copy report</span>
                  <span class="export-card__desc">Formatted summary to clipboard</span>
                </button>
              </div>
              <section class="ev-exttools v11-section export-exttools">
                <div class="v11-eyebrow"><h3>External tools</h3><span class="v11-eyebrow__ln" aria-hidden="true"></span></div>
                ${_buildExternalToolsHTML(evidence)}
              </section>
            </div>
          </div>
        </div>
      `;

      const cardEl = document.getElementById('result-header');
      if (cardEl) cardEl.classList.add('load-fade-in');

      hvChatReset();
      hvChatSyncFab();
    }

    function renderStage2Criteria(criteria, gene, hgvs_c, evidence, variantId) {
      const DB_LINKS = _buildDbLinks(gene, hgvs_c, evidence, variantId);

      _hvInitOverrides(
        criteria, gene, _stripHgvsPrefix(hgvs_c || ''),
        (lastVariant && lastVariant.inheritance_input) || '',
        lastResult ? lastResult.points_total : null,
        lastResult ? lastResult.classification : '',
        DB_LINKS,
      );

      const area = document.getElementById('criteria-area');
      if (area) {
        area.innerHTML = _hvCriteriaSectionsHTML();

        area.classList.remove('gen-reveal');

        void area.offsetWidth;
        area.classList.add('gen-reveal');
      }
      markCriteriaTabReady();
    }

    function _hvCriteriaSectionsHTML() {
      const criteria = HV_OVR.effective;
      const DB_LINKS = HV_OVR.dbLinks;
      const met   = criteria.filter(c =>
        c.status === 'met' && (c.direction === 'pathogenic' || c.direction === 'benign')
      );
      const insuf = criteria.filter(c => c.status === 'insufficient_data');

      const notAssessed = criteria.filter(c => c.status === 'not_assessed');
      const na    = criteria.filter(c => c.status === 'na' || c.status === 'not_met');

      return _hvCriteriaVerdictHTML() +
        _rsec('Criteria met', met, DB_LINKS, { kind: 'met' }) +
        _rsec('Insufficient data to assess', insuf, DB_LINKS, { kind: 'insuf' }) +
        _rsec('Not assessed', notAssessed, DB_LINKS, { kind: 'na' }) +
        _rsec('Not applicable / not met', na, DB_LINKS, { kind: 'na' });
    }

    function _hvTierTone(tier) {
      if (tier === 'Pathogenic' || tier === 'Likely pathogenic') return 'warn';
      if (tier === 'Likely benign' || tier === 'Benign') return 'good';
      if (!tier || tier === '—') return 'muted';
      return 'info';
    }

    function _hvCriteriaVerdictHTML() {
      const n = _hvOverrideCount();
      const adj = HV_OVR.adjusted;
      const isAdj = !!(n && adj);
      const pending = !!n && !adj;
      const baseTier = HV_OVR.baseTier || '—';
      const basePts = Number.isFinite(HV_OVR.basePoints) ? formatCritPoints(HV_OVR.basePoints) : '—';
      const metN = HV_OVR.base.filter(c =>
        c.status === 'met' && (c.direction === 'pathogenic' || c.direction === 'benign')
      ).length;

      const tile = (k, v, sub, tone) =>
        `<div class="ev-snap__tile" data-tone="${tone}">
          <div class="ev-snap__k">${k}</div>
          <div class="ev-snap__v is-text">${v}</div>
          ${sub ? `<div class="ev-snap__sub">${sub}</div>` : ''}
        </div>`;

      const tiles = [tile(
        'HeartVar classification', _hvEsc(baseTier),
        `${basePts} pts <span class="nb">·</span> ${metN} criteria met`,
        _hvTierTone(baseTier),
      )];
      if (pending) {
        tiles.push(tile('Curator-adjusted', '<span class="ev-snap__pend">recomputing…</span>', '', 'muted'));
      } else if (isAdj) {
        const adjPts = formatCritPoints(adj.points_total);
        const same = adj.classification === baseTier;
        tiles.push(tile(
          'Curator-adjusted', _hvEsc(adj.classification),
          `${adjPts} pts <span class="nb">·</span> ${n} ${n === 1 ? 'change' : 'changes'}`
            + (same ? ' <span class="nb">·</span> same tier' : '')
            + ` <button type="button" class="crit-ovr-reset" onclick="hvResetOverrides()">Reset</button>`,
          _hvTierTone(adj.classification),
        ));
      }

      const ident = [HV_OVR.gene, HV_OVR.hgvsC].filter(Boolean).map(_hvEsc);
      return `<div class="tab-overview">
        <div class="tab-overview__head">
          <div class="tab-overview__title">${ident.length
            ? `<span class="tab-overview__gene">${ident[0]}</span>${ident[1] ? ' ' + ident[1] : ''}`
            : 'ACMG/AMP criteria'}</div>
          <div class="tab-overview__sub">ACMG/AMP criteria — Richards 2015 definitions, Tavtigian 2020 points</div>
        </div>
        <div class="ev-snap ev-tilerow">${tiles.join('')}</div>
      </div>`;
    }
