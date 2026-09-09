# Data-cache build scripts

These scripts build the local data caches HeartVar reads at runtime. They are
deploy-host tooling, not part of the running app: `scripts/` is excluded
from the Docker image (`.dockerignore`), so run them on the host that prepares
the mounted data volume, after `pip install -r requirements.txt`.

Outputs land in two places:
- **`data/`** (gitignored, large, regenerable) — the SQLite DBs, the SpliceAI
  slice, the cardiac-panel BED, the AlphaFold structures.
- **`backend/data/`** (tracked, small curated reference) — a few JSON files
  (e.g. `cvd_gene_panel.json`, `clingen_gene_validity.json`, `gene_mechanism.json`)
  ship in the repo; the scripts only regenerate them.

## Build order

Most scripts are independent and can run in any order/in parallel. Three
dependencies matter:

1. **`build_uniprot_db.py` → `build_alphafold_structures.py`.** AlphaFold
   structure selection looks up UniProt accessions in `data/uniprot.db`. Build
   UniProt first.
2. **`build_cvd_gene_panel.py` (only if regenerating the panel) → `build_alphafold_structures.py`,
   `build_gene_mechanism.py`, `build_clingen_gene_validity.py`.** These three
   filter to `backend/data/cvd_gene_panel.json`. That file is **tracked**, so
   they work out of the box; only if you regenerate the panel must you do it
   before rebuilding the three consumers.
3. **`data/cardiac_panel.grch38.bed` → `build_gnomad_freq_db.py`, `build_spliceai_db.py`.**
   Both slice to the cardiac-panel intervals. The BED is auto-built on first
   use (resolved from Ensembl REST by `_cardiac_panel.py`), so simply running
   either slicer builds it; it is then reused. The runtime gnomAD and SpliceAI
   clients also read this BED for panel membership.

Everything else (`build_clinvar_db`, `build_biogrid_db`, `build_fetal_heart_db`,
`build_gtex_db`, `build_medgen_db`, `build_mgi_db`, `build_hpo_labels`,
`build_gnomad_constraint_db`, `build_panelapp_snapshot`, `build_hgnc_alias_db`,
`build_erepo_dump`) has no inter-script dependency.

**AlphaMissense** is not built by a script: download the precomputed table
(`AlphaMissense_hg38.tsv.gz` + `.tbi`, ~600 MB, from Zenodo record 8208688) and
point `ALPHAMISSENSE_PATH` at it.

## Data-source licences

Most sources are open (CC0/CC-BY/public domain): ClinVar, ClinGen, HGNC,
HPO, MGI, GTEx, UniProt, GenCC. Three carry use restrictions the institute
should note:

- **SpliceAI slice** (`build_spliceai_db.py`) — source VCF is **CC BY-NC 4.0
  (non-commercial)**.
- **AlphaMissense** — **CC BY-NC-SA 4.0 (non-commercial)**.
- **PanelApp** — no explicit redistribution licence.
