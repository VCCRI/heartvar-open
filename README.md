# HeartVar

HeartVar is a research tool for preliminary ACMG/AMP classification of variants
found in cardiovascular disease cohorts. Given a gene,
variant, and clinical context, it queries around 20 public databases in
parallel and assembles the evidence into a per-criterion view. Classification is hybrid. Nineteen criteria are decided
programmatically in Python from structured data; seven that require interpretation are decided by a single LLM call, which is optional and off by
default. 

## What's in here

| Entry | What it is |
|---|---|
| `backend/` | The FastAPI service, the ACMG engine (`acmg/`), the per-source clients (`clients/`), shipped reference JSON (`data/`) and the test suite (`tests/`). |
| `static/` | Frontend CSS and JavaScript, the vendored 3Dmol viewer, and `acmg_constants.json`. |
| `index.html` | The whole single-page frontend. |
| `scripts/` | Deploy-host tooling that builds the local data caches. |
| `Dockerfile` | Runtime image. `Dockerfile.builder` is the data-mirror build job. |
| `.env.example` | Every environment variable. |

## Licence and data sources

The code is MIT, see [LICENSE](LICENSE). The data the build scripts fetch is
not covered by that licence and several sources carry their own terms. The
full per-source list is in [scripts/README.md](scripts/README.md).

## Citing this work

> HeartVar: An LLM-Assisted Tool for Clinical Classification of Variants in Cardiovascular Disease Cohorts.
> Jamie-Lee Thompson, Debjani Das, Sally L Dunwoodie, Eleni Giannoulatou.
> bioRxiv 2026.09.10.750569; doi: https://doi.org/10.64898/2026.09.10.750569

## Acknowledgements

Supported by Anthropic's AI for Science program. Built on public resources from
Ensembl/EMBL-EBI, the Broad Institute (gnomAD, SpliceAI), Google DeepMind
(AlphaMissense), NCBI (ClinVar, PubMed, PubTator3, MedGen), UniProt, EBI ProtVar,
GTEx, PanelApp Australia, GenCC, the Alliance of Genome Resources/MGI, BioGRID,
Open Targets, ClinGen (eRepo) and AlphaFold. The CHDgene gene list is maintained
by Victor Chang Cardiac Research Institute.
