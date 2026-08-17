# hEDS Biomarker Discovery

A literature-based discovery (LBD) pipeline that mines PubMed for genes indirectly
implicated in hypermobile Ehlers-Danlos Syndrome (hEDS) biology but never directly
studied in hEDS itself — surfacing a small, ranked shortlist of hypotheses for further
research, not a diagnostic test.

## Why this exists

hEDS is the one Ehlers-Danlos subtype with no clinically validated diagnostic
biomarker. Diagnosis today is purely clinical (the 2017 international criteria:
Beighton score, systemic manifestation checklist, family history, exclusion of other
conditions). Real genetic hits have started to emerge — see [Prior work](#prior-work) —
but nothing yet resolves into a test.

Given that, the honest and interesting angle isn't pretending to reverse-engineer a
biomarker from a symptom checklist (that's just re-deriving the existing clinical
criteria). It's **literature-based discovery**: apply Swanson's classic ABC model —
if A (hEDS) is linked to B (some phenotype/pathway) in the literature, and B is linked
to C (a gene) in the literature, but A and C are never linked directly, C is a
literature-implied candidate worth a second look — using only free public data.

## Method

1. **Fetch** (`fetch_literature.py`) — Query PubMed for hEDS/HSD-specific terminology
   (deliberately excludes bare "Ehlers-Danlos syndrome" and other-subtype terms, so
   classical/vascular-EDS genes don't pollute the corpus), then pull entity annotations
   (genes, diseases, chemicals) for each abstract from NCBI's **PubTator3** API.
2. **Graph** (`build_graph.py`) — Build a co-occurrence knowledge graph: nodes are
   bioconcepts, edges are weighted by how often two concepts appear in the same
   abstract, with a PMI score for association strength.
3. **Discover** (`discover_candidates.py`) — ABC-model link prediction: genes directly
   co-mentioned with hEDS/HSD are "already studied" and excluded; genes reachable only
   through a shared bridge concept (2 hops, never 1) are ranked as candidates via the
   Adamic-Adar index, which down-weights generic/promiscuous bridge terms (e.g. "pain")
   more than specific ones — a known false-positive source in co-occurrence-based LBD.
4. **Validate** (`validate_candidates.py`) — Cross-check the candidate shortlist against
   independent public evidence: STRING (protein-interaction connectivity to known
   hEDS-relevant genes), GTEx (expression in connective-tissue-proxy tissues,
   principally cultured fibroblasts), and ClinVar (existing variant evidence).
5. **Output**: `results/candidates.csv` — a ranked shortlist (~10-20 genes) with
   evidence columns and a plain-language rationale per row.

All data sources are free, public, no-auth APIs. No patient data is used anywhere in
this pipeline.

## Prior work

This project doesn't operate in a vacuum — 2024-2025 brought the first real hEDS
genetic findings, and the reference gene sets in `validate_candidates.py` are built
from them, not from other EDS subtypes' genes:

- **KLK gene family cluster** (esp. *KLK15*) — recurrent variant found via whole-exome
  sequencing of 200 hEDS patients, reproduces connective-tissue defects in mouse
  knock-ins. [Norris Lab, 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11213194/)
- **ACKR3, SLC39A13** — the first hEDS GWAS meta-analysis (~1,800 cases, ~5,000
  controls) found genome-wide-significant signal at these loci, pointing to a
  neuroimmune/stromal model rather than a single collagen gene.
  [Petrucci-Nelson et al., medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.09.19.25336146v1)
- **MIA3** — a separately proposed 2024-2025 candidate, still unresolved.
- **TNXB** (partial/haploinsufficient) — the longest-studied lead, but explains only
  ~1% of cases, and serum tenascin-X failed as a screening test in follow-up work.
  Kept as a weak legacy positive control, not a strong reference.
- **Plasma ECM-fragmentation signature** (fibronectin/collagen-I/tenascin fragments) —
  the closest thing to an actual proposed biomarker in the literature, still
  unvalidated. [Ritelli et al.](https://pubmed.ncbi.nlm.nih.gov/39225014/)
- **HEDGE Study** — the Ehlers-Danlos Society's ongoing population-scale sequencing
  effort (1,000 hEDS patients); results still emerging.

No prior work applying literature-based discovery, knowledge-graph mining, or the
Swanson ABC model to hEDS specifically was found — this project is, as far as this
review could tell, the first application of this method to this disease.

### Angles deliberately not pursued, and why

- **Symptom/comorbidity subtype clustering** — already done, repeatedly, at real scale
  on real patient cohorts (Mayo Clinic, K-means/UMAP on ~2,100-2,700 patients,
  [Petrucci et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38779137/)). Re-running
  clustering on synthetic or self-collected data would replicate published results with
  worse data, not extend them.
- **Digital/wearable biomarkers** (AI-scored Beighton assessment from video, wearable
  HRV/autonomic monitoring) — real, active, and promising research directions, but
  every dataset behind them is access-gated ("available upon reasonable request"), not
  publicly downloadable. Confirms literature-mining is the feasible public-data lane
  for a project like this, not a workaround.
- **SemMedDB** (semantic-predicate literature relations, a higher-precision alternative
  to raw co-occurrence) — considered, but it was deprecated December 2024, is now a
  frozen archive current only through May 2024 (would miss the 2025 GWAS and 2024
  KLK15 findings above), and requires a UMLS/UTS license. Not worth the tradeoff for a
  lean, current pipeline; PubTator3 co-occurrence with a hub-node degree cutoff is used
  instead.

## Usage

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt

python pipeline.py --top-n 20
```

Runtime is dominated by API rate-limit pacing (~1 req/sec), not compute — expect
roughly 10-15 minutes end to end for the full ~1,300-abstract corpus. Results are
cached under `cache/` so re-runs are fast unless `--skip-cache` is passed.

Flags: `--retmax` (PubMed result cap), `--min-weight` (minimum co-occurring abstracts
to count an edge), `--max-bridge-degree` (hub-node exclusion threshold), `--top-n`
(shortlist size).

## Sample results (this run)

On the corpus available at the time of writing (1,301 annotated abstracts, 1,345
entity nodes, 20,238 co-occurrence edges), 2 clean hEDS/HSD-specific seed nodes were
found. 18 genes were correctly bucketed as "already studied" — including **SLC39A13**
and **TNXB**, both genuinely-studied hEDS-relevant genes, confirming the pipeline's
direct-neighbor exclusion works as intended.

12 candidates survived to the shortlist (`results/candidates.csv`). The top-ranked
candidate by literature bridge score, **C1R (complement C1r)**, is notable: it wasn't
put there by design, but it independently lines up with a 2025/2026 hEDS proteomics
paper reporting complement-cascade proteins differentially expressed in patient plasma
— an encouraging (if small) convergence between this literature-mining approach and
independent wet-lab evidence. Most of the remaining candidates (COL3A1, COL6A3, PLOD1,
B4GALT7, DSE, FKBP14, AEBP1) are genes that define *other* rare connective-tissue/EDS-
spectrum disorders — biologically coherent as a "look at this related gene family"
hypothesis, if not individually surprising.

**Honest limitation**: the very newest 2024-2025 hEDS findings are not fully
represented yet. ACKR3 (the 2025 GWAS hit) is present in the corpus (graph degree 19)
but doesn't clear the co-occurrence weight threshold to register as "already studied" —
consistent with only a handful of very recent papers discussing it in an hEDS context
so far. KLK15 and MIA3 don't appear in the corpus at all, most likely because PubTator3
hasn't finished indexing their (2024-2025) source papers yet. This is a real
corpus-freshness constraint of the underlying data source, not a pipeline bug, and
should resolve as PubTator3's index catches up and the pipeline is re-run over time.

## Verification

- `tests/test_build_graph.py` — a synthetic-graph smoke test proving the core
  discovery logic actually distinguishes "directly studied" from "literature-implied
  but never directly studied" candidates, with no network dependency. Run with
  `venv\Scripts\python tests\test_build_graph.py`.
- **Sanity check on real output**: known hEDS-relevant genes that appear in the real
  corpus (KLK15, ACKR3, SLC39A13, TNXB) should land in the excluded "already studied"
  bucket, not the final shortlist. If any of them surface as "novel," that's either a
  corpus-date-cutoff issue or a genuinely interesting observation worth a closer look.
- **Required human spot-check**: for the top candidates in `results/candidates.csv`,
  read the rationale and the underlying bridge concept, and confirm it reflects real
  biology rather than an annotation artifact. This is a hypothesis-generation tool —
  every output needs that check before it means anything.

## Limitations and ethics

This is a **hypothesis-generating research tool, not a diagnostic instrument**. No
confirmed diagnostic biomarker for hEDS exists in the published literature today. Every
candidate this pipeline surfaces is, at best, a starting point for wet-lab or clinical
follow-up — never something to act on directly, and never something that should be
used to include, exclude, or influence an actual diagnosis.
