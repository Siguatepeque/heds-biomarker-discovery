# hEDS Biomarker Discovery

![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-2c7a64)
[![Full write-up: live](https://img.shields.io/badge/full%20write--up-live-8c611f)](https://siguatepeque.github.io/heds-biomarker-discovery/)

This pipeline mines the research literature on hypermobile Ehlers-Danlos Syndrome
(hEDS) for genes that show up near the disease in the writing without ever being
studied in it directly. What comes out is a short, ranked list of genes worth a
second look, which is a starting point for research and not a diagnostic test of
any kind.

**Full write-up:** [siguatepeque.github.io/heds-biomarker-discovery](https://siguatepeque.github.io/heds-biomarker-discovery/)
reads more like a normal article. This README is the developer-facing version.

## Abstract

There is no lab test for hypermobile Ehlers-Danlos Syndrome (hEDS); diagnosis is
still checklist-based under the 2017 criteria. This project applies literature-based
discovery (Swanson's ABC model) to PubMed to surface candidate genes that co-occur
with hEDS only indirectly in the research literature, then validates each one against
STRING, GTEx, and ClinVar. Backtested on literature published only through 2021, more
than three and a half years before the first hEDS GWAS meta-analysis, the pipeline
independently flagged SLC39A13 as a candidate; that same gene later reached
genome-wide significance in the 2025 meta-analysis [[1]](#references). This is a
hypothesis-generating tool, not a diagnostic one.

## Key results

> [!IMPORTANT]
> **The pipeline has a track record, and it says two more genes are worth checking
> right now.** Backtested on pre-2022 literature only, it flagged **SLC39A13** years
> before an independent GWAS confirmed it. **VWF** and **C1R** are on today's
> candidate list through that same kind of indirect evidence: unconfirmed, but in
> the exact position SLC39A13 was in before anyone checked it.

| Gene | Status | Evidence |
|---|---|---|
| **SLC39A13** | Backtest-confirmed | Flagged from pre-2022 literature alone; reached genome-wide significance in the 2025 GWAS meta-analysis [[1]](#references) |
| **VWF** | Live candidate, unconfirmed | Real bridge via a documented bleeding/bruising symptom link; not yet directly studied in hEDS |
| **C1R** | Live candidate, unconfirmed | Matches an independent 2025-26 proteomics finding the pipeline was never pointed at |

> [!NOTE]
> This is a backtest, not a prediction, and it only held up after two rounds of
> finding and fixing real bugs in my own code, see [The real test](#the-real-test).
> Full candidate list and methodology in [What we found](#what-we-found), and the
> weak points in [What to be skeptical of](#what-to-be-skeptical-of).

## Contents

1. [Abstract](#abstract)
2. [Key results](#key-results)
3. [The problem, in plain terms](#the-problem-in-plain-terms)
4. [The theory](#the-theory)
5. [How the robot works](#how-the-robot-works)
6. [The real test](#the-real-test)
7. [Negative control](#negative-control)
8. [What we found](#what-we-found)
9. [What's new here](#whats-new-here)
10. [What backs this up](#what-backs-this-up)
11. [What I decided not to do, and why](#what-i-decided-not-to-do-and-why)
12. [What to be skeptical of](#what-to-be-skeptical-of)
13. [What this is not](#what-this-is-not)
14. [Running it](#running-it)
15. [Does it actually work?](#does-it-actually-work)
16. [References](#references)

## The problem, in plain terms

hEDS is the one Ehlers-Danlos subtype nobody has a lab test for. You get diagnosed by
checklist: how far your joints bend, which of a long list of other symptoms you have,
whether it runs in your family, and whether a doctor has ruled out everything else it
could be instead. That is the 2017 diagnostic standard, and it is still the standard.
Real genetic leads have started showing up in just the last couple of years (see
[What backs this up](#what-backs-this-up)), but nothing has turned into an actual test
yet.

## The theory

I did not want to build the obvious thing: a classifier trained on symptom checklists,
which just re-learns the diagnostic criteria and calls it a "biomarker." Instead I went
after something nobody seems to have tried for this disease: literature-based
discovery.

It is an old idea. Don Swanson used it in the 1980s to connect fish oil and Raynaud's
disease before anyone had run the clinical trial, just by noticing that both terms kept
showing up near a third concept in completely separate sets of papers. The logic
carries over cleanly here. If hEDS is linked to some phenotype B in the literature, and
B is linked to gene C in a totally different set of papers, but nobody has ever written
a paper connecting hEDS to C directly, then C is worth a look. All of this runs on free,
public data. No lab access required, no patient data touched at any point.

## How the robot works

![The pipeline: fetch, graph, discover, validate](docs/assets/pipeline-diagram.svg)

1. **Fetch** (`fetch_literature.py`). Pulls hEDS and HSD specific papers from PubMed.
   I deliberately left out bare "Ehlers-Danlos syndrome" and other-subtype terms, so
   genes belonging to classical or vascular EDS do not sneak into the corpus. Each
   abstract's genes, diseases, and chemicals get tagged by NCBI's PubTator3 service,
   so I did not have to write my own entity recognizer.
2. **Graph** (`build_graph.py`). Turns all of that into a co-occurrence graph. Two
   concepts get an edge if they show up in the same abstract, weighted by how often,
   plus a PMI score for how surprising that co-occurrence is.
3. **Discover** (`discover_candidates.py`). The actual ABC-model step. Anything
   directly co-mentioned with hEDS in the same abstract is filed as "already studied"
   and set aside. What is left, genes reachable only through a shared bridge concept
   and never directly, gets ranked by an Adamic-Adar score. That score is common
   neighbor math that discounts generic bridge terms (things like "pain," which
   connect to almost everything) more than specific ones. This distinction matters a
   lot. Raw co-occurrence counting alone is a known source of false positives in this
   kind of work, and it bit me twice during development (see [What to be skeptical
   of](#what-to-be-skeptical-of)).
4. **Validate** (`validate_candidates.py`). Before anything makes the final list, it
   gets checked against three independent public sources: STRING (does it physically
   interact with a known hEDS-relevant protein), GTEx (is it actually expressed in
   connective tissue), and ClinVar (does it have any documented variant record at all).
5. Out comes `results/candidates.csv`, usually 10 to 20 genes, each with its evidence
   and a plain language reason it made the list.

## The real test

Before trusting any candidate this pipeline produces, it's worth asking whether its
shortlists mean anything at all. The only honest way to check is to go back in time
and see whether an old shortlist would have held up. `backtest.py` reruns the exact
same pipeline using only literature published before a chosen cutoff year (the cached
documents each carry their own publication date, so this needs no new network calls),
and checks whether genes that real 2024 to 2025 hEDS genetics later confirmed were
already showing up as candidates before that confirmation existed.

```
python backtest.py --cutoffs 2019 2020 2021 2022 2023 2024 2025 2026
```

![Timeline: SLC39A13 flagged from pre-2022 literature, more than three and a half years before a 2025 GWAS found genome-wide-significant signal at the same gene](docs/assets/backtest-timeline.svg)

Using only papers published through 2021, more than three and a half years before the
first hEDS GWAS meta-analysis went up on medRxiv (19 September 2025) [[1]](#references), the pipeline
flags SLC39A13 as a candidate through a purely indirect literature path. In September
2025, an independent GWAS meta-analysis (about 1,800 cases, 5,000 controls, a
completely different method: genotyping, not literature) found genome-wide significant
signal at that exact gene.

I did not trust that number on the first pass, and I do not think anyone should trust
a result like this without trying to break it. Getting here took two rounds of finding
real bugs in my own code. The full story, including exactly which papers, which bridge
concepts, and hand-verified arithmetic, is in `results/backtest.txt` and told in full
in the [research paper](https://siguatepeque.github.io/heds-biomarker-discovery/).

One more honest thing, checked just now. SLC39A13 is still classified as a candidate
today, not yet "already studied," even eleven months after the GWAS. Exactly one
document in the whole corpus directly co-mentions hEDS and SLC39A13: the GWAS preprint
itself. No second paper has echoed the connection yet. That is a separate finding
about the corpus rather than about the method: the literature index can lag a real
discovery by the better part of a year after it is published, and that lag is the
gap this kind of tool sits in.

That is the positive control: a real gene, correctly flagged early. A method that only
ever says yes is not proof of anything, so it needs a negative control too.

## Negative control

The real test above shows the method can get ahead of the field on hEDS. It says
nothing about whether the method would do that for *any* disease, whether or not a real
gene is there to find. So I picked fibromyalgia and ran the identical pipeline
(`control_disease.py`, `controls.py`, `control_backtest.py`) against it, unchanged
except for swapping the PubMed query and seed terms.

Fibromyalgia is a deliberate, harder case: a literature almost four times the size of
hEDS's (4,938 documents fetched, retmax=5000), dominated by generic hub terms (pain,
fatigue, cytokines) that are exactly what the Adamic-Adar weighting is built to
discount, and, as of when this project's design was decided, no gene that had reached
genome-wide significance. That stopped being true partway through this project:
Kerrebijn et al., preprinted in September 2025 and published in *Nature Medicine* in
July 2026 [[5]](#references), found the first genome-wide-significant fibromyalgia
loci: HTT, GPR52, DCC, DRD2, NCAM1, MDGA2, CELF4, and CAMKV. Rather than treat that as a problem, I used it: it turns
fibromyalgia into a second, independent check of the same shape as the SLC39A13 story,
just run in the opposite direction. Does the pipeline also backfill these genes as
early candidates for fibromyalgia, the way it did for SLC39A13 in hEDS?

Across all 8 cutoff years, including today's full corpus, the answer is no. None of the
8 real fibromyalgia GWAS genes is ever flagged as a candidate. Seven sit as "present in
the corpus but below the bridge floor" even at the 2026 cutoff; the eighth, DCC, moves
straight to "already directly studied" once a 2026 document (almost certainly the GWAS
paper itself) mentions it alongside fibromyalgia by name, same pattern as SLC39A13's own
GWAS preprint in the real backtest. Zero false "years-early" flags on the disease this
control was built to stress-test.

The live pipeline run (`results/control_candidates.csv`) still surfaces 6 candidates for
fibromyalgia today: Nfe2l2, Bax, S100A7, PPARGC1A, Myog, IL23A. None of them match any
of the 8 real GWAS genes, and none match the older, weaker fibromyalgia candidate-gene
literature either (GCH1, COMT, OPRM1). What they have in common is thematic, not
genetic: Nrf2, Bax, and PGC-1alpha are recurring names in fibromyalgia's
oxidative-stress/mitochondrial-dysfunction mechanistic literature, a popular framing
that is not itself a genetic finding. That is a real, named limitation, not a clean
pass: bridge degree alone doesn't distinguish a specific-and-wrong bridge concept from a
specific-and-right one, and a large enough non-genetic literature can still clear the
hub-degree filter. Full cutoff-by-cutoff numbers, citations, and the established-gene
check are in `results/control_backtest.txt`.

## What we found

On the corpus I had while writing this (1,301 abstracts, 1,345 entities, 20,238
co-occurrence edges), the pipeline correctly filed 24 genes as "already studied,"
including COL3A1, COL6A3, SMAD3, and TNXB, and surfaced 10 candidates.

![Composite score per candidate gene, PLOD1 highest at 7.69 down to BCL2A1 at 2.78](docs/assets/candidates-chart.svg)

| Gene | Score | Why it's here |
|---|---:|---|
| PLOD1 | 7.69 | Defines kyphoscoliotic EDS |
| DSE | 5.40 | Defines musculocontractural EDS |
| B4GALT7 | 5.40 | Defines spondylodysplastic EDS |
| AEBP1 | 5.29 | Defines classical-like EDS |
| TGFB1 | 5.22 | The TGF-beta ligand itself |
| **SLC39A13** | 4.68 | The gene from the retrospective backtest above |
| FKBP14 | 3.81 | Defines kyphoscoliotic EDS, type 2 |
| VWF | 2.86 | Real bridge via "bleedings," a documented hEDS symptom |
| C1R | 2.83 | Matches an independent 2025-26 proteomics finding |
| BCL2A1 | 2.78 | Weak. Low expression, likely noise |

Topping the list is **PLOD1**, which defines kyphoscoliotic EDS. Most of the rest
(AEBP1, DSE, B4GALT7, FKBP14) are similarly genes that define other EDS-spectrum
disorders. That is a reasonable "maybe this whole gene family matters here too"
hypothesis, not individually shocking.

Three entries are worth a closer look on their own merits.

#### SLC39A13

This is the gene from the retrospective backtest above, years-early evidence later
confirmed by GWAS, still classified as a candidate today rather than "already
studied," even now.

#### C1R (complement C1r)

Independently lines up with a 2025 to 2026 proteomics paper that found
complement-cascade proteins showing up differently in hEDS patients' blood. I did not
point the pipeline at it. It came out of the graph structure on its own, and it
happens to match a completely separate wet-lab finding.

#### VWF (von Willebrand factor)

Reaches the list through a real, coherent bridge concept ("bleedings," a symptom term
appearing in 18 real hEDS papers), and bleeding or bruising tendency is a documented
hEDS comorbidity. It has not been directly studied in hEDS by name yet, but the path
here is clean, the same shape as the SLC39A13 case before anyone confirmed it.

Not everything on the list deserves that kind of confidence. **BCL2A1** is here too,
and I think it is mostly noise: low connective-tissue expression (3.5 TPM), reached
only through a chain of classical (not hypermobile) EDS papers, and no clear
biological reason to expect it. I am leaving it in the output because hiding a weak
result is worse than showing one, but it should be read as weak.

ACKR3 and MIA3, two of the newest 2024 to 2025 hEDS genes, do not show up in the
corpus at all yet. That is a data-freshness limit of the underlying literature index,
not a bug in the pipeline.

## What's new here

Before writing any code, I checked whether this had been done. As far as I could find,
nobody has pointed literature-based discovery or knowledge-graph mining at hEDS before.
Not even at neighboring diseases with a similar "diagnosis of exclusion" problem, like
long COVID or fibromyalgia. This angle, at least, seems to be new.

## What backs this up

I did not want to build this in a vacuum, so before writing any code I spent a while
figuring out what had already been tried. Turns out 2024 to 2025 was a genuinely big
couple of years for hEDS genetics, and the reference gene sets baked into
`validate_candidates.py` come from that work, not from other EDS subtypes' better
known genes.

- **KLK gene family**, especially KLK15. A recurrent variant turned up in whole exome
  sequencing of 200 hEDS patients, and it reproduces connective tissue defects when
  reintroduced in mice. [[2]](#references)
- **ACKR3 and SLC39A13**, the first hEDS GWAS meta-analysis, pointing toward a
  neuroimmune and stromal story rather than a single collagen gene. [[1]](#references)
- **MIA3**, another 2024 to 2025 candidate, still unresolved.
- **TNXB** (partial deficiency), the oldest lead in this space, though it only
  explains about 1 percent of cases, and serum tenascin-X testing failed as a
  screening tool when people actually tried it. Kept as a weak legacy check.
- **A plasma ECM-fragmentation signature** (fibronectin, collagen-I, and tenascin
  fragments), probably the closest thing to an actual proposed biomarker that exists
  right now, still unvalidated. [[3]](#references)
- **The HEDGE Study**, the Ehlers-Danlos Society's own large scale sequencing effort
  covering 1,000 patients, still producing results.

## What I decided not to do, and why

<details>
<summary>Three things worth naming, and why each one was skipped</summary>

**Clustering patients into subtypes by symptoms.** Already done, and done well. Mayo
Clinic ran K-means and UMAP on more than 2,100 real patients and published real
subtypes. [[4]](#references) I do not have that kind of data, and running the same
analysis on something worse would not add anything.

**Wearables and video-based biomarkers.** Genuinely exciting research: AI-scored
hypermobility from video, wearable heart-rate-variability tracking. Every dataset
behind it is "available on request," not public. That confirmed literature mining was
the right lane for a project built entirely on public data.

**Swapping in SemMedDB** instead of raw co-occurrence. SemMedDB gives semantic
relationships instead of "these two things appeared near each other," which sounds
strictly better. It was deprecated in December 2024 and has not been updated since, so
it would miss both the 2025 GWAS paper and the 2024 KLK15 finding above, and it needs a
UMLS license to even access. Not worth the tradeoff for something meant to stay current
and easy to run.

</details>

## What to be skeptical of

I would rather name the weak points myself than have you find them.

- **VWF and C1R are today's live candidates, and they deserve harder scrutiny than
  a single pass.** Both reach the list the way SLC39A13 did before its GWAS
  confirmation: through indirect evidence, not direct study. I traced VWF's bridge
  concept and it looks real, but it has not gone through the same two rounds of
  adversarial re-checking that SLC39A13's result did. C1R's match to the 2025-26
  proteomics finding is suggestive, not something I have independently re-verified.
  Run `backtest.py` and the graph queries in `results/backtest.txt` yourself before
  treating either as more than a lead.
- **BCL2A1 is flagged weak on purpose.** I left it in the shortlist instead of quietly
  dropping it. Read it as "the pipeline's noise floor," not as a finding.
- **The seed-matching patterns in `discover_candidates.py` are a hand-built list.** I
  found and fixed two false positives in that exact list while building this project.
  There could be a third I have not found. If you add a new disease term to the search
  scope, check it against a real, messy abstract before trusting the output.
- **The "already studied" bucket depends entirely on PubTator3's entity tagging being
  right.** It mostly is, but automated tagging makes mistakes, and I have not manually
  audited every one of the 24 genes in that bucket.
- **Run the tests before you trust a change.** `tests/test_build_graph.py` and
  `tests/test_backtest.py` cover the logic that has actually broken before. If you
  touch `discover_candidates.py`, run them first.

## What this is not

> [!WARNING]
> This is a hypothesis-generating research tool, not a diagnostic instrument. There
> is no confirmed diagnostic biomarker for hEDS in the published literature: not from
> this project, not from anyone, not yet. Nothing here should be used to include,
> exclude, or otherwise influence an actual diagnosis. At best, a candidate on this
> list is a starting point for someone with a lab to look closer.

## Running it

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt

python pipeline.py --top-n 20
```

Most of the runtime is just being polite to free APIs (about 1 request per second),
not actual computation. Figure 10 to 15 minutes for the full corpus. Everything gets
cached under `cache/`, so re-runs are quick unless you pass `--skip-cache`.

Other flags: `--retmax` (how many PubMed results to pull), `--min-weight` (how many
co-occurring abstracts before an edge counts), `--max-bridge-degree` (the hub-node
cutoff), `--top-n` (shortlist size).

Tested on Python 3.11.6; dependency versions in `requirements.txt` were pinned to
the latest release of each still compatible with Python 3.9+ as of 2026-08-18.

## Does it actually work?

- `tests/test_build_graph.py` is a small synthetic-graph test that checks the core
  logic actually separates "directly studied" from "implied but never studied," and
  covers the two seed-matching regressions found during development. No network
  required. Run it with `venv\Scripts\python tests\test_build_graph.py`.
- `tests/test_backtest.py` checks the date-filtering logic behind the retrospective
  validation above.
- `tests/test_controls.py` checks that the control disease's seed patterns
  (`controls.py`) don't cross-match hEDS text or vice versa, and that swapping the
  pipeline over to the control config and back is fully reversible.
- On real output, the sanity check is calibrated to what is actually established in
  the literature, not just what is true biologically: TNXB, COL3A1, and COL6A3 should
  land in "already studied" because they are repeatedly, directly co-mentioned with
  hEDS. SLC39A13 is expected to stay a candidate until a second independent paper
  echoes the GWAS finding directly.
- And the one I cannot automate away: read the top candidates yourself. Look at the
  bridge concept connecting them to hEDS and ask whether it is real biology or an
  annotation glitch. Nothing in the output means much until someone has done that.

## References

1. Petrucci-Nelson T, Guilhaumou S, Berrandou T-E, et al. Complex Genetics and
   Regulatory Drivers of Hypermobile Ehlers-Danlos Syndrome: Insights from
   Genome-Wide Association Study Meta-analysis. medRxiv. 2025 Sep 21. Preprint.
   doi:[10.1101/2025.09.19.25336146](https://doi.org/10.1101/2025.09.19.25336146)
2. Gensemer C, Petrucci T, Beck T, et al. KLK15 alters connective tissues in
   hypermobile Ehlers-Danlos syndrome. iScience. 2025;28(9):113343.
   doi:[10.1016/j.isci.2025.113343](https://doi.org/10.1016/j.isci.2025.113343)
3. Ritelli M, Chiarelli N, Cinquina V, et al. Bridging the Diagnostic Gap for
   Hypermobile Ehlers-Danlos Syndrome and Hypermobility Spectrum Disorders: Evidence
   of a Common Extracellular Matrix Fragmentation Pattern in Patient Plasma as a
   Potential Biomarker. Am J Med Genet A. 2025 Jan.
   doi:[10.1002/ajmg.a.63857](https://doi.org/10.1002/ajmg.a.63857)
4. Petrucci T, Barclay SJ, Gensemer C, et al. Phenotypic Clusters and Multimorbidity
   in Hypermobile Ehlers-Danlos Syndrome. Mayo Clin Proc Innov Qual Outcomes. 2024 Jun.
   doi:[10.1016/j.mayocpiqo.2024.04.001](https://doi.org/10.1016/j.mayocpiqo.2024.04.001)
5. Kerrebijn I, et al. The genetic architecture of fibromyalgia across 2.5 million
   individuals. medRxiv. 2025 Sep 19. Preprint.
   doi:[10.1101/2025.09.18.25335914](https://doi.org/10.1101/2025.09.18.25335914).
   Published: Nat Med. 2026 Jul 28.
   doi:[10.1038/s41591-026-04492-6](https://doi.org/10.1038/s41591-026-04492-6)
