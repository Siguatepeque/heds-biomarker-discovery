# hEDS Biomarker Discovery

A pipeline that reads PubMed for a living — specifically, it mines the literature on
hypermobile Ehlers-Danlos Syndrome (hEDS) looking for genes that show up *near* the
disease in the research, without ever being studied *in* it directly. The output is a
short, ranked list of hypotheses worth a second look. Not a diagnosis, not a test —
a lead list.

## Why I built this

hEDS is the one Ehlers-Danlos subtype nobody has a lab test for. You get diagnosed by
checklist: how far your joints bend (the Beighton score), which of a long list of
other symptoms you have, whether it runs in your family, and whether a doctor has
ruled out everything else it could be instead. That's it. That's the 2017 diagnostic
standard, and it's still the standard. Real genetic leads have started showing up in
just the last couple of years (see [Prior work](#prior-work) below), but nothing has
turned into an actual test yet.

Given that, I didn't want to build the obvious thing — a classifier trained on symptom
checklists, which just re-learns the diagnostic criteria and calls it a "biomarker."
Instead I went after something nobody seems to have tried for this disease: literature-
based discovery. It's an old idea (Don Swanson used it in the 1980s to link fish oil
and Raynaud's disease before anyone had run the trial). The logic: if hEDS is linked to
some phenotype B in the literature, and B is linked to gene C in a totally different
set of papers, but nobody has ever written a paper connecting hEDS to C directly — C is
worth a look. All from free, public data, no lab access required.

## How it works

1. **Fetch** (`fetch_literature.py`) — pulls hEDS/HSD-specific papers from PubMed.
   I deliberately left out bare "Ehlers-Danlos syndrome" and other-subtype terms, so
   genes belonging to classical or vascular EDS don't sneak into the corpus. Each
   abstract's genes, diseases, and chemicals get tagged by NCBI's **PubTator3**, so I
   didn't have to write my own entity recognizer.
2. **Graph** (`build_graph.py`) — turns all that into a co-occurrence graph. Two
   concepts get an edge if they show up in the same abstract, weighted by how often,
   plus a PMI score for how surprising that co-occurrence is.
3. **Discover** (`discover_candidates.py`) — the actual ABC-model step. Anything
   directly co-mentioned with hEDS is filed as "already studied" and set aside. What's
   left — genes reachable only through a shared bridge concept, never directly — gets
   ranked by Adamic-Adar score, which is just common-neighbors math that discounts
   generic bridge terms (like "pain," which connects to everything) more than specific
   ones. That distinction matters a lot; raw co-occurrence counting alone is a known
   source of false positives in this kind of work.
4. **Validate** (`validate_candidates.py`) — before anything makes the final list, it
   gets checked against three independent public sources: STRING (does it physically
   interact with a known hEDS-relevant protein?), GTEx (is it actually expressed in
   connective tissue?), and ClinVar (does it have any documented variant record at
   all?).
5. Out comes `results/candidates.csv` — usually 10-20 genes, each with its evidence
   and a plain-English reason it's on the list.

Every data source here is free and public. No patient data touches this pipeline at
any point.

## Prior work

I didn't want to build this in a vacuum, so before writing any code I spent a while
figuring out what's already been tried. Turns out 2024-2025 was a genuinely big couple
of years for hEDS genetics — the reference gene sets baked into `validate_candidates.py`
come from that work, not from the other EDS subtypes' better-known genes:

- **KLK gene family** (especially *KLK15*) — a recurrent variant turned up in
  whole-exome sequencing of 200 hEDS patients, and it reproduces connective-tissue
  defects when reintroduced in mice.
  [Norris Lab, 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11213194/)
- **ACKR3 and SLC39A13** — the first hEDS GWAS meta-analysis (roughly 1,800 cases
  against 5,000 controls) found genome-wide-significant signal here, pointing toward a
  neuroimmune/stromal story rather than a single collagen gene.
  [Petrucci-Nelson et al., medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.09.19.25336146v1)
- **MIA3** — another 2024-2025 candidate, still unresolved.
- **TNXB** (partial deficiency) — the oldest lead in this space, but it only explains
  about 1% of cases, and serum tenascin-X testing flopped as a screening tool when
  people actually tried it. I kept it in as a weak legacy check, not a strong signal.
- **A plasma ECM-fragmentation signature** (fibronectin/collagen-I/tenascin fragments)
  — probably the closest thing to an actual proposed biomarker that exists right now,
  still unvalidated. [Ritelli et al.](https://pubmed.ncbi.nlm.nih.gov/39225014/)
- **The HEDGE Study** — the Ehlers-Danlos Society's own large-scale sequencing effort
  (1,000 patients), still producing results.

As far as I could find — and I looked, specifically, more than once — nobody has
pointed literature-based discovery or knowledge-graph mining at hEDS before. Not even
at neighboring diseases with a similar "diagnosis of exclusion" problem, like long
COVID or fibromyalgia. So this angle, at least, seems to be new.

### What I decided not to do, and why

- **Clustering patients into subtypes by symptoms.** Already done, and done well —
  Mayo Clinic ran K-means and UMAP on 2,100+ real patients and published real subtypes.
  [Petrucci et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38779137/) I don't have that
  data, and re-running the same analysis on something worse wouldn't add anything.
- **Wearables and video-based biomarkers.** This is genuinely exciting research —
  AI-scored hypermobility from video, wearable HRV tracking — but every dataset behind
  it is "available on request," not public. Which, honestly, just confirmed that
  literature mining was the right call for a project built entirely on public data.
- **Swapping in SemMedDB** instead of raw co-occurrence. SemMedDB gives you semantic
  relationships instead of "these two things appeared near each other," which sounds
  strictly better. But it was deprecated in December 2024 and hasn't been updated
  since — it would miss both the 2025 GWAS paper and the 2024 KLK15 finding above —
  and it needs a UMLS license to even access. Not worth it for something meant to stay
  current and easy to run.

## Running it

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt

python pipeline.py --top-n 20
```

Most of the runtime is just being polite to free APIs (~1 request/sec), not actual
computation — figure 10-15 minutes for the full ~1,300-abstract corpus. Everything gets
cached under `cache/`, so re-runs are quick unless you pass `--skip-cache`.

Other flags: `--retmax` (how many PubMed results to pull), `--min-weight` (how many
co-occurring abstracts before an edge counts), `--max-bridge-degree` (the hub-node
cutoff), `--top-n` (shortlist size).

## The actual test: would this have gotten ahead of the field?

A shortlist is a claim about the future — "these are worth a second look." The only
honest way to test that claim is to go back in time. `backtest.py` reruns the exact
same pipeline using only literature published before a chosen cutoff year (the cached
PubTator3 documents each carry their own publication date, so this needs no new
network calls), and checks whether genes that real 2024-2025 hEDS genetics later
confirmed were already showing up as candidates — before that confirmation existed.

```
python backtest.py --cutoffs 2019 2020 2021 2022 2023 2024 2025 2026
```

**The result**: using only papers published through 2021 — more than three and a half
years before the first hEDS GWAS meta-analysis went up on medRxiv (19 Sept 2025) —
the pipeline flags **SLC39A13** as a candidate, through a purely indirect literature
path. In September 2025, an independent GWAS meta-analysis (~1,800 cases, ~5,000
controls, a completely different method — genotyping, not literature) found
genome-wide-significant signal at that exact gene.

I don't just want to report that number — I want to say how I got to trust it, because
I didn't at first. This result went through two rounds of me trying to break it:

**Round 1.** An unrelated exploratory run (loosening the co-occurrence threshold to see
what a noisier pass surfaces) turned up KLK15 — another 2024 hEDS gene — as a weak
signal. Tracing it back led to an actual paper, "KLK15 alters connective tissues in
hypermobile Ehlers-Danlos syndrome" (iScience, Aug 2025), which explicitly studies
KLK15 in hEDS by name. But my pipeline's seed detection had been matching against the
graph's single stored display name per entity ID — and PubTator had normalized this
paper's very specific "hypermobile Ehlers-Danlos syndrome" mention down to the generic
Ehlers-Danlos syndrome ID shared by every subtype, so the match never fired. Fix: seed
detection now scans each document's own raw annotation text directly, not an
aggregated per-ID name that happens to be whichever text was recorded first.

**Round 2.** Re-running the SLC39A13 backtest after that fix, it initially showed as
"already directly studied" by 2022 — which would have meant no prediction case at all.
Digging into why turned up two bugs in my own seed-pattern list: my "type 3" pattern
(meant to catch the pre-2017 name for hEDS) was a naive substring match that also fired
on "Ehlers-Danlos syndrome, **spondylodysplastic** form type 3" — a completely
different subtype, the one SLC39A13 is actually already known to cause. And separately,
I'd included "generalized joint hypermobility" as a seed term, but it's a clinical sign,
not a diagnosis — I found a real paper studying it as a population explicitly distinct
from EDS. Fixed both; added regression tests for both.

Only after those fixes did I trust the result enough to publish it. I checked it
directly: no document before 2022 co-mentions an hEDS/HSD seed term and SLC39A13 in the
same abstract, and neither of the two papers that discuss SLC39A13 at all (2020, 2021 —
both about the unrelated rare subtype it was already known to cause) is itself a seed
document. The link is built from several independently weak literature bridges, most
substantially a shared "juvenile connective tissue diseases" concept appearing in 54
real hEDS/HSD papers — I hand-verified the score's arithmetic against the raw graph,
not just the code path. Full mechanism and sources are in `results/backtest.txt`.

**One more honest thing, checked just now**: SLC39A13 is *still* classified as a
candidate today, not yet "already studied" — even eleven months after the GWAS.
Exactly one document in the whole corpus directly co-mentions hEDS and SLC39A13: the
GWAS preprint itself. No second paper has echoed the connection yet, so it hasn't
cleared this method's own 2-paper bar for "established." That's not a weakness in the
method — it's a real, separate observation: the literature index can lag a genuine
discovery by the better part of a year even after publication, which is exactly the
kind of gap a tool like this is positioned to notice.

Worth being precise about the boundary here too. SLC39A13 wasn't a total unknown going
in — it's the established cause of a different, rare EDS subtype, which is exactly why
it had *some* indirect literature trail to follow. ACKR3, the GWAS's other major locus,
has no connective-tissue literature footprint at all, at any cutoff I tested including
today, and the pipeline correctly has nothing to say about it. KLK15 has exactly one
paper with a matched disease mention in the entire corpus — one paper isn't enough
independent evidence, and correctly doesn't count. This method gets ahead of the field
when a candidate has *some* textual trail, however thin; it has nothing to offer for a
genuinely de novo finding with no trace at all, and it's honest about the difference.

## What it actually found (the current shortlist)

On the corpus I had when writing this (1,301 abstracts, 1,345 entities, 20,238
co-occurrence edges), the pipeline correctly filed 24 genes as "already studied" —
including COL3A1, COL6A3, SMAD3, and TNXB — and surfaced 10 candidates. Topping the
list is **PLOD1**, which defines kyphoscoliotic EDS; most of the rest (AEBP1, DSE,
B4GALT7, FKBP14) are similarly genes that define *other* EDS-spectrum disorders — a
reasonable "maybe this whole gene family matters here too" hypothesis, not individually
shocking. **SLC39A13** sits in the middle of the list, exactly where the backtest above
says it should be. Further down, **C1R (complement C1r)** independently lines up with a
2025/2026 proteomics paper that found complement-cascade proteins in hEDS patients'
blood — a second, smaller case of two unrelated methods landing on the same idea.

ACKR3 and MIA3 don't show up in the corpus at all yet, most likely because PubTator3
hasn't finished indexing their source papers or the field hasn't cited them enough yet
for co-occurrence evidence to build up. That's a data-freshness limit of the underlying
literature index, not a bug in the pipeline, and it should resolve as the field
catches up and this gets re-run down the line.

## Does it actually work?

- `tests/test_build_graph.py` is a small synthetic-graph test that checks the core
  logic actually separates "directly studied" from "implied but never studied" —
  no network required. Run it with `venv\Scripts\python tests\test_build_graph.py`.
- `tests/test_backtest.py` checks the date-filtering logic behind the retrospective
  validation above — the one piece of that script not already covered by the core
  discovery test.
- On real output, the sanity check is calibrated to what's actually established in the
  literature, not just what's true biologically: TNXB, COL3A1, and COL6A3 should land
  in "already studied" (they're repeatedly, directly co-mentioned with hEDS). SLC39A13
  is expected to stay a *candidate* until a second independent paper echoes the GWAS
  finding directly — see "The actual test" above for why that's the correct behavior,
  not a bug.
- And the one I can't automate away: read the top candidates yourself. Look at the
  bridge concept connecting them to hEDS and ask whether it's real biology or an
  annotation glitch. This is a hypothesis-generator, not an oracle — that check is
  what makes the output mean anything.

## The honest disclaimer

This is a hypothesis-generating research tool, not a diagnostic instrument. There is no
confirmed diagnostic biomarker for hEDS in the published literature — not from this
project, not from anyone, not yet. Nothing here should be used to include, exclude, or
otherwise influence an actual diagnosis. At best, a candidate on this list is a
starting point for someone with a lab to look closer.
