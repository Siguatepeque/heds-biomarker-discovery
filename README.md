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
python backtest.py --cutoffs 2020 2022 2023 2024 2025
```

**The result**: using only papers published through 2024 — strictly before the first
hEDS GWAS meta-analysis went up on medRxiv (19 Sept 2025) — the pipeline flags
**SLC39A13** as a candidate. No paper in that pre-2025 corpus ever mentions hEDS/HSD and
SLC39A13 in the same abstract; I checked directly. The connection is entirely indirect,
through a shared "juvenile connective tissue diseases" bridge concept whose link to the
hEDS/HSD seed term crossed the co-occurrence threshold specifically in papers dated
2024. Then, in September 2025, an independent GWAS meta-analysis (~1,800 cases, ~5,000
controls, a completely different method — genotyping, not literature) found
genome-wide-significant signal at that exact gene.

Worth being precise about what this does and doesn't show. SLC39A13 wasn't a total
unknown going in — it's the established cause of a separate, rare, recessive EDS
subtype (spondylocheirodysplastic EDS), and that's exactly why it had enough of an
indirect literature trail for the method to reach it. ACKR3, the GWAS's other major
locus, has essentially no connective-tissue literature footprint before 2025 — it's
absent from the corpus at every cutoff I tested — and the pipeline correctly has
nothing to say about it. That's the honest boundary here: this method can get ahead of
the field when a candidate has *any* indirect textual trail to follow, and it has
nothing to offer for a genuinely de novo finding with no prior trace at all. SLC39A13
is the case where the trail existed, and following it landed on the right answer
months early. Full mechanism, exact edge weights, and sources are in
`results/backtest.txt`.

## What it actually found (the full shortlist)

On the corpus I had when writing this (1,301 abstracts, 1,345 entities, 20,238
co-occurrence edges — i.e. everything through 2026, not the pre-2025 backtest slice
above), the pipeline found 2 clean hEDS/HSD seed terms and correctly filed 18 genes as
"already studied" — including **SLC39A13** and **TNXB**. Notice that's a change from
the backtest: SLC39A13 has moved from "indirect candidate" to "directly studied" as
2025's papers actually caught up to it. That's not a contradiction, it's the point —
the field closed the gap this pipeline flagged.

12 candidates made the final shortlist, topped by **COL3A1** — the gene that defines
vascular EDS. That's not a surprise, biologically; it's a close nomenclature-and-
literature neighbor of hEDS, so it makes sense it'd surface this way. Most of the top
of the list is like that: COL6A3, PLOD1, B4GALT7, DSE, FKBP14, AEBP1 all define *other*
EDS-spectrum disorders — a reasonable "maybe this whole gene family matters here too"
hypothesis, not individually shocking.

The one that actually stopped me, further down the list at #11, is **C1R (complement
C1r)**. I didn't point the pipeline at it — it came out of the graph structure on its
own — and it happens to line up with a 2025/2026 proteomics paper that found
complement-cascade proteins showing up differently in hEDS patients' blood. It's not
the strongest score in the shortlist, but it's the one place where this cheap
literature-mining method and an unrelated wet-lab study landed on the same idea
independently. That's worth more to me than a higher rank would be.

Worth being upfront about: the newest 2024-2025 findings aren't well represented yet.
ACKR3 shows up in the graph but hasn't cleared the co-occurrence threshold to count as
"studied" — probably because only a handful of very recent papers mention it in an hEDS
context so far. KLK15 and MIA3 don't show up in the corpus at all, most likely because
PubTator3 hasn't finished indexing their source papers. That's a data-freshness problem
with the underlying literature index, not a bug in the pipeline, and it should sort
itself out as PubTator3 catches up and this gets re-run down the line.

## Does it actually work?

- `tests/test_build_graph.py` is a small synthetic-graph test that checks the core
  logic actually separates "directly studied" from "implied but never studied" —
  no network required. Run it with `venv\Scripts\python tests\test_build_graph.py`.
- `tests/test_backtest.py` checks the date-filtering logic behind the retrospective
  validation above — the one piece of that script not already covered by the core
  discovery test.
- On real output, the sanity check is: known hEDS-relevant genes (KLK15, ACKR3,
  SLC39A13, TNXB) should land in "already studied," not the shortlist. If one of them
  ever shows up as a "candidate," that's either a corpus-timing issue or something
  worth actually looking into.
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
