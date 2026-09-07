# hEDS Biomarker Discovery

Pipeline that mines hEDS-related literature for genes linked to hEDS/HSD only
through shared bridge concepts, ranks them, and checks each against STRING,
GTEx, and ClinVar. Output is a ranked shortlist for follow-up study. This is
hypothesis generation, not a diagnostic test. No candidate listed here is a
validated biomarker.

Read the [research notes](https://siguatepeque.github.io/heds-biomarker-discovery/)
for the evidence review. This file documents the code.

## Archive status

`results/candidates.csv` was refreshed on 7 September 2026 with the corrected
code and live STRING/GTEx/ClinVar annotation: 1,301 unique documents, 1,345
entities, 20,238 edges, 132 directly mentioned genes excluded, 4 indirect
candidates. SLC39A13, C1R, and VWF are excluded because they have direct
mentions. Annotation values are live lookups from that date, not fixed constants;
a rerun can return different STRING/GTEx/ClinVar values.

`results/graph.graphml`, `results/control_candidates.csv`, and
`results/control_graph.graphml` remain snapshots from before the September 2026
code corrections: left as committed, not refreshed. `results/backtest.txt` and
`results/control_backtest.txt` are legacy reports. Their cutoff blocks and numeric
tables are preserved as run; the surrounding interpretation has been revised.

Where behavior differs from the original method, it is described in
[Method](#method).

## Method

1. Fetch (`fetch_literature.py`). Pulls papers matching an hEDS/HSD query from
   PubMed: hEDS-related search results only, not PubMed as a whole, with no
   separate search for papers connecting bridge concepts to genes. PubTator3
   supplies the entity annotations. Graph and discovery steps count each PMID
   once. Missing annotations leave gaps in the graph, not evidence against a gene.
2. Graph (`build_graph.py`). Builds a co-occurrence graph from abstracts.
   Edges carry a co-occurrence weight and a stored PMI value. PMI is stored
   only. Ranking does not use it.
3. Discover (`discover_candidates.py`). Applies the ABC pattern: genes that
   share a bridge concept with hEDS/HSD but are never directly co-mentioned
   with it in the same abstract become candidates. Matching uses each
   annotation's own surface text against the alias list in that file.
4. Score. Candidates are ranked by Adamic-Adar over bridge nodes: a ranking
   sum, not a probability.
5. Annotate (`validate_candidates.py`). Checks each candidate against STRING
   functional associations (can include text mining, not guaranteed physical
   independent interactions), GTEx normal-tissue expression in fibroblast,
   skin, and aorta proxies (not a case-control result), and ClinVar keyword
   match over up to 5 records (a weak signal). Adds a composite score. Failed or
   missing lookups remain unknown and add no points; an unresolved gene symbol
   no longer drops the candidate. Scores depend on which sources are available.

Definitions used here:

* `min_weight` (default 2) is bridge evidence only. A bridge concept needs
  co-occurrence in at least that many seed documents, and a gene-bridge edge
  needs at least that weight. It does not set a novelty threshold.
* Direct-mention exclusion is independent of `min_weight`. Any document that co-mentions
  hEDS/HSD seed text and a gene in the same abstract files that gene as
  directly mentioned. A co-mention does not establish that a paper investigated
  the gene in hEDS. The archived snapshot required at least two such documents,
  which left some previously studied genes on its candidate list.
* Once an eligible bridge admits a candidate, the Adamic-Adar score also counts
  single-document edges to other eligible bridge nodes. That scoring rule is
  unchanged; the two-document floor is an admission rule, not a floor for every
  score contribution.

## Current output

Full corpus: 1,301 abstracts, 1,345 entities, 20,238 co-occurrence edges.
132 genes directly mentioned alongside hEDS/HSD, 4 indirect candidates.
Scores below are composite scores from the refreshed `results/candidates.csv`
(7 September 2026, live annotation).

| Gene | Composite | Note |
|---|---|---|
| PLOD1 | 7.96 | Defines kyphoscoliotic EDS |
| AEBP1 | 6.29 | Defines classical-like EDS |
| FKBP14 | 3.81 | Defines kyphoscoliotic EDS type 2 |
| BCL2A1 | 2.78 | Low expression, little support |

The original 10-candidate snapshot (with DSE, B4GALT7, TGFB1, SLC39A13, VWF,
C1R) is preserved in git history. Those six left the list under the corrected
direct-mention rule. What remains is a gene-family overlap pattern, not a
per-gene finding. Context for the three excluded entries of interest:

* SLC39A13. The 2025 GWAS meta-analysis (Petrucci-Nelson et al., medRxiv
  2025.09.19.25336146, 1,815 cases and 5,008 controls) reports a
  multi-candidate locus including SLC39A13 and PSMC3: an association, not
  established causality. SLC39A13 was not new to connective-tissue literature:
  it defines a rare recessive EDS subtype and appears near generalized joint
  hypermobility background. See `results/backtest.txt`.
* VWF. The bridge runs through bleeding/bruising terms, a documented hEDS
  comorbidity. The corpus already holds a 2013 bleeding-clinic paper on
  symptomatic hypermobility (PMID 23030528) and a 2019 review using the phrase
  hypermobile-EDS (PMID 31329366) that the old matcher missed.
* C1R. A 2016 serum study (Watanabe et al., PMID 26709396) directly assayed
  C1R, with western follow-up in JHS/EDS-HT under pre-2017 diagnoses, not the
  2017 hEDS criteria. That paper is already corpus input.

BCL2A1 has little supporting evidence here. Its 3.5 TPM annotation is not enough
to decide biological relevance. ACKR3 appears in the cached GWAS abstract but
was not an indirect candidate; MIA3 was absent from the cached corpus.
AEBP1's ClinVar record matched on the September 2026 lookup; that value can
change between runs because it is a live keyword search, not a curated label.

## Backtest (archived, original method)

`backtest.py` reruns discovery on papers dated before each cutoff year. Full
cutoffs and bug history are preserved in `results/backtest.txt`.

Archived headline numbers under the original method, papers before 2022: 733
documents, 8 candidates, SLC39A13 at rank 8 with Adamic-Adar 0.860365.
The corrected offline run uses the same 733 dated documents and retrieves
SLC39A13 at rank 3 of 3, still scoring 0.860365. It has no matched direct mention
in that subset. The smaller list reflects stricter exclusions, not new evidence.

The code and target list were developed after the GWAS was known. This checks
retrospective retrieval, not prospective prediction. The method has not yet
been compared with a shortlist of known EDS genes or publication-frequency baselines.

The GWAS preprint itself (PMID 41001447) is the single direct co-mention of
hEDS/HSD seed text with SLC39A13 in the full corpus. The original method kept
SLC39A13 as a candidate because direct filing then required two documents:
a threshold effect, not an indexing-delay finding. No claim is made of first
use of literature-based discovery in neighboring diseases.

## Fibromyalgia comparison (archived, original method)

`control_disease.py`, `controls.py`, and `control_backtest.py` run the same
pipeline for fibromyalgia with swapped query and seed terms. Details and the
full cutoff table are preserved in `results/control_backtest.txt`.

The 8 fibromyalgia GWAS targets (Kerrebijn et al., medRxiv 2025.09.18.25335914,
published Nat Med 2026) are never candidates at any cutoff. That is a
coverage and recall limit, not proof of specificity. Before 2020 the control
corpus holds 1 document, so early cutoffs test almost nothing. The archived
control run lists 6 candidates (Nfe2l2, Bax, S100A7, PPARGC1A, Myog, IL23A),
none matching the 8 GWAS targets or the older set (GCH1, COMT, OPRM1). The
oxidative-stress theme across some of them shows the failure mode: a large
non-genetic literature can clear the bridge filter and look identical in
structure to a genetic bridge.

## Limits

* Seed matching depends on an alias list and PubTator3 annotations. A paper
  comparing several EDS subtypes can co-mention genes from another subtype.
* STRING, GTEx, and ClinVar supply annotations rather than required pass/fail
  filters. Code tests do not establish biological or diagnostic validity.
* The backtest filters publication dates in today's cache, not a historical
  snapshot of PubMed. Search limits, missing papers, and changed annotations
  can affect which connections are available at each cutoff.

## Prospective sources, not analyzed

These datasets provide a next step beyond literature ranking. Their measurements
have not been reanalyzed in this repository.

* Griggs et al. 2025, 29 hEDS plus 29 matched controls, mass spectrometry,
  [PRIDE PXD062941](https://www.ebi.ac.uk/pride/archive/projects/PXD062941), expanded
  ELISA 41 plus 38 (not automatically independent).
  No replication of the disease-specific ECM-fragment pattern of Ritelli et
  al. doi 10.1002/ajmg.a.63857. doi 10.1093/immhor/vlaf044, PMC12448790.
* Cinquina et al. 2026, 88 hEDS plus 88 HSD plus 176 controls, Olink, public
  NPX supplements. Reported no significant hEDS/HSD difference (not
  equivalence). C1R and VWF not measured. Some US patients, all controls
  Italian. doi 10.1186/s12014-026-09588-2, PMC13081554. Start with a within-Italy
  comparison to assess collection-site effects; check participant overlap before
  calling agreement between studies independent replication.

## Running it

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt

python pipeline.py --top-n 20
```

The scripts pause between API calls. Runtime depends on corpus size and service
availability. PubTator responses are cached under `cache/`; search and candidate
annotation still use the network. Pass `--skip-cache` to refetch PubTator data.

After fetching a corpus, run the graph-only backtest without network calls:

```
python backtest.py --cutoffs 2019 2020 2021 2022 2023 2024 2025 2026
```

Flags: `--retmax` (PubMed result count), `--min-weight` (bridge evidence
floor), `--max-bridge-degree` (hub-node cutoff), `--top-n` (shortlist size).

The September review ran on the existing Python 3.12.6 environment.
`requirements.txt` retains the original pins; that pinned environment was not
reinstalled during this review.

## Tests

Each check runs directly without network access:

```
venv\Scripts\python tests\test_build_graph.py
venv\Scripts\python tests\test_backtest.py
venv\Scripts\python tests\test_controls.py
venv\Scripts\python tests\test_validate_candidates.py
venv\Scripts\python tests\test_site.py
```

* `test_build_graph.py` checks single-mention exclusion, seed aliases, subtype
  exclusions, and duplicate-document invariance.
* `test_backtest.py` checks cutoff date filtering, including that undated
  documents never pass.
* `test_controls.py` checks that hEDS and control seed patterns do not
  cross-match and that the control swap restores hEDS config.

* `test_validate_candidates.py` distinguishes missing API data from negative
  results and checks that unresolved candidates remain in the output.
* `test_site.py` checks page anchors, table headers, and the static HTML structure.

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
5. Kerrebijn I, et al. The genetic architecture of fibromyalgia across 2.5
   million individuals. medRxiv. 2025 Sep 19. Preprint.
   doi:[10.1101/2025.09.18.25335914](https://doi.org/10.1101/2025.09.18.25335914).
   Published: Nat Med. 2026 Jul 28. doi:[10.1038/s41591-026-04492-6](https://doi.org/10.1038/s41591-026-04492-6)
6. Watanabe et al. 2016. Serum study that directly assayed C1R, with western
   follow-up in JHS/EDS-HT under pre-2017 diagnoses. PMID 26709396.
   doi:[10.3892/ijmm.2015.2437](https://doi.org/10.3892/ijmm.2015.2437)
7. Bleeding-clinic study of symptomatic hypermobility. 2013. PMID 23030528.
   doi:[10.1111/hae.12020](https://doi.org/10.1111/hae.12020)
8. Review using the phrase hypermobile-EDS. 2019. PMID 31329366.
   doi:[10.1111/hae.13800](https://doi.org/10.1111/hae.13800)
9. Griggs et al. 2025. 29 hEDS and 29 matched controls, mass spectrometry,
   PRIDE PXD062941, expanded ELISA 41 and 38.
   doi:[10.1093/immhor/vlaf044](https://doi.org/10.1093/immhor/vlaf044).
   PMC12448790.
10. Cinquina et al. 2026. 88 hEDS, 88 HSD, 176 controls, Olink.
    doi:[10.1186/s12014-026-09588-2](https://doi.org/10.1186/s12014-026-09588-2).
    PMC13081554.
