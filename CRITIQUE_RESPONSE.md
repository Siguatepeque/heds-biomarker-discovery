# Response to the supplied evidence briefing

The supplied PDF is titled *hEDS biomarker discovery by literature mining: why
the concept does not hold*, dated 8 September 2026. This response addresses its
eight arguments against the repository. It is not a claim that software changes
validate a biomarker, and it is not an independent verification of every numerical
detail or review citation in the PDF.

**The central methodological criticism is justified:** abstract co-occurrence
alone cannot establish novelty, causal relevance, or diagnostic performance.
The defensible output here is an auditable map of literature connections. Whether
it is useful for selecting experiments remains unproved.

## 1. A complex trait is not a missing-single-gene retrieval problem

Agreed about the mismatch between the available evidence and a single causal-gene
or diagnostic claim. The README and discovery docstring now state a narrower
purpose. Reference genes are research anchors, not confirmed hEDS causes or a
complete set of positive labels. The audit treats them as retrospective retrieval
targets only.

The PDF's stronger inference does not follow: an unknown or complex genetic basis
does **not** rule out measurable biomarkers, a multianalyte panel, modifiers, or
genetically distinct subsets. Nor is “no single gene” true by definition for all
future explanations of hEDS. GeneReviews describes an unknown etiology and likely
genetic heterogeneity, alongside familial inheritance. A susceptibility locus,
an effector gene, a circulating protein, and a diagnostic test are different
objects; this pipeline does not establish the links between them.

**Still required:** specify the intended biological hypothesis and clinical use
case before claiming discovery. No single-gene or diagnostic endpoint is measured
by the present code.

## 2. Coverage and diagnostic terminology undermine novelty

Agreed. A gene without a detected annotated co-mention is not necessarily absent
from hEDS literature. The fetch is disease-query scoped, not an independent search
for every B-C relation. An abstract count alone neither proves adequacy nor gives
a universal minimum sample size for discovery.

The offline audit now records exact document IDs, date gaps, input-file hashes,
effective-corpus hashes and active seed patterns. It rejects conflicting versions
of a document instead of silently choosing one. The local full-cache run has
1,301 documents and 891 seed documents; four input documents lack usable dates.
The old single-mention and alias corrections remain in place.

**Not fixed:** missing abstracts, full-text coverage, omitted annotations, species
normalization, or equivalence of hEDS, HSD and pre-2017 JHS/EDS-HT diagnoses.
The cache is a directory union; its filenames do not prove historical query
membership. The audit explicitly warns against mixing disease caches. Its
fingerprints are provenance, not a substitute for the uncommitted input cache.

## 3. A bridge score cannot recover biological context

Agreed. The score is unchanged; no claim is made that degree filtering, PMI, or
Adamic-Adar separates mechanism from comorbidity, or measurement from speculation.
The precision percentages quoted from other ABC evaluations are not measured
performance of this repository and are not imported as its expected accuracy.

`evidence.json` now exposes **both** legs of every scored bridge, the source
document IDs and original passages/annotations, bridge degrees and individual
score contributions. It distinguishes admission edges from weaker edges counted
after admission. Tests verify that contributions sum to the actual discovery
score. Every entry is marked `unreviewed`.

**Not fixed:** context extraction and expert adjudication. Before using a lead,
review both legs for diagnosis/subtype, organism, assay, study design, negation,
direction of effect and whether the claimed relation was actually investigated.
Exporting passages makes that work possible; it does not mean it has been done.

## 4. The shortlist mostly recovers other EDS subtypes

Agreed with that description of the current output. A frozen, cited
`eds_genes.py` mapping now flags genes in GeneReviews' EDS differential table by
human Entrez ID. This includes the table's brittle cornea syndrome entries.
PLOD1, AEBP1 and FKBP14 are flagged in the fresh output. Unflagged is not equivalent
to hEDS-specific, and a subtype flag is not grounds for silently deleting a gene:
other-disease involvement does not logically exclude a modifier role in hEDS.

The remaining `Gene:597` has the PubTator display text `alpha1`; the archived
Entrez lookup resolved it to BCL2A1. The offline audit preserves the identifier
and source display text rather than pretending to have re-resolved the symbol.
Its expression annotation cannot establish either relevance or irrelevance.

**Still required:** evidence beyond family membership for any proposed hEDS role.
No candidate is promoted to a biomarker or a confirmed false positive here.

## 5. Retrospective retrieval needs baselines and a future test

Agreed. The new audit compares ABC against document frequency over all graph
genes not directly mentioned with the seed, and against that same ordering
restricted to the frozen differential set. The baselines are not restricted to
ABC survivors. The targets and differential set are present-day, post hoc choices,
not historical snapshots.

The actual pre-2022 run uses 733 documents:

| Method | SLC39A13 rank | Target groups retrieved at k=3 |
|---|---:|---:|
| ABC / Adamic-Adar | 3 | 1 / 4 |
| Document frequency | 4 | 0 / 4 |
| Known-EDS document frequency | 3 | 1 / 4 |

The denominator groups KLK15's alternate records once and includes absent targets.
These are association-target retrieval counts, not recall of established causal
genes. The comparison does **not** show improvement over known-EDS shortlisting.
The cutoff and k are illustrative, not preregistered. In the full-cache audit,
all methods retrieve 0/4 target groups, with zero eligible groups: KLK15, ACKR3
and SLC39A13 are directly mentioned, while MIA3 is absent.

**Still required:** freeze the corpus, rules, target endpoint and evaluation plan
before future evidence is available; compare at a fixed review budget across
multiple held-out tasks. A new retrospective report cannot manufacture a
prospective prediction.

## 6. Fibromyalgia is not a specificity test

Agreed. The archived 0/8 result and oxidative-stress theme are evidence of a
retrieval limitation and possible theme bias, not proof of specificity. The
targets are not a comprehensive set of causal genes, and non-target candidates
are not automatically demonstrated false positives. Neither sensitivity nor
specificity for a diagnostic use case can be inferred from these lists.

**Not rerun:** the corrected fibromyalgia comparison. The local cache audited
here reproduces the hEDS counts; it is not a frozen control corpus. The control
script still searches PubMed. The archived control output is explicitly labeled
as such, not passed off as a new offline result. Both live pipeline entry points
use annotation-only ranking and no longer impose a hidden 50-candidate long-list
cap before applying the requested shortlist size.

**Still required:** a separately frozen control corpus with usable historical
coverage, identical scoring rules, baseline comparisons, and reviewed relevance
labels if a precision or specificity estimate is desired.

## 7. Annotation must not masquerade as independent validation

Agreed. The arbitrary composite is removed from current code and new output.
STRING, GTEx and ClinVar cannot reorder a candidate. Tests cover a high-AA gene
with unknown annotations outranking a lower-AA gene with all positive annotations,
deterministic ties, and retention of unresolved symbols and subtype flags.

STRING remains functional association context that can include text mining;
changing its label does not make it independent. GTEx remains normal proxy-tissue
expression. ClinVar remains a keyword lookup among at most five records, not
variant adjudication. Unknown responses stay distinct from successful negatives.
Live rows now label their lookup run time, STRING version, GTEx dataset, and
`live_unfrozen` status.

The offline audit produces matched graph/CSV/evidence/report artifacts and records
their hashes, source-file hashes and runtime versions. Existing annotated and
control snapshots are preserved and clearly separated from current graph-only
results. Historical composite values are not silently rewritten as new findings.

**Not fixed:** reproducible live API response replay, channel-specific independent
STRING evidence, disease-specific expression analysis, or ClinVar variant review.
A timestamp is not a frozen response archive. Removing annotation points prevents
one circular use of the data; it does not remove literature bias from ABC itself.

## 8. Direct measurement is essential, but null results have a scope

Agreed that no abstract shortlist can substitute for a biomarker assay. Griggs
reports complement differences and failure to reproduce the proposed disease-
specific ECM-fragment pattern. Its western-blot methods explicitly describe plasma
from the same individuals and blood draws as the earlier Ritelli work: that
comparison should not be called independent-cohort replication. An expanded ELISA
cohort likewise should not automatically be counted as independent replication.

Cinquina reports 54 hEDS-control and 49 HSD-control DEPs, no significant hEDS-HSD
DEPs, and 69 for the combined patient comparison. No significant difference is
not proof of equivalence, and a shared hEDS/HSD signal would not by itself fail
an explicitly combined-population research endpoint. It would not establish
hEDS-specific diagnosis either.

The PDF's “site effect unaddressed” wording needs correction. The paper reports
plate randomization balanced across groups and countries and no DEPs in Italy-US
patient comparisons. Those checks do not eliminate residual site/handling
confounding: all controls were Italian, and sample processing timelines differed.
Neither paper demonstrates clinical performance for this repository's shortlist.

**Not performed:** proteomics reanalysis. A measurement study needs a prespecified
endpoint (hEDS vs HSD is different from combined patients vs controls), assay
coverage checks, within-Italy sensitivity analyses, participant-overlap and
relatedness checks, confounder adjustment and multiple-testing control. Predictive
work needs patient/cohort-separated validation with feature selection and tuning
inside training folds, followed by genuinely independent replication. Do not
reuse literature or measurements used to select candidates as independent validation.

## Sources checked and evidence left available

* Repository code, archived outputs, and newly executed offline full-cache and
  pre-2022 audits: `results/audit/`, `results/audit_pre2022/`.
* [GeneReviews: Hypermobile Ehlers-Danlos Syndrome](https://www.ncbi.nlm.nih.gov/books/NBK1279/),
  chapter updated 22 February 2024: diagnostic criteria, unknown etiology and
  differential-diagnosis Table 1. The 20 curated Entrez ID/symbol pairs were also
  checked against NCBI Gene esummary; no mismatches were returned.
* [Griggs et al., ImmunoHorizons 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12448790/),
  especially cohort recruitment, western analyses, results and replication limits.
* [Cinquina et al., Clinical Proteomics 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13081554/),
  especially participant selection, storage/handling, statistical analysis and
  Italy-US comparisons.
* Genetics, earlier serum studies and the archived fibromyalgia comparison retain
  their citations in the README. This response does not independently endorse
  every effect-size figure, percentage or review-specific precision estimate in
  the supplied PDF.

The practical outcome is an inspectable literature experiment with a baseline
that matches its highlighted retrospective hit. Its added value for experimental
prioritization, and any biomarker claim, remain to be demonstrated.
