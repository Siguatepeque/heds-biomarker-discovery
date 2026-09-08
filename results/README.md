# Saved results

* `audit/`: current offline full-cache graph, candidates, bridge provenance and
  fingerprinted report (1,301 documents, 4 candidates). No live annotation or composite.
* `audit_pre2022/`: current offline pre-2022 check (733 documents, k=3).
  SLC39A13 is rank 3 by ABC, 4 by frequency, and 3 by known-EDS frequency.
  See the root README for reproducible commands and cache limitations.
* `candidates.csv`: archived 7 September 2026 live STRING/GTEx/ClinVar annotation
  (4 indirect candidates). Its composite column is historical, no longer emitted
  by current code. This snapshot has not been relabeled as a new annotation run.
* `graph.graphml`, `control_candidates.csv`, `control_graph.graphml`:
  snapshots from before the September 2026 code corrections. Left as
  committed, not refreshed.
* `backtest.txt`, `control_backtest.txt`: legacy reports. Cutoff blocks and
  numeric tables preserved as run; surrounding interpretation revised. The
  root README records the corrected hEDS offline run. The fibromyalgia run
  has not been repeated.
