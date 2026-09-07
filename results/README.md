# Saved results

The numeric files retain the original results from before the September 2026
code corrections. They have not been regenerated in this review.

* `candidates.csv`, `graph.graphml`, `control_candidates.csv`,
  `control_graph.graphml`: tracked snapshot from before the code corrections.
  Left as committed, not refreshed, not current outputs. The pipeline writes
  into `results/` on a deliberate rerun, so this is a snapshot, not a
  permanent archive; a rerun after code corrections overwrites these files
  and will differ.
* `backtest.txt`, `control_backtest.txt`: legacy reports. Cutoff blocks and
  numeric tables preserved as run. Surrounding interpretation revised to
  correct their interpretation. The root README records the corrected hEDS
  offline run. The fibromyalgia run has not been repeated.
