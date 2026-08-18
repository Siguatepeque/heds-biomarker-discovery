"""Negative-control backtest: rerun backtest.py's retrospective-validation machinery
on the fibromyalgia control corpus instead of hEDS.

Fibromyalgia has no gene that has reached genome-wide significance and replicated the
way SLC39A13 did for hEDS (see controls.py), so the honest expected outcome is that no
control candidate should show that "flagged years early, later confirmed" pattern.
This script doesn't hardcode a target-gene list to check against, the way backtest.py
checks KLK15/ACKR3/SLC39A13/MIA3, because there is no established fibromyalgia gene to
check against - that's the point of the control. Read results/control_backtest.txt for
the actual established-gene check, done by hand against OMIM/GWAS Catalog.

Reuses backtest.filter_before and backtest.run_cutoff completely unchanged. Only the
corpus (control_disease.fetch_control_corpus) and seed patterns
(controls.CONTROL_SEED_PATTERNS) differ, swapped in the same way control_disease.py
swaps them for the main control pipeline.
"""
import argparse

import discover_candidates
from backtest import run_cutoff
from control_disease import fetch_control_corpus
from controls import CONTROL_EXCLUDE_PATTERNS, CONTROL_SEED_PATTERNS

# Fibromyalgia's actual established genes, as of this writing: the first genome-wide-
# significant/prioritized loci for fibromyalgia, from Kerrebijn et al., "The genetic
# architecture of fibromyalgia across 2.5 million individuals" (N=2.5M, 54,629 cases),
# preprinted on medRxiv 2025-09-19 (doi:10.1101/2025.09.18.25335914) and published in
# Nature Medicine 2026-07-28 (doi:10.1038/s41591-026-04492-6) - one study, not two.
# Before this paper, the GWAS Catalog had zero fibromyalgia associations - this control
# disease was only a clean "no known genetic basis" case at the time this project's
# design was decided, not necessarily by the time the corpus was actually fetched.
# Checking against these real genes, the same way backtest.py checks hEDS candidates
# against SLC39A13 etc., is the honest version of this control - see
# results/control_backtest.txt.
TARGET_GENES = {
    "Gene:3064": "HTT", "Gene:9293": "GPR52", "Gene:1630": "DCC", "Gene:1813": "DRD2",
    "Gene:4684": "NCAM1", "Gene:161357": "MDGA2", "Gene:56853": "CELF4", "Gene:79012": "CAMKV",
}


def report_control(cutoff_year, subset, graph, long_list, already_studied):
    print(f"\n=== Cutoff: papers before {cutoff_year} ===")
    print(f"  {len(subset)} documents, {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    if long_list is None:
        print("  No fibromyalgia seed mentions in this subset - too little corpus yet.")
        return
    print(f"  {len(long_list)} candidates, {len(already_studied)} already-studied genes")
    for _, row in long_list.head(5).iterrows():
        print(f"  top candidate: {row['name']:20s} (score {row['adamic_adar_score']:.2f})")

    candidate_ids = set(long_list["node_id"]) if len(long_list) else set()
    for node_id, label in TARGET_GENES.items():
        if node_id in already_studied:
            print(f"  {label:8s} -> ALREADY DIRECTLY STUDIED by {cutoff_year}")
        elif node_id in candidate_ids:
            score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
            print(f"  {label:8s} -> CANDIDATE before {cutoff_year} (score {score:.2f}) - would have been flagged")
        elif node_id in graph.nodes:
            print(f"  {label:8s} -> present in corpus but not a candidate (degree {graph.degree(node_id)})")
        else:
            print(f"  {label:8s} -> not present in corpus before {cutoff_year}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs", type=int, nargs="+",
                         default=[2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026])
    parser.add_argument("--skip-cache", action="store_true")
    args = parser.parse_args()

    documents = fetch_control_corpus(skip_cache=args.skip_cache)
    print(f"Fetched {len(documents)} control (fibromyalgia) documents total")

    original_seeds = discover_candidates.SEED_PATTERNS
    original_excludes = discover_candidates.EXCLUDE_PATTERNS
    discover_candidates.SEED_PATTERNS = CONTROL_SEED_PATTERNS
    discover_candidates.EXCLUDE_PATTERNS = CONTROL_EXCLUDE_PATTERNS
    try:
        for year in args.cutoffs:
            subset, graph, long_list, already_studied = run_cutoff(documents, year)
            report_control(year, subset, graph, long_list, already_studied)
    finally:
        discover_candidates.SEED_PATTERNS = original_seeds
        discover_candidates.EXCLUDE_PATTERNS = original_excludes
