"""Plausibility-control backtest: rerun backtest.py's retrospective-retrieval machinery
on the fibromyalgia control corpus instead of hEDS.

Control candidate yield is context for how the method behaves on another disease,
not a specificity proof for any hEDS candidate: fibromyalgia now has
genome-wide-significant loci of its own (see below), so presence or absence of a
"retrieved early" pattern here neither validates nor invalidates the hEDS run.

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


def report_control(cutoff_year, subset, graph, long_list, directly_mentioned):
    print(f"\n=== Cutoff: papers before {cutoff_year} ===")
    print(f"  {len(subset)} documents, {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    if long_list is None:
        print("  No fibromyalgia seed mentions in this subset - too little corpus yet.")
        return
    total = len(long_list)
    print(f"  {total} indirect candidates (genes never directly mentioned alongside fibromyalgia in this subset)")
    for _, row in long_list.head(5).iterrows():
        print(f"  top candidate: {row['name']:20s} (score {row['adamic_adar_score']:.2f})")

    rank = {row["node_id"]: i + 1 for i, row in long_list.iterrows()}
    candidate_ids = set(long_list["node_id"]) if total else set()
    for node_id, label in TARGET_GENES.items():
        if node_id in directly_mentioned:
            print(f"  {label:8s} -> DIRECTLY MENTIONED by {cutoff_year} (not an indirect-candidate case)")
        elif node_id in candidate_ids:
            score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
            print(f"  {label:8s} -> retrieved as candidate rank {rank[node_id]}/{total} before {cutoff_year} (score {score:.2f})")
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
            subset, graph, long_list, directly_mentioned = run_cutoff(documents, year)
            report_control(year, subset, graph, long_list, directly_mentioned)
    finally:
        discover_candidates.SEED_PATTERNS = original_seeds
        discover_candidates.EXCLUDE_PATTERNS = original_excludes
