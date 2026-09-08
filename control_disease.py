"""CLI: run the same literature-based discovery pipeline as pipeline.py, but on the
fibromyalgia plausibility control instead of hEDS.

fetch (PubMed + PubTator3, control query) -> graph -> discover (ABC-model, control
seeds) -> annotate (STRING/GTEx/ClinVar plausibility signals, unchanged) ->
results/control_candidates.csv

Reuses fetch_corpus / build_graph / discover / validate completely unchanged. The
only difference from pipeline.py is which query and seed patterns are active while
those functions run - swapped in by temporarily overriding the relevant module-level
constants in fetch_literature and discover_candidates, then restoring them. Those
two modules' own code, and everything on the hEDS path, is untouched. Control
candidate yield is context for the method's behavior, not a specificity proof for
any hEDS candidate.
"""
import argparse
from pathlib import Path

import pandas as pd

import discover_candidates
import fetch_literature
from build_graph import build_graph, save_graph
from controls import CONTROL_EXCLUDE_PATTERNS, CONTROL_QUERY, CONTROL_SEED_PATTERNS
from discover_candidates import DEFAULT_MAX_BRIDGE_DEGREE, DEFAULT_MIN_WEIGHT, discover
from validate_candidates import validate

RESULTS_DIR = Path(__file__).parent / "results"


def fetch_control_corpus(retmax=5000, skip_cache=False):
    """fetch_literature.fetch_corpus(), pointed at CONTROL_QUERY instead of QUERY."""
    original_query = fetch_literature.QUERY
    fetch_literature.QUERY = CONTROL_QUERY
    try:
        return fetch_literature.fetch_corpus(retmax=retmax, skip_cache=skip_cache)
    finally:
        fetch_literature.QUERY = original_query


def discover_control(graph, documents, min_weight=DEFAULT_MIN_WEIGHT,
                      max_bridge_degree=DEFAULT_MAX_BRIDGE_DEGREE, top_n=50):
    """discover_candidates.discover(), pointed at the control seed patterns instead."""
    original_seeds = discover_candidates.SEED_PATTERNS
    original_excludes = discover_candidates.EXCLUDE_PATTERNS
    discover_candidates.SEED_PATTERNS = CONTROL_SEED_PATTERNS
    discover_candidates.EXCLUDE_PATTERNS = CONTROL_EXCLUDE_PATTERNS
    try:
        return discover(graph, documents, min_weight=min_weight,
                         max_bridge_degree=max_bridge_degree, top_n=top_n)
    finally:
        discover_candidates.SEED_PATTERNS = original_seeds
        discover_candidates.EXCLUDE_PATTERNS = original_excludes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=5000)
    parser.add_argument("--min-weight", type=int, default=DEFAULT_MIN_WEIGHT)
    parser.add_argument("--max-bridge-degree", type=int, default=DEFAULT_MAX_BRIDGE_DEGREE)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--skip-cache", action="store_true")
    args = parser.parse_args()

    print("Fetching control (fibromyalgia) literature...")
    documents = fetch_control_corpus(retmax=args.retmax, skip_cache=args.skip_cache)
    print(f"  {len(documents)} annotated documents")

    print("Building co-occurrence graph...")
    graph = build_graph(documents)
    print(f"  {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    RESULTS_DIR.mkdir(exist_ok=True)
    save_graph(graph, RESULTS_DIR / "control_graph.graphml")

    print("Running ABC-model discovery...")
    long_list, directly_mentioned = discover_control(
        graph, documents, min_weight=args.min_weight, max_bridge_degree=args.max_bridge_degree,
        top_n=len(graph)
    )
    print(f"  {len(long_list)} candidates, {len(directly_mentioned)} directly-mentioned genes excluded")

    print("Annotating candidates with STRING/GTEx/ClinVar plausibility signals...")
    shortlist = validate(long_list, top_n=args.top_n)
    out_path = RESULTS_DIR / "control_candidates.csv"
    pd.DataFrame(shortlist).to_csv(out_path, index=False)
    print(f"Wrote {len(shortlist)} candidates to {out_path}")


if __name__ == "__main__":
    main()
