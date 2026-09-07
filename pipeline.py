"""CLI: run the full hEDS literature-based discovery pipeline end to end.

fetch (PubMed + PubTator3) -> graph (co-occurrence) -> discover (ABC-model) ->
validate (STRING/GTEx/ClinVar) -> results/candidates.csv
"""
import argparse
from pathlib import Path

import pandas as pd

from fetch_literature import fetch_corpus
from build_graph import build_graph, save_graph
from discover_candidates import discover, DEFAULT_MIN_WEIGHT, DEFAULT_MAX_BRIDGE_DEGREE
from validate_candidates import validate

RESULTS_DIR = Path(__file__).parent / "results"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=5000)
    parser.add_argument("--min-weight", type=int, default=DEFAULT_MIN_WEIGHT)
    parser.add_argument("--max-bridge-degree", type=int, default=DEFAULT_MAX_BRIDGE_DEGREE)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--skip-cache", action="store_true")
    args = parser.parse_args()

    print("Fetching literature...")
    documents = fetch_corpus(retmax=args.retmax, skip_cache=args.skip_cache)
    print(f"  {len(documents)} annotated documents")

    print("Building co-occurrence graph...")
    graph = build_graph(documents)
    print(f"  {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    RESULTS_DIR.mkdir(exist_ok=True)
    save_graph(graph, RESULTS_DIR / "graph.graphml")

    print("Running ABC-model discovery...")
    long_list, directly_mentioned = discover(
        graph, documents, min_weight=args.min_weight, max_bridge_degree=args.max_bridge_degree
    )
    print(f"  {len(long_list)} candidates, {len(directly_mentioned)} directly-mentioned genes excluded")

    print("Annotating candidates with STRING/GTEx/ClinVar plausibility signals...")
    shortlist = validate(long_list, top_n=args.top_n)
    out_path = RESULTS_DIR / "candidates.csv"
    pd.DataFrame(shortlist).to_csv(out_path, index=False)
    print(f"Wrote {len(shortlist)} candidates to {out_path}")


if __name__ == "__main__":
    main()
