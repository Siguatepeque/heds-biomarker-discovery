"""Retrospective validation: rerun discovery using only literature published before a
cutoff date, and check whether genes later confirmed by real 2024-2025 hEDS genetics
(KLK15, ACKR3, SLC39A13, MIA3) already surface as candidates - i.e. would this pipeline
have flagged them before the field did?

Reuses the already-cached PubTator3 documents (no new network calls); each cached
document carries its own publication date.
"""
import argparse
import glob
import json

from build_graph import build_graph
from discover_candidates import discover, seed_documents

# The genes real 2024-2025 hEDS research confirmed. node_id uses PubTator's own
# Entrez Gene ID, resolved once via NCBI esearch (db=gene) - see README for citations.
TARGET_GENES = {
    "Gene:55554": "KLK15",
    "Gene:58197": "KLK15 (alt record)",
    "Gene:57007": "ACKR3",
    "Gene:91252": "SLC39A13",
    "Gene:375056": "MIA3",
}

# Same-family/pathway members: if the exact gene doesn't show up, a relative showing up
# is still a meaningful "pointed in the right direction early" signal.
FAMILY_PREFIXES = {"kallikrein": "KLK", "atypical chemokine receptor": "ACKR"}


def load_all_documents():
    documents = []
    for path in glob.glob("cache/*.json"):
        documents.extend(json.loads(open(path, encoding="utf-8").read()))
    return documents


def filter_before(documents, cutoff_year):
    kept = []
    for doc in documents:
        date = doc.get("date")
        if date and int(date[:4]) < cutoff_year:
            kept.append(doc)
    return kept


def run_cutoff(documents, cutoff_year, min_weight=2, max_bridge_degree=500):
    subset = filter_before(documents, cutoff_year)
    graph = build_graph(subset)
    if not seed_documents(subset):
        return subset, graph, None, None

    long_list, already_studied = discover(
        graph, subset, min_weight=min_weight, max_bridge_degree=max_bridge_degree, top_n=len(graph)
    )
    return subset, graph, long_list, already_studied


def report(cutoff_year, subset, graph, long_list, already_studied):
    print(f"\n=== Cutoff: papers before {cutoff_year} ===")
    print(f"  {len(subset)} documents, {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    if long_list is None:
        print("  No hEDS/HSD seed mentions in this subset - too little corpus yet.")
        return

    candidate_ids = set(long_list["node_id"]) if len(long_list) else set()
    for node_id, label in TARGET_GENES.items():
        if node_id in already_studied:
            print(f"  {label:24s} -> ALREADY DIRECTLY STUDIED by {cutoff_year} (not a prediction case)")
        elif node_id in candidate_ids:
            score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
            print(f"  {label:24s} -> CANDIDATE before {cutoff_year} (score {score:.2f}) - would have been flagged")
        elif node_id in graph.nodes:
            print(f"  {label:24s} -> present in corpus but not a candidate (degree {graph.degree(node_id)})")
        else:
            print(f"  {label:24s} -> not present in corpus before {cutoff_year}")

    # Broader family/pathway check across all Gene nodes in this subset's candidate list
    for node_id in candidate_ids:
        name = graph.nodes[node_id].get("name", "")
        for family_desc, prefix in FAMILY_PREFIXES.items():
            if name.upper().startswith(prefix) and node_id not in TARGET_GENES:
                score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
                print(f"  [family hit] {name} ({family_desc}) -> candidate before {cutoff_year} (score {score:.2f})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs", type=int, nargs="+", default=[2020, 2022, 2023, 2024])
    args = parser.parse_args()

    documents = load_all_documents()
    print(f"Loaded {len(documents)} cached documents total")

    for year in args.cutoffs:
        subset, graph, long_list, already_studied = run_cutoff(documents, year)
        report(year, subset, graph, long_list, already_studied)
