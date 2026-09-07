"""Retrospective retrieval check: rerun discovery using only literature published
before a cutoff date, and check whether retrospective association targets from
2024-2025 hEDS genetics (KLK15, ACKR3, SLC39A13, MIA3) were already retrievable as
indirect candidates - i.e. would this pipeline have surfaced them before those
association reports? These are retrieval/association targets, not confirmed causal
genes, and a pre-cutoff ranking is not a validated prediction of future findings.

Reuses the already-cached PubTator3 documents (no new network calls); each cached
document carries its own publication date.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

from build_graph import build_graph, dedupe_documents, doc_identity
from discover_candidates import discover, seed_documents

# Retrospective retrieval/association targets from 2024-2025 hEDS research (not
# confirmed causal genes). node_id uses PubTator's own Entrez Gene ID, resolved
# once via NCBI esearch (db=gene) - see README for citations.
TARGET_GENES = {
    "Gene:55554": "KLK15",
    "Gene:58197": "KLK15 (alt record)",
    "Gene:57007": "ACKR3",
    "Gene:91252": "SLC39A13",
    "Gene:375056": "MIA3",
}

# Same-family/pathway members: if the exact gene is not retrieved, a relative showing
# up is still a "pointed in the same neighborhood early" retrieval signal, not proof
# of a causal prediction.
FAMILY_PREFIXES = {"kallikrein": "KLK", "atypical chemokine receptor": "ACKR"}


def load_all_documents():
    cache_dir = Path(__file__).parent / "cache"
    documents = []
    for path in sorted(cache_dir.glob("*.json")):
        documents.extend(json.loads(path.read_text(encoding="utf-8")))
    documents = dedupe_documents(documents)
    documents.sort(key=lambda d: doc_identity(d) or "")
    return documents


def doc_year(date):
    """Calendar year of an ISO date/datetime string, or None when missing/malformed.

    Only real calendar dates count: invalid months, days, times, and non-positive
    years are rejected, so malformed dates can never count as early evidence.
    """
    if not isinstance(date, str):
        return None
    text = date.strip()
    for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            parsed = datetime.strptime(text, fmt)
        except ValueError:
            continue
        return parsed.year if parsed.year > 0 else None
    return None


def filter_before(documents, cutoff_year):
    kept = []
    for doc in dedupe_documents(documents):
        year = doc_year(doc.get("date"))
        if year is not None and year < cutoff_year:
            kept.append(doc)
    return kept


def run_cutoff(documents, cutoff_year, min_weight=2, max_bridge_degree=500):
    subset = filter_before(documents, cutoff_year)
    graph = build_graph(subset)
    if not seed_documents(subset):
        return subset, graph, None, None

    long_list, directly_mentioned = discover(
        graph, subset, min_weight=min_weight, max_bridge_degree=max_bridge_degree, top_n=len(graph)
    )
    return subset, graph, long_list, directly_mentioned


def report(cutoff_year, subset, graph, long_list, directly_mentioned):
    print(f"\n=== Cutoff: papers before {cutoff_year} ===")
    print(f"  {len(subset)} documents, {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    if long_list is None:
        print("  No hEDS/HSD seed mentions in this subset - too little corpus yet.")
        return

    total = len(long_list)
    print(f"  {total} indirect candidates (genes never directly mentioned alongside hEDS/HSD in this subset)")
    rank = {row["node_id"]: i + 1 for i, row in long_list.iterrows()}
    candidate_ids = set(long_list["node_id"]) if total else set()
    for node_id, label in TARGET_GENES.items():
        if node_id in directly_mentioned:
            print(f"  {label:24s} -> DIRECTLY MENTIONED by {cutoff_year} (not an indirect-candidate case)")
        elif node_id in candidate_ids:
            score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
            print(f"  {label:24s} -> retrieved as candidate rank {rank[node_id]}/{total} before {cutoff_year} (score {score:.2f})")
        elif node_id in graph.nodes:
            print(f"  {label:24s} -> present in corpus but not a candidate (degree {graph.degree(node_id)})")
        else:
            print(f"  {label:24s} -> not present in corpus before {cutoff_year}")

    # Broader family/pathway check across all Gene nodes in this subset's candidate list
    for node_id in sorted(candidate_ids):
        name = graph.nodes[node_id].get("name", "")
        for family_desc, prefix in FAMILY_PREFIXES.items():
            if name.upper().startswith(prefix) and node_id not in TARGET_GENES:
                score = long_list.loc[long_list["node_id"] == node_id, "adamic_adar_score"].iloc[0]
                print(f"  [family hit] {name} ({family_desc}) -> retrieved as candidate rank {rank[node_id]}/{total} before {cutoff_year} (score {score:.2f})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs", type=int, nargs="+", default=[2020, 2022, 2023, 2024])
    args = parser.parse_args()

    documents = load_all_documents()
    print(f"Loaded {len(documents)} cached documents total")

    for year in args.cutoffs:
        subset, graph, long_list, directly_mentioned = run_cutoff(documents, year)
        report(year, subset, graph, long_list, directly_mentioned)
