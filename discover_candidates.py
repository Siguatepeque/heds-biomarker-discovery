"""Swanson ABC-model literature-based discovery.

Surfaces genes that are indirectly linked to hEDS/HSD via shared literature context
(a bridge concept both co-occur with) but are never directly co-mentioned with it -
i.e. candidates implied by the literature but not yet directly studied in hEDS.
"""
import pandas as pd
import networkx as nx

# hEDS/HSD seed phenotype terms, matched case-insensitively as substrings of PubTator's
# Disease-entity display names. Substring matching on surface text (rather than a single
# hardcoded MeSH/MONDO ID) is deliberate: PubTator normalizes the same real-world concept
# to different database IDs across papers, and text matching catches all of them cheaply.
# Deliberately specific (mirrors fetch_literature.QUERY): bare "hypermobil" or
# "ehlers-danlos" would also match generic "joint hypermobility" (a huge non-specific
# phenotype hub) and OTHER EDS subtypes (classic, vascular, type VII, ...), which
# defeats the whole point of scoping seeds to hEDS/HSD specifically.
SEED_PATTERNS = [
    "hypermobile ehlers-danlos",
    "ehlers-danlos syndrome hypermobility type",
    "ehlers-danlos syndrome type iii",
    "ehlers-danlos syndrome type 3",
    "eds-ht",
    "joint hypermobility syndrome",
    "hypermobility spectrum disorder",
    "hypermobile spectrum disorder",
    "generalized joint hypermobility",
]
# Legacy/combined MeSH headings that would otherwise slip through the patterns above and
# drag in OTHER subtypes' well-known genes - e.g. the old combined "type III and type IV"
# heading mixes hEDS with vascular EDS (defined by COL3A1), and "Marfanoid joint
# hypermobility syndrome" is a distinct, FBN1-related syndrome, not hEDS.
EXCLUDE_PATTERNS = ["and type iv", "marfanoid"]

DEFAULT_MIN_WEIGHT = 2  # require >=2 co-occurring abstracts - filters one-off mentions
# Hub-node cutoff: literature-based-discovery research documents that generic, high-degree
# bridge concepts (e.g. "pain", "fatigue") inflate false positives beyond what Adamic-Adar's
# implicit degree penalty alone corrects for. Bridge candidates above this structural degree
# are excluded outright rather than just down-weighted.
DEFAULT_MAX_BRIDGE_DEGREE = 500


def seed_nodes(graph):
    seeds = set()
    for n, data in graph.nodes(data=True):
        if data.get("type") != "Disease":
            continue
        name = data.get("name", "").lower().replace(",", "")
        if any(p in name for p in SEED_PATTERNS) and not any(p in name for p in EXCLUDE_PATTERNS):
            seeds.add(n)
    return seeds


def discover(graph, min_weight=DEFAULT_MIN_WEIGHT, max_bridge_degree=DEFAULT_MAX_BRIDGE_DEGREE, top_n=50):
    """Returns (long_list_df, already_studied_genes). long_list_df has columns
    node_id, name, adamic_adar_score, ranked descending."""
    seeds = seed_nodes(graph)
    if not seeds:
        raise ValueError("No hEDS/HSD seed nodes found in graph - check corpus/query scope")

    direct_neighbor_genes = {
        n for s in seeds for n in graph.neighbors(s)
        if graph.nodes[n].get("type") == "Gene" and graph[s][n]["weight"] >= min_weight
    }

    bridge_nodes = {
        n for s in seeds for n in graph.neighbors(s)
        if graph[s][n]["weight"] >= min_weight and graph.degree(n) <= max_bridge_degree
    }

    candidates = set()
    for b in bridge_nodes:
        for c in graph.neighbors(b):
            if (graph.nodes[c].get("type") == "Gene"
                    and c not in direct_neighbor_genes
                    and c not in seeds
                    and graph[b][c]["weight"] >= min_weight):
                candidates.add(c)

    pairs = [(s, c) for s in seeds for c in candidates]
    scores = {}
    for s, c, aa in nx.adamic_adar_index(graph, pairs):
        scores[c] = max(scores.get(c, 0.0), aa)

    rows = [
        {"node_id": c, "name": graph.nodes[c].get("name", c), "adamic_adar_score": scores.get(c, 0.0)}
        for c in candidates
    ]
    rows.sort(key=lambda r: r["adamic_adar_score"], reverse=True)

    long_list = pd.DataFrame(rows[:top_n], columns=["node_id", "name", "adamic_adar_score"])
    return long_list, direct_neighbor_genes
