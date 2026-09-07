"""Swanson ABC-model literature-based discovery.

Surfaces genes that are indirectly linked to hEDS/HSD via shared literature context
(a bridge concept both co-occur with) but are never directly co-mentioned with it -
i.e. candidates implied by the literature but not directly mentioned alongside hEDS.

Seed detection works off each mention's raw annotation text, document by document -
not off the graph's single stored display name per entity ID. PubTator sometimes
normalizes text that literally says "hypermobile Ehlers-Danlos syndrome" down to the
generic Ehlers-Danlos syndrome ID (shared with every other subtype), which would make
it invisible to ID-based matching even though the actual sentence is unambiguous.
Concrete case this fixed: iScience, PMID 40949095, "KLK15 alters connective tissues in
hypermobile Ehlers-Danlos syndrome" - PubTator tags the disease mention with the
generic ID, so graph-node-name matching missed it; text-level matching does not.
"""
import math
import re

import pandas as pd

from build_graph import dedupe_documents, document_entities

# hEDS/HSD seed phenotype terms, matched against a normalized form of each
# annotation's own surface text (lowercased, unicode dashes -> space, punctuation
# [,.;:/(){}] -> space, whitespace collapsed). Normalization is what makes
# "Ehlers-Danlos syndrome, hypermobility type", "Ehlers–Danlos syndrome-
# hypermobility type" and "hypermobile-EDS" match their canonical patterns.
# Deliberately specific (mirrors fetch_literature.QUERY):
# bare "hypermobil" or "ehlers-danlos" would also match generic "joint hypermobility"
# (a huge non-specific phenotype hub) and OTHER EDS subtypes (classic, vascular, type
# VII, ...), which defeats the whole point of scoping seeds to hEDS/HSD specifically.
# "heds" is matched as a whole token only (word boundaries), so arbitrary substrings
# cannot match; "hypermobile eds" covers the "hypermobile-EDS"/"Hypermobile EDS"
# abbreviation spellings.
SEED_PATTERNS = [
    "hypermobile ehlers-danlos",
    "hypermobile eds",
    "heds",
    "ehlers-danlos syndrome hypermobility type",
    "ehlers-danlos syndrome type iii",
    "ehlers-danlos syndrome type 3",
    "ehlers-danlos type iii",
    "ehlers-danlos type 3",
    "eds-ht",
    "joint hypermobility syndrome",
    "hypermobility spectrum disorder",
    "hypermobility spectrum disorders",
    "hypermobile spectrum disorder",
    "hypermobile spectrum disorders",
]
# Deliberately NOT a seed pattern: "generalized joint hypermobility" (GJH). It's part
# of the diagnostic criteria for hEDS, but it's a clinical sign, not a diagnosis - real
# papers study GJH as its own population explicitly distinct from EDS (e.g. PMID
# 31594391: "crucial to distinguish the genetic basis of GJH from... EDS"). Treating it
# as a seed term produced a false "directly mentioned" hit during development; still
# fine to keep in fetch_literature.QUERY, which only decides corpus *scope*, not seed status.
# Legacy/combined MeSH headings and name collisions that would otherwise slip through
# the patterns above and drag in OTHER subtypes' well-known genes - e.g. the old
# combined "type III and type IV" heading mixes hEDS with vascular EDS (defined by
# COL3A1), "Marfanoid joint hypermobility syndrome" is a distinct, FBN1-related
# syndrome, and "spondylodysplastic Ehlers-Danlos syndrome type 3" is a *different*,
# SLC39A13-defined subtype whose name happens to contain "Ehlers-Danlos syndrome ...
# type 3" as a literal substring - naive matching mistook it for hEDS.
EXCLUDE_PATTERNS = ["and type iv", "marfanoid", "spondylodysplastic", "spondylo-dysplastic"]

DEFAULT_MIN_WEIGHT = 2  # require >=2 seed documents / co-occurring abstracts
# Hub-node cutoff: literature-based-discovery research documents that generic, high-degree
# bridge concepts (e.g. "pain", "fatigue") inflate false positives beyond what Adamic-Adar's
# implicit degree penalty alone corrects for. Bridge candidates above this structural degree
# are excluded outright rather than just down-weighted.
DEFAULT_MAX_BRIDGE_DEGREE = 500


def _normalize(text):
    t = text.lower()
    for dash in ("\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2015", "\u2212", "\u00ad"):
        t = t.replace(dash, "-")
    t = t.replace("-", " ")
    for ch in ".,;:/()[]{}":
        t = t.replace(ch, " ")
    return re.sub(r"\s+", " ", t).strip()


def _is_seed_mention(text):
    t = _normalize(text)
    if not t:
        return False
    for p in EXCLUDE_PATTERNS:
        if _normalize(p) and _normalize(p) in t:
            return False
    for p in SEED_PATTERNS:
        np = _normalize(p)
        # Whole-phrase boundary match, so a longer token cannot match through a
        # shared prefix (e.g. "type iiia" must not match via "type iii"). Known
        # plurals are listed explicitly in SEED_PATTERNS instead.
        if np and re.search(r"\b" + re.escape(np) + r"\b", t):
            return True
    return False


def seed_documents(documents):
    """Documents containing at least one Disease annotation whose own surface text
    names hEDS/HSD specifically - checked per mention, not via the graph's aggregated
    per-ID display name (see module docstring for why that distinction matters)."""
    documents = dedupe_documents(documents)
    hits = []
    for doc in documents:
        for passage in doc.get("passages", []):
            for ann in passage.get("annotations", []):
                if ann.get("infons", {}).get("type") != "Disease":
                    continue
                if _is_seed_mention(ann.get("text", "")):
                    hits.append(doc)
                    break
            else:
                continue
            break
    return hits


def _adamic_adar(graph, bridge_nodes, candidates):
    scores = {c: 0.0 for c in candidates}
    for b in sorted(bridge_nodes):
        degree = graph.degree(b)
        if degree <= 1:
            continue  # log(1) = 0 -> undefined contribution; standard to skip
        weight = 1.0 / math.log(degree)
        for c in graph.neighbors(b):
            if c in scores:
                scores[c] += weight
    return scores


def discover(graph, documents, min_weight=DEFAULT_MIN_WEIGHT,
             max_bridge_degree=DEFAULT_MAX_BRIDGE_DEGREE, top_n=50):
    """Returns (long_list_df, directly_mentioned_genes). long_list_df has columns
    node_id, name, adamic_adar_score, ranked descending (ties broken by node_id
    for determinism). A gene co-mentioned in ANY seed document counts as directly
    mentioned - even a single co-mention - independently of the min_weight bridge
    floor, which only governs bridge/candidate edge strength."""
    documents = dedupe_documents(documents)
    seeds = seed_documents(documents)
    if not seeds:
        raise ValueError("No hEDS/HSD seed mentions found in documents - check corpus/query scope")

    # Count, across seed documents only, how many times each co-occurring entity appears.
    counts = {}
    for doc in seeds:
        for node_id in document_entities(doc):
            counts[node_id] = counts.get(node_id, 0) + 1

    directly_mentioned_genes = {
        n for n in counts
        if n in graph.nodes and graph.nodes[n].get("type") == "Gene"
    }

    bridge_nodes = {
        n for n, count in counts.items()
        if count >= min_weight and n in graph.nodes and graph.degree(n) <= max_bridge_degree
    }

    candidates = set()
    for b in bridge_nodes:
        for c in graph.neighbors(b):
            if (graph.nodes[c].get("type") == "Gene"
                    and c not in directly_mentioned_genes
                    and graph[b][c]["weight"] >= min_weight):
                candidates.add(c)

    scores = _adamic_adar(graph, bridge_nodes, candidates)

    rows = [
        {"node_id": c, "name": graph.nodes[c].get("name", c), "adamic_adar_score": scores.get(c, 0.0)}
        for c in candidates
    ]
    rows.sort(key=lambda r: (-r["adamic_adar_score"], r["node_id"]))

    long_list = pd.DataFrame(rows[:top_n], columns=["node_id", "name", "adamic_adar_score"])
    return long_list, directly_mentioned_genes
