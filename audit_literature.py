"""Offline, retrospective hEDS literature audit; no annotation API calls.

Writes candidates.csv (all ABC candidates), graph.graphml, evidence.json (both
bridge legs plus original passages for manual review), and report.json (coverage,
rankings, target retrieval, input/code fingerprints). --top-n sets metric k only.
Cache input is a directory union, NOT a reconstructed historical PubMed query.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import networkx as nx
import pandas as pd

import discover_candidates as dc
from backtest import TARGET_GENES, doc_year, filter_before
from build_graph import build_graph, doc_identity, document_entities, save_graph
from eds_genes import OTHER_EDS_GENES

ROOT = Path(__file__).resolve().parent
COLUMNS = ["node_id", "name", "adamic_adar_score", "document_frequency",
           "frequency_rank", "known_eds_rank", "other_eds_subtype"]


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def fingerprint(value):
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def unique_documents(documents):
    """Do not silently choose between conflicting annotation versions of a PMID."""
    unique = {}
    for doc in documents:
        if not isinstance(doc, dict) or not doc_identity(doc):
            raise ValueError("Every document must be an object with a usable PMID/id")
        identifier = doc_identity(doc)
        passages = doc.get("passages", [])
        if not isinstance(passages, list):
            raise ValueError(f"Document {identifier}: passages must be a list")
        for passage in passages:
            if not isinstance(passage, dict) or not isinstance(passage.get("annotations", []), list):
                raise ValueError(f"Document {identifier}: invalid passage/annotations")
            for annotation in passage.get("annotations", []):
                if (not isinstance(annotation, dict) or not isinstance(annotation.get("infons", {}), dict)
                        or not isinstance(annotation.get("text", ""), str)):
                    raise ValueError(f"Document {identifier}: invalid annotation")
                infons = annotation.get("infons", {})
                if any(infons.get(key) is not None and not isinstance(infons[key], str)
                       for key in ("type", "identifier")):
                    raise ValueError(f"Document {identifier}: annotation type/identifier must be strings")
        if identifier in unique and canonical(unique[identifier]) != canonical(doc):
            raise ValueError(f"Conflicting cached documents for {identifier}; select one frozen corpus")
        unique[identifier] = doc
    return [unique[key] for key in sorted(unique)]


def load_documents(cache_dir):
    documents, hashes = [], {}
    for path in sorted(Path(cache_dir).glob("pubtator_*.json")):
        raw = path.read_bytes()
        try:
            batch = json.loads(raw)
        except (ValueError, UnicodeError) as exc:
            raise ValueError(f"Invalid JSON in {path.name}: {exc}") from exc
        if not isinstance(batch, list):
            raise ValueError(f"{path.name}: expected a PubTator document list")
        hashes[path.name] = hashlib.sha256(raw).hexdigest()
        documents.extend(batch)
    if not documents:
        raise ValueError(f"No documents in {cache_dir}; supply a populated, disease-specific cache directory")
    return unique_documents(documents), hashes


def target_retrieval(rankings, graph, directly_mentioned, k, testable):
    """One denominator entry per target label, not per alternate Entrez record."""
    groups = {}
    for node, label in TARGET_GENES.items():
        groups.setdefault(label.split(" (")[0], []).append(node)
    ranks = {method: {node: i + 1 for i, node in enumerate(nodes)}
             for method, nodes in rankings.items()}
    targets = []
    for label, nodes in groups.items():
        statuses = {node: ("absent" if node not in graph else
                          "directly_mentioned" if node in directly_mentioned else
                          "candidate" if node in ranks["abc"] else "present_not_candidate")
                    for node in nodes}
        target_ranks = {method: min((table[node] for node in nodes if node in table), default=None)
                        for method, table in ranks.items()}
        targets.append({"target": label, "node_status": statuses, "ranks": target_ranks})
    metrics = {}
    for method in rankings:
        hits = sum(t["ranks"][method] is not None and t["ranks"][method] <= k for t in targets)
        metrics[method] = {"k": k, "hits": hits if testable else None,
                           "target_denominator": len(groups),
                           "recall_at_k": hits / len(groups) if testable and groups else None}
    return {"targets": targets, "metrics": metrics,
            "graph_covered_target_groups": sum(any(n in graph for n in nodes) for nodes in groups.values()),
            "eligible_target_groups": sum(any(n in ranks["frequency"] for n in nodes)
                                          for nodes in groups.values()) if testable else None}


def audit(documents, cutoff=None, min_weight=2, max_bridge_degree=500, top_n=20):
    for name, value in (("cutoff", cutoff), ("min_weight", min_weight),
                        ("max_bridge_degree", max_bridge_degree), ("top_n", top_n)):
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
            raise ValueError(f"{name} must be a positive integer")
    original = unique_documents(documents)
    if not original:
        raise ValueError("Cannot audit an empty input corpus")
    documents = filter_before(original, cutoff) if cutoff is not None else original
    graph = build_graph(documents)
    seeds = dc.seed_documents(documents)
    seed_ids = {doc_identity(doc) for doc in seeds}
    entities = {doc_identity(doc): document_entities(doc) for doc in documents}
    occurrences = {}
    for identifier, nodes in entities.items():
        for node in nodes:
            occurrences.setdefault(node, set()).add(identifier)
    if seeds:
        long_list, direct = dc.discover(graph, documents, min_weight=min_weight,
                                        max_bridge_degree=max_bridge_degree, top_n=len(graph))
        candidates = long_list.to_dict("records")
    else:
        candidates, direct = [], set()
    frequency = sorted((n for n in graph if graph.nodes[n]["type"] == "Gene" and n not in direct),
                       key=lambda n: (-len(occurrences[n]), n))
    known_eds = [n for n in frequency if n in OTHER_EDS_GENES]
    rankings = {"abc": [r["node_id"] for r in candidates], "frequency": frequency, "known_eds": known_eds}
    frequency_ranks = {n: i + 1 for i, n in enumerate(frequency)}
    known_ranks = {n: i + 1 for i, n in enumerate(known_eds)}
    for row in candidates:
        node = row["node_id"]
        row.update(document_frequency=len(occurrences[node]), frequency_rank=frequency_ranks[node],
                   known_eds_rank=known_ranks.get(node),
                   other_eds_subtype=OTHER_EDS_GENES.get(node, {}).get("subtype"))

    # Mirror discover's bridge rule; the test checks contribution sums against its scores.
    bridges = {n for n, ids in occurrences.items()
               if len(ids & seed_ids) >= min_weight and graph.degree(n) <= max_bridge_degree}
    evidence, referenced = [], set()
    for row in candidates:
        node, legs = row["node_id"], []
        for bridge in sorted(bridges & set(graph.neighbors(node))):
            degree = graph.degree(bridge)
            if degree <= 1:
                continue
            ab = sorted(occurrences[bridge] & seed_ids)
            bc = sorted(occurrences[bridge] & occurrences[node])
            referenced.update(ab + bc)
            legs.append({"bridge_id": bridge, "bridge_name": graph.nodes[bridge]["name"],
                         "degree": degree, "adamic_adar_contribution": 1 / math.log(degree),
                         "seed_document_ids": ab, "gene_bridge_document_ids": bc,
                         "gene_bridge_weight": graph[bridge][node]["weight"],
                         "admission_edge": graph[bridge][node]["weight"] >= min_weight})
        evidence.append({"node_id": node, "review_status": "unreviewed", "bridges": legs})
    evidence = {"candidates": evidence, "documents": {
        doc_identity(doc): {"date": doc.get("date"), "passages": doc.get("passages", [])}
        for doc in documents if doc_identity(doc) in referenced},
        "review_instructions": "Review both legs for species, diagnosis/subtype, negation, study design, "
                               "comorbidity vs mechanism and measurement vs hypothesis. Co-occurrence "
                               "does not establish a relation. No context review has been performed."}
    report = {
        "status": "testable_retrospective_retrieval" if seeds else "untestable_no_seed_mentions",
        "scope": "Cached-directory union; historical PubMed query membership unknown. "
                 "Not a novelty, causal, diagnostic or prospective validation claim.",
        "baseline_definition": "Frequency: all graph Gene nodes minus directly mentioned IDs, ordered "
                               "by unique document count descending then node_id. Known-EDS: same "
                               "ordering restricted to the frozen present-day differential set. "
                               "Neither baseline is restricted to ABC-admitted candidates.",
        "target_definition": "Present-day retrospective association targets; alternate IDs share one "
                             "denominator entry. Not confirmed causal genes. No precision/specificity estimate.",
        "parameters": {"cutoff_before_year": cutoff, "min_weight": min_weight,
                       "max_bridge_degree": max_bridge_degree, "metric_k": top_n},
        "seed_patterns": dc.SEED_PATTERNS, "exclude_patterns": dc.EXCLUDE_PATTERNS,
        "counts": {"input_documents": len(original), "documents": len(documents),
                   "input_missing_or_invalid_dates": sum(doc_year(d.get("date")) is None for d in original),
                   "seed_documents": len(seeds), "nodes": len(graph), "edges": graph.number_of_edges(),
                   "documents_with_kept_entities": sum(bool(nodes) for nodes in entities.values()),
                   "frequency_baseline_genes": len(frequency), "known_eds_baseline_genes": len(known_eds),
                   "directly_mentioned_genes": len(direct), "candidates": len(candidates)},
        "input_corpus_sha256": fingerprint(original), "effective_corpus_sha256": fingerprint(documents),
        "document_ids": [doc_identity(d) for d in documents], "directly_mentioned_ids": sorted(direct),
        "rankings": rankings, "retrieval": target_retrieval(rankings, graph, direct, top_n, bool(seeds)),
        "versions": {"python": platform.python_version(), "pandas": pd.__version__, "networkx": nx.__version__},
        "code_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in
                        ("audit_literature.py", "build_graph.py", "discover_candidates.py", "backtest.py", "eds_genes.py")},
    }
    return graph, candidates, evidence, report


def write_outputs(output_dir, graph, candidates, evidence, report):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_graph(graph, output_dir / "graph.graphml")
    with (output_dir / "candidates.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(candidates)
    (output_dir / "evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    report = dict(report, output_sha256={
        name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
        for name in ("graph.graphml", "candidates.csv", "evidence.json")})
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "cache")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "audit")
    parser.add_argument("--cutoff", type=positive_int)
    parser.add_argument("--min-weight", type=positive_int, default=2)
    parser.add_argument("--max-bridge-degree", type=positive_int, default=500)
    parser.add_argument("--top-n", type=positive_int, default=20)
    args = parser.parse_args()
    try:
        documents, cache_hashes = load_documents(args.cache_dir)
        graph, candidates, evidence, report = audit(documents, args.cutoff, args.min_weight,
                                                    args.max_bridge_degree, args.top_n)
        report["cache_file_sha256"] = cache_hashes
        write_outputs(args.output_dir, graph, candidates, evidence, report)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Audit failed: {exc}\n")
    print(f"{report['status']}: {report['counts']}")
    print("WARNING: cached-directory union, unknown historical query membership; context unreviewed.")
    print(f"Wrote offline audit to {args.output_dir}")


if __name__ == "__main__":
    main()
