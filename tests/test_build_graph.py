"""Smoke test for the ABC-model discovery logic: does it correctly separate
'directly studied' from 'literature-implied but never directly studied' candidates?
Plain assert-based, no network, no test framework.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from build_graph import build_graph
from discover_candidates import discover


def _doc(pmid, entities):
    """entities: list of (type, identifier, text)"""
    return {
        "pmid": pmid,
        "passages": [{
            "text": "synthetic",
            "annotations": [
                {"text": text, "infons": {"type": etype, "identifier": ident}}
                for etype, ident, text in entities
            ],
        }],
    }


def build_fixture():
    seed = ("Disease", "D001", "joint hypermobility syndrome")
    bridge = ("Disease", "D002", "bridge phenotype")
    gene_x = ("Gene", "100", "GeneX")   # directly co-occurs with the seed -> "already studied"
    gene_y = ("Gene", "200", "GeneY")   # only reachable via the bridge -> candidate
    gene_z = ("Gene", "300", "GeneZ")   # isolated, no path to seed at all -> absent

    return [
        _doc("1", [seed, gene_x]),
        _doc("2", [seed, gene_x]),   # weight 2: clears the min-weight floor
        _doc("3", [seed, bridge]),
        _doc("4", [seed, bridge]),   # weight 2
        _doc("5", [bridge, gene_y]),
        _doc("6", [bridge, gene_y]),  # weight 2
        _doc("7", [gene_z]),
    ]


def main():
    graph = build_graph(build_fixture())
    long_list, already_studied = discover(graph, min_weight=2, max_bridge_degree=500, top_n=50)

    candidate_ids = set(long_list["node_id"])

    assert "Gene:100" in already_studied, "directly co-occurring gene should be 'already studied'"
    assert "Gene:100" not in candidate_ids, "already-studied gene must not appear as a candidate"

    assert "Gene:200" in candidate_ids, "bridge-connected gene should surface as a candidate"
    y_score = long_list.loc[long_list["node_id"] == "Gene:200", "adamic_adar_score"].iloc[0]
    assert y_score > 0, "candidate must have a positive Adamic-Adar score"

    assert "Gene:300" not in candidate_ids, "unconnected gene must not appear as a candidate"
    assert "Gene:300" not in already_studied, "unconnected gene must not appear as already-studied"

    print("All checks passed.")


if __name__ == "__main__":
    main()
