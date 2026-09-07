"""Smoke test for the ABC-model discovery logic: does it correctly separate
'directly mentioned' from 'literature-implied but never directly mentioned' candidates?
Plain assert-based, no network, no test framework.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from build_graph import build_graph, dedupe_documents
from discover_candidates import discover, _is_seed_mention


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
    documents = build_fixture()
    graph = build_graph(documents)
    long_list, already_studied = discover(graph, documents, min_weight=2, max_bridge_degree=500, top_n=50)

    candidate_ids = set(long_list["node_id"])

    assert "Gene:100" in already_studied, "directly co-occurring gene should be 'already studied'"
    assert "Gene:100" not in candidate_ids, "already-studied gene must not appear as a candidate"

    assert "Gene:200" in candidate_ids, "bridge-connected gene should surface as a candidate"
    y_score = long_list.loc[long_list["node_id"] == "Gene:200", "adamic_adar_score"].iloc[0]
    assert y_score > 0, "candidate must have a positive Adamic-Adar score"

    assert "Gene:300" not in candidate_ids, "unconnected gene must not appear as a candidate"
    assert "Gene:300" not in already_studied, "unconnected gene must not appear as already-studied"

    # Single-paper direct mention excludes at every min_weight floor: GeneY shares
    # one seed document, so it is directly mentioned even where min_weight > 1.
    _seed = ("Disease", "D001", "joint hypermobility syndrome")
    one_hit = build_fixture() + [_doc("8", [_seed, ("Gene", "200", "GeneY")])]
    for floor in (1, 2, 3):
        g1 = build_graph(one_hit)
        ll1, dm1 = discover(g1, one_hit, min_weight=floor, max_bridge_degree=500, top_n=50)
        assert "Gene:200" in dm1, f"single seed co-mention must count at min_weight={floor}"
        assert "Gene:200" not in set(ll1["node_id"]), f"single-hit gene must not be a candidate at min_weight={floor}"

    # Duplicate documents (same PMID twice, incl. int/str PMID forms) must not
    # change the candidate set or scores.
    dup = [dict(d, pmid=str(d["pmid"])) for d in documents]
    dup = dup + [dict(d) for d in dup[:3]] + [{"pmid": 5, "passages": dup[4]["passages"]}]
    g2 = build_graph(dup)
    ll2, _ = discover(g2, dup, min_weight=2, max_bridge_degree=500, top_n=50)
    assert set(ll2["node_id"]) == candidate_ids, "duplicate PMIDs must not change candidates"
    for _, r in long_list.iterrows():
        s2 = ll2.loc[ll2["node_id"] == r["node_id"], "adamic_adar_score"]
        assert len(s2) == 1 and abs(s2.iloc[0] - r["adamic_adar_score"]) < 1e-9, \
            "duplicate PMIDs must not change scores"

    # BioC 'id' fallback: distinct copies sharing {'id': '123'} (int/str forms)
    # collapse, while distinct anonymous docs (no pmid/id) never collapse.
    bioc_a = {"id": "123", "passages": documents[0]["passages"]}
    bioc_b = {"id": 123, "passages": documents[1]["passages"]}
    assert len(dedupe_documents([bioc_a, bioc_b])) == 1, \
        "copied BioC docs with the same id must collapse"
    anon_a = {"passages": documents[0]["passages"]}
    anon_b = {"passages": documents[0]["passages"]}
    assert len(dedupe_documents([anon_a, anon_b])) == 2, \
        "distinct anonymous docs must not collapse"

    # Regression cases: two real false positives found during development, where a
    # DIFFERENT EDS-spectrum condition's name contains an hEDS pattern as a literal
    # substring, or a shared clinical sign gets mistaken for a seed diagnosis.
    assert not _is_seed_mention("Ehlers-Danlos syndrome, spondylodysplastic form type 3"), \
        "a different EDS subtype must not match just because its name contains 'type 3'"
    assert not _is_seed_mention("generalized joint hypermobility"), \
        "a clinical sign shared with non-hEDS populations must not count as a seed diagnosis"
    assert _is_seed_mention("hypermobile Ehlers-Danlos syndrome"), \
        "the actual condition must still match"
    assert _is_seed_mention("hypermobile-EDS"), "hyphenated abbreviation must match"
    assert _is_seed_mention("hEDS"), "hEDS abbreviation must match as a whole token"
    assert not _is_seed_mention("sheds"), "heds must not match inside arbitrary substrings"
    assert _is_seed_mention("Ehlers-Danlos syndrome, hypermobility type"), \
        "comma variant must match"
    assert _is_seed_mention("Ehlers\u2013Danlos syndrome hypermobility type"), \
        "unicode-dash variant must match"
    assert not _is_seed_mention("Vascular Ehlers-Danlos Syndrome"), \
        "other subtypes must not match"
    assert not _is_seed_mention("Classical Ehlers-Danlos syndrome"), \
        "other subtypes must not match"
    assert not _is_seed_mention("Ehlers-Danlos syndrome type III and type IV"), \
        "combined III/IV heading must stay excluded"
    assert not _is_seed_mention("Ehlers-Danlos syndrome type iiia"), \
        "a longer token must not match through a shared 'type iii' prefix"
    assert _is_seed_mention("hypermobility spectrum disorders"), \
        "known plural must match via its explicit plural pattern"
    assert _is_seed_mention("hypermobile spectrum disorders"), \
        "known plural must match via its explicit plural pattern"

    print("All checks passed.")


if __name__ == "__main__":
    main()
