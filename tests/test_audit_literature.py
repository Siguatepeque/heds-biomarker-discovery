"""Run directly: python tests/test_audit_literature.py (no network)."""
import json
import hashlib
import math
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import audit_literature as al


def document(identifier, nodes, date="2019-01-01"):
    return {"id": identifier, "date": date, "passages": [{
        "text": "Synthetic source text for manual review.",
        "annotations": [{"text": name, "infons": {"type": node.split(":")[0],
                                                   "identifier": node.split(":")[1]}}
                        for node, name in nodes]}]}


def must_fail(call, message):
    try:
        call()
    except ValueError as exc:
        assert message in str(exc), str(exc)
    else:
        raise AssertionError("Expected ValueError")


def main():
    seed = ("Disease:h", "hypermobile Ehlers-Danlos syndrome")
    bridge = ("Disease:b", "bridge")
    weak = ("Chemical:w", "second bridge")
    gene = ("Gene:91252", "SLC39A13")
    plod = ("Gene:5351", "PLOD1")
    direct = ("Gene:55554", "KLK15")
    background = ("Gene:999", "noncandidate frequent gene")
    docs = [document("s1", [seed, bridge, weak, direct]),
            document("s2", [seed, bridge, weak]),
            document("bc1", [bridge, gene, plod]), document("bc2", [bridge, gene, plod]),
            document("weak", [weak, gene])]
    docs += [document(f"freq{i}", [background]) for i in range(5)]
    docs += [document("undated", [background], None)]
    graph, rows, evidence, report = al.audit(docs, top_n=1)
    assert report["counts"]["documents"] == 11
    assert report["counts"]["input_missing_or_invalid_dates"] == 1
    assert len(rows) == 2, "metric k must not truncate the candidate list"
    assert rows[0]["node_id"] == gene[0] and rows[0]["frequency_rank"] == 2
    assert rows[0]["known_eds_rank"] == 1
    assert rows[0]["other_eds_subtype"] == "spondylodysplastic EDS"
    assert report["rankings"]["frequency"][0] == background[0]
    assert direct[0] not in report["rankings"]["frequency"]
    assert direct[0] in report["directly_mentioned_ids"]
    for row, candidate in zip(rows, evidence["candidates"]):
        assert candidate["review_status"] == "unreviewed"
        assert math.isclose(sum(b["adamic_adar_contribution"] for b in candidate["bridges"]),
                            row["adamic_adar_score"])
        for leg in candidate["bridges"]:
            assert leg["seed_document_ids"] == ["s1", "s2"]
            assert leg["gene_bridge_weight"] == len(leg["gene_bridge_document_ids"])
    gene_legs = evidence["candidates"][0]["bridges"]
    weak_leg = next(b for b in gene_legs if b["bridge_id"] == weak[0])
    assert weak_leg["gene_bridge_document_ids"] == ["weak"] and not weak_leg["admission_edge"]
    assert evidence["documents"]["bc1"]["passages"][0]["text"].startswith("Synthetic")
    metrics = report["retrieval"]["metrics"]
    assert metrics["abc"]["target_denominator"] == 4, "KLK15 aliases share one denominator entry"
    assert metrics["abc"]["hits"] == 1 and metrics["frequency"]["hits"] == 0
    assert report["retrieval"]["graph_covered_target_groups"] == 2
    assert report["retrieval"]["eligible_target_groups"] == 1
    assert al.audit(list(reversed(docs)) + [docs[0]], top_n=1)[1:] == (rows, evidence, report)
    changed = dict(docs[0], date="2020-01-01")
    must_fail(lambda: al.audit(docs + [changed]), "Conflicting")
    must_fail(lambda: al.audit([{"passages": []}]), "usable PMID/id")
    must_fail(lambda: al.audit([{"id": "x", "passages": None}]), "passages must be a list")
    must_fail(lambda: al.audit([{"id": "x", "passages": [{"annotations": [None]}]}]), "invalid annotation")
    bad_annotation = document("bad", [("Gene:1", "x")])
    bad_annotation["passages"][0]["annotations"][0]["infons"]["type"] = []
    must_fail(lambda: al.audit([bad_annotation]), "must be strings")
    must_fail(lambda: al.audit([]), "empty input")
    must_fail(lambda: al.audit(docs, top_n=0), "positive integer")
    cutoff = al.audit(docs, cutoff=2020)[3]
    assert cutoff["counts"]["documents"] == 10 and "undated" not in cutoff["document_ids"]
    empty = al.audit(docs, cutoff=2019)
    assert not empty[1] and empty[3]["status"] == "untestable_no_seed_mentions"
    assert empty[3]["retrieval"]["metrics"]["abc"]["recall_at_k"] is None
    seed_only = al.audit(docs[:2])
    assert seed_only[3]["status"] == "testable_retrospective_retrieval" and not seed_only[1]
    assert seed_only[3]["retrieval"]["metrics"]["abc"]["recall_at_k"] == 0
    with patch.dict(al.TARGET_GENES, {"Gene:a": "X", "Gene:b": "X (alt record)"}, clear=True):
        grouped = al.target_retrieval({"abc": ["Gene:a", "Gene:b"], "frequency": [], "known_eds": []},
                                      graph, set(), 2, True)
        assert grouped["metrics"]["abc"]["hits"] == 1
        assert grouped["metrics"]["abc"]["target_denominator"] == 1
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        tmp = Path(tmp)
        must_fail(lambda: al.load_documents(tmp), "No documents")
        cache = tmp / "pubtator_test.json"
        cache.write_text(json.dumps(docs + [docs[0]]), encoding="utf-8")
        loaded, hashes = al.load_documents(tmp)
        assert len(loaded) == len(docs) and len(hashes[cache.name]) == 64
        cache.write_text("{}", encoding="utf-8")
        must_fail(lambda: al.load_documents(tmp), "expected a PubTator document list")
        cache.write_text("{bad", encoding="utf-8")
        must_fail(lambda: al.load_documents(tmp), "Invalid JSON")
        out = tmp / "outputs"
        al.write_outputs(out, graph, rows, evidence, report)
        first = {p.name: p.read_bytes() for p in out.iterdir()}
        manifest = json.loads(first["report.json"])
        for name, digest in manifest["output_sha256"].items():
            assert hashlib.sha256(first[name]).hexdigest() == digest
        al.write_outputs(out, *al.audit(list(reversed(docs)), top_n=1))
        assert first == {p.name: p.read_bytes() for p in out.iterdir()}
        al.write_outputs(out, *empty)
        assert (out / "candidates.csv").read_text().strip() == ",".join(al.COLUMNS)
    print("All offline audit checks passed.")


if __name__ == "__main__":
    with patch("requests.get", side_effect=AssertionError("Offline audit attempted a network call")):
        main()
