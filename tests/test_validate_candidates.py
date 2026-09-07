"""Network-free checks for validate_candidates.py unknown handling.

HTTP failures and malformed/null/empty API bodies must yield None (unknown,
unscored) - never a crash and never silent biological False/0. Only a returned
measurement of median 0 is a genuine zero. Uses unittest.mock, no network.
"""
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import requests

import validate_candidates as vc


class _Resp:
    def __init__(self, status=200, payload=None, raises=None):
        self.status_code = status
        self._payload = payload
        self._raises = raises

    def json(self):
        if self._raises:
            raise self._raises
        return self._payload


def main():
    # HTTP failure (exception) is unknown, not a negative.
    with patch.object(vc.requests, "get", side_effect=requests.RequestException("down")):
        assert vc.string_connected_to_reference("X") is None
        assert vc.gtex_expression("X") is None
        assert vc.clinvar_relevant("X") is None
        assert vc.gene_symbol("1") is None
    # Non-200 is unknown, not a negative.
    with patch.object(vc.requests, "get", return_value=_Resp(status=500, payload={})):
        assert vc.string_connected_to_reference("X") is None
        assert vc.gtex_expression("X") is None
        assert vc.clinvar_relevant("X") is None
    # Malformed/null/wrong-shaped bodies never crash and never become negatives.
    for bad in (None, "oops", 42, {"unexpected": 1}):
        with patch.object(vc.requests, "get", return_value=_Resp(payload=bad)):
            assert vc.string_connected_to_reference("X") is None
            assert vc.gtex_expression("X") is None
            assert vc.clinvar_relevant("X") is None
            assert vc.gene_symbol("1") is None
    # Malformed STRING rows are unknown, not negatives: a skipped malformed row
    # could hide a reference hit. Only a clean empty list or clean named
    # non-reference rows are genuine negatives.
    for bad_rows in (["not", "a", "dict"], [{"no_name": 1}], [{"preferredName_B": ""}],
                     [{"preferredName_B": None}], [{"preferredName_B": 42}],
                     [{"preferredName_B": "KLK15"}, "junk"]):
        with patch.object(vc.requests, "get", return_value=_Resp(payload=bad_rows)):
            assert vc.string_connected_to_reference("X") is None, \
                f"malformed STRING rows {bad_rows!r} must be unknown"
    with patch.object(vc.requests, "get", return_value=_Resp(payload=[])):
        assert vc.string_connected_to_reference("X") is False
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload=[{"preferredName_B": "TP53"}])):
        assert vc.string_connected_to_reference("X") is False
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload=[{"preferredName_B": "KLK15"}])):
        assert vc.string_connected_to_reference("X") is True
    with patch.object(vc.requests, "get", return_value=_Resp(raises=ValueError("bad json"))):
        assert vc.string_connected_to_reference("X") is None
        assert vc.gtex_expression("X") is None
        assert vc.clinvar_relevant("X") is None
    # GTEx: empty successful expression result is NO DATA (None); a returned
    # median of 0 is a genuine zero.
    seq = [_Resp(payload={"data": [{"gencodeId": "G1"}]}), _Resp(payload={"data": []})]
    with patch.object(vc.requests, "get", side_effect=list(seq)):
        assert vc.gtex_expression("X") is None, "empty expression rows must be NO DATA, not 0.0"
    seq = [_Resp(payload={"data": [{"gencodeId": "G1"}]}),
           _Resp(payload={"data": [{"median": 0}, {"median": 0.0}]})]
    with patch.object(vc.requests, "get", side_effect=list(seq)):
        assert vc.gtex_expression("X") == 0, "returned median 0 must stay a genuine zero"
    # GTEx rows without a numeric median are unusable -> unknown.
    seq = [_Resp(payload={"data": [{"gencodeId": "G1"}]}),
           _Resp(payload={"data": [{"median": "high"}, {}]})]
    with patch.object(vc.requests, "get", side_effect=list(seq)):
        assert vc.gtex_expression("X") is None
    # Bool/non-finite/negative medians are rejected as unknown, and one bad row
    # poisons the lookup instead of silently narrowing the claimed max.
    for bad_median in (True, float("inf"), float("nan"), -1.0):
        seq = [_Resp(payload={"data": [{"gencodeId": "G1"}]}),
               _Resp(payload={"data": [{"median": 9.0}, {"median": bad_median}]})]
        with patch.object(vc.requests, "get", side_effect=list(seq)):
            assert vc.gtex_expression("X") is None, \
                f"median {bad_median!r} must make GTEx unknown, not narrow the max"

    # True negatives stay distinguishable from unknowns: a present-but-empty
    # ClinVar idlist is a genuine False; a missing idlist is unknown.
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload={"esearchresult": {"idlist": []}})):
        assert vc.clinvar_relevant("X") is False
    with patch.object(vc.requests, "get", return_value=_Resp(payload={"esearchresult": {}})):
        assert vc.clinvar_relevant("X") is None, "missing idlist must be unknown"
    with patch.object(vc.requests, "get", return_value=_Resp(payload={"esearchresult": None})):
        assert vc.clinvar_relevant("X") is None
    # Summary needs actual records: missing/null/non-dict/empty results are
    # unknown, and an error string carrying a keyword must never read as biology.
    for bad_result in (None, "joint error", 42, {}, {"result": None},
                       {"result": "clinvar joint error"}, {"result": {}},
                       {"result": {"uids": ["1"]}},
                       {"result": {"1": {"error": "joint query failed"}}}):
        with patch.object(vc.requests, "get", side_effect=[
                _Resp(payload={"esearchresult": {"idlist": [1]}}),
                _Resp(payload=bad_result)]):
            assert vc.clinvar_relevant("X") is None, \
                f"summary {bad_result!r} must be unknown"
    for condition, expected in (("Ehlers-Danlos syndrome", True), ("unrelated condition", False)):
        with patch.object(vc.requests, "get", side_effect=[
                _Resp(payload={"esearchresult": {"idlist": ["1"]}}),
                _Resp(payload={"result": {"uids": ["1"], "1": {"condition": condition}}})]):
            assert vc.clinvar_relevant("X") is expected
    # gene_symbol returns only a valid non-empty string.
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload={"result": {"1": {"name": 42}}})):
        assert vc.gene_symbol("1") is None
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload={"result": {"1": {"name": "  "}}})):
        assert vc.gene_symbol("1") is None
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload={"result": {"1": {"name": "GENE1"}}})):
        assert vc.gene_symbol("1") == "GENE1"
    with patch.object(vc.requests, "get",
                       return_value=_Resp(payload={"esearchresult": {"idlist": [1]}})):
        # second call fails -> unknown, not False
        calls = {"n": 0}

        def _two(url, **kw):
            calls["n"] += 1
            if calls["n"] == 1:
                return _Resp(payload={"esearchresult": {"idlist": [1]}})
            raise requests.RequestException("down")

        with patch.object(vc.requests, "get", side_effect=_two):
            assert vc.clinvar_relevant("X") is None

    # validate(): unknowns add nothing and are marked unscored/unknown, while
    # genuine negatives score the same but read as False/0.0 in the rationale.
    import pandas as pd
    ll = pd.DataFrame([{"node_id": "Gene:1", "name": "G", "adamic_adar_score": 2.0}])
    with patch.object(vc, "gene_symbol", return_value="GENE1"), \
         patch.object(vc, "string_connected_to_reference", return_value=None), \
         patch.object(vc, "gtex_expression", return_value=None), \
         patch.object(vc, "clinvar_relevant", return_value=None), \
         patch.object(vc.time, "sleep", return_value=None):
        rows = vc.validate(ll, top_n=5)
    assert rows[0]["composite_score"] == 2.0, "unknown sources must add nothing"
    assert "unknown" in rows[0]["rationale"] and "unscored" in rows[0]["rationale"]
    assert rows[0]["string_connected"] is None and rows[0]["gtex_max_median_tpm"] is None
    with patch.object(vc, "gene_symbol", return_value="GENE1"), \
         patch.object(vc, "string_connected_to_reference", return_value=False), \
         patch.object(vc, "gtex_expression", return_value=0.0), \
         patch.object(vc, "clinvar_relevant", return_value=False), \
         patch.object(vc.time, "sleep", return_value=None):
        neg = vc.validate(ll, top_n=5)
    assert neg[0]["composite_score"] == 2.0
    assert "False" in neg[0]["rationale"] and "0.0" in neg[0]["rationale"]
    assert "unknown" not in neg[0]["rationale"], "true negatives must not read as unknown"

    # Symbol lookup failure must not drop the candidate: the row is preserved
    # with a null symbol and unscored sources, and no protein API is queried.
    with patch.object(vc, "gene_symbol", return_value=None), \
         patch.object(vc, "string_connected_to_reference",
                      side_effect=AssertionError("must not query STRING without a symbol")), \
         patch.object(vc, "gtex_expression",
                      side_effect=AssertionError("must not query GTEx without a symbol")), \
         patch.object(vc, "clinvar_relevant",
                      side_effect=AssertionError("must not query ClinVar without a symbol")), \
         patch.object(vc.time, "sleep", return_value=None):
        kept = vc.validate(ll, top_n=5)
    assert len(kept) == 1, "unresolved candidates must not disappear"
    assert kept[0]["gene_symbol"] is None
    assert kept[0]["node_id"] == "Gene:1" and "Gene:1" in kept[0]["rationale"]
    assert kept[0]["composite_score"] == 2.0, "unresolved rows keep only the bridge score"
    assert kept[0]["string_connected"] is None and kept[0]["gtex_max_median_tpm"] is None \
        and kept[0]["clinvar_relevant"] is None

    print("All checks passed.")


if __name__ == "__main__":
    main()
