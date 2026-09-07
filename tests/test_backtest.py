"""Smoke test for backtest.py's date filtering - the one piece of new, untested logic
the retrospective-validation script adds on top of the already-tested discovery core.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from backtest import filter_before


def main():
    docs = [
        {"pmid": "1", "date": "2019-06-01T00:00:00Z"},
        {"pmid": "2", "date": "2023-01-15T00:00:00Z"},
        {"pmid": "3", "date": "2024-12-31T00:00:00Z"},
        {"pmid": "4", "date": None},  # no date - must never pass the filter
        {"pmid": "5", "date": ""},
        {"pmid": "6", "date": "not-a-date"},
        {"pmid": "7", "date": "202x-01-01"},
        {"pmid": 8, "date": 20200101},
        {"pmid": "2", "date": "2023-01-15T00:00:00Z"},  # duplicate PMID - counted once
        {"pmid": "9", "date": "2020-garbage"},
        {"pmid": "10", "date": "0000-00-00"},
        {"pmid": "11", "date": "2021-13-01"},
        {"pmid": "12", "date": "2021-02-30"},
        {"pmid": "13", "date": "2021-01-01T25:00:00Z"},
        {"pmid": "14", "date": "-001-01-01"},
        {"pmid": "15", "date": "2020-02-29"},  # real leap day - kept
    ]

    kept = filter_before(docs, 2024)
    kept_pmids = {str(d["pmid"]) for d in kept}

    assert kept_pmids == {"1", "2", "15"}, f"expected pmids 1, 2 and 15 only, got {kept_pmids}"
    assert "3" not in kept_pmids, "2024-dated doc must be excluded by a 2024 cutoff (strictly before)"
    assert "4" not in kept_pmids, "undated docs must never be treated as pre-cutoff"
    for bad in ("5", "6", "7", "8", "9", "10", "11", "12", "13", "14"):
        assert bad not in kept_pmids, f"malformed date pmid {bad} must never count as early evidence"
    assert len(kept) == 3, f"duplicate PMIDs must be counted once, got {len(kept)} docs"

    print("All checks passed.")


if __name__ == "__main__":
    main()
