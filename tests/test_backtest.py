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
    ]

    kept = filter_before(docs, 2024)
    kept_pmids = {d["pmid"] for d in kept}

    assert kept_pmids == {"1", "2"}, f"expected pmids 1 and 2 only, got {kept_pmids}"
    assert "3" not in kept_pmids, "2024-dated doc must be excluded by a 2024 cutoff (strictly before)"
    assert "4" not in kept_pmids, "undated docs must never be treated as pre-cutoff"

    print("All checks passed.")


if __name__ == "__main__":
    main()
