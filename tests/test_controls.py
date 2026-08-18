"""Smoke test for the control-disease seed patterns (controls.py): the control and
hEDS seed patterns must not cross-match each other's texts, and the monkeypatch swap
control_disease.py / control_backtest.py use to point discover_candidates at the
control config must restore the hEDS config afterward.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import discover_candidates
from controls import CONTROL_EXCLUDE_PATTERNS, CONTROL_SEED_PATTERNS


def main():
    original_seeds = list(discover_candidates.SEED_PATTERNS)
    original_excludes = list(discover_candidates.EXCLUDE_PATTERNS)

    # hEDS patterns (module default) must not match fibromyalgia text, and must not
    # be affected by anything in controls.py just by importing it.
    assert not discover_candidates._is_seed_mention("fibromyalgia"), \
        "hEDS seed patterns must not match the control disease's own name"
    assert not discover_candidates._is_seed_mention("fibromyalgia syndrome"), \
        "hEDS seed patterns must not match the control disease's own name"

    # Swap in the control config, same mechanism control_disease.py uses.
    discover_candidates.SEED_PATTERNS = CONTROL_SEED_PATTERNS
    discover_candidates.EXCLUDE_PATTERNS = CONTROL_EXCLUDE_PATTERNS
    try:
        assert discover_candidates._is_seed_mention("fibromyalgia"), \
            "control seed patterns must match the control disease's own name"
        assert discover_candidates._is_seed_mention("fibrositis"), \
            "control seed patterns must match fibromyalgia's pre-1990s name"
        assert not discover_candidates._is_seed_mention("hypermobile Ehlers-Danlos syndrome"), \
            "control seed patterns must not match hEDS text"
        assert not discover_candidates._is_seed_mention("joint hypermobility syndrome"), \
            "control seed patterns must not match hEDS text"
        assert not discover_candidates._is_seed_mention("fibromyalgia-like symptoms in long COVID"), \
            "control exclude patterns must filter out comparator usage, not a real diagnosis"
    finally:
        discover_candidates.SEED_PATTERNS = original_seeds
        discover_candidates.EXCLUDE_PATTERNS = original_excludes

    # Swap must be fully reversible: hEDS behavior restored afterward.
    assert discover_candidates._is_seed_mention("hypermobile Ehlers-Danlos syndrome"), \
        "hEDS seed patterns must be restored after the control swap"
    assert not discover_candidates._is_seed_mention("fibromyalgia"), \
        "hEDS seed patterns must not match fibromyalgia after the swap is undone"
    assert discover_candidates.SEED_PATTERNS == original_seeds
    assert discover_candidates.EXCLUDE_PATTERNS == original_excludes

    print("All checks passed.")


if __name__ == "__main__":
    main()
