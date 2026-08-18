"""Negative-control disease configuration: fibromyalgia.

Mirrors fetch_literature.QUERY and discover_candidates.SEED_PATTERNS /
EXCLUDE_PATTERNS, but scoped to fibromyalgia instead of hEDS. Used by
control_disease.py and control_backtest.py to rerun the exact same pipeline code on
a disease with no reproducible monogenic cause, as a check on whether the method
just emits candidates for anything it's pointed at.

Fibromyalgia is the deliberate choice: a huge literature (15,000+ PubMed hits),
dominated by generic hub terms (pain, fatigue, cytokines) that are exactly the kind
of thing discover_candidates.py's Adamic-Adar / hub-degree logic is designed to
discount, and - as of 2026 - no gene that has reached genome-wide significance and
replicated. If the pipeline is meaningful rather than a co-occurrence-count
rubber stamp, it should not reproduce an SLC39A13-shaped result here.
"""

# Mirrors fetch_literature.QUERY's structure: PubMed title/abstract field search.
# "fibrositis" is fibromyalgia's own pre-1990s name (same condition, different era
# of terminology), not a different disease - kept for the same reason
# fetch_literature.QUERY includes hEDS's own older names ("EDS-HT", "type III").
CONTROL_QUERY = (
    '("fibromyalgia"[tiab] OR "fibromyalgia syndrome"[tiab] OR "fibrositis"[tiab])'
)

# Matched case-insensitively as substrings of each annotation's own surface text,
# same mechanism as discover_candidates.SEED_PATTERNS.
CONTROL_SEED_PATTERNS = [
    "fibromyalgia",
    "fibrositis",
]

# Hygiene guard, same spirit as discover_candidates.py's GJH exclusion: "fibromyalgia"
# also gets invoked as a loose symptom-pattern comparator in papers about a DIFFERENT
# primary condition (e.g. "long COVID patients reported fibromyalgia-like widespread
# pain"), which is a comparison, not a fibromyalgia diagnosis being studied. Left
# unguarded, that usage would misclassify someone else's paper as a fibromyalgia seed
# document.
CONTROL_EXCLUDE_PATTERNS = ["fibromyalgia-like", "fibromyalgia like"]
