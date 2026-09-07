"""Negative-control disease configuration: fibromyalgia.

Mirrors fetch_literature.QUERY and discover_candidates.SEED_PATTERNS /
EXCLUDE_PATTERNS, but scoped to fibromyalgia instead of hEDS. Used by
control_disease.py and control_backtest.py to rerun the same pipeline code on a
different disease as a plausibility/context check on candidate yield - not as a
specificity proof. A control yielding fewer (or more) candidates does not by itself
prove the hEDS candidates are true or false; fibromyalgia itself now has
genome-wide-significant loci (see control_backtest.py), so the "no known genetic
basis" framing is outdated and must not be used to claim specificity from absence.

Fibromyalgia is kept as the control because it is a large literature (15,000+
PubMed hits) dominated by generic hub terms (pain, fatigue, cytokines) that are
exactly the kind of thing discover_candidates.py's Adamic-Adar / hub-degree logic
is designed to discount.
"""

# Mirrors fetch_literature.QUERY's structure: PubMed title/abstract field search.
# "fibrositis" is fibromyalgia's own pre-1990s name (same condition, different era
# of terminology), not a different disease - kept for the same reason
# fetch_literature.QUERY includes hEDS's own older names ("EDS-HT", "type III").
CONTROL_QUERY = (
    '("fibromyalgia"[tiab] OR "fibromyalgia syndrome"[tiab] OR "fibrositis"[tiab])'
)

# Matched against the same normalized surface text as discover_candidates.SEED_PATTERNS.
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
