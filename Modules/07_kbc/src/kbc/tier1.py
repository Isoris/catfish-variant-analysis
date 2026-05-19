"""Tier 1 damaging-variant classifier.

Tier 1 is the manuscript headline (SPEC_KBC.md §1.8): high-confidence
loss-of-function variants whose LoF claim is mechanistically defensible
from SnpEff annotation alone, no model layer required.

Members of Tier 1:
- stop_gained
- frameshift_variant
- start_lost
- splice_donor_variant       (SnpEff calls these at canonical GT donor sites)
- splice_acceptor_variant    (SnpEff calls these at canonical AG acceptor sites)
"""

from __future__ import annotations

TIER1_LOF_CONSEQUENCES: frozenset[str] = frozenset({
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "splice_donor_variant",
    "splice_acceptor_variant",
})


def is_tier1_lof(consequence: str | None) -> bool:
    """True if `consequence` (a SnpEff consequence string) is a Tier 1 LoF.

    SnpEff sometimes emits multiple consequences per variant joined by
    '&', '|', ',' or ';'. We treat any token matching a Tier 1 label as
    a Tier 1 hit — the most-severe convention used elsewhere in the
    pipeline matches this.
    """
    if consequence is None:
        return False
    tokens = (
        consequence.replace("&", ",")
        .replace("|", ",")
        .replace(";", ",")
        .split(",")
    )
    return any(t.strip() in TIER1_LOF_CONSEQUENCES for t in tokens)
