"""Tier 1 LoF classifier tests (SPEC_KBC.md §1.8)."""

from __future__ import annotations

from kbc.tier1 import is_tier1_lof


def test_canonical_tier1_labels():
    for label in (
        "stop_gained",
        "frameshift_variant",
        "start_lost",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ):
        assert is_tier1_lof(label), f"expected Tier 1: {label}"


def test_non_tier1_labels():
    for label in (
        "missense_variant",
        "synonymous_variant",
        "splice_region_variant",
        "intron_variant",
        "5_prime_UTR_variant",
        "stop_lost",  # stop_lost is NOT Tier 1 per §1.8
        "",
    ):
        assert not is_tier1_lof(label), f"expected NOT Tier 1: {label!r}"


def test_none_is_not_tier1():
    assert is_tier1_lof(None) is False


def test_multi_consequence_string_is_tier1_if_any_token_matches():
    # SnpEff sometimes emits compound consequences joined by '&' or '|'.
    assert is_tier1_lof("missense_variant&stop_gained")
    assert is_tier1_lof("splice_region_variant|splice_donor_variant")
    assert is_tier1_lof("frameshift_variant,downstream_gene_variant")


def test_multi_consequence_all_non_tier1():
    assert not is_tier1_lof("missense_variant&splice_region_variant")
