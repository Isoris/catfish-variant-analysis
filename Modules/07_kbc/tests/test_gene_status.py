"""Per-gene status classifier tests (SPEC_KBC.md §5)."""

from __future__ import annotations

import pandas as pd

from kbc.gene_status import (
    STATUS_COMPOUND_HET_UNKNOWN,
    STATUS_HET_MASKED,
    STATUS_HOM_EXPOSED,
    STATUS_REFERENCE_LIKE,
    GeneVariants,
    classify_per_gene_status,
)
from kbc.io import GenotypeFrame


def _gf(rows: list[tuple[str, str, str]]) -> GenotypeFrame:
    df = pd.DataFrame(rows, columns=["sample_id", "variant_id", "gt"])
    return GenotypeFrame(df=df)


def _gene(*variant_ids: str) -> GeneVariants:
    return GeneVariants(
        gene_id="G",
        variant_ids=variant_ids,
        chrom="chr1",
        gene_start=1000,
        gene_end=2000,
    )


def test_case1_all_ref_returns_reference_like():
    gts = _gf([("S1", "v1", "0/0"), ("S1", "v2", "0/0")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_REFERENCE_LIKE


def test_case2_any_hom_alt_returns_hom_exposed():
    gts = _gf([("S1", "v1", "1/1"), ("S1", "v2", "0/0")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_HOM_EXPOSED


def test_case2_hom_alt_plus_het_returns_hom_exposed():
    """1/1 always wins regardless of other genotypes (spec §5 case 2)."""
    gts = _gf([("S1", "v1", "1/1"), ("S1", "v2", "0/1")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_HOM_EXPOSED


def test_case3_single_het_no_roh_returns_het_masked():
    gts = _gf([("S1", "v1", "0/1"), ("S1", "v2", "0/0")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_HET_MASKED


def test_case3_single_het_with_roh_promoted_to_hom_exposed():
    gts = _gf([("S1", "v1", "0/1"), ("S1", "v2", "0/0")])
    roh = pd.DataFrame([{"chrom": "chr1", "start": 800, "end": 2500}])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts, roh_df=roh) == STATUS_HOM_EXPOSED


def test_case3_single_het_with_non_overlapping_roh_stays_het_masked():
    gts = _gf([("S1", "v1", "0/1"), ("S1", "v2", "0/0")])
    roh = pd.DataFrame([{"chrom": "chr1", "start": 5000, "end": 6000}])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts, roh_df=roh) == STATUS_HET_MASKED


def test_case3_roh_on_different_chrom_does_not_promote():
    gts = _gf([("S1", "v1", "0/1")])
    roh = pd.DataFrame([{"chrom": "chr2", "start": 800, "end": 2500}])
    assert classify_per_gene_status("S1", _gene("v1"), gts, roh_df=roh) == STATUS_HET_MASKED


def test_case4_multiple_hets_returns_compound_het_unknown():
    gts = _gf([("S1", "v1", "0/1"), ("S1", "v2", "0/1")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_COMPOUND_HET_UNKNOWN


def test_missing_genotypes_drop_from_observed():
    """Missing GT carries no information; gene with all-missing is reference_like."""
    gts = _gf([("S1", "v1", "./."), ("S1", "v2", "./.")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_REFERENCE_LIKE


def test_partial_missing_uses_observed_only():
    gts = _gf([("S1", "v1", "./."), ("S1", "v2", "0/1")])
    assert classify_per_gene_status("S1", _gene("v1", "v2"), gts) == STATUS_HET_MASKED
