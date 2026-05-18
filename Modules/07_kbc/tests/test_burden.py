"""End-to-end MVP 1 burden test on the mini fixture.

Validates SPEC_KBC.md §3.A behaviour: per (sample, inversion, segment)
hom-exposed / het-masked / compound-het-unknown counts at Tier 1 only.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from kbc.burden import build_table_a
from kbc.io import (
    load_genotypes_tsv,
    load_inversion_intervals,
    load_karyotypes,
    load_roh_intervals,
    load_sample_metadata,
    load_variant_master,
)


@pytest.fixture
def table_a(variant_master_path, genotypes_path, karyotypes_path, inversions_path, samples_path, roh_dir):
    vm = load_variant_master(variant_master_path)
    gt = load_genotypes_tsv(genotypes_path)
    kary = load_karyotypes(karyotypes_path)
    invs = load_inversion_intervals(inversions_path)
    samples = load_sample_metadata(samples_path)
    roh = load_roh_intervals(roh_dir)
    return build_table_a(
        variant_master=vm,
        genotypes=gt,
        karyotypes=kary,
        inversion_intervals=invs,
        sample_metadata=samples,
        scoring_weights=None,
        roh_intervals=roh,
    )


def _row(df: pd.DataFrame, sample_id: str, inversion_id: str, pod_segment: str) -> pd.Series:
    sub = df[
        (df["sample_id"] == sample_id)
        & (df["inversion_id"] == inversion_id)
        & (df["pod_segment"] == pod_segment)
    ]
    assert len(sub) == 1, f"expected exactly one row for {sample_id}/{inversion_id}/{pod_segment}, got {len(sub)}"
    return sub.iloc[0]


def test_row_count_is_samples_x_segments(table_a):
    # 6 samples × 4 (inversion, segment) cells = 24 rows
    assert len(table_a) == 6 * 4


def test_karyotype_attached(table_a):
    assert _row(table_a, "S001", "inv1", "whole")["karyotype_class"] == "AA"
    assert _row(table_a, "S003", "inv2", "L")["karyotype_class"] == "AB"
    assert _row(table_a, "S005", "inv1", "whole")["karyotype_class"] == "BB"


def test_is_pruned_unrelated_flagged(table_a):
    assert _row(table_a, "S001", "inv1", "whole")["is_pruned_unrelated"] is True or \
           _row(table_a, "S001", "inv1", "whole")["is_pruned_unrelated"] == True  # noqa
    assert not _row(table_a, "S002", "inv1", "whole")["is_pruned_unrelated"]


def test_s001_inv1_whole_AA_hom_at_G1_het_at_G3(table_a):
    # S001: v1=1/1 (G1 hom), v2=0/0, v3 excluded (missense), v4=0/1 (G3 het).
    r = _row(table_a, "S001", "inv1", "whole")
    assert r["n_damaging_variants"] == 2          # v1 and v4
    assert r["hom_exposed_gene_count"] == 1        # G1
    assert r["het_masked_gene_count"] == 1         # G3
    assert r["compound_het_unknown_gene_count"] == 0
    assert r["burden_score_sum"] == 3.0            # 2 (hom v1) + 1 (het v4)


def test_s002_inv1_whole_double_hom(table_a):
    # S002: v1=1/1, v4=1/1.
    r = _row(table_a, "S002", "inv1", "whole")
    assert r["n_damaging_variants"] == 2
    assert r["hom_exposed_gene_count"] == 2        # G1, G3
    assert r["het_masked_gene_count"] == 0
    assert r["compound_het_unknown_gene_count"] == 0
    assert r["burden_score_sum"] == 4.0


def test_s003_inv1_whole_AB_compound_het_at_G1(table_a):
    # S003 AB: v1=0/1, v2=0/1 in G1 → compound_het_unknown; v4=0/1 in G3 → het_masked.
    r = _row(table_a, "S003", "inv1", "whole")
    assert r["n_damaging_variants"] == 3           # v1, v2, v4
    assert r["hom_exposed_gene_count"] == 0
    assert r["het_masked_gene_count"] == 1         # G3
    assert r["compound_het_unknown_gene_count"] == 1   # G1
    assert r["burden_score_sum"] == 3.0


def test_s004_inv1_whole_AB_two_het_masked(table_a):
    # S004 AB: v1=0/1 (G1 single het, no other), v4=0/1 (G3 het).
    r = _row(table_a, "S004", "inv1", "whole")
    assert r["n_damaging_variants"] == 2
    assert r["hom_exposed_gene_count"] == 0
    assert r["het_masked_gene_count"] == 2         # G1 single-het, G3 single-het
    assert r["compound_het_unknown_gene_count"] == 0
    assert r["burden_score_sum"] == 2.0


def test_s005_inv1_whole_BB_hom_at_G1(table_a):
    # S005 BB: v2=1/1 (G1 hom), v4=0/1 (G3 het).
    r = _row(table_a, "S005", "inv1", "whole")
    assert r["n_damaging_variants"] == 2
    assert r["hom_exposed_gene_count"] == 1
    assert r["het_masked_gene_count"] == 1
    assert r["burden_score_sum"] == 3.0


def test_s006_inv1_whole_BB_roh_promotion(table_a):
    # S006 BB has a single het at v2 in G1. ROH (chr1:1000-3000) covers G1 → hom_exposed by promotion.
    r = _row(table_a, "S006", "inv1", "whole")
    assert r["n_damaging_variants"] == 1           # only v2
    assert r["hom_exposed_gene_count"] == 1        # ROH-promoted
    assert r["het_masked_gene_count"] == 0
    assert r["compound_het_unknown_gene_count"] == 0


def test_inv2_segment_resolution_L_M_R(table_a):
    # S003 has v5=1/1 in inv2.L only.
    r = _row(table_a, "S003", "inv2", "L")
    assert r["hom_exposed_gene_count"] == 1
    assert _row(table_a, "S003", "inv2", "M")["hom_exposed_gene_count"] == 0
    assert _row(table_a, "S003", "inv2", "R")["hom_exposed_gene_count"] == 0

    # S004 has v6=1/1 in inv2.M only.
    assert _row(table_a, "S004", "inv2", "L")["hom_exposed_gene_count"] == 0
    r = _row(table_a, "S004", "inv2", "M")
    assert r["hom_exposed_gene_count"] == 1
    assert _row(table_a, "S004", "inv2", "R")["hom_exposed_gene_count"] == 0

    # S005 has v7=1/1 in inv2.R only.
    assert _row(table_a, "S005", "inv2", "L")["hom_exposed_gene_count"] == 0
    assert _row(table_a, "S005", "inv2", "M")["hom_exposed_gene_count"] == 0
    r = _row(table_a, "S005", "inv2", "R")
    assert r["hom_exposed_gene_count"] == 1


def test_missense_variant_excluded_from_n_damaging(table_a):
    # v3 (chr1:3000 missense) is not Tier 1; it should never contribute to n_damaging_variants.
    # S003 has v3=0/0 anyway, but the strongest check: no sample should ever see G2 (missense-only gene) in any gene count.
    # If v3 had leaked through, S003's G1+G3+G2 (compound_het via v3=0/0+v1=0/1+v2=0/1) story would differ.
    # We assert by aggregate: the sum of (hom+het+comp) across S003.inv1.whole equals exactly 2 (G1, G3).
    r = _row(table_a, "S003", "inv1", "whole")
    assert r["hom_exposed_gene_count"] + r["het_masked_gene_count"] + r["compound_het_unknown_gene_count"] == 2
