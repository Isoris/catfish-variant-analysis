"""Schema-compliance test: Table A output validates against its JSON schema."""

from __future__ import annotations

import pandas as pd

from kbc.burden import build_table_a
from kbc.io import (
    load_genotypes_tsv,
    load_inversion_intervals,
    load_karyotypes,
    load_roh_intervals,
    load_sample_metadata,
    load_variant_master,
)
from kbc.schema import validate_rows


def test_table_a_validates_against_schema(
    variant_master_path, genotypes_path, karyotypes_path, inversions_path, samples_path, roh_dir
):
    vm = load_variant_master(variant_master_path)
    gt = load_genotypes_tsv(genotypes_path)
    kary = load_karyotypes(karyotypes_path)
    invs = load_inversion_intervals(inversions_path)
    samples = load_sample_metadata(samples_path)
    roh = load_roh_intervals(roh_dir)
    table_a = build_table_a(
        variant_master=vm,
        genotypes=gt,
        karyotypes=kary,
        inversion_intervals=invs,
        sample_metadata=samples,
        roh_intervals=roh,
    )
    validate_rows(table_a, "A_kbc_sample_inversion_burden.schema.json")
