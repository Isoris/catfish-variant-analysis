#!/usr/bin/env python3
"""KBC MVP 1 driver — produces Table A.

Usage:
    python 03_per_sample_per_segment.py \
        --variant-master      variant_master_scored.tsv \
        --joint-vcf           joint.vcf.gz \
        --karyotypes          karyotypes.tsv \
        --inversion-intervals intervals.tsv \
        --sample-metadata     samples.tsv \
        --out                 A_kbc_sample_inversion_burden.tsv \
        [--scoring-weights    scoring_weights.tsv] \
        [--roh-dir            roh_per_sample/] \
        [--genotypes-tsv      genotypes.tsv]    # synthetic / test path

Cohort guard: KBC operates only on the 226-sample pure C. gariepinus
hatchery cohort. The sample metadata's `cohort` column (when present)
must be exclusively the project's hatchery tag. The script aborts if
samples from other cohorts are present.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow `python scripts/03_*.py` from the module root without installing.
MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT / "src"))

from kbc.burden import build_table_a  # noqa: E402
from kbc.io import (  # noqa: E402
    load_gene_spans,
    load_genotypes_tsv,
    load_genotypes_vcf,
    load_inversion_intervals,
    load_karyotypes,
    load_roh_intervals,
    load_sample_metadata,
    load_scoring_weights,
    load_variant_master,
)
from kbc.schema import write_table_a  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-master", required=True, type=Path)
    parser.add_argument("--joint-vcf", type=Path, default=None,
                        help="Joint multisample VCF (mutually exclusive with --genotypes-tsv)")
    parser.add_argument("--genotypes-tsv", type=Path, default=None,
                        help="Long-form genotypes TSV (sample_id, variant_id, gt). For tests / no-cyvcf2 environments.")
    parser.add_argument("--karyotypes", required=True, type=Path)
    parser.add_argument("--inversion-intervals", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--scoring-weights", type=Path, default=None)
    parser.add_argument("--roh-dir", type=Path, default=None)
    parser.add_argument("--gff3", type=Path, default=None,
                        help="GFF3 annotation. Used for proper gene spans in the ROH-promotion check (§5 case 3). "
                             "When absent, gene span falls back to the bounding box of Tier-1 damaging variants.")
    parser.add_argument("--cohort-tag", type=str, default=None,
                        help="If set, asserts every sample's `cohort` column equals this tag.")
    args = parser.parse_args(argv)

    if (args.joint_vcf is None) == (args.genotypes_tsv is None):
        parser.error("exactly one of --joint-vcf and --genotypes-tsv must be provided")

    variant_master = load_variant_master(args.variant_master)
    karyotypes = load_karyotypes(args.karyotypes)
    intervals = load_inversion_intervals(args.inversion_intervals)
    metadata = load_sample_metadata(args.sample_metadata)
    weights = load_scoring_weights(args.scoring_weights)
    roh = load_roh_intervals(args.roh_dir)
    gene_spans = load_gene_spans(args.gff3)

    if args.cohort_tag is not None and "cohort" in metadata.columns:
        bad = metadata[metadata["cohort"] != args.cohort_tag]
        if not bad.empty:
            print(
                f"ERROR: cohort guard violated. {len(bad)} samples not in cohort='{args.cohort_tag}'",
                file=sys.stderr,
            )
            return 2

    if args.joint_vcf is not None:
        gts = load_genotypes_vcf(args.joint_vcf, restrict_variant_ids=set(variant_master["variant_id"].astype(str)))
    else:
        gts = load_genotypes_tsv(args.genotypes_tsv)

    table_a = build_table_a(
        variant_master=variant_master,
        genotypes=gts,
        karyotypes=karyotypes,
        inversion_intervals=intervals,
        sample_metadata=metadata,
        scoring_weights=weights,
        roh_intervals=roh,
        gene_spans=gene_spans,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_table_a(table_a, args.out)
    print(f"wrote {len(table_a)} rows → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
