"""Per-(sample, gene, segment) status classification.

Implements SPEC_KBC.md §5 case 1–4:

    case 1: all genotypes 0/0                          → reference_like
    case 2: at least one variant 1/1                   → hom_exposed
    case 3: exactly one variant 0/1, no 1/1
            if ROH(sample) covers the gene             → hom_exposed   (ROH-promotion)
            else                                       → het_masked
    case 4: ≥ 2 variants 0/1, no 1/1                   → compound_het_unknown

Missing genotypes (./.) carry no information and are dropped from the
per-gene classification; a gene with all-missing genotypes for a sample
is reported as reference_like (no observation = no exposure).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd

from .io import GT_HET, GT_HOM_ALT, GT_HOM_REF, GT_MISSING, GenotypeFrame

PerGeneStatus = str  # one of: reference_like / hom_exposed / het_masked / compound_het_unknown

STATUS_REFERENCE_LIKE = "reference_like"
STATUS_HOM_EXPOSED = "hom_exposed"
STATUS_HET_MASKED = "het_masked"
STATUS_COMPOUND_HET_UNKNOWN = "compound_het_unknown"


@dataclass
class GeneVariants:
    """The damaging variants in one gene × segment, in the format the classifier needs."""

    gene_id: str
    variant_ids: tuple[str, ...]
    chrom: str
    # ROH overlap is queried by chrom + gene span; gene span is approximated as
    # the min/max position of damaging variants in the gene for MVP 1 (a gene's
    # full coordinate span would need GFF; not required for ROH-promotion since
    # an ROH covering any damaging variant in the gene implies hom-coverage
    # for that variant locus).
    gene_min_pos: int
    gene_max_pos: int


def gene_overlaps_any_roh(gene: GeneVariants, roh_df: pd.DataFrame | None) -> bool:
    """True if any ROH interval for this sample overlaps the gene span."""
    if roh_df is None or roh_df.empty:
        return False
    hits = roh_df[
        (roh_df["chrom"] == gene.chrom)
        & (roh_df["start"] < gene.gene_max_pos)
        & (roh_df["end"] > gene.gene_min_pos)
    ]
    return not hits.empty


def classify_per_gene_status(
    sample_id: str,
    gene: GeneVariants,
    genotypes: GenotypeFrame,
    roh_df: pd.DataFrame | None = None,
) -> PerGeneStatus:
    """Apply §5 case-by-case classification for one (sample, gene, segment)."""
    gts: list[str] = [genotypes.lookup(sample_id, v) for v in gene.variant_ids]
    observed = [g for g in gts if g != GT_MISSING]

    if not observed:
        return STATUS_REFERENCE_LIKE

    if any(g == GT_HOM_ALT for g in observed):
        return STATUS_HOM_EXPOSED

    n_het = sum(1 for g in observed if g == GT_HET)
    if n_het == 0:
        return STATUS_REFERENCE_LIKE

    if n_het == 1:
        if gene_overlaps_any_roh(gene, roh_df):
            return STATUS_HOM_EXPOSED  # ROH-promotion
        return STATUS_HET_MASKED

    return STATUS_COMPOUND_HET_UNKNOWN
