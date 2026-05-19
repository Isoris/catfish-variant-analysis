"""Per-(sample × inversion × segment) aggregation → Table A.

Implements SPEC_KBC.md §3.A and §4 Step 2:

For each (sample, inversion, segment):
    - count Tier-1 damaging variants this sample carries  (n_damaging_variants)
    - per-gene classification → het_masked / hom_exposed / compound_het_unknown counts
    - burden_score_sum = sum(weight × allele_dose) over Tier-1 variants

The burden_score_sum weighting uses per-consequence weights from
scoring_weights.tsv when provided; otherwise weight = 1.0 (so the
score collapses to allele dose). The spec ("reused as-is; not
re-derived") leaves the weights table's schema underspecified — see
HANDOFF.md if a different shape is needed.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Mapping

import pandas as pd

from .gene_status import (
    STATUS_COMPOUND_HET_UNKNOWN,
    STATUS_HET_MASKED,
    STATUS_HOM_EXPOSED,
    STATUS_REFERENCE_LIKE,
    GeneVariants,
    classify_per_gene_status,
)
from .io import GT_HET, GT_HOM_ALT, GT_MISSING, GenotypeFrame
from .tier1 import is_tier1_lof


@dataclass
class TableARow:
    sample_id: str
    inversion_id: str
    pod_segment: str
    karyotype_class: str
    n_damaging_variants: int
    het_masked_gene_count: int
    hom_exposed_gene_count: int
    compound_het_unknown_gene_count: int
    burden_score_sum: float
    is_pruned_unrelated: bool
    roh_overlap_frac: float | None = None
    f_roh: float | None = None
    ancestry_group: str | None = None

    def as_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "inversion_id": self.inversion_id,
            "pod_segment": self.pod_segment,
            "karyotype_class": self.karyotype_class,
            "n_damaging_variants": self.n_damaging_variants,
            "het_masked_gene_count": self.het_masked_gene_count,
            "hom_exposed_gene_count": self.hom_exposed_gene_count,
            "compound_het_unknown_gene_count": self.compound_het_unknown_gene_count,
            "burden_score_sum": self.burden_score_sum,
            "is_pruned_unrelated": self.is_pruned_unrelated,
            "roh_overlap_frac": self.roh_overlap_frac,
            "f_roh": self.f_roh,
            "ancestry_group": self.ancestry_group,
        }


def _allele_dose(gt: str) -> int:
    if gt == GT_HOM_ALT:
        return 2
    if gt == GT_HET:
        return 1
    return 0


def _annotate_variants_with_segment(
    variant_master: pd.DataFrame,
    inversion_intervals: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join variants to (inversion_id, pod_segment) by (chrom, pos).

    A variant outside every inversion interval is dropped (out of scope
    for KBC). A variant inside multiple inversions or multiple segments
    of the same inversion is emitted once per overlap.
    """
    vm = variant_master.copy()
    vm["pos"] = vm["pos"].astype(int)
    annot_rows: list[dict] = []
    for inv in inversion_intervals.itertuples(index=False):
        hits = vm[
            (vm["chrom"] == inv.chrom)
            & (vm["pos"] >= inv.start)
            & (vm["pos"] < inv.end)
        ]
        for v in hits.itertuples(index=False):
            annot_rows.append({
                "variant_id": v.variant_id,
                "chrom": v.chrom,
                "pos": int(v.pos),
                "gene_id": v.gene_id,
                "consequence": v.consequence,
                "impact": getattr(v, "impact", None),
                "inversion_id": inv.inversion_id,
                "pod_segment": inv.pod_segment,
            })
    return pd.DataFrame(annot_rows)


def _genes_in_segment(
    annotated: pd.DataFrame,
    gene_spans: Mapping[str, tuple[str, int, int]] | None = None,
) -> dict[tuple[str, str, str], GeneVariants]:
    """Group Tier-1 damaging variants by (inversion_id, pod_segment, gene_id).

    When a GFF-derived gene span is available for the gene, use it; otherwise
    fall back to the bounding box of the gene's damaging variants in this
    segment (the MVP 1 approximation).
    """
    tier1 = annotated[annotated["consequence"].map(is_tier1_lof)]
    grouped: dict[tuple[str, str, str], GeneVariants] = {}
    for (inv, seg, gene), sub in tier1.groupby(["inversion_id", "pod_segment", "gene_id"]):
        if pd.isna(gene) or gene == "":
            continue
        gene_id = str(gene)
        chrom = str(sub["chrom"].iloc[0])
        gff = gene_spans.get(gene_id) if gene_spans else None
        if gff is not None and gff[0] == chrom:
            _, start, end = gff
            span_source = "gff"
        else:
            start = int(sub["pos"].min())
            end = int(sub["pos"].max()) + 1   # half-open from inclusive max-position
            span_source = "variant_bbox"
        grouped[(inv, seg, gene_id)] = GeneVariants(
            gene_id=gene_id,
            variant_ids=tuple(sub["variant_id"].astype(str).tolist()),
            chrom=chrom,
            gene_start=start,
            gene_end=end,
            span_source=span_source,
        )
    return grouped


def build_table_a(
    variant_master: pd.DataFrame,
    genotypes: GenotypeFrame,
    karyotypes: pd.DataFrame,
    inversion_intervals: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    scoring_weights: Mapping[str, float] | None = None,
    roh_intervals: Mapping[str, pd.DataFrame] | None = None,
    gene_spans: Mapping[str, tuple[str, int, int]] | None = None,
) -> pd.DataFrame:
    """Build Table A: one row per (sample, inversion, segment). Tier 1 only.

    `gene_spans` (gene_id → (chrom, start, end), 0-based half-open) is consulted
    for the ROH-overlap test in §5 case 3. When absent or when a gene_id is not
    in the map, the gene span falls back to the bounding box of Tier-1 damaging
    variants in that gene (see Modules/07_kbc/SPEC_KBC.md §5 note).
    """
    scoring_weights = scoring_weights or {}
    roh_intervals = roh_intervals or {}

    annotated = _annotate_variants_with_segment(variant_master, inversion_intervals)
    genes_by_seg = _genes_in_segment(annotated, gene_spans=gene_spans)

    # Tier 1 variants × their (inversion, segment) annotations
    tier1_var_df = annotated[annotated["consequence"].map(is_tier1_lof)]
    tier1_var_index: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    # (inversion_id, pod_segment) → list of (variant_id, consequence)
    for v in tier1_var_df.itertuples(index=False):
        tier1_var_index[(v.inversion_id, v.pod_segment)].append((v.variant_id, v.consequence))

    karyo_lookup = {
        (str(r.sample_id), str(r.inversion_id)): str(r.karyotype_class)
        for r in karyotypes.itertuples(index=False)
    }

    meta_lookup = {str(r.sample_id): r for r in sample_metadata.itertuples(index=False)}

    rows: list[TableARow] = []
    samples = sample_metadata["sample_id"].astype(str).tolist()
    inversions_segments = list({(inv, seg) for (inv, seg, _g) in genes_by_seg.keys()})
    # Also include (inv, seg) pairs that have no Tier-1 variants but still need a row per sample,
    # so the table is rectangular over the candidate set.
    inversions_segments = sorted(
        set(inversions_segments)
        | {(str(r.inversion_id), str(r.pod_segment)) for r in inversion_intervals.itertuples(index=False)}
    )

    for sample_id in samples:
        meta = meta_lookup.get(sample_id)
        is_pruned = bool(getattr(meta, "is_pruned_unrelated", False)) if meta is not None else False
        f_roh = getattr(meta, "f_roh", None) if meta is not None else None
        if pd.isna(f_roh):
            f_roh = None
        ancestry_group = getattr(meta, "ancestry_group", None) if meta is not None else None
        if isinstance(ancestry_group, float) and pd.isna(ancestry_group):
            ancestry_group = None

        sample_roh = roh_intervals.get(sample_id)

        for inv, seg in inversions_segments:
            karyo = karyo_lookup.get((sample_id, inv), "NA")

            # n_damaging_variants and burden_score_sum from Tier-1 variants
            n_dmg = 0
            burden = 0.0
            for v_id, conseq in tier1_var_index.get((inv, seg), []):
                gt = genotypes.lookup(sample_id, v_id)
                if gt == GT_MISSING:
                    continue
                dose = _allele_dose(gt)
                if dose == 0:
                    continue
                n_dmg += 1  # count the variant once per sample regardless of dose
                w = scoring_weights.get(str(conseq), 1.0)
                burden += w * dose

            # per-gene status counts
            het_masked = 0
            hom_exp = 0
            comp_het = 0
            for (inv2, seg2, _gene_id), gene in genes_by_seg.items():
                if inv2 != inv or seg2 != seg:
                    continue
                status = classify_per_gene_status(sample_id, gene, genotypes, sample_roh)
                if status == STATUS_HOM_EXPOSED:
                    hom_exp += 1
                elif status == STATUS_HET_MASKED:
                    het_masked += 1
                elif status == STATUS_COMPOUND_HET_UNKNOWN:
                    comp_het += 1
                # STATUS_REFERENCE_LIKE → no count

            rows.append(TableARow(
                sample_id=sample_id,
                inversion_id=inv,
                pod_segment=seg,
                karyotype_class=karyo,
                n_damaging_variants=n_dmg,
                het_masked_gene_count=het_masked,
                hom_exposed_gene_count=hom_exp,
                compound_het_unknown_gene_count=comp_het,
                burden_score_sum=burden,
                is_pruned_unrelated=is_pruned,
                roh_overlap_frac=None,    # MVP 1: not computed; needs ROH-vs-segment overlap calculation
                f_roh=f_roh if isinstance(f_roh, (int, float)) else None,
                ancestry_group=ancestry_group,
            ))

    return pd.DataFrame([r.as_dict() for r in rows])
