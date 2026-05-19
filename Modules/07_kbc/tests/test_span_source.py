"""Verify the GFF vs variant_bbox fallback paths in burden._genes_in_segment."""

from __future__ import annotations

import pandas as pd

from kbc.burden import _annotate_variants_with_segment, _genes_in_segment
from kbc.io import load_gene_spans, load_inversion_intervals, load_variant_master


def test_span_source_is_gff_when_spans_provided(variant_master_path, inversions_path, gff_path):
    vm = load_variant_master(variant_master_path)
    invs = load_inversion_intervals(inversions_path)
    spans = load_gene_spans(gff_path)
    annotated = _annotate_variants_with_segment(vm, invs)
    genes = _genes_in_segment(annotated, gene_spans=spans)

    g1 = next(g for (inv, seg, gid), g in genes.items() if gid == "G1")
    assert g1.span_source == "gff"
    assert (g1.gene_start, g1.gene_end) == (1000, 3000)


def test_span_source_falls_back_to_variant_bbox_when_no_spans(variant_master_path, inversions_path):
    vm = load_variant_master(variant_master_path)
    invs = load_inversion_intervals(inversions_path)
    annotated = _annotate_variants_with_segment(vm, invs)
    genes = _genes_in_segment(annotated, gene_spans=None)

    g1 = next(g for (inv, seg, gid), g in genes.items() if gid == "G1")
    assert g1.span_source == "variant_bbox"
    # G1 has v1@1500 and v2@2500 → bbox [1500, 2501) (half-open from max+1)
    assert (g1.gene_start, g1.gene_end) == (1500, 2501)


def test_span_source_falls_back_when_gene_not_in_gff(variant_master_path, inversions_path):
    """A gene_id missing from the GFF dict falls back to bbox even if other genes are GFF."""
    vm = load_variant_master(variant_master_path)
    invs = load_inversion_intervals(inversions_path)
    annotated = _annotate_variants_with_segment(vm, invs)
    partial_spans = {"G1": ("chr1", 999, 3001)}    # only G1
    genes = _genes_in_segment(annotated, gene_spans=partial_spans)

    g1 = next(g for (inv, seg, gid), g in genes.items() if gid == "G1")
    assert g1.span_source == "gff"
    g3 = next(g for (inv, seg, gid), g in genes.items() if gid == "G3")
    assert g3.span_source == "variant_bbox"


def test_span_source_falls_back_when_gff_chrom_mismatches(variant_master_path, inversions_path):
    """A genuine paralog/duplicate gene_id on the wrong chromosome must not be applied."""
    vm = load_variant_master(variant_master_path)
    invs = load_inversion_intervals(inversions_path)
    annotated = _annotate_variants_with_segment(vm, invs)
    spans = {"G1": ("chr99", 0, 1000)}    # wrong chromosome
    genes = _genes_in_segment(annotated, gene_spans=spans)

    g1 = next(g for (inv, seg, gid), g in genes.items() if gid == "G1")
    assert g1.span_source == "variant_bbox"
