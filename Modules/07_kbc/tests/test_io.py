"""Loader tests, with a focus on the GFF3 gene-span loader (B5)."""

from __future__ import annotations

from kbc.io import load_gene_spans


def test_gff3_loads_six_genes(gff_path):
    spans = load_gene_spans(gff_path)
    assert set(spans.keys()) == {"G1", "G2", "G3", "G4", "G5", "G6"}


def test_gff3_coordinates_are_0_based_half_open(gff_path):
    spans = load_gene_spans(gff_path)
    # GFF row "chr1 mini gene 1001 3000" → 0-based half-open [1000, 3000)
    assert spans["G1"] == ("chr1", 1000, 3000)
    # G3 row: "chr1 mini gene 3801 4200" → [3800, 4200)
    assert spans["G3"] == ("chr1", 3800, 4200)


def test_gff3_loader_returns_empty_for_missing_path(tmp_path):
    spans = load_gene_spans(tmp_path / "does_not_exist.gff3")
    assert spans == {}


def test_gff3_loader_returns_empty_for_none():
    spans = load_gene_spans(None)
    assert spans == {}


def test_gff3_loader_skips_comment_and_non_gene_rows(tmp_path):
    p = tmp_path / "mixed.gff3"
    p.write_text(
        "##gff-version 3\n"
        "# arbitrary header\n"
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=A\n"
        "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=mRNA-A\n"
        "chr1\tx\texon\t1\t50\t.\t+\t.\tID=exon-A1\n"
        "chr1\tx\tgene\t200\t300\t.\t-\t.\tID=B\n"
    )
    spans = load_gene_spans(p)
    assert spans == {"A": ("chr1", 0, 100), "B": ("chr1", 199, 300)}


def test_gff3_loader_strips_gene_prefix(tmp_path):
    p = tmp_path / "ens_style.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=gene-MYGENE;Name=MYGENE\n"
    )
    spans = load_gene_spans(p)
    assert spans == {"MYGENE": ("chr1", 0, 100)}
