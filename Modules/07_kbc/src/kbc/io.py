"""I/O loaders for KBC MVP 1.

Reads:
- variant_master_scored.tsv  (subset of columns relevant for Tier 1)
- per-sample genotypes  (from a joint VCF via cyvcf2, or from a TSV used by tests)
- inversion karyotype table  (sample × inversion → AA/AB/BB/NA)
- inversion interval table   (inversion × POD segment → genomic interval)
- sample metadata
- scoring weights  (optional; consequence → weight)
- ROH BEDs per sample  (optional; used by gene_status case 3)
- GFF3 gene spans          (optional; supplies proper gene coordinates
                            for the ROH-overlap test in gene_status case 3.
                            When absent, gene_status falls back to the
                            bounding box of Tier-1 damaging variants in
                            the gene — a documented MVP 1 approximation.)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

GT_HOM_REF = "0/0"
GT_HET = "0/1"
GT_HOM_ALT = "1/1"
GT_MISSING = "./."


# --- variant master --------------------------------------------------------

VARIANT_MASTER_MIN_COLUMNS = (
    "variant_id", "chrom", "pos", "gene_id", "consequence", "impact",
)


def load_variant_master(path: str | Path) -> pd.DataFrame:
    """Load variant_master_scored.tsv. Only columns relevant to MVP 1 are required."""
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "variant_id": str, "gene_id": str})
    missing = [c for c in VARIANT_MASTER_MIN_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"variant_master missing required columns: {missing}")
    return df


# --- genotypes -------------------------------------------------------------

@dataclass
class GenotypeFrame:
    """Long-form genotype table: one row per (sample_id, variant_id, gt).

    `gt` is a string in {'0/0', '0/1', '1/1', './.'}. Multiallelic sites
    must be split upstream (asserted on load).
    """

    df: pd.DataFrame

    def __post_init__(self) -> None:
        for col in ("sample_id", "variant_id", "gt"):
            if col not in self.df.columns:
                raise ValueError(f"GenotypeFrame missing column: {col}")
        bad = set(self.df["gt"].unique()) - {GT_HOM_REF, GT_HET, GT_HOM_ALT, GT_MISSING}
        if bad:
            raise ValueError(f"unrecognised GT strings: {sorted(bad)} (multiallelics must be split upstream)")

    def lookup(self, sample_id: str, variant_id: str) -> str:
        rows = self.df[(self.df["sample_id"] == sample_id) & (self.df["variant_id"] == variant_id)]
        if rows.empty:
            return GT_MISSING
        return str(rows["gt"].iloc[0])

    def samples(self) -> list[str]:
        return sorted(self.df["sample_id"].unique().tolist())


def load_genotypes_tsv(path: str | Path) -> GenotypeFrame:
    """Load a long-form TSV with columns sample_id, variant_id, gt. Used by tests."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    return GenotypeFrame(df=df)


def load_genotypes_vcf(path: str | Path, restrict_variant_ids: Iterable[str] | None = None) -> GenotypeFrame:
    """Load genotypes from a (joint) VCF/BCF via cyvcf2.

    cyvcf2 is a heavy dependency; we import it lazily so test environments
    that don't ship it can still exercise the rest of the pipeline through
    `load_genotypes_tsv`.
    """
    try:
        from cyvcf2 import VCF  # type: ignore
    except ImportError as exc:  # pragma: no cover — exercised in production env only
        raise RuntimeError(
            "cyvcf2 is required to load VCF input. Install via envs/kbc.yaml or "
            "use load_genotypes_tsv for synthetic data."
        ) from exc

    keep = set(restrict_variant_ids) if restrict_variant_ids is not None else None
    vcf = VCF(str(path), gts012=False)
    samples = list(vcf.samples)
    rows: list[dict[str, str]] = []
    for rec in vcf:
        vid = f"{rec.CHROM}:{rec.POS}:{rec.REF}:{','.join(rec.ALT)}"
        if keep is not None and vid not in keep:
            continue
        if len(rec.ALT) != 1:
            raise ValueError(f"multiallelic site at {vid}; split upstream before loading")
        for i, sample in enumerate(samples):
            gt = rec.genotypes[i]  # [a1, a2, phased]
            if gt[0] < 0 or gt[1] < 0:
                gtstr = GT_MISSING
            else:
                a1, a2 = sorted((gt[0], gt[1]))
                gtstr = f"{a1}/{a2}"
            rows.append({"sample_id": sample, "variant_id": vid, "gt": gtstr})
    return GenotypeFrame(df=pd.DataFrame(rows))


# --- inversion karyotype table --------------------------------------------

KARYOTYPE_VALUES = ("AA", "AB", "BB", "NA")


def load_karyotypes(path: str | Path) -> pd.DataFrame:
    """Long-form: sample_id, inversion_id, karyotype_class."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ("sample_id", "inversion_id", "karyotype_class"):
        if col not in df.columns:
            raise ValueError(f"karyotype table missing column: {col}")
    bad = set(df["karyotype_class"].unique()) - set(KARYOTYPE_VALUES)
    if bad:
        raise ValueError(f"karyotype_class contains unexpected values: {sorted(bad)}")
    return df


# --- inversion interval table ---------------------------------------------

POD_SEGMENTS = ("L", "M", "R", "whole")


def load_inversion_intervals(path: str | Path) -> pd.DataFrame:
    """Long-form: inversion_id, pod_segment, chrom, start, end.

    Single-crossover candidates emit one row with pod_segment='whole'.
    Double-crossover candidates emit three rows L/M/R.
    """
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "inversion_id": str, "pod_segment": str})
    required = ("inversion_id", "pod_segment", "chrom", "start", "end")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"inversion intervals missing columns: {missing}")
    bad = set(df["pod_segment"].unique()) - set(POD_SEGMENTS)
    if bad:
        raise ValueError(f"pod_segment contains unexpected values: {sorted(bad)}")
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


# --- sample metadata -------------------------------------------------------

def load_sample_metadata(path: str | Path) -> pd.DataFrame:
    """sample_id, cohort, batch, is_pruned_unrelated (and optionally f_roh, ancestry_group)."""
    df = pd.read_csv(path, sep="\t", dtype={"sample_id": str})
    if "sample_id" not in df.columns:
        raise ValueError("sample metadata missing sample_id")
    if "is_pruned_unrelated" in df.columns:
        df["is_pruned_unrelated"] = df["is_pruned_unrelated"].astype(bool)
    return df


# --- scoring weights -------------------------------------------------------

def load_scoring_weights(path: str | Path | None) -> dict[str, float]:
    """Map consequence → weight. Returns empty dict if path is None or file missing.

    Expected schema: TSV with columns `consequence` and `weight`. Anything
    else will be filed as 'unsupported schema' in MVP 1 and the loader
    will return an empty dict (i.e. burden_score_sum falls back to
    allele-dose count).
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    df = pd.read_csv(p, sep="\t")
    if not {"consequence", "weight"}.issubset(df.columns):
        return {}
    return {str(r.consequence): float(r.weight) for r in df.itertuples(index=False)}


# --- GFF3 gene spans -------------------------------------------------------

# Attribute keys (in order of preference) used to resolve the gene_id from a
# GFF3 row's 9th column. The first key that yields a non-empty value wins.
# When the value is prefixed with the GFF Sequence-Ontology prefix `gene-`
# (Ensembl/NCBI convention), the prefix is stripped to recover the bare ID.
GFF3_GENE_ID_KEYS = ("gene_id", "ID", "Name", "gene_name")


def _parse_gff3_attributes(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for tok in field.strip().split(";"):
        tok = tok.strip()
        if not tok or "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _resolve_gene_id(attrs: dict[str, str]) -> str | None:
    for k in GFF3_GENE_ID_KEYS:
        v = attrs.get(k)
        if v:
            return v[len("gene-"):] if v.startswith("gene-") else v
    return None


def load_gene_spans(path: str | Path | None) -> dict[str, tuple[str, int, int]]:
    """Parse a GFF3 file and return gene_id → (chrom, start, end).

    Looks at lines where the type (column 3) is exactly `gene`. Both 1-based
    GFF coordinates are converted to 0-based half-open (start-1, end).

    Returns an empty dict if `path` is None or the file is missing.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    spans: dict[str, tuple[str, int, int]] = {}
    with p.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "gene":
                continue
            chrom = cols[0]
            try:
                start1 = int(cols[3])
                end1 = int(cols[4])
            except ValueError:
                continue
            attrs = _parse_gff3_attributes(cols[8])
            gene_id = _resolve_gene_id(attrs)
            if gene_id is None:
                continue
            spans[gene_id] = (chrom, start1 - 1, end1)
    return spans


# --- ROH BEDs --------------------------------------------------------------

def load_roh_intervals(roh_dir: str | Path | None) -> dict[str, pd.DataFrame]:
    """Map sample_id → DataFrame[chrom, start, end] of ROH intervals.

    Looks for {roh_dir}/{sample_id}.bed (3-column BED). Returns empty
    dict if roh_dir is None or doesn't exist.
    """
    if roh_dir is None:
        return {}
    d = Path(roh_dir)
    if not d.exists():
        return {}
    out: dict[str, pd.DataFrame] = {}
    for bed in d.glob("*.bed"):
        sample_id = bed.stem
        df = pd.read_csv(bed, sep="\t", header=None, names=["chrom", "start", "end"], dtype={"chrom": str})
        out[sample_id] = df
    return out
