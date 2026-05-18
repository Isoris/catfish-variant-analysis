# 07_kbc — Karyotype Burden Contrast

**Status:** SPEC ONLY — awaiting audit before implementation.
**Manuscript role:** headline POD-compatible result for the v20 inversions
manuscript (`MS_Inversions_North_african_catfish`).

KBC is a cohort-level burden contrast across the three karyotype classes
(AA / AB / BB) of each candidate inversion, segment by segment
(L / M / R for double-crossover PODs, whole for single-crossover). The
headline statistic is the per-sample hom-exposed gene count for
high-confidence loss-of-function variants (Tier 1: stop-gained,
frameshift, start-lost, canonical splice donor/acceptor); AB is the
masking reference.

## Documents in this module

| File | Purpose |
|---|---|
| `SPEC_KBC.md` | the full spec (v0.8) |
| `HANDOFF.md` | implementation handoff, including the four open questions to resolve with the project owner before coding |
| `MANUSCRIPT_PARAGRAPH.md` | manuscript-ready prose (methods + results + abstract versions) |
| `UMBRELLA_README.md` | how KBC sits beside its sibling specs (HPP, HAPS) |

## Layout

```text
Modules/07_kbc/
├── SPEC_KBC.md
├── HANDOFF.md
├── MANUSCRIPT_PARAGRAPH.md
├── UMBRELLA_README.md
├── config/                empty — populated during MVP 1
├── schemas/                JSON-Schema 2020-12 contracts for tables A–F
├── scripts/                empty — populated during MVP 1
├── src/kbc/                empty — populated during MVP 1
├── tests/fixtures/         empty — populated during MVP 1
└── outputs/                local dev only; canonical results land under
                            ${BASE}/results/catfish-variant-analysis/07_kbc/
```

## Hard rules (excerpted from `HANDOFF.md`)

1. Three-cohort rule absolute: 226-sample pure *C. gariepinus* hatchery only.
2. Interpretation is POD-*compatible*, never POD-found. No claim of
   balancing selection, true overdominance, or fitness.
3. HWE deviation is never used as evidence.
4. Each inversion is one independent test unit. No pooled raw-count test
   across inversions.
5. Synonymous variants are a parallel negative control, not a denominator.
6. Manuscript headline uses Tier 1 (high-confidence LoF) only.
7. Tier 2 stays disabled until the splice-validation six-check gate
   passes (`SPEC_KBC.md` §1.9).
8. Family hubs are confounds, not units of analysis.
9. Per-variant het count inside inversions is not a burden statistic;
   the headline is per-gene hom-exposed count.

See `HANDOFF.md` for the four open questions to resolve before MVP 1.

## Sibling specs

- **HPP** (Haplotype Projection from Pedigree): individual-level
  inheritance projection from ngsPedigree Stage 3. Lives in the
  `ngsPedigree` repo. Blocking dependency: Stage 3 readiness.
- **HAPS** (Haplotype-Aware Protein Scoring): bcftools/csq + VESM
  rescore on haplotype-reconstructed proteins. Lives in
  `MODULE_CONSERVATION`. Deferred until after the manuscript.
