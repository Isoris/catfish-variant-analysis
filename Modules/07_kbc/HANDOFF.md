# HANDOFF — KBC implementation into `catfish-variant-analysis`

**For:** Claude Code session in the `catfish-variant-analysis` repo.
**You are receiving:** `KBC_SPEC.md` (the spec), `MANUSCRIPT_PARAGRAPH.md` (manuscript-ready prose), `UMBRELLA_README.md` (where KBC sits in the four-sibling stack).
**Status:** SPEC ONLY — **awaiting audit before implementation.**
**Manuscript role:** **headline POD result for the v20 inversions manuscript.**

---

## What KBC is, in one paragraph

KBC (Karyotype Burden Contrast) is the cohort-level burden test for the inversions manuscript. For each candidate inversion, it stratifies the 226-sample pure *C. gariepinus* hatchery cohort by karyotype class (AA / AB / BB) and POD segment (L / M / R for double-crossover candidates; whole interval for single-crossover), then tests whether AA and BB homozygotes carry more hom-exposed deleterious genes than AB heterozygotes. AB is the masking reference. Synonymous variants run in parallel as a negative control. Each inversion is one independent test unit — there is no pooled raw-count test across inversions. The headline statistic uses **Tier 1 high-confidence loss-of-function variants only** (stop-gained, frameshift, start-lost, canonical splice donor/acceptor); Tier 2 (validated splice Class A) is a robustness check; Tier 3 (strong missense) is exploratory.

## Where KBC lives in the repo

Proposed location:

```
catfish-variant-analysis/Modules/NN_kbc/
  SPEC_KBC.md
  config/
  schemas/
  scripts/
  src/kbc/
  tests/
  outputs/             ← .gitignored
```

Numbering: pick the next free `NN_` slot in the repo's `Modules/` directory. KBC is a new module; it does not replace anything.

Results land under `${BASE}/results/catfish-variant-analysis/NN_kbc/` per the repo's results-layout convention.

## Inputs KBC consumes (and where they are)

All required inputs already exist in the project:

| Input | Source | Where |
|---|---|---|
| `variant_master_scored.tsv` | MODULE_CONSERVATION STEP 16 | already produced; canonical CSQ + SIFT + VESM + splice + scoring |
| Joint VCF | MODULE_CONSERVATION STEP 03 | already produced |
| Inversion karyotype table | PCAngsd K=3 + Hungarian arrangement-label matching | inversion atlas |
| Inversion intervals + L/M/R POD partition | inversion atlas | already defined |
| Sample metadata (incl. NAToRA-pruned flag) | project sample sheet | existing |
| Scoring weights | MODULE_CONSERVATION `scoring_weights.tsv` | existing |
| ngsRelate kinship matrix | existing 23-col `.res` | for the mixed-model secondary mode |

Optional inputs: ancestry Q matrix (NGSadmix), ROH BED per sample (MODULE_3), local ancestry (Engine B), callable mask.

**KBC does not depend on:** GERP, phastCons, phyloP, orthology-tier, CAFE, family-evolution scores. The EGO pipeline is decoupled.

## Hard rules — DO NOT violate

1. **Three-cohort rule absolute.** KBC operates only on the 226-sample pure *C. gariepinus* hatchery cohort. F1 hybrid and *C. macrocephalus* wild cohorts are not inputs.
2. **Headline interpretation is POD-*compatible*, never POD-found.** No claim of balancing selection, true overdominance, or fitness.
3. **HWE deviation is not used as evidence** for anything (mechanically inevitable in a hatchery cohort).
4. **Each inversion is one test unit.** No pooled raw-count test across inversions. Across-inversion summaries are descriptive (signal-class counts) or meta-analytic (inverse-variance weighted), never sum-of-raw-counts. See §1.6 of the spec.
5. **Synonymous variants are a parallel test, not a denominator.** No dN/dS-style ratio as the headline statistic. See §1.7.
6. **Headline uses Tier 1 (high-confidence LoF) only.** Tier 2 and Tier 3 are robustness/exploratory diagnostics. See §1.8.
7. **Tier 2 disabled until splice validation passes** the six-check gate in §1.9.
8. **Family hubs are confounds, not units of analysis.** Relatedness is handled by NAToRA-pruning (primary) + kinship-corrected GLM (secondary). No `family_transmission_summary` table — that's HPP's job (the sibling spec).
9. **Per-variant het count inside inversions is NOT a burden statistic.** AB carries more het variants by construction; this is karyotype, not load. The headline is **per-gene hom-exposed count**. See §1.3 and the §6.2 explanation in the spec.

## The four open questions to settle BEFORE coding

1. **The `failed_synonymous_control` cutoff** (currently `|syn delta| > 0.5 × |damaging delta|`). Calibrated against what? Likely needs an empirical distribution from a small pilot run. The spec flags this for audit.
2. **Kinship matrix source** — use the ngsRelate kinship matrix directly, or the NAToRA-derived structure? Project precedent should win.
3. **Splice validation checklist (§1.9) status** — has any of the six checks been run against the existing splice module? If not, Tier 2 stays disabled. Manuscript can ship with Tier 1 + Tier 3 only.
4. **PCAngsd K=3 posterior cutoff** for `karyotype_class = NA`. The spec drops low-posterior samples to NA but the threshold is project-determined.

## Staged MVP

Build in this order, **stop after each stage and audit before continuing**:

- **MVP 1** — Table A (per-sample, per-inversion, per-segment hom-exposed and het-masked counts at Tier 1 only). Sanity check against existing per-sample burden tables.
- **MVP 2** — Table B (variant-to-arrangement assignment from karyotype distribution). Cross-check against the inversion atlas's per-arrangement variant frequency tracks for consistency.
- **MVP 3** — Table C + Table D + Table E (per-karyotype summary, synonymous control, per-inversion report). Run for Tier 1 first. This is the manuscript headline.
- **MVP 4** — Add Tier 2 if and only if splice validation has passed. Add Tier 3 as exploratory.
- **MVP 5** — Mixed-model robustness + L vs M secondary test for double-CO PODs.
- **MVP 6** — Diagnostic figures + Table F across-inversion summary.

## What to deliver back

After MVP 3, the manuscript paragraph in `MANUSCRIPT_PARAGRAPH.md` should drop in directly. The bracketed placeholders ([N], [X], [Y], etc.) map to specific columns in Table E and Table F. The paragraph writes itself once the numbers exist.

## What NOT to do

- Do not extend KBC to multi-cohort comparisons.
- Do not start with HPP integration — HPP is a sibling spec, ships separately (it lives in the `ngsPedigree` repo, see `UMBRELLA_README.md`).
- Do not start with HAPS — HAPS is deferred until after the manuscript.
- Do not add features beyond what the spec specifies without a separate audit chat.
- Do not invent variant-impact thresholds beyond what `variant_master_scored.tsv` already provides — those decisions are MODULE_CONSERVATION's, not KBC's.

## First message I would send

> Read `KBC_SPEC.md` in full before writing any code. Then read the four open questions in this handoff. Resolve them with the project owner (Quentin) by asking — do not guess defaults. After that, build MVP 1 (Table A only, Tier 1 only) and stop for review.
