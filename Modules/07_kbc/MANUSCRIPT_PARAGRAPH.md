# Closing statistics paragraph — POD-compatible burden test

For the final manuscript paragraph (methods / results boundary). Two
versions: a methods-section version and a shorter results-section
version. Choose one based on where it fits best in the manuscript
structure.

---

## Version A — methods-section paragraph (recommended, ~320 words)

Each candidate inversion was treated as an independent analysis unit.
Within each inversion, samples were stratified by karyotype class
(AA, AB, BB) derived from PCAngsd K=3 analysis with Hungarian
arrangement-label matching, and partitioned by POD segment (L, M, R
for double-crossover candidates; whole interval for single-crossover
candidates). The headline burden statistic was the per-sample
hom-exposed gene count for **high-confidence loss-of-function (LoF)
variants — stop-gained, frameshift, start-lost, and canonical splice
donor / acceptor variants (Tier 1)** — for which the LoF claim is
mechanistically defensible from sequence annotation alone. Two
robustness tiers were run in parallel: Tier 2 added Class A splice
subclasses from the splice annotation module (enabled only after the
module passed a six-check validation against canonical site
identification, GFF coordinate consistency, strand handling,
boundary classification, SnpEff concordance, and isoform stability);
Tier 3 added SIFT4G-deleterious and VESM-strong missense as
exploratory variants. For each (inversion × segment × tier) cell,
per-sample hom-exposed gene counts were compared between karyotype
classes with AB as the masking-reference class, motivated by the
gene-level complementation expected in heterozygotes that carry one
copy of each arrangement. The headline contrasts were AA vs AB and
BB vs AB on hom-exposed gene count, tested by Mann-Whitney U on the
NAToRA-pruned subset of 81 unrelated individuals as the primary mode,
with a secondary negative-binomial generalized linear mixed model on
the full 226-sample cohort using the ngsRelate kinship matrix as a
random effect to confirm robustness to cryptic relatedness. The
identical contrast structure was run in parallel on synonymous
variants as a within-cohort negative control for arrangement
divergence, callable-region effects, and variant-density bias;
inversions in which the synonymous-variant contrast tracked the
damaging contrast in direction and magnitude were classified as
failed-synonymous-control and reported but not used as evidence. No
ratio-form normalisation (such as dN/dS) was applied as a headline
statistic; synonymous variants entered the analysis as a parallel
test rather than as a denominator. Per-inversion p-values were
Benjamini-Hochberg-adjusted across the candidate set at q = 0.05.
**The manuscript-headline POD-compatibility classification of each
inversion was made on Tier 1 (high-confidence LoF) only**; Tier 2 and
Tier 3 results are reported as concordance diagnostics — the
tier-stability column in Supplementary Table SX records whether the
three tiers agree, with the strongest claim reserved for inversions
that are POD-compatible at all three. For double-crossover inversions,
a secondary contrast tested whether the central recombinant island
(M) showed a modified burden architecture relative to the flanking
inverted segments (L, R) within each karyotype class. Across-inversion
summaries report the number of candidate inversions falling in each
signal class (POD-compatible-both, AA-exposed-only, BB-exposed-only,
no-signal, failed-synonymous-control, inconclusive-small-n); a
random-effects meta-analytic pooled effect size on the well-behaved
subset of inversions is reported as a descriptive summary, with
inverse-variance weighting and never as a sum of raw counts. The
reported pattern is structural and mutational — consistent with
pseudo-overdominance through gene-level complementation of
arrangement-private recessive load — and, given a single-timepoint
hatchery cohort without fitness measurement, does not demonstrate
balancing selection, true overdominance, or fitness consequences.

---

## Version B — results-section closing paragraph (~190 words)

Across the [N] candidate inversions tested at the high-confidence
LoF tier, [X] showed the POD-compatible pattern on both arrangement
backgrounds (AA and BB homozygotes both exposed elevated hom-exposed
LoF-gene counts relative to AB heterozygotes, with the
synonymous-variant negative control flat); [Y] showed an asymmetric
pattern in which only one arrangement carried significantly elevated
hom-exposed LoF load; [Z] failed the synonymous-control gate and were
not used as evidence; [W] showed no significant signal; and [V] were
classified as inconclusive due to insufficient karyotype-class sample
sizes. The random-effects pooled effect across the well-behaved
subset of inversions was [Δ AA-vs-AB = X.X (95% CI x.x – x.x)] and
[Δ BB-vs-AB = Y.Y (95% CI y.y – y.y)] on hom-exposed LoF gene count,
with between-inversion heterogeneity [I² = Z%]. The same contrast
on the broader Tier 3 set (LoF + splice + strong missense) agreed in
direction at [k]/[N] inversions; [m] inversions reached significance
only at Tier 3 and are reported with a model-dependence caveat. All
headline results were confirmed by the full-cohort kinship-corrected
mixed-model robustness analysis. For double-crossover inversions,
[j of n] showed a significantly different burden architecture in the
central recombinant island relative to the flanking inverted
segments, consistent with the expected loss of arrangement-haplotype
identity through double-crossover events. These patterns are
POD-compatible structural-mutational signatures and are not, in the
absence of direct fitness measurement, evidence of balancing
selection or true overdominance.

---

## Version C — single tight sentence for the abstract (~55 words)

In a karyotype-stratified burden contrast on the 226-sample hatchery
cohort, AA and BB homozygotes carried significantly more
hom-exposed high-confidence loss-of-function genes than AB
heterozygotes at [X] of [N] candidate inversions, with the same
contrast on synonymous variants flat — a POD-compatible
structural-mutational pattern that does not demonstrate balancing
selection or fitness.

---

## Notes for filling in the bracketed numbers

The bracketed placeholders in Version B and C are the only numbers
that need to be filled in after running KBC. Specifically:

- `N` = total candidate inversions (table E row count)
- `X`, `Y`, `Z`, `W`, `V` = counts in each `kbc_signal_class`
  (table F sheet 1)
- pooled effect sizes + 95% CI + I² come from table F sheet 2
- `k`, `n` for the L-vs-M contrast come from `dco_segment_test_p`
  filtered at FDR q = 0.05

All five of these can be auto-generated from KBC's output tables once
the pipeline runs; the manuscript paragraph is a template waiting for
the numbers, not a manual writing exercise.

---

## What I would NOT put in this paragraph

- Anything about overdominance proven, balancing selection
  established, fitness measured, or POD-found. Always
  POD-*compatible*, never POD-found.
- Anything about HWE deviation as evidence (mechanically inevitable
  in a structured hatchery cohort; carries no selective information).
- dN/dS, ω, or substitution-rate language. Synonymous variants are a
  parallel test, not a rate.
- Family-hub-aggregated results. Hubs are a relatedness confound,
  not a unit of summary.
- Heterozygous variant counts inside inversions as a burden
  statistic. AB carries more het variants by construction; the
  gene-level hom-exposed count is the headline.
- Cross-cohort comparisons. The three-cohort rule applies: this is
  the 226-sample pure *C. gariepinus* hatchery cohort only.

---

*Drafted as part of KBC SPEC v0.7. Numbers to be filled in after the
pipeline runs.*
