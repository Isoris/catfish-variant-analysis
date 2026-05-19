#!/usr/bin/env bash
# Modules/07_kbc/config/kbc_config.sh
#
# Configuration for the KBC (Karyotype Burden Contrast) module.
# All paths are resolved relative to ${BASE} (the project root on LANTA).
# Override any of these by exporting before sourcing.
#
# Status: MVP 1 — Table A (per-sample × inversion × POD-segment), Tier 1 only.
# Later MVPs will populate the remaining tables (B–F) and the secondary
# kinship-corrected GLM mode.

set -euo pipefail

# --- Project root -----------------------------------------------------------
: "${BASE:?BASE must be set to the project root}"

# --- KBC module location ---------------------------------------------------
: "${KBC_DIR:=${BASE}/Modules/07_kbc}"
: "${KBC_RESULTS:=${BASE}/results/catfish-variant-analysis/07_kbc}"

mkdir -p "${KBC_RESULTS}"

# --- Inputs (canonical project locations) ----------------------------------
: "${VARIANT_MASTER:=${BASE}/results/MODULE_CONSERVATION/step16/variant_master_scored.tsv}"
: "${JOINT_VCF:=${BASE}/results/MODULE_CONSERVATION/step03/joint.vcf.gz}"
: "${GFF3:=${BASE}/00-samples/fClaHyb_Gar_LG.gff3}"
: "${INVERSION_KARYOTYPES:=${BASE}/results/inversion_atlas/karyotypes_PCAngsd_K3_hungarian.tsv}"
: "${INVERSION_INTERVALS:=${BASE}/results/inversion_atlas/intervals_with_pod_partition.tsv}"
: "${SAMPLE_METADATA:=${BASE}/01_inputs_check/samples.tsv}"
: "${SCORING_WEIGHTS:=${BASE}/results/MODULE_CONSERVATION/scoring_weights.tsv}"

# --- Optional inputs --------------------------------------------------------
: "${ROH_BED_DIR:=${BASE}/results/MODULE_3/roh_per_sample}"   # one BED per sample
: "${KINSHIP_MATRIX:=${BASE}/results/ngsRelate/kinship.res}"   # MVP 5+

# --- Outputs ---------------------------------------------------------------
: "${TABLE_A_OUT:=${KBC_RESULTS}/A_kbc_sample_inversion_burden.tsv}"
: "${TABLE_B_OUT:=${KBC_RESULTS}/B_kbc_variant_arrangement_assignments.tsv}"
: "${TABLE_C_OUT:=${KBC_RESULTS}/C_kbc_karyotype_burden_summary.tsv}"
: "${TABLE_D_OUT:=${KBC_RESULTS}/D_kbc_synonymous_control_summary.tsv}"
: "${TABLE_E_OUT:=${KBC_RESULTS}/E_kbc_per_inversion_report.tsv}"
: "${TABLE_F_SHEET1_OUT:=${KBC_RESULTS}/F_kbc_across_inversion_class_counts.tsv}"
: "${TABLE_F_SHEET2_OUT:=${KBC_RESULTS}/F_kbc_across_inversion_meta_analysis.tsv}"

# --- Cohort guard ----------------------------------------------------------
# KBC is locked to the 226-sample pure C. gariepinus hatchery cohort.
# Any code path that ingests samples must assert this.
: "${KBC_COHORT_TAG:=gariepinus_hatchery_226}"

# --- Tier knobs ------------------------------------------------------------
# Tier 1 is the manuscript headline; always enabled.
# Tier 2 stays disabled until the splice-validation six-check gate (§1.9)
# passes. The flag here is the operational kill-switch.
: "${KBC_TIER1_ENABLED:=1}"
: "${KBC_TIER2_ENABLED:=0}"
: "${KBC_TIER3_ENABLED:=0}"

# --- Statistical knobs (placeholder defaults; awaiting audit) --------------
# These are NOT calibrated — they exist so MVP 1 has named values to refer
# back to. Each corresponds to an open question in HANDOFF.md.
: "${KBC_SYN_CONTROL_FRACTION_CUTOFF:=0.5}"       # see HANDOFF.md open Q1
: "${KBC_MIN_KARYOTYPE_CELL_SIZE:=5}"              # signal-classification gate
: "${KBC_PCANGSD_POSTERIOR_CUTOFF:=}"             # see HANDOFF.md open Q4 — empty = inherit upstream call
