#!/usr/bin/env bash
# =============================================================================
# 00_config.sh — catfish-variant-analysis
# =============================================================================
# Single source of truth for paths, parameters, and SLURM defaults used by
# every module in this repo. Source from any script:
#
#   source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"
#       — or, from a Module subdir —
#   source "$(dirname "${BASH_SOURCE[0]}")/../../00_config.sh"
#
# This repo is for the 226-sample pure *C. gariepinus* hatchery cohort.
# Do not use for other cohorts without auditing every parameter.
# =============================================================================

set -euo pipefail

# ── Project root ─────────────────────────────────────────────────────────────
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── This repo's location ────────────────────────────────────────────────────
export REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MODULES_DIR="${REPO_DIR}/Modules"

# ── Sibling repos ───────────────────────────────────────────────────────────
export DIVERSITY_REPO="${DIVERSITY_REPO:-${BASE}/catfish-diversity-analysis}"
export POPULATION_REPO="${POPULATION_REPO:-${BASE}/catfish-population-analysis}"

# ── Reference / cohort-stable inputs ────────────────────────────────────────
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"
export BAM_DIR="${BAM_DIR:-${BASE}/01-bams}"
export SAMPLES_IND="${BASE}/het_roh/01_inputs_check/samples.ind"
export BAMLIST="${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt"
export SAMPLE_LIST="${BASE}/01_inputs_check/samples_226_pure_gariepinus.txt"

# Linkage groups
export CHROM_LIST=(C_gar_LG01 C_gar_LG02 C_gar_LG03 C_gar_LG04 C_gar_LG05
                   C_gar_LG06 C_gar_LG07 C_gar_LG08 C_gar_LG09 C_gar_LG10
                   C_gar_LG11 C_gar_LG12 C_gar_LG13 C_gar_LG14 C_gar_LG15
                   C_gar_LG16 C_gar_LG17 C_gar_LG18 C_gar_LG19 C_gar_LG20
                   C_gar_LG21 C_gar_LG22 C_gar_LG23 C_gar_LG24 C_gar_LG25
                   C_gar_LG26 C_gar_LG27 C_gar_LG28)
export N_CHROM=${#CHROM_LIST[@]}

# ── Outputs ─────────────────────────────────────────────────────────────────
export OUTROOT="${OUTROOT:-${BASE}/results/catfish-variant-analysis}"
export OUT_CALLABLE="${OUTROOT}/01_callable_mask"
export OUT_BISNP="${OUTROOT}/02_biSNP"
export OUT_CLAIR3="${OUTROOT}/03_clair3"
export OUT_DELLY="${OUTROOT}/04_sv_delly"
export OUT_MANTA="${OUTROOT}/05_sv_manta"
export OUT_DUAL="${OUTROOT}/06_sv_dual_validation"
export LOG_DIR="${OUTROOT}/logs"

# Canonical paths downstream consumers can rely on (referenced from sibling
# repos by these well-known names, NOT by digging into module subfolders):
export CALLABLE_SITES="${OUT_CALLABLE}/callable_sites.bed.gz"
export BEAGLE_DIR="${OUT_BISNP}/beagle_byRF_majmin"

# ── biSNP discovery parameters ──────────────────────────────────────────────
export MAF_MIN=0.05
export SNP_PVAL=1e-6
export MIN_DEPTH=4
export MAX_DEPTH=99
export MIN_MAPQ=30
export MIN_BASEQ=20

# ── clair3 parameters ───────────────────────────────────────────────────────
export CLAIR3_MODEL_DIR="${CLAIR3_MODEL_DIR:-/path/to/clair3/models}"
export CLAIR3_PLATFORM="${CLAIR3_PLATFORM:-ilmn}"

# ── SV calling parameters ───────────────────────────────────────────────────
export DELLY_SV_TYPES=(DEL DUP INV BND TRA INS)
export MANTA_RUNTIME_MEM_GB=16

# ── Compute / SLURM ─────────────────────────────────────────────────────────
export RSCRIPT_BIN="${RSCRIPT_BIN:-/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript}"
export MAMBA_ENV="${MAMBA_ENV:-assembly}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
export SLURM_PARTITION="${SLURM_PARTITION:-compute}"

# ── Helpers ─────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%F %T')] [var] $*"; }
die() { echo "[$(date '+%F %T')] [var] [ERROR] $*" >&2; exit 1; }

init_dirs() {
  mkdir -p "${OUTROOT}" \
           "${OUT_CALLABLE}" "${OUT_BISNP}" "${OUT_CLAIR3}" \
           "${OUT_DELLY}" "${OUT_MANTA}" "${OUT_DUAL}" \
           "${LOG_DIR}"
}

config_print() {
  cat <<EOF
=== catfish-variant-analysis ===
BASE             : ${BASE}
REPO_DIR         : ${REPO_DIR}
OUTROOT          : ${OUTROOT}
N_CHROM          : ${N_CHROM}
MAF_MIN          : ${MAF_MIN}
RSCRIPT_BIN      : ${RSCRIPT_BIN}
================================
EOF
}

if [[ "${1:-}" == "print" ]]; then config_print; fi
