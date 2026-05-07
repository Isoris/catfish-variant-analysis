# catfish-variant-analysis

Variant discovery and callable-region characterization for the
*Clarias gariepinus* hatchery cohort (n = 226):

- biallelic SNPs (ANGSD genotype likelihoods, Beagle GL format)
- callable mask
- small variants (clair3, per-sample VCF)
- structural variants (DELLY + Manta, dual-caller validation)

This is one of three sibling **catfish-{population,diversity,variant}-analysis**
repos that together produce the population-genetic primitives consumed by
[`inversion-atlas`](https://github.com/Isoris/inversion-atlas) and any future
papers on this cohort.

## What this repo produces

The variant data layer that everything else is built on:

| Output | Description | Consumed by |
|---|---|---|
| Callable-site mask (BED) | per-position pass/fail across the cohort | diversity, population, this repo |
| biSNP set (Beagle GL) | MAF ≥ 0.05, biallelic, thinned variants | population (PCA, NGSadmix), diversity (ROH) |
| SNP majmin assignments by RF | major/minor allele orientation per RF | population (Beagle subsets) |
| clair3 small-variant VCFs | per-sample, hard-genotyped | inversion-atlas (GHSL haplotype matrix) |
| DELLY SV calls | DEL / DUP / INV / BND / TRA / INS, per-sample | inversion-atlas (SV evidence layer) |
| Manta SV calls | dual-caller corroboration | inversion-atlas (SV evidence layer) |
| Dual-caller breakpoint-validated SV set | DELLY ∩ Manta | inversion-atlas |

## Inputs

| Input | From |
|---|---|
| BAMs + BAI | upstream read prep (currently `${BASE}/01-bams/` on LANTA) |
| Reference FASTA | `${BASE}/00-samples/fClaHyb_Gar_LG.fa` |
| Sample manifest | `${BASE}/01_inputs_check/` |

## Engines used

- ANGSD (system) — biSNP discovery via genotype likelihoods
- [`angsd_fixed_HWE`](https://github.com/Isoris/angsd_fixed_HWE) — patched
  ANGSD with fixed-F EM (used where standard HWE-based MAF estimation gives
  biased calls in the family-structured hatchery cohort)
- clair3 (system) — small-variant calling
- DELLY2 (system) — SV calling, primary
- Manta (system) — SV calling, secondary / corroborator
- bcftools, samtools, bedtools (system)

## Layout

```text
catfish-variant-analysis/
├── 00_config.sh                  root config
├── Modules/
│   ├── 01_callable_mask/         per-cohort callable-site BED
│   ├── 02_biSNP_discovery/       ANGSD biallelic SNP set + Beagle GLs
│   ├── 03_clair3/                small-variant calling per sample
│   ├── 04_sv_delly/              DELLY DEL/DUP/INV/BND/TRA/INS
│   ├── 05_sv_manta/              Manta SV calls
│   └── 06_sv_dual_validation/    DELLY ∩ Manta breakpoint validation
├── envs/
├── docs/
│   ├── module_contracts/
│   └── methods/
├── tests/
└── README.md
```

## Status

Scaffold. Pipelines exist on LANTA but live outside any git repo today
(spread across `${BASE}/popstruct_thin/04_beagle_byRF_majmin/`,
`${BASE}/03-variant-calls/`, etc). They will be migrated into `Modules/`
over time, one module at a time.

## Citation

Project umbrella DOI: TBD (Zenodo, will be issued at v1.0 tag).

## Cohort note

**This repo is for the 226-sample pure *Clarias gariepinus* hatchery cohort
only.** Do not use it for the F₁ hybrid (*C. gariepinus* × *C. macrocephalus*)
genome assembly cohort or any future *C. macrocephalus* wild cohort — those
are separate manuscripts and may need different parameter choices.
