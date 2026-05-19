# atlas_schemas

JSON Schemas (Draft 2020-12) for the four data-layer pipelines the
inversion atlas depends on. Each pipeline has one IN schema (manifest
of paths + parameters) and one OUT schema (one schema = one row of the
emitted table).

| Pipeline | IN | OUT | Project location |
|---|---|---|---|
| Clair3 small-variant calling | `clair3.in.schema.json` | `clair3.out.schema.json` | `Modules/03_clair3/` |
| VESM missense scoring | `vesm.in.schema.json` | `vesm.out.schema.json` | MODULE_CONSERVATION step 14 |
| Deleterious / conservation annotation | `deleterious_conservation.in.schema.json` | `deleterious_conservation.out.schema.json` | MODULE_CONSERVATION steps 12–16, emits `variant_master_scored.tsv` |
| Sample groups | `sample_groups.in.schema.json` | `sample_groups.out.schema.json` | derived from PCAngsd / NGSadmix / NAToRA / ngsRelate / ROH outputs |

## Flow

```
                       ┌─ Clair3 OUT (per-sample VCF rows)
BAMs ──► Clair3 IN ───►│
                       └─ joint VCF (bcftools merge of all samples)
                                       │
                                       ▼
                       ┌──────────────────────────────┐
GFF3 ────────────────► │ Deleterious/Conservation IN  │
reference FASTA ─────► │                              │
SnpEff / SIFT4G / VESM │   per-variant annotation     │
splice module / scoring│                              │
weights                └─────────────┬────────────────┘
                                     │
        ┌────────────────────────────┘
        │
        ▼   (for missense rows only)
VESM IN ──► VESM OUT ──► merged back as the vesm_llr column
                                     │
                                     ▼
                  variant_master_scored.tsv  (= Deleterious/Conservation OUT)
                                     │
                                     ▼
                                  KBC, HAPS, downstream

──────────────────────────────────────────────────────────────

PCAngsd Q matrix    ┐
NGSadmix Q matrix   ├─► Sample groups IN ──► Sample groups OUT
NAToRA keep-list    │                          (one row per sample,
ngsRelate .res      │                           plus per-inversion
ROH BEDs            │                           karyotype sub-objects)
family hubs TSV     ┘
```

## Conventions

- **Row vs manifest.** OUT schemas describe one row of a tabular output
  (the consumer iterates rows). IN schemas describe one invocation /
  manifest object — paths + parameters.
- **Strictness.** All schemas set `additionalProperties: false`. Extend
  via PR rather than smuggling extra columns.
- **Coordinates.** VCF/GFF row schemas use VCF's native 1-based pos.
  The KBC internal `GeneVariants.gene_start/gene_end` is 0-based
  half-open — that is a KBC-internal representation and not part of
  this atlas-schema vocabulary.
- **Variant identity.** Every per-variant table joins on
  `variant_id = '{chrom}:{pos}:{ref}:{alt}'` after left-normalization
  and biallelic split.
- **Nulls.** Optional columns use `["<type>", "null"]` rather than
  omitting the key, so consumers see one stable shape.

## Validation

Schemas are not currently wired into a runtime validator outside of
KBC's own `test_schema_compliance.py`. Suggested next step: a small
script that takes a TSV + a schema and validates row-by-row using
`jsonschema`, runnable in CI.
