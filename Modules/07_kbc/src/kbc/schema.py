"""JSON Schema validation helpers.

The schemas under Modules/07_kbc/schemas/ are the contract for the
per-table outputs. Each output writer validates rows against the
appropriate schema before emitting the TSV.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
from jsonschema import Draft202012Validator

SCHEMA_DIR = Path(__file__).resolve().parents[2] / "schemas"


def load_schema(name: str) -> dict:
    """Load a JSON schema by filename (e.g. 'A_kbc_sample_inversion_burden.schema.json')."""
    p = SCHEMA_DIR / name
    return json.loads(p.read_text())


def validate_rows(df: pd.DataFrame, schema_name: str) -> None:
    """Validate every row of `df` against the named schema. Raises on first failure."""
    schema = load_schema(schema_name)
    validator = Draft202012Validator(schema)
    for i, row in enumerate(df.to_dict(orient="records")):
        # pandas turns NaN into float('nan'); the schema wants null for nullable cols
        clean = {k: (None if (isinstance(v, float) and v != v) else v) for k, v in row.items()}
        errors = list(validator.iter_errors(clean))
        if errors:
            msgs = "; ".join(f"{list(e.absolute_path)}: {e.message}" for e in errors)
            raise ValueError(f"row {i} fails {schema_name}: {msgs}")


def write_table_a(df: pd.DataFrame, path: str | Path) -> None:
    """Validate against A_kbc_sample_inversion_burden.schema.json, then write TSV."""
    validate_rows(df, "A_kbc_sample_inversion_burden.schema.json")
    df.to_csv(path, sep="\t", index=False, na_rep="")
