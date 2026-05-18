"""pytest fixtures for KBC MVP 1 tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add Modules/07_kbc/src to sys.path so `import kbc` works in tests.
MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT / "src"))

FIXTURES = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def variant_master_path() -> Path:
    return FIXTURES / "mini_variant_master.tsv"


@pytest.fixture
def genotypes_path() -> Path:
    return FIXTURES / "mini_genotypes.tsv"


@pytest.fixture
def karyotypes_path() -> Path:
    return FIXTURES / "mini_karyotype.tsv"


@pytest.fixture
def inversions_path() -> Path:
    return FIXTURES / "mini_inversions.tsv"


@pytest.fixture
def samples_path() -> Path:
    return FIXTURES / "mini_samples.tsv"


@pytest.fixture
def roh_dir() -> Path:
    return FIXTURES / "roh"
