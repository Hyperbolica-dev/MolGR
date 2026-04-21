# pyright: reportMissingImports=false

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from molgr.fallback import xyz2omol
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


_MOLFILE_CASES_SPEC = importlib.util.spec_from_file_location(
    "molgr_cases_molfile",
    Path("scripts/molgr_cases_molfile.py").resolve(),
)
_MOLFILE_CASES_MODULE = importlib.util.module_from_spec(_MOLFILE_CASES_SPEC)
assert _MOLFILE_CASES_SPEC.loader is not None
_MOLFILE_CASES_SPEC.loader.exec_module(_MOLFILE_CASES_MODULE)
load_molfile_cases = _MOLFILE_CASES_MODULE.load_molfile_cases


def test_fallback_monnmo_molfile_regression() -> None:
    case = load_molfile_cases(Path("tests/data/sdf/MoNNMo.sdf"), limit=1)[0]

    result = xyz2omol(
        case["xyz_block"],
        total_charge=case["total_charge"],
        total_radical_electrons=case["total_radical_electrons"],
    )

    assert result is not None
    equivalent, info = check_equivalence(case["ground_truth_rdmol"], pybel_to_rdmol(result))
    assert equivalent, info.reason
