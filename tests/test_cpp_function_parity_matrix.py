# pyright: reportMissingImports=false

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from molgr import _core
from molgr.fallback.pipeline.reconstruct_with_metals import (
    _build_metal_states as py_build_metal_states,
)
from molgr.fallback.pipeline.reconstruct_with_metals import (
    combine_metal_with_omol as py_combine_metal_with_omol,
)
from molgr.fallback.pipeline.reconstruct_with_metals import (
    get_possible_metal_radicals as py_get_possible_metal_radicals,
)
from molgr.fallback.pipeline.reconstruct_with_metals import xyz2omol as py_xyz2omol
from molgr.fallback.pipeline.reconstruct_without_metals import (
    xyz_to_omol_no_metal as py_xyz_to_omol_no_metal,
)
from molgr.fallback.pipeline.resonance import (
    get_radical_resonances as py_get_radical_resonances,
)
from molgr.fallback.pipeline.resonance import process_resonance as py_process_resonance
from molgr.fallback.stages.break_bond import break_deformed_ene as py_break_deformed_ene
from molgr.fallback.stages.break_bond import break_one_bond as py_break_one_bond
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated as py_clean_carbene_neighbor_unsaturated,
)
from molgr.fallback.stages.clean import clean_neighbor_radicals as py_clean_neighbor_radicals
from molgr.fallback.stages.clean import clean_resonances as py_clean_resonances
from molgr.fallback.stages.eliminate import eliminate_1_3_dipole as py_eliminate_1_3_dipole
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom as py_eliminate_carbene_neighbor_heteroatom,
)
from molgr.fallback.stages.eliminate import eliminate_carboxyl as py_eliminate_carboxyl
from molgr.fallback.stages.eliminate import (
    eliminate_charge_spliting as py_eliminate_charge_spliting,
)
from molgr.fallback.stages.eliminate import eliminate_CN_in_doubt as py_eliminate_cn_in_doubt
from molgr.fallback.stages.eliminate import (
    eliminate_high_positive_charge_atoms as py_eliminate_high_positive_charge_atoms,
)
from molgr.fallback.stages.eliminate import (
    eliminate_negative_charges as py_eliminate_negative_charges,
)
from molgr.fallback.stages.eliminate import eliminate_NNN as py_eliminate_nnn
from molgr.fallback.stages.eliminate import (
    eliminate_positive_charges as py_eliminate_positive_charges,
)
from molgr.fallback.stages.fresh import (
    assign_charge_radical_for_atom as py_assign_charge_radical_for_atom,
)
from molgr.fallback.stages.fresh import assign_radical_dots as py_assign_radical_dots
from molgr.fallback.stages.fresh import fresh_omol_charge_radical as py_fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections as py_make_connections
from molgr.fallback.stages.preprocess import pre_clean as py_pre_clean
from molgr.fallback.stages.preprocess import validate_omol as py_validate_omol
from molgr.interface import mol_data_to_rdkit, pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


_pipeline: Any = _core.pipeline
_with_metals: Any = _pipeline.reconstruct_with_metals
_stages: Any = _core.stages
_BENCHMARK_HARD_CASE_INDICES: tuple[int, ...] = (7, 20, 24, 33, 35, 47, 49, 52)


@dataclass(frozen=True)
class FunctionParityRow:
    function_name: str
    case_name: str
    case_input: str
    cpp_source_path: str
    run_pair: Callable[[], tuple[Any, Any]]
    normalize: Callable[[Any], Any]


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _dump(value: Any) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=True)


def _diff_payload(fallback_value: Any, cpp_value: Any) -> dict[str, Any]:
    if fallback_value == cpp_value:
        return {"equal": True}
    return {
        "equal": False,
        "fallback_type": type(fallback_value).__name__,
        "cpp_type": type(cpp_value).__name__,
        "fallback": fallback_value,
        "cpp": cpp_value,
    }


def _smiles_token(mol: pybel.Molecule) -> str:
    smi = mol.write("smi")
    assert smi is not None
    return smi.split()[0]


def _smiles_from_mol_data(mol_data: _core.utils.MoleculeData) -> str:
    rdmol = mol_data_to_rdkit(mol_data)
    return Chem.MolToSmiles(rdmol, canonical=True)


def _molecule_signature(mol: pybel.Molecule) -> dict[str, Any]:
    atoms: list[dict[str, int]] = []
    for idx in range(1, mol.OBMol.NumAtoms() + 1):
        atom = mol.OBMol.GetAtom(idx)
        atoms.append(
            {
                "idx": idx,
                "atomic_num": atom.GetAtomicNum(),
                "formal_charge": atom.GetFormalCharge(),
                "spin": atom.GetSpinMultiplicity(),
            }
        )

    bonds: list[dict[str, int]] = []
    for bond in ob.OBMolBondIter(mol.OBMol):
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bonds.append(
            {
                "begin": begin_idx,
                "end": end_idx,
                "order": bond.GetBondOrder(),
            }
        )
    bonds.sort(key=lambda item: (item["begin"], item["end"], item["order"]))

    return {
        "smiles": _smiles_token(mol),
        "atoms": atoms,
        "bonds": bonds,
    }


def _charge_and_radicals(mol: pybel.Molecule) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.atoms:
        charge += atom.OBAtom.GetFormalCharge()
        radicals += atom.OBAtom.GetSpinMultiplicity()
    return charge, radicals


def _break_bond_signature(mol: pybel.Molecule) -> dict[str, Any]:
    total_charge, total_radicals = _charge_and_radicals(mol)
    return {
        "smiles": _smiles_token(mol),
        "total_charge": total_charge,
        "total_radicals": total_radicals,
        "atoms": [
            {
                "idx": atom.idx,
                "formal_charge": atom.OBAtom.GetFormalCharge(),
                "spin": atom.OBAtom.GetSpinMultiplicity(),
            }
            for atom in mol.atoms
        ],
        "bonds": sorted(
            [
                (
                    min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                    max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                    bond.GetBondOrder(),
                )
                for bond in ob.OBMolBondIter(mol.OBMol)
            ]
        ),
    }


def _set_twisted_ene_coordinates(mol: pybel.Molecule) -> None:
    coordinates = {
        1: (0.0, 0.0, 0.0),
        2: (1.2, 0.2, 0.0),
        3: (2.3, 1.0, 0.2),
        4: (3.0, 1.5, 1.4),
    }
    for idx, (x, y, z) in coordinates.items():
        atom = mol.OBMol.GetAtom(idx)
        if atom is not None:
            atom.SetVector(x, y, z)


def _seed_xyz_case_from_smiles(smiles: str) -> tuple[str, int, int]:
    seed = pybel.readstring("smi", smiles)
    xyz_block = str(seed.write("xyz"))
    total_charge, total_radical_electrons = _charge_and_radicals(seed)
    return xyz_block, total_charge, total_radical_electrons


def _load_test_case_smiles() -> list[str]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    lines = [
        line.strip() for line in csv_path.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    assert lines and lines[0] == "smiles"
    return lines[1:]


def _has_carbon_double_bond(smiles: str) -> bool:
    mol = pybel.readstring("smi", smiles)
    for bond in ob.OBMolBondIter(mol.OBMol):
        begin_atom = mol.OBMol.GetAtom(bond.GetBeginAtomIdx())
        end_atom = mol.OBMol.GetAtom(bond.GetEndAtomIdx())
        if (
            bond.GetBondOrder() == 2
            and begin_atom.GetAtomicNum() == 6
            and end_atom.GetAtomicNum() == 6
        ):
            return True
    return False


def _radical_atom_indices(mol: pybel.Molecule) -> tuple[int, ...]:
    return tuple(atom.idx for atom in mol.atoms if atom.OBAtom.GetSpinMultiplicity() > 0)


def _load_benchmark_case_map() -> dict[int, dict[str, Any]]:
    case_map: dict[int, dict[str, Any]] = {}
    smiles_rows = _load_test_case_smiles()
    for case_idx in _BENCHMARK_HARD_CASE_INDICES:
        input_smiles = smiles_rows[case_idx - 1]
        mol = Chem.MolFromSmiles(input_smiles)
        assert mol is not None
        mol_h = Chem.AddHs(mol)
        embed_code = rdDistGeom.EmbedMolecule(mol_h)
        assert int(embed_code) == 0
        total_charge = sum(
            int(mol_h.GetAtomWithIdx(atom_idx).GetFormalCharge())
            for atom_idx in range(mol_h.GetNumAtoms())
        )
        total_radical_electrons = sum(
            int(mol_h.GetAtomWithIdx(atom_idx).GetNumRadicalElectrons())
            for atom_idx in range(mol_h.GetNumAtoms())
        )
        case_map[case_idx] = {
            "case_idx": case_idx,
            "input_smiles": input_smiles,
            "xyz_block": Chem.MolToXYZBlock(mol_h),
            "total_charge": total_charge,
            "total_radical_electrons": total_radical_electrons,
        }
    assert len(case_map) == len(_BENCHMARK_HARD_CASE_INDICES)
    return case_map


def _run_benchmark_semantic_no_metal_case(
    case_idx: int,
    input_smiles: str,
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    fallback_status = "not_run"
    fallback_conversion = "not_run"
    fallback_error: str | None = None
    cpp_status = "not_run"
    cpp_conversion = "not_run"
    cpp_error: str | None = None
    equivalence_reason = "not_run"
    equivalent: bool | None = None

    fallback_rdmol: Chem.Mol | None = None
    cpp_rdmol: Chem.Mol | None = None

    try:
        py_output = py_xyz_to_omol_no_metal(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        if py_output is None:
            fallback_status = "reconstruct_none"
            fallback_conversion = "skipped"
        else:
            fallback_status = "reconstruct_ok"
            fallback_rdmol = pybel_to_rdmol(py_output)
            fallback_conversion = "rdkit_ok"
    except Exception as exc:  # noqa: BLE001
        fallback_status = "error"
        fallback_conversion = "rdkit_failed"
        fallback_error = f"{type(exc).__name__}: {exc}"

    try:
        cpp_output = _pipeline.reconstruct_without_metals.xyz_to_omol_no_metal(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        if cpp_output is None:
            cpp_status = "reconstruct_none"
            cpp_conversion = "skipped"
        else:
            cpp_status = "reconstruct_ok"
            cpp_rdmol = mol_data_to_rdkit(cpp_output)
            cpp_conversion = "rdkit_ok"
    except Exception as exc:  # noqa: BLE001
        cpp_status = "error"
        cpp_conversion = "rdkit_failed"
        cpp_error = f"{type(exc).__name__}: {exc}"

    if fallback_rdmol is not None and cpp_rdmol is not None:
        try:
            equivalent, info = check_equivalence(
                fallback_rdmol,
                cpp_rdmol,
                use_chirality=True,
                max_resonance=100,
            )
            equivalence_reason = info.reason or "unknown"
        except Exception as exc:  # noqa: BLE001
            equivalent = False
            equivalence_reason = f"equivalence_check_failed: {type(exc).__name__}: {exc}"

    benchmark_match = (
        fallback_conversion == "rdkit_ok" and cpp_conversion == "rdkit_ok" and equivalent is True
    )

    payload = {
        "benchmark_match": benchmark_match,
        "case_idx": case_idx,
        "input_smiles": input_smiles,
        "fallback_status": fallback_status,
        "fallback_conversion": fallback_conversion,
        "fallback_error": fallback_error,
        "cpp_status": cpp_status,
        "cpp_conversion": cpp_conversion,
        "cpp_error": cpp_error,
        "equivalent": equivalent,
        "equivalence_reason": equivalence_reason,
    }
    expected = dict(payload)
    return payload, expected


def _build_seed_obmol() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    b = obmol.NewAtom()
    b.SetAtomicNum(5)
    b.SetFormalCharge(0)
    b.SetSpinMultiplicity(0)

    for _ in range(4):
        h = obmol.NewAtom()
        h.SetAtomicNum(1)

    ne = obmol.NewAtom()
    ne.SetAtomicNum(10)
    ne.SetFormalCharge(7)
    ne.SetSpinMultiplicity(0)

    o = obmol.NewAtom()
    o.SetAtomicNum(8)
    o.SetFormalCharge(1)
    o.SetSpinMultiplicity(1)

    c = obmol.NewAtom()
    c.SetAtomicNum(6)
    c.SetFormalCharge(0)
    c.SetSpinMultiplicity(0)

    obmol.AddBond(1, 2, 1)
    obmol.AddBond(1, 3, 1)
    obmol.AddBond(1, 4, 1)
    obmol.AddBond(1, 5, 1)
    obmol.AddBond(8, 7, 2)

    obmol.EndModify()
    return obmol


def _clear_all_bonds(obmol: ob.OBMol) -> None:
    obmol.BeginModify()
    bonds = list(ob.OBMolBondIter(obmol))
    for bond in bonds:
        obmol.DeleteBond(bond)
    obmol.EndModify()


def _build_multi_rule_seed_a() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    atom1 = obmol.NewAtom()
    atom1.SetAtomicNum(6)
    atom1.SetFormalCharge(-1)

    atom2 = obmol.NewAtom()
    atom2.SetAtomicNum(7)
    atom2.SetFormalCharge(1)

    atom3 = obmol.NewAtom()
    atom3.SetAtomicNum(6)
    atom3.SetFormalCharge(0)

    atom4 = obmol.NewAtom()
    atom4.SetAtomicNum(6)
    atom4.SetFormalCharge(0)

    atom5 = obmol.NewAtom()
    atom5.SetAtomicNum(6)
    atom5.SetFormalCharge(0)

    atom6 = obmol.NewAtom()
    atom6.SetAtomicNum(6)
    atom6.SetFormalCharge(-1)

    obmol.AddBond(1, 2, 2)
    obmol.AddBond(2, 3, 2)
    obmol.AddBond(4, 5, 2)
    obmol.AddBond(5, 6, 2)

    obmol.EndModify()
    return obmol


def _make_seed(smiles: str, radical_atom_indices: Sequence[int]) -> pybel.Molecule:
    mol = pybel.readstring("smi", smiles)
    for idx in radical_atom_indices:
        mol.OBMol.GetAtom(idx).SetSpinMultiplicity(1)
    return mol


def _require_api(module_name: str, attr_name: str) -> bool:
    module = getattr(_stages, module_name, None)
    if module is None:
        return False
    return hasattr(module, attr_name)


def _matrix_rows() -> list[FunctionParityRow]:
    rows: list[FunctionParityRow] = []

    xyz_li_co = """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
"""
    xyz_no_metal = """2
NO
N 0.000 0.000 0.000
O 1.200 0.000 0.000
"""

    def run_xyz2omol_case(
        xyz_block: str,
        total_charge: int | None = None,
        total_radical_electrons: int | None = None,
    ) -> tuple[Any, Any]:
        if total_charge is None or total_radical_electrons is None:
            seed = pybel.readstring("xyz", xyz_block)
            charge, radicals = _charge_and_radicals(seed)
        else:
            charge = total_charge
            radicals = total_radical_electrons
        py_result = py_xyz2omol(xyz_block, charge, radicals)
        cpp_result = _with_metals.xyz2omol(xyz_block, charge, radicals)
        if py_result is None:
            py_norm = None
        else:
            py_norm = _smiles_from_mol_data(
                _core.utils.extract_molecule_data(_get_ptr(py_result.OBMol))
            )
        cpp_norm = None if cpp_result is None else _smiles_from_mol_data(cpp_result)
        return py_norm, cpp_norm

    rows.append(
        FunctionParityRow(
            function_name="pipeline.reconstruct_with_metals.xyz2omol",
            case_name="synthetic_li_co",
            case_input="xyz=LiCO,total_charge=0,total_radical=0",
            cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
            run_pair=lambda: run_xyz2omol_case(xyz_li_co),
            normalize=lambda value: value,
        )
    )

    hard_case_smiles = _load_test_case_smiles()
    benchmark_case_map = _load_benchmark_case_map()

    def run_break_one_bond_hard_case(
        smiles: str,
        given_charge: int,
        given_radical: int,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        py_mol = pybel.readstring("smi", smiles)
        cpp_mol = pybel.readstring("smi", smiles)
        py_after, py_charge = py_break_one_bond(
            py_mol,
            given_charge=given_charge,
            given_radical=given_radical,
        )
        cpp_charge = _stages.break_bond.break_one_bond_ptr(
            _get_ptr(cpp_mol.OBMol),
            given_charge,
            given_radical,
        )
        return (
            {
                "charge": int(py_charge),
                "molecule": _break_bond_signature(py_after),
            },
            {
                "charge": int(cpp_charge),
                "molecule": _break_bond_signature(cpp_mol),
            },
        )

    def run_break_deformed_ene_hard_case(
        smiles: str,
        given_charge: int,
        given_radical: int,
    ) -> tuple[pybel.Molecule, pybel.Molecule]:
        py_mol = pybel.readstring("smi", smiles)
        cpp_mol = pybel.readstring("smi", smiles)
        py_after = py_break_deformed_ene(
            py_mol,
            given_charge=given_charge,
            given_radical=given_radical,
            tolerance=5.0,
        )
        _stages.break_bond.break_deformed_ene_ptr(
            _get_ptr(cpp_mol.OBMol),
            given_charge,
            given_radical,
            5.0,
        )
        return py_after, cpp_mol

    for case_index in _BENCHMARK_HARD_CASE_INDICES:
        smiles = hard_case_smiles[case_index - 1]
        xyz_block, total_charge, total_radical_electrons = _seed_xyz_case_from_smiles(smiles)
        case_name = f"benchmark_case_{case_index}"
        base_case_input = (
            f"benchmark_case_index={case_index},"
            f"smiles={smiles},"
            f"total_charge={total_charge},"
            f"total_radical={total_radical_electrons}"
        )

        rows.append(
            FunctionParityRow(
                function_name="pipeline.reconstruct_without_metals.xyz_to_omol_no_metal",
                case_name=case_name,
                case_input=base_case_input,
                cpp_source_path="src/cpp/src/pipeline/reconstruct_without_metals.cpp",
                run_pair=lambda x=xyz_block, c=total_charge, r=total_radical_electrons: (
                    run_xyz_to_omol_no_metal_case(
                        x,
                        c,
                        r,
                    )
                ),
                normalize=lambda value: value,
            )
        )

        rows.append(
            FunctionParityRow(
                function_name="stages.break_bond.break_one_bond_ptr",
                case_name=case_name,
                case_input=base_case_input,
                cpp_source_path="src/cpp/src/stages/break_bond.cpp",
                run_pair=lambda s=smiles, c=total_charge, r=total_radical_electrons: (
                    run_break_one_bond_hard_case(
                        s,
                        c,
                        r,
                    )
                ),
                normalize=lambda value: value,
            )
        )

        if _has_carbon_double_bond(smiles):
            rows.append(
                FunctionParityRow(
                    function_name="stages.break_bond.break_deformed_ene_ptr",
                    case_name=case_name,
                    case_input=f"{base_case_input},tolerance=5.0",
                    cpp_source_path="src/cpp/src/stages/break_bond.cpp",
                    run_pair=lambda s=smiles, c=total_charge, r=total_radical_electrons: (
                        run_break_deformed_ene_hard_case(
                            s,
                            c,
                            r,
                        )
                    ),
                    normalize=_break_bond_signature,
                )
            )

        py_intermediate = py_xyz_to_omol_no_metal(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        if py_intermediate is None:
            continue
        radical_atom_indices = _radical_atom_indices(py_intermediate)
        if not radical_atom_indices:
            continue

        rows.append(
            FunctionParityRow(
                function_name="pipeline.resonance.process_resonance_ptr",
                case_name=case_name,
                case_input=(
                    f"benchmark_case_index={case_index},"
                    f"intermediate_from_xyz_to_omol_no_metal,"
                    f"smiles={smiles},"
                    f"radical_atom_indices={list(radical_atom_indices)},"
                    f"charge={total_charge}"
                ),
                cpp_source_path="src/cpp/src/pipeline/resonance.cpp",
                run_pair=lambda x=xyz_block, c=total_charge, r=total_radical_electrons: (
                    _run_process_resonance_from_xyz_intermediate_case(
                        x,
                        c,
                        r,
                        c,
                    )
                ),
                normalize=lambda value: value,
            )
        )

    for case_index in _BENCHMARK_HARD_CASE_INDICES:
        benchmark_case = benchmark_case_map[case_index]
        input_smiles = str(benchmark_case["input_smiles"])
        xyz_block = benchmark_case["xyz_block"]
        total_charge = benchmark_case["total_charge"]
        total_radical_electrons = benchmark_case["total_radical_electrons"]

        assert isinstance(xyz_block, str)
        assert isinstance(total_charge, int)
        assert isinstance(total_radical_electrons, int)

        rows.append(
            FunctionParityRow(
                function_name="benchmark.reconstruct_without_metals.xyz_to_omol_no_metal_equivalence",
                case_name=f"benchmark_case_{case_index}",
                case_input=(
                    f"case_idx={case_index},"
                    f"input_smiles={input_smiles},"
                    f"provider=rdkit_embed_xyz,"
                    f"total_charge={total_charge},"
                    f"total_radical={total_radical_electrons},"
                    f"use_chirality=True,max_resonance=100"
                ),
                cpp_source_path="src/cpp/src/pipeline/reconstruct_without_metals.cpp",
                run_pair=lambda idx=case_index, smiles=input_smiles, x=xyz_block, c=total_charge, r=total_radical_electrons: (
                    _run_benchmark_semantic_no_metal_case(
                        idx,
                        smiles,
                        x,
                        c,
                        r,
                    )
                ),
                normalize=lambda value: value,
            )
        )

    for case_name, xyz_block, total_charge, total_radical_electrons in [
        (
            "synthetic_na_no",
            """3
NaNO
Na 0.0 0.0 0.0
N 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
            0,
            0,
        ),
        (
            "synthetic_li_na_co",
            """4
LiNaCO
Li 0.0 0.0 0.0
Na 0.0 1.5 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
            0,
            0,
        ),
        (
            "metal_only_li",
            """1
Li
Li 0.0 0.0 0.0
""",
            0,
            0,
        ),
    ]:
        rows.append(
            FunctionParityRow(
                function_name="pipeline.reconstruct_with_metals.xyz2omol",
                case_name=case_name,
                case_input=(
                    f"xyz={case_name},"
                    f"total_charge={total_charge},"
                    f"total_radical={total_radical_electrons}"
                ),
                cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
                run_pair=lambda x=xyz_block, c=total_charge, r=total_radical_electrons: (
                    run_xyz2omol_case(x, c, r)
                ),
                normalize=lambda value: value,
            )
        )

    def run_xyz_to_omol_no_metal_case(
        xyz_block: str,
        total_charge: int,
        total_radical_electrons: int,
    ) -> tuple[Any, Any]:
        py_result = py_xyz_to_omol_no_metal(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        cpp_result = _pipeline.reconstruct_without_metals.xyz_to_omol_no_metal(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        if py_result is None:
            py_norm = None
        else:
            py_norm = _smiles_from_mol_data(
                _core.utils.extract_molecule_data(_get_ptr(py_result.OBMol))
            )
        cpp_norm = None if cpp_result is None else _smiles_from_mol_data(cpp_result)
        return py_norm, cpp_norm

    for idx, smiles in enumerate(
        [
            "[H]C#Cc1c([C-]=O)nc([H])n1[H]",
            "O=C=C1C=CNC([O-])=C1",
            "N#C[C]([O-])N=C1C[NH2+]C1",
            "[C]#Cc1nnc(C#C)o1",
        ],
        start=1,
    ):
        xyz_block, total_charge, total_radical_electrons = _seed_xyz_case_from_smiles(smiles)
        rows.append(
            FunctionParityRow(
                function_name="pipeline.reconstruct_without_metals.xyz_to_omol_no_metal",
                case_name=f"seeded_case_{idx}",
                case_input=(
                    f"smiles={smiles},"
                    f"total_charge={total_charge},"
                    f"total_radical={total_radical_electrons}"
                ),
                cpp_source_path="src/cpp/src/pipeline/reconstruct_without_metals.cpp",
                run_pair=lambda x=xyz_block, c=total_charge, r=total_radical_electrons: (
                    run_xyz_to_omol_no_metal_case(x, c, r)
                ),
                normalize=lambda value: value,
            )
        )
    rows.append(
        FunctionParityRow(
            function_name="pipeline.reconstruct_with_metals.xyz2omol",
            case_name="no_metal_edge",
            case_input="xyz=NO,total_charge=0,total_radical=0",
            cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
            run_pair=lambda: run_xyz2omol_case(xyz_no_metal),
            normalize=lambda value: value,
        )
    )

    for metal, valence in [("Fe", 2), ("Cu", 2), ("Li", 10), ("UUnobtainium", 1)]:
        rows.append(
            FunctionParityRow(
                function_name="pipeline.reconstruct_with_metals.get_possible_metal_radicals",
                case_name=f"{metal}_{valence}",
                case_input=f"metal={metal},valence={valence}",
                cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
                run_pair=lambda m=metal, v=valence: (
                    py_get_possible_metal_radicals(m, v),
                    _with_metals.get_possible_metal_radicals(m, v),
                ),
                normalize=lambda value: sorted(value),
            )
        )

    rows.append(
        FunctionParityRow(
            function_name="pipeline.reconstruct_with_metals.build_metal_states_ptr",
            case_name="li_co_atom_1",
            case_input="xyz=LiCO,atom_idx=1",
            cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
            run_pair=_run_build_metal_states_case,
            normalize=lambda states: [
                {
                    "idx": state.idx,
                    "symbol": state.symbol,
                    "element_idx": state.element_idx,
                    "valence": state.valence,
                    "radical_num": state.radical_num,
                    "x": state.position_x,
                    "y": state.position_y,
                    "z": state.position_z,
                }
                for state in states
            ],
        )
    )

    def run_combine_metal_case() -> tuple[Any, Any]:
        metal_seed_py = pybel.readstring("xyz", xyz_li_co)
        metal_seed_cpp = pybel.readstring("xyz", xyz_li_co)
        py_state = py_build_metal_states(metal_seed_py.OBMol.GetAtom(1))[0]
        cpp_state = _with_metals.build_metal_states_ptr(_get_ptr(metal_seed_cpp.OBMol), 1)[0]

        organic_xyz = """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
"""
        organic_py = pybel.readstring("xyz", organic_xyz)
        organic_cpp = pybel.readstring("xyz", organic_xyz)
        py_after = py_combine_metal_with_omol(organic_py, [py_state])
        _with_metals.combine_metal_with_omol_ptr(_get_ptr(organic_cpp.OBMol), [cpp_state])
        return py_after, organic_cpp

    rows.append(
        FunctionParityRow(
            function_name="pipeline.reconstruct_with_metals.combine_metal_with_omol_ptr",
            case_name="single_metal_insert",
            case_input="organic=CO,metal=Li",
            cpp_source_path="src/cpp/src/pipeline/reconstruct_with_metals.cpp",
            run_pair=run_combine_metal_case,
            normalize=_molecule_signature,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.fresh.fresh_omol_charge_radical_ptr",
            case_name="seed_obmol",
            case_input="manual_seed_with_b_h4_ne_o_c",
            cpp_source_path="src/cpp/src/stages/fresh.cpp",
            run_pair=lambda: (
                py_fresh_omol_charge_radical(pybel.Molecule(_build_seed_obmol())),
                _run_fresh_cpp(),
            ),
            normalize=_molecule_signature,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.fresh.assign_radical_dots_ptr",
            case_name="atom_6_neon",
            case_input="manual_seed_atom_idx=6",
            cpp_source_path="src/cpp/src/stages/fresh.cpp",
            run_pair=_run_assign_radical_dots_case,
            normalize=lambda value: value,
        )
    )

    def run_assign_charge_case() -> tuple[Any, Any]:
        py_mol = pybel.Molecule(_build_seed_obmol())
        cpp_mol = pybel.Molecule(_build_seed_obmol())
        py_atom = py_mol.OBMol.GetAtom(8)
        py_assign_charge_radical_for_atom(py_atom)
        _stages.fresh.assign_charge_radical_for_atom_ptr(_get_ptr(cpp_mol.OBMol), 8)
        return py_mol, cpp_mol

    rows.append(
        FunctionParityRow(
            function_name="stages.fresh.assign_charge_radical_for_atom_ptr",
            case_name="atom_8_oxygen",
            case_input="manual_seed_atom_idx=8",
            cpp_source_path="src/cpp/src/stages/fresh.cpp",
            run_pair=run_assign_charge_case,
            normalize=_molecule_signature,
        )
    )

    def run_make_connections_case() -> tuple[Any, Any]:
        xyz_block = """2
no_bond
N 0.000 0.000 0.000
O 1.200 0.000 0.000
"""
        py_mol = pybel.readstring("xyz", xyz_block)
        cpp_mol = pybel.readstring("xyz", xyz_block)
        _clear_all_bonds(py_mol.OBMol)
        _clear_all_bonds(cpp_mol.OBMol)
        py_make_connections(py_mol, factor=1.4)
        _stages.preprocess.make_connections_ptr(_get_ptr(cpp_mol.OBMol), factor=1.4)
        return py_mol, cpp_mol

    rows.append(
        FunctionParityRow(
            function_name="stages.preprocess.make_connections_ptr",
            case_name="no_bond_no_edge",
            case_input="xyz=NO,factor=1.4",
            cpp_source_path="src/cpp/src/stages/preprocess.cpp",
            run_pair=run_make_connections_case,
            normalize=_molecule_signature,
        )
    )

    def run_pre_clean_case() -> tuple[Any, Any]:
        xyz_block = """6
sif5_setup
Si 0.000 0.000 0.000
F 1.600 0.000 0.000
F -1.600 0.000 0.000
F 0.000 1.600 0.000
F 0.000 -1.600 0.000
F 0.000 0.000 1.600
"""
        py_mol = pybel.readstring("xyz", xyz_block)
        cpp_mol = pybel.readstring("xyz", xyz_block)
        for mol in (py_mol, cpp_mol):
            _clear_all_bonds(mol.OBMol)
            mol.OBMol.BeginModify()
            for idx in range(2, 7):
                mol.OBMol.AddBond(1, idx, 1)
            mol.OBMol.EndModify()
        py_pre_clean(py_mol)
        _stages.preprocess.pre_clean_ptr(_get_ptr(cpp_mol.OBMol))
        return py_mol, cpp_mol

    rows.append(
        FunctionParityRow(
            function_name="stages.preprocess.pre_clean_ptr",
            case_name="sif5_setup",
            case_input="xyz=SiF5,bond_order=single",
            cpp_source_path="src/cpp/src/stages/preprocess.cpp",
            run_pair=run_pre_clean_case,
            normalize=_molecule_signature,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.preprocess.validate_omol_ptr",
            case_name="carbene_spin_mod2",
            case_input="xyz=C,total_charge=0,total_radical=1",
            cpp_source_path="src/cpp/src/stages/preprocess.cpp",
            run_pair=_run_validate_case,
            normalize=lambda value: value,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_1_3_dipole_ptr",
            case_name="o_minus_n_double_c",
            case_input="smiles=[O-]N=C,given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_13_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_positive_charges_ptr",
            case_name="single_carbon_radical",
            case_input="smiles=C(atom1.spin=2),given_charge=1",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_positive_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_negative_charges_ptr",
            case_name="single_oxygen_radical",
            case_input="smiles=O(atom1.spin=2),given_charge=-1",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_negative_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_nnn_ptr",
            case_name="n_n_n_single_bond_chain",
            case_input="manual_seed=N-N-N,given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_nnn_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_high_positive_charge_atoms_ptr",
            case_name="cation_single_bond_oxygen",
            case_input="manual_seed=C+1-O(spin=1),given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_high_positive_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_cn_in_doubt_ptr",
            case_name="c_equals_n_plus",
            case_input="smiles=[CH2]=[NH2+],given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_cn_in_doubt_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_carboxyl_ptr",
            case_name="neutral_carboxyl",
            case_input="smiles=OC=O,given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_carboxyl_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_carbene_neighbor_heteroatom_ptr",
            case_name="carbene_next_to_oxygen",
            case_input="smiles=CO(atom1.spin=2),given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_carbene_neighbor_heteroatom_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.eliminate.eliminate_charge_spliting_ptr",
            case_name="three_radicals_cco",
            case_input="smiles=CCO(atom1,2,3.spin=1),given_charge=0",
            cpp_source_path="src/cpp/src/stages/eliminate.cpp",
            run_pair=_run_eliminate_charge_spliting_case,
            normalize=lambda value: (
                _molecule_signature(value[0]),
                value[1],
            ),
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.clean.clean_neighbor_radicals_ptr",
            case_name="adjacent_radicals_cc",
            case_input="smiles=CC(atom1.spin=1,atom2.spin=1)",
            cpp_source_path="src/cpp/src/stages/clean.cpp",
            run_pair=_run_clean_neighbor_radicals_case,
            normalize=_molecule_signature,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.clean.clean_carbene_neighbor_unsaturated_ptr",
            case_name="carbene_allyl_shift",
            case_input="smiles=CC=C(atom1.spin=2)",
            cpp_source_path="src/cpp/src/stages/clean.cpp",
            run_pair=_run_clean_carbene_neighbor_unsaturated_case,
            normalize=_molecule_signature,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="stages.clean.clean_resonances_ptr",
            case_name="multi_rule_seed_a",
            case_input="manual_seed_multi_rule_a",
            cpp_source_path="src/cpp/src/stages/clean.cpp",
            run_pair=lambda: (
                py_clean_resonances(pybel.Molecule(_build_multi_rule_seed_a())),
                _run_clean_cpp(),
            ),
            normalize=_molecule_signature,
        )
    )

    def run_break_deformed_ene_case() -> tuple[Any, Any]:
        py_mol = pybel.readstring("smi", "CC(=C)C")
        cpp_mol = pybel.readstring("smi", "CC(=C)C")
        _set_twisted_ene_coordinates(py_mol)
        _set_twisted_ene_coordinates(cpp_mol)
        py_after = py_break_deformed_ene(py_mol, given_charge=0, given_radical=2, tolerance=5.0)
        _stages.break_bond.break_deformed_ene_ptr(_get_ptr(cpp_mol.OBMol), 0, 2, 5.0)
        return py_after, cpp_mol

    rows.append(
        FunctionParityRow(
            function_name="stages.break_bond.break_deformed_ene_ptr",
            case_name="twisted_substituted_ene",
            case_input="smiles=CC(=C)C,given_charge=0,given_radical=2,tolerance=5.0",
            cpp_source_path="src/cpp/src/stages/break_bond.cpp",
            run_pair=run_break_deformed_ene_case,
            normalize=_break_bond_signature,
        )
    )

    def run_break_one_bond_charge_case() -> tuple[Any, Any]:
        py_mol = pybel.readstring("smi", "[NH+]=C")
        cpp_mol = pybel.readstring("smi", "[NH+]=C")
        py_after, py_charge = py_break_one_bond(py_mol, given_charge=0, given_radical=1)
        cpp_charge = _stages.break_bond.break_one_bond_ptr(_get_ptr(cpp_mol.OBMol), 0, 1)
        return (
            {
                "charge": int(py_charge),
                "molecule": _break_bond_signature(py_after),
            },
            {
                "charge": int(cpp_charge),
                "molecule": _break_bond_signature(cpp_mol),
            },
        )

    rows.append(
        FunctionParityRow(
            function_name="stages.break_bond.break_one_bond_ptr",
            case_name="n_plus_double_bond_charge_transfer",
            case_input="smiles=[NH+]=C,given_charge=0,given_radical=1",
            cpp_source_path="src/cpp/src/stages/break_bond.cpp",
            run_pair=run_break_one_bond_charge_case,
            normalize=lambda value: value,
        )
    )

    def run_break_one_bond_single_delete_case() -> tuple[Any, Any]:
        py_mol = pybel.readstring("smi", "CC")
        cpp_mol = pybel.readstring("smi", "CC")
        py_after, py_charge = py_break_one_bond(py_mol, given_charge=0, given_radical=3)
        cpp_charge = _stages.break_bond.break_one_bond_ptr(_get_ptr(cpp_mol.OBMol), 0, 3)
        return (
            {
                "charge": int(py_charge),
                "molecule": _break_bond_signature(py_after),
            },
            {
                "charge": int(cpp_charge),
                "molecule": _break_bond_signature(cpp_mol),
            },
        )

    rows.append(
        FunctionParityRow(
            function_name="stages.break_bond.break_one_bond_ptr",
            case_name="single_bond_delete_fallback",
            case_input="smiles=CC,given_charge=0,given_radical=3",
            cpp_source_path="src/cpp/src/stages/break_bond.cpp",
            run_pair=run_break_one_bond_single_delete_case,
            normalize=lambda value: value,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="pipeline.resonance.get_radical_resonances_ptr",
            case_name="cc_equals_c_atom1",
            case_input="smiles=CC=C,radical_atom_indices=[1]",
            cpp_source_path="src/cpp/src/pipeline/resonance.cpp",
            run_pair=lambda: _run_get_radical_resonances_case("CC=C", (1,)),
            normalize=lambda value: value,
        )
    )

    rows.append(
        FunctionParityRow(
            function_name="pipeline.resonance.process_resonance_ptr",
            case_name="c_equals_cc_equals_c_charge0",
            case_input="smiles=C=CC=C,radical_atom_indices=[2],charge=0",
            cpp_source_path="src/cpp/src/pipeline/resonance.cpp",
            run_pair=lambda: _run_process_resonance_case("C=CC=C", (2,), 0),
            normalize=lambda value: value,
        )
    )

    return rows


def _run_fresh_cpp() -> pybel.Molecule:
    mol = pybel.Molecule(_build_seed_obmol())
    _stages.fresh.fresh_omol_charge_radical_ptr(_get_ptr(mol.OBMol))
    return mol


def _run_build_metal_states_case() -> tuple[Any, Any]:
    py_mol = pybel.readstring(
        "xyz",
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    cpp_mol = pybel.readstring(
        "xyz",
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    return (
        py_build_metal_states(py_mol.OBMol.GetAtom(1)),
        _with_metals.build_metal_states_ptr(_get_ptr(cpp_mol.OBMol), 1),
    )


def _run_assign_radical_dots_case() -> tuple[int, int]:
    py_mol = pybel.Molecule(_build_seed_obmol())
    cpp_mol = pybel.Molecule(_build_seed_obmol())
    return (
        py_assign_radical_dots(py_mol.OBMol.GetAtom(6)),
        _stages.fresh.assign_radical_dots_ptr(_get_ptr(cpp_mol.OBMol), 6),
    )


def _run_clean_cpp() -> pybel.Molecule:
    mol = pybel.Molecule(_build_multi_rule_seed_a())
    _stages.clean.clean_resonances_ptr(_get_ptr(mol.OBMol))
    return mol


def _run_validate_case() -> tuple[bool, bool]:
    xyz_block = """1
carbene
C 0.000 0.000 0.000
"""
    py_mol = pybel.readstring("xyz", xyz_block)
    cpp_mol = pybel.readstring("xyz", xyz_block)

    py_atom = py_mol.atoms[0].OBAtom
    cpp_atom = cpp_mol.atoms[0].OBAtom
    py_atom.SetFormalCharge(0)
    cpp_atom.SetFormalCharge(0)
    py_atom.SetSpinMultiplicity(1)
    cpp_atom.SetSpinMultiplicity(1)

    return py_validate_omol(py_mol, 0, 1), _stages.preprocess.validate_omol_ptr(
        _get_ptr(cpp_mol.OBMol), 0, 1
    )


def _run_eliminate_13_case() -> tuple[tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]]:
    py_mol = pybel.readstring("smi", "[O-]N=C")
    cpp_mol = pybel.readstring("smi", "[O-]N=C")
    py_out = py_eliminate_1_3_dipole(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_1_3_dipole_ptr(_get_ptr(cpp_mol.OBMol), 0)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_positive_case() -> tuple[tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]]:
    py_mol = pybel.readstring("smi", "C")
    cpp_mol = pybel.readstring("smi", "C")
    py_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    cpp_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    py_out = py_eliminate_positive_charges(py_mol, 1)
    cpp_charge = _stages.eliminate.eliminate_positive_charges_ptr(_get_ptr(cpp_mol.OBMol), 1)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_negative_case() -> tuple[tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]]:
    py_mol = pybel.readstring("smi", "O")
    cpp_mol = pybel.readstring("smi", "O")
    py_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    cpp_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    py_out = py_eliminate_negative_charges(py_mol, -1)
    cpp_charge = _stages.eliminate.eliminate_negative_charges_ptr(_get_ptr(cpp_mol.OBMol), -1)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_nnn_case() -> tuple[tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]]:
    py_mol = pybel.readstring("smi", "N=N=N")
    cpp_mol = pybel.readstring("smi", "N=N=N")
    py_out = py_eliminate_nnn(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_nnn_ptr(_get_ptr(cpp_mol.OBMol), 0)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_high_positive_case() -> tuple[
    tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]
]:
    py_mol = pybel.readstring("smi", "CO")
    cpp_mol = pybel.readstring("smi", "CO")
    py_mol.OBMol.GetAtom(1).SetFormalCharge(1)
    cpp_mol.OBMol.GetAtom(1).SetFormalCharge(1)
    py_mol.OBMol.GetAtom(2).SetSpinMultiplicity(1)
    cpp_mol.OBMol.GetAtom(2).SetSpinMultiplicity(1)
    py_out = py_eliminate_high_positive_charge_atoms(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_high_positive_charge_atoms_ptr(
        _get_ptr(cpp_mol.OBMol), 0
    )
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_cn_in_doubt_case() -> tuple[
    tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]
]:
    py_mol = pybel.readstring("smi", "[CH2]=[NH2+]")
    cpp_mol = pybel.readstring("smi", "[CH2]=[NH2+]")
    py_out = py_eliminate_cn_in_doubt(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_cn_in_doubt_ptr(_get_ptr(cpp_mol.OBMol), 0)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_carboxyl_case() -> tuple[tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]]:
    py_mol = pybel.readstring("smi", "OC=O")
    cpp_mol = pybel.readstring("smi", "OC=O")
    py_out = py_eliminate_carboxyl(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_carboxyl_ptr(_get_ptr(cpp_mol.OBMol), 0)
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_carbene_neighbor_heteroatom_case() -> tuple[
    tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]
]:
    py_mol = pybel.readstring("smi", "CO")
    cpp_mol = pybel.readstring("smi", "CO")
    py_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    cpp_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    py_out = py_eliminate_carbene_neighbor_heteroatom(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_carbene_neighbor_heteroatom_ptr(
        _get_ptr(cpp_mol.OBMol),
        0,
    )
    return py_out, (cpp_mol, cpp_charge)


def _run_eliminate_charge_spliting_case() -> tuple[
    tuple[pybel.Molecule, int], tuple[pybel.Molecule, int]
]:
    py_mol = pybel.readstring("smi", "CCO")
    cpp_mol = pybel.readstring("smi", "CCO")
    for idx in (1, 2, 3):
        py_mol.OBMol.GetAtom(idx).SetSpinMultiplicity(1)
        cpp_mol.OBMol.GetAtom(idx).SetSpinMultiplicity(1)
    py_out = py_eliminate_charge_spliting(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_charge_spliting_ptr(_get_ptr(cpp_mol.OBMol), 0)
    return py_out, (cpp_mol, cpp_charge)


def _run_clean_neighbor_radicals_case() -> tuple[pybel.Molecule, pybel.Molecule]:
    py_mol = pybel.readstring("smi", "CC")
    cpp_mol = pybel.readstring("smi", "CC")
    py_mol.OBMol.GetAtom(1).SetSpinMultiplicity(1)
    py_mol.OBMol.GetAtom(2).SetSpinMultiplicity(1)
    cpp_mol.OBMol.GetAtom(1).SetSpinMultiplicity(1)
    cpp_mol.OBMol.GetAtom(2).SetSpinMultiplicity(1)
    py_after = py_clean_neighbor_radicals(py_mol)
    _stages.clean.clean_neighbor_radicals_ptr(_get_ptr(cpp_mol.OBMol))
    return py_after, cpp_mol


def _run_clean_carbene_neighbor_unsaturated_case() -> tuple[pybel.Molecule, pybel.Molecule]:
    py_mol = pybel.readstring("smi", "CC=C")
    cpp_mol = pybel.readstring("smi", "CC=C")
    py_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    cpp_mol.OBMol.GetAtom(1).SetSpinMultiplicity(2)
    py_after = py_clean_carbene_neighbor_unsaturated(py_mol)
    _stages.clean.clean_carbene_neighbor_unsaturated_ptr(_get_ptr(cpp_mol.OBMol))
    return py_after, cpp_mol


def _run_get_radical_resonances_case(
    smiles: str,
    radical_atom_indices: Sequence[int],
) -> tuple[list[str], list[str]]:
    seed = _make_seed(smiles, radical_atom_indices)
    py_seed = seed.clone
    cpp_seed = seed.clone
    py_tokens = [_smiles_token(mol) for mol in py_get_radical_resonances(py_seed)]
    cpp_ptrs = [int(ptr) for ptr in _pipeline.get_radical_resonances_ptr(_get_ptr(cpp_seed.OBMol))]
    try:
        cpp_tokens = [_pipeline.smiles_token_ptr(ptr) for ptr in cpp_ptrs]
    finally:
        for ptr in cpp_ptrs:
            _core.free_obmol_ptr(ptr)
    return py_tokens, cpp_tokens


def _run_process_resonance_case(
    smiles: str,
    radical_atom_indices: Sequence[int],
    charge: int,
) -> tuple[list[tuple[str, int]], list[tuple[str, int]]]:
    seed = _make_seed(smiles, radical_atom_indices)
    py_seed = seed.clone
    cpp_seed = seed.clone
    py_resonances = py_get_radical_resonances(py_seed)
    cpp_ptrs = [int(ptr) for ptr in _pipeline.get_radical_resonances_ptr(_get_ptr(cpp_seed.OBMol))]
    try:
        py_out: list[tuple[str, int]] = []
        cpp_out: list[tuple[str, int]] = []
        for py_resonance, cpp_ptr in zip(py_resonances, cpp_ptrs):
            py_processed, py_charge = py_process_resonance(py_resonance.clone, charge)
            cpp_processed_ptr, cpp_charge = _pipeline.process_resonance_ptr(cpp_ptr, charge)
            cpp_processed_ptr = int(cpp_processed_ptr)
            try:
                py_out.append((_smiles_token(py_processed), int(py_charge)))
                cpp_out.append((_pipeline.smiles_token_ptr(cpp_processed_ptr), int(cpp_charge)))
            finally:
                _core.free_obmol_ptr(cpp_processed_ptr)
    finally:
        for ptr in cpp_ptrs:
            _core.free_obmol_ptr(ptr)
    return py_out, cpp_out


def _run_process_resonance_from_xyz_intermediate_case(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    charge: int,
) -> tuple[list[tuple[str, int]], list[tuple[str, int]]]:
    py_seed = py_xyz_to_omol_no_metal(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    cpp_seed_data = _pipeline.reconstruct_without_metals.xyz_to_omol_no_metal(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    py_out: list[tuple[str, int]] = []
    if py_seed is not None:
        for py_resonance in py_get_radical_resonances(py_seed.clone):
            py_processed, py_charge = py_process_resonance(py_resonance.clone, charge)
            py_out.append((_smiles_token(py_processed), int(py_charge)))

    cpp_out: list[tuple[str, int]] = []
    if cpp_seed_data is None:
        return py_out, cpp_out

    cpp_seed_ptr = _core.utils.molecule_data_to_obmol_ptr(cpp_seed_data)
    try:
        cpp_ptrs = [int(ptr) for ptr in _pipeline.get_radical_resonances_ptr(cpp_seed_ptr)]
        try:
            for cpp_ptr in cpp_ptrs:
                cpp_processed_ptr, cpp_charge = _pipeline.process_resonance_ptr(cpp_ptr, charge)
                cpp_processed_ptr = int(cpp_processed_ptr)
                try:
                    cpp_out.append((_pipeline.smiles_token_ptr(cpp_processed_ptr), int(cpp_charge)))
                finally:
                    _core.free_obmol_ptr(cpp_processed_ptr)
        finally:
            for ptr in cpp_ptrs:
                _core.free_obmol_ptr(ptr)
    finally:
        _core.free_obmol_ptr(cpp_seed_ptr)
    return py_out, cpp_out


MATRIX_ROWS = _matrix_rows()


@pytest.mark.parametrize(
    "row",
    MATRIX_ROWS,
    ids=lambda row: f"{row.function_name}::{row.case_name}",
)
def test_cpp_parity_function_matrix(row: FunctionParityRow) -> None:
    py_output, cpp_output = row.run_pair()
    py_norm = row.normalize(py_output)
    cpp_norm = row.normalize(cpp_output)

    assert py_norm == cpp_norm, "\n".join(
        [
            f"FAIL_ID={row.function_name}::{row.case_name}",
            f"input={row.case_input}",
            f"fallback={_dump(py_norm)}",
            f"cpp={_dump(cpp_norm)}",
            f"diff_payload={_dump(_diff_payload(py_norm, cpp_norm))}",
            f"likely_cpp_source={row.cpp_source_path}",
        ]
    )
