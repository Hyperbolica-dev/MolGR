from __future__ import annotations

from importlib import import_module

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod


METHOD_IDS: tuple[str, ...] = (
    "rdkit_determine_bonds",
    "openbabel_read_xyz",
    "cell2mol_v2",
    "molgr_fallback",
    "molgr_cpp",
)


def get_method_registry() -> list[BenchmarkMethod]:
    rdkit_method_cls = import_module(
        "benchmarks.smiles_xyz_benchmark.methods.rdkit_determine_bonds"
    ).RDKitDetermineBondsMethod
    openbabel_method_cls = import_module(
        "benchmarks.smiles_xyz_benchmark.methods.openbabel_read_xyz"
    ).OpenBabelReadXYZMethod
    cell2mol_method_cls = import_module(
        "benchmarks.smiles_xyz_benchmark.methods.cell2mol_v2"
    ).Cell2MolV2Method
    molgr_fallback_method_cls = import_module(
        "benchmarks.smiles_xyz_benchmark.methods.molgr_fallback"
    ).MolGRFallbackMethod
    molgr_cpp_method_cls = import_module(
        "benchmarks.smiles_xyz_benchmark.methods.molgr_cpp"
    ).MolGRCppMethod
    return [
        rdkit_method_cls(),
        openbabel_method_cls(),
        cell2mol_method_cls(),
        molgr_fallback_method_cls(),
        molgr_cpp_method_cls(),
    ]


__all__ = ["METHOD_IDS", "get_method_registry", "BenchmarkMethod"]
