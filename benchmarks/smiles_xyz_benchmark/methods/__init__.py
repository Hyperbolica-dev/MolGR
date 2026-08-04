from __future__ import annotations

from importlib import import_module

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod


METHOD_IDS: tuple[str, ...] = (
    "rdkit_determine_bonds",
    "openbabel_read_xyz",
    "cell2mol_v2",
    "molgr_fallback",
    "molgr_cpp",
    "xyzgraph_cheminf_full",
)

_METHOD_CLASS_IMPORTS: dict[str, tuple[str, str]] = {
    "rdkit_determine_bonds": (
        "benchmarks.smiles_xyz_benchmark.methods.rdkit_determine_bonds",
        "RDKitDetermineBondsMethod",
    ),
    "openbabel_read_xyz": (
        "benchmarks.smiles_xyz_benchmark.methods.openbabel_read_xyz",
        "OpenBabelReadXYZMethod",
    ),
    "cell2mol_v2": (
        "benchmarks.smiles_xyz_benchmark.methods.cell2mol_v2",
        "Cell2MolV2Method",
    ),
    "molgr_fallback": (
        "benchmarks.smiles_xyz_benchmark.methods.molgr_fallback",
        "MolGRFallbackMethod",
    ),
    "molgr_cpp": (
        "benchmarks.smiles_xyz_benchmark.methods.molgr_cpp",
        "MolGRCppMethod",
    ),
    "xyzgraph_cheminf_full": (
        "benchmarks.smiles_xyz_benchmark.methods.xyzgraph_cheminf_full",
        "XYZGraphCheminfFullMethod",
    ),
}


def get_method_registry(method_ids: tuple[str, ...] | None = None) -> list[BenchmarkMethod]:
    selected_ids = METHOD_IDS if method_ids is None else method_ids
    registry: list[BenchmarkMethod] = []
    for method_id in selected_ids:
        if method_id not in _METHOD_CLASS_IMPORTS:
            continue
        module_name, class_name = _METHOD_CLASS_IMPORTS[method_id]
        method_cls = getattr(import_module(module_name), class_name)
        registry.append(method_cls())
    return registry


__all__ = ["METHOD_IDS", "get_method_registry", "BenchmarkMethod"]
