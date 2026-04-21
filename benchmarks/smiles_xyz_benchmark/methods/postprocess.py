from __future__ import annotations

from importlib import import_module

from rdkit import Chem


def finalize_rdmol_with_dative_bonds(rdmol: Chem.Mol) -> tuple[Chem.Mol, str]:
    make_dative_bond = import_module("molgr.utils.post_process").make_dative_bond
    mol_with_dative = make_dative_bond(rdmol)
    mol_no_h = Chem.RemoveHs(mol_with_dative)
    predicted_smiles = Chem.MolToSmiles(mol_no_h, canonical=True, isomericSmiles=True)
    return mol_no_h, predicted_smiles


__all__ = ["finalize_rdmol_with_dative_bonds"]
