"""
Scoring functions
"""
from __future__ import annotations
import typing
__all__: list[str] = ['omol_score_from_ptr', 'test_deviation_score', 'test_physchem_penalty', 'test_symmetry_penalty', 'test_total_score']
def omol_score_from_ptr(mol_ptr: typing.SupportsInt) -> float:
    """
            Calculate total OMolScore using a memory pointer to an OpenBabel::OBMol.
            This allows compatibility with SWIG-wrapped OpenBabel objects.
            
            Parameters:
                mol_ptr (int): The memory address of the OBMol object (use `int(mol.this)` in Python).
                
            Returns:
                float: Total score.
    """
def test_deviation_score(xyz_block: str, atom_idx: typing.SupportsInt) -> float:
    """
    Calculate geometry deviation for atom (1-based index) from XYZ (For Testing)
    """
def test_physchem_penalty(smiles: str) -> float:
    """
    Calculate PhysChem penalty from SMILES (For Testing)
    """
def test_symmetry_penalty(smiles: str) -> float:
    """
    Calculate symmetry penalty from SMILES (For Testing)
    """
def test_total_score(xyz_block: str) -> float:
    """
    Calculate total OMolScore from XYZ block (For Testing)
    """
