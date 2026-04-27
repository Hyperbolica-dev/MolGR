from typing import List, cast

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molgr import _core as core


OB_RDKIT_BOND_ORDER_MAPPING = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
    5: Chem.BondType.AROMATIC,
}


def mol_data_to_pybel(mol_data: core.utils.MoleculeData) -> pybel.Molecule:
    """
    Convert MoleculeData to Pybel Molecule.
    """
    obmol = ob.OBMol()
    obmol.BeginModify()
    for atom in mol_data.atoms:
        obatom: ob.OBAtom = obmol.NewAtom()
        obatom.SetAtomicNum(atom.atomic_num)
        obatom.SetFormalCharge(atom.formal_charge)
        obatom.SetSpinMultiplicity(atom.radical_num)
        obatom.SetVector(atom.x, atom.y, atom.z)
    for bond in mol_data.bonds:
        obmol.AddBond(bond.begin_atom_idx, bond.end_atom_idx, bond.order)
    obmol.EndModify()
    return pybel.Molecule(obmol)


def mol_data_to_rdkit(mol_data: core.utils.MoleculeData, sanitize: bool = True) -> Chem.Mol:
    """
    Convert MoleculeData to RDKit Mol.
    """
    # Keep the C++ and Python backend boundary identical by routing MoleculeData
    # through the same OpenBabel -> RDKit conversion path used by the fallback.
    return pybel_to_rdmol(mol_data_to_pybel(mol_data), sanitize=sanitize)


def pybel_to_rdmol(omol: pybel.Molecule, sanitize: bool = True) -> Chem.Mol:
    """
    Convert Pybel Molecule to RDKit Mol.
    """
    bonds = [
        (
            cast(ob.OBBond, bond).GetBeginAtomIdx() - 1,
            cast(ob.OBBond, bond).GetEndAtomIdx() - 1,
            cast(ob.OBBond, bond).GetBondOrder(),
        )
        for bond in ob.OBMolBondIter(omol.OBMol)
    ]
    formal_charges: List[int] = [
        cast(ob.OBAtom, atom).GetFormalCharge() for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    formal_radicals: List[int] = [
        cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(omol.write("xyz")))
    for bond in bonds:
        rwmol.AddBond(
            bond[0], bond[1], OB_RDKIT_BOND_ORDER_MAPPING.get(bond[2], Chem.BondType.ZERO)
        )
    for atom_id, (charge, radical) in enumerate(zip(formal_charges, formal_radicals)):
        atom = rwmol.GetAtomWithIdx(atom_id)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        atom.SetNumRadicalElectrons(radical)
    rdmol = Chem.MolFromMolBlock(Chem.MolToMolBlock(rwmol), removeHs=False)
    rdmol.UpdatePropertyCache(strict=False)
    if sanitize:
        Chem.SanitizeMol(rdmol)
    else:
        rdmol.UpdatePropertyCache(strict=False)
    rdmol.UpdatePropertyCache(strict=False)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.AssignStereochemistryFrom3D(rdmol)
    Chem.AssignCIPLabels(rdmol)
    Chem.Kekulize(rdmol)
    return rdmol
