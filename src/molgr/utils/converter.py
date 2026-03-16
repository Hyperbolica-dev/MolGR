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
        obmol.NewBond(bond.begin_atom_idx, bond.end_atom_idx, bond.order)
    obmol.EndModify()
    return pybel.Molecule(obmol)


def mol_data_to_rdkit(mol_data: core.utils.MoleculeData, sanitize: bool = True) -> Chem.Mol:
    """
    Convert MoleculeData to RDKit Mol.
    """
    rwmol = Chem.RWMol()
    for atom in mol_data.atoms:
        atom_idx = rwmol.AddAtom(Chem.Atom(atom.atomic_num))
        rdatom = rwmol.GetAtomWithIdx(atom_idx)
        rdatom.SetNoImplicit(True)
        rdatom.SetFormalCharge(atom.formal_charge)
        rdatom.SetNumRadicalElectrons(atom.radical_num)
    for bond_data in mol_data.bonds:
        rwmol.AddBond(
            bond_data.begin_atom_idx - 1,
            bond_data.end_atom_idx - 1,
            OB_RDKIT_BOND_ORDER_MAPPING[bond_data.order],
        )

    rdmol = Chem.MolFromMolBlock(
        Chem.MolToMolBlock(rwmol),
        sanitize=False,
        removeHs=False,
    )
    if rdmol is None:
        raise ValueError("MolFromMolBlock failed")

    conf = Chem.Conformer(rdmol.GetNumAtoms())
    for atom_idx, atom in enumerate(mol_data.atoms):
        conf.SetAtomPosition(atom_idx, (atom.x, atom.y, atom.z))
    rdmol.RemoveAllConformers()
    rdmol.AddConformer(conf)

    if sanitize:
        Chem.SanitizeMol(rdmol)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.AssignStereochemistryFrom3D(rdmol)
    Chem.AssignCIPLabels(rdmol)
    mol = Chem.Mol(rdmol)
    for bond_idx in range(mol.GetNumBonds()):
        rd_bond = mol.GetBondWithIdx(bond_idx)
        if rd_bond.GetStereo() == Chem.BondStereo.STEREONONE:
            rd_bond.SetBondDir(Chem.BondDir.NONE)
    Chem.Kekulize(mol)
    return mol


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
        rwmol.AddBond(bond[0], bond[1], OB_RDKIT_BOND_ORDER_MAPPING.get(bond[2], Chem.BondType.ZERO))
    for atom_id, (charge, radical) in enumerate(zip(formal_charges, formal_radicals)):
        atom = rwmol.GetAtomWithIdx(atom_id)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        atom.SetNumRadicalElectrons(radical)
    rdmol = Chem.MolFromMolBlock(Chem.MolToMolBlock(rwmol), removeHs=False)
    if sanitize:
        Chem.SanitizeMol(rdmol)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.AssignStereochemistryFrom3D(rdmol)
    Chem.AssignCIPLabels(rdmol)
    Chem.Kekulize(rdmol)
    return rdmol
