from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from . import _core as core


OB_RDKIT_BOND_ORDER_MAPPING = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
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
    rdmol = Chem.RWMol()
    for atom in mol_data.atoms:
        atom_idx = rdmol.AddAtom(Chem.Atom(atom.atomic_num))
        rdatom = rdmol.GetAtomWithIdx(atom_idx)
        rdatom.SetFormalCharge(atom.formal_charge)
        rdatom.SetNumRadicalElectrons(atom.radical_num)
    conf = Chem.Conformer(rdmol.GetNumAtoms())
    for atom in mol_data.atoms:
        conf.SetAtomPosition(atom.atomic_num - 1, (atom.x, atom.y, atom.z))
    rdmol.AddConformer(conf)
    for bond in mol_data.bonds:
        rdmol.AddBond(
            bond.begin_atom_idx - 1, bond.end_atom_idx - 1, OB_RDKIT_BOND_ORDER_MAPPING[bond.order]
        )

    if sanitize:
        Chem.SanitizeMol(rdmol)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.AssignStereochemistryFrom3D(rdmol)
    Chem.AssignCIPLabels(rdmol)
    Chem.Kekulize(rdmol)
    return rdmol.GetMol()
