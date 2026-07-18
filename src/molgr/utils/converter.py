from typing import List, cast

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molgr import _core as core
from molgr.fallback.utils.consts import NON_METAL_DICT


OB_RDKIT_BOND_ORDER_MAPPING = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
    5: Chem.BondType.AROMATIC,
}

METAL_UNPAIRED_ELECTRONS_PROP = "MOLGR_METAL_UNPAIRED_ELECTRONS"


def _is_metal_atomic_num(atomic_num: int) -> bool:
    return atomic_num not in NON_METAL_DICT


def _set_rdkit_unpaired_electrons(atom: Chem.Atom, unpaired_electrons: int) -> None:
    if _is_metal_atomic_num(int(atom.GetAtomicNum())):
        atom.SetIntProp(METAL_UNPAIRED_ELECTRONS_PROP, int(unpaired_electrons))
        atom.SetNumRadicalElectrons(0)
        return
    atom.SetNumRadicalElectrons(int(unpaired_electrons))


def _finalize_metal_unpaired_electrons(rdmol: Chem.Mol) -> None:
    has_metal = False
    for atom in rdmol.GetAtoms():  # pyright: ignore[reportCallIssue]
        if not _is_metal_atomic_num(int(atom.GetAtomicNum())):
            continue
        has_metal = True
        atom.SetNumRadicalElectrons(0)
    if has_metal:
        Chem.CreateAtomIntPropertyList(rdmol, METAL_UNPAIRED_ELECTRONS_PROP)
    rdmol.UpdatePropertyCache(strict=False)


def get_atom_unpaired_electrons(atom: Chem.Atom) -> int:
    """Return MolGR's unpaired-electron count for an RDKit atom."""

    if _is_metal_atomic_num(int(atom.GetAtomicNum())) and atom.HasProp(
        METAL_UNPAIRED_ELECTRONS_PROP
    ):
        try:
            return int(atom.GetIntProp(METAL_UNPAIRED_ELECTRONS_PROP))
        except (RuntimeError, TypeError, ValueError):
            return int(atom.GetProp(METAL_UNPAIRED_ELECTRONS_PROP))
    return int(atom.GetNumRadicalElectrons())


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


def mol_data_to_rdkit(
    mol_data: core.utils.MoleculeData,
    sanitize: bool = True,
    kekulize: bool = True,
) -> Chem.Mol:
    """
    Convert MoleculeData to RDKit Mol.
    """
    rwmol = Chem.RWMol()
    for atom in mol_data.atoms:
        atom_id = rwmol.AddAtom(Chem.Atom(atom.atomic_num))
        rd_atom = rwmol.GetAtomWithIdx(atom_id)
        rd_atom.SetNoImplicit(True)
        rd_atom.SetFormalCharge(atom.formal_charge)
        _set_rdkit_unpaired_electrons(rd_atom, atom.radical_num)
    for bond_data in mol_data.bonds:
        rwmol.AddBond(
            bond_data.begin_atom_idx - 1,
            bond_data.end_atom_idx - 1,
            OB_RDKIT_BOND_ORDER_MAPPING.get(bond_data.order, Chem.BondType.ZERO),
        )
    rdmol = rwmol.GetMol()
    conf = Chem.Conformer(rdmol.GetNumAtoms())
    for atom_idx, atom in enumerate(mol_data.atoms):
        conf.SetAtomPosition(atom_idx, (atom.x, atom.y, atom.z))
    rdmol.RemoveAllConformers()
    rdmol.AddConformer(conf)

    if sanitize:
        Chem.SanitizeMol(rdmol)
    if kekulize:
        Chem.Kekulize(rdmol)
    _finalize_metal_unpaired_electrons(rdmol)
    return rdmol


def pybel_to_rdmol(
    omol: pybel.Molecule,
    sanitize: bool = True,
    kekulize: bool = True,
) -> Chem.Mol:
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
        _set_rdkit_unpaired_electrons(atom, radical)
    rdmol = rwmol.GetMol()
    rdmol.UpdatePropertyCache(strict=False)
    if sanitize:
        Chem.SanitizeMol(rdmol)
    if kekulize:
        Chem.Kekulize(rdmol)
    _finalize_metal_unpaired_electrons(rdmol)
    return rdmol


__all__ = [
    "METAL_UNPAIRED_ELECTRONS_PROP",
    "get_atom_unpaired_electrons",
    "mol_data_to_pybel",
    "mol_data_to_rdkit",
    "pybel_to_rdmol",
]
