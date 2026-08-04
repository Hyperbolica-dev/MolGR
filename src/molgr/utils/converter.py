from typing import List, cast

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molgr import _core as core
from molgr.fallback.utils.consts import NON_METAL_DICT
from molgr.fallback.utils.electrons import (
    LONE_PAIR_COUNT_PROP,
    UNRESOLVED_TWO_ELECTRON_CENTER_PROP,
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)


OB_RDKIT_BOND_ORDER_MAPPING = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
    5: Chem.BondType.AROMATIC,
}

METAL_UNPAIRED_ELECTRONS_PROP = "MOLGR_METAL_UNPAIRED_ELECTRONS"


def _sanitize_preserving_explicit_bonds(rdmol: Chem.Mol, *, kekulize: bool) -> None:
    """Sanitize a graph without rejecting a valid explicit-bond representation.

    Reconstruction can leave an aromaticity assignment that RDKit cannot
    Kekulize after charge/bond edits. In that case retain the explicit bond
    orders supplied by Open Babel and skip only the failing Kekule step.
    """

    if not kekulize:
        Chem.SanitizeMol(rdmol)
        return
    try:
        Chem.SanitizeMol(rdmol)
        Chem.Kekulize(rdmol)
        return
    except Chem.KekulizeException:
        rdmol.ClearComputedProps()
        sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        Chem.SanitizeMol(rdmol, sanitizeOps=sanitize_ops)
        Chem.KekulizeIfPossible(rdmol, clearAromaticFlags=False)


def _is_metal_atomic_num(atomic_num: int) -> bool:
    return atomic_num not in NON_METAL_DICT


def _set_rdkit_unpaired_electrons(atom: Chem.Atom, unpaired_electrons: int) -> None:
    """Map real unpaired electrons to RDKit, using a side property for metals."""

    if _is_metal_atomic_num(int(atom.GetAtomicNum())):
        atom.SetIntProp(METAL_UNPAIRED_ELECTRONS_PROP, int(unpaired_electrons))
        atom.SetNumRadicalElectrons(0)
        return
    atom.SetNumRadicalElectrons(int(unpaired_electrons))


def _set_rdkit_lone_pair_count(atom: Chem.Atom, lone_pair_count: int) -> None:
    """Persist active lone pairs as an RDKit atom property without changing spin."""

    if lone_pair_count:
        atom.SetIntProp(LONE_PAIR_COUNT_PROP, int(lone_pair_count))
    elif atom.HasProp(LONE_PAIR_COUNT_PROP):
        atom.ClearProp(LONE_PAIR_COUNT_PROP)


def _set_rdkit_unresolved_two_electron_center(atom: Chem.Atom, value: bool) -> None:
    """Persist deferred two-electron occupancy without selecting singlet/triplet."""

    if value:
        atom.SetBoolProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP, True)
    elif atom.HasProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP):
        atom.ClearProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP)


def _finalize_metal_unpaired_electrons(rdmol: Chem.Mol) -> None:
    has_metal = False
    for atom in rdmol.GetAtoms():  # pyright: ignore[reportCallIssue]
        if not _is_metal_atomic_num(int(atom.GetAtomicNum())):
            continue
        has_metal = True
        atom.SetNumRadicalElectrons(0)
    if has_metal:
        Chem.CreateAtomIntPropertyList(rdmol, METAL_UNPAIRED_ELECTRONS_PROP)
    if any(
        atom.HasProp(LONE_PAIR_COUNT_PROP)
        for atom in rdmol.GetAtoms()  # pyright: ignore[reportCallIssue]
    ):
        Chem.CreateAtomIntPropertyList(rdmol, LONE_PAIR_COUNT_PROP)
    if any(
        atom.HasProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP)
        for atom in rdmol.GetAtoms()  # pyright: ignore[reportCallIssue]
    ):
        Chem.CreateAtomBoolPropertyList(rdmol, UNRESOLVED_TWO_ELECTRON_CENTER_PROP)
    rdmol.UpdatePropertyCache(strict=False)


def _restore_rdkit_unpaired_electrons(
    rdmol: Chem.Mol,
    unpaired_electrons: List[int],
) -> None:
    """Restore MolGR electron assignments overwritten by RDKit sanitization."""

    for atom, count in zip(
        rdmol.GetAtoms(),  # pyright: ignore[reportCallIssue]
        unpaired_electrons,
    ):
        _set_rdkit_unpaired_electrons(atom, count)


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


def get_atom_lone_pair_count(atom: Chem.Atom) -> int:
    """Return MolGR's explicit lone-pair count for an RDKit atom."""

    if not atom.HasProp(LONE_PAIR_COUNT_PROP):
        return 0
    try:
        return int(atom.GetIntProp(LONE_PAIR_COUNT_PROP))
    except (RuntimeError, TypeError, ValueError):
        return int(atom.GetProp(LONE_PAIR_COUNT_PROP))


def has_atom_unresolved_two_electron_center(atom: Chem.Atom) -> bool:
    """Return whether an RDKit atom carries MolGR's unresolved-center marker."""

    if not atom.HasProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP):
        return False
    try:
        return bool(atom.GetBoolProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP))
    except (RuntimeError, TypeError, ValueError):
        return atom.GetProp(UNRESOLVED_TWO_ELECTRON_CENTER_PROP).lower() in {
            "1",
            "true",
        }


def mol_data_to_pybel(mol_data: core.utils.MoleculeData) -> pybel.Molecule:
    """Copy all three electron classifications from MoleculeData to Open Babel.

    ``radical_num`` remains real unpaired electrons; lone-pair and unresolved
    fields use independent generic data and are not inferred from one another.
    """
    obmol = ob.OBMol()
    obmol.BeginModify()
    for atom in mol_data.atoms:
        obatom: ob.OBAtom = obmol.NewAtom()
        obatom.SetAtomicNum(atom.atomic_num)
        obatom.SetFormalCharge(atom.formal_charge)
        set_unpaired_electron_count(obatom, atom.radical_num)
        set_lone_pair_count(obatom, atom.lone_pair_count)
        set_unresolved_two_electron_center(
            obatom,
            atom.unresolved_two_electron_center,
        )
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
    """Convert MoleculeData while preserving independent electron classifications.

    Nonmetal unpaired electrons use RDKit radicals, metal unpaired electrons use
    ``MOLGR_METAL_UNPAIRED_ELECTRONS``, and active lone-pair/unresolved fields use
    separate atom properties. Sanitization is not allowed to redefine these
    stored MolGR assignments.
    """
    unpaired_electrons = [int(atom.radical_num) for atom in mol_data.atoms]
    rwmol = Chem.RWMol()
    for atom in mol_data.atoms:
        atom_id = rwmol.AddAtom(Chem.Atom(atom.atomic_num))
        rd_atom = rwmol.GetAtomWithIdx(atom_id)
        rd_atom.SetNoImplicit(True)
        rd_atom.SetFormalCharge(atom.formal_charge)
        _set_rdkit_unpaired_electrons(rd_atom, atom.radical_num)
        _set_rdkit_lone_pair_count(rd_atom, atom.lone_pair_count)
        _set_rdkit_unresolved_two_electron_center(
            rd_atom,
            atom.unresolved_two_electron_center,
        )
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
        _sanitize_preserving_explicit_bonds(rdmol, kekulize=kekulize)
    elif kekulize:
        Chem.KekulizeIfPossible(rdmol)
    _restore_rdkit_unpaired_electrons(rdmol, unpaired_electrons)
    _finalize_metal_unpaired_electrons(rdmol)
    return rdmol


def pybel_to_rdmol(
    omol: pybel.Molecule,
    sanitize: bool = True,
    kekulize: bool = True,
) -> Chem.Mol:
    """Convert Open Babel topology while preserving all MolGR electron fields.

    Real unpaired electrons, active lone pairs, and unresolved centers are read
    before RDKit sanitization and restored independently afterwards.
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
        get_unpaired_electron_count(cast(ob.OBAtom, atom)) for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    lone_pair_counts: List[int] = [
        get_lone_pair_count(cast(ob.OBAtom, atom)) for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    unresolved_two_electron_centers: List[bool] = [
        has_unresolved_two_electron_center(cast(ob.OBAtom, atom))
        for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(omol.write("xyz")))
    for bond in bonds:
        rwmol.AddBond(
            bond[0], bond[1], OB_RDKIT_BOND_ORDER_MAPPING.get(bond[2], Chem.BondType.ZERO)
        )
    for atom_id, (charge, radical, lone_pair_count, unresolved_center) in enumerate(
        zip(
            formal_charges,
            formal_radicals,
            lone_pair_counts,
            unresolved_two_electron_centers,
        )
    ):
        atom = rwmol.GetAtomWithIdx(atom_id)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        _set_rdkit_unpaired_electrons(atom, radical)
        _set_rdkit_lone_pair_count(atom, lone_pair_count)
        _set_rdkit_unresolved_two_electron_center(atom, unresolved_center)
    rdmol = rwmol.GetMol()
    rdmol.UpdatePropertyCache(strict=False)
    if sanitize:
        _sanitize_preserving_explicit_bonds(rdmol, kekulize=kekulize)
    elif kekulize:
        Chem.KekulizeIfPossible(rdmol)
    _restore_rdkit_unpaired_electrons(rdmol, formal_radicals)
    _finalize_metal_unpaired_electrons(rdmol)
    return rdmol


__all__ = [
    "METAL_UNPAIRED_ELECTRONS_PROP",
    "get_atom_lone_pair_count",
    "get_atom_unpaired_electrons",
    "has_atom_unresolved_two_electron_center",
    "mol_data_to_pybel",
    "mol_data_to_rdkit",
    "pybel_to_rdmol",
]
