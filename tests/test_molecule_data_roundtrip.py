# pyright: reportMissingImports=false

import pytest


pytest.importorskip("openbabel")


def _get_ptr(obmol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def test_molecule_data_roundtrip_obmol_preserves_fields() -> None:
    from openbabel import openbabel as ob

    from molgr import _core  # type: ignore

    obmol = ob.OBMol()
    obmol.BeginModify()

    atom1 = obmol.NewAtom()
    atom1.SetAtomicNum(6)
    atom1.SetFormalCharge(-1)
    atom1.SetSpinMultiplicity(2)
    atom1.SetVector(0.1, 0.2, 0.3)

    atom2 = obmol.NewAtom()
    atom2.SetAtomicNum(8)
    atom2.SetFormalCharge(1)
    atom2.SetSpinMultiplicity(1)
    atom2.SetVector(1.1, 1.2, 1.3)

    atom3 = obmol.NewAtom()
    atom3.SetAtomicNum(7)
    atom3.SetFormalCharge(0)
    atom3.SetSpinMultiplicity(3)
    atom3.SetVector(-2.1, -2.2, -2.3)

    obmol.AddBond(1, 2, 2)
    obmol.AddBond(2, 3, 1)
    obmol.EndModify()

    md_in = _core.utils.extract_molecule_data(_get_ptr(obmol))
    roundtrip_ptr = _core.utils.molecule_data_to_obmol_ptr(md_in)

    try:
        md_out = _core.utils.extract_molecule_data(roundtrip_ptr)
    finally:
        _core.free_obmol_ptr(roundtrip_ptr)

    assert len(md_out.atoms) == len(md_in.atoms)
    assert len(md_out.bonds) == len(md_in.bonds)
    assert md_out.total_charge == md_in.total_charge
    assert md_out.total_radical_num == md_in.total_radical_num

    for atom_in, atom_out in zip(md_in.atoms, md_out.atoms):
        assert atom_out.atomic_num == atom_in.atomic_num
        assert atom_out.formal_charge == atom_in.formal_charge
        assert atom_out.radical_num == atom_in.radical_num
        assert atom_out.x == pytest.approx(atom_in.x)
        assert atom_out.y == pytest.approx(atom_in.y)
        assert atom_out.z == pytest.approx(atom_in.z)

    for bond_in, bond_out in zip(md_in.bonds, md_out.bonds):
        assert bond_out.begin_atom_idx == bond_in.begin_atom_idx
        assert bond_out.end_atom_idx == bond_in.end_atom_idx
        assert bond_out.order == bond_in.order
        assert bond_out.begin_atom_idx >= 1
        assert bond_out.end_atom_idx >= 1
