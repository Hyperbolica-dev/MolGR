#include "molgr/utils/conversions.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

namespace molgr
{
    namespace utils
    {
        OpenBabel::OBMol MolFromMoleculeData(const MoleculeData &data)
        {
            OpenBabel::OBMol mol;
            mol.BeginModify();

            for (const AtomData &atom_data : data.atoms)
            {
                OpenBabel::OBAtom *atom = mol.NewAtom();
                atom->SetAtomicNum(atom_data.atomic_num);
                atom->SetFormalCharge(atom_data.formal_charge);
                atom->SetSpinMultiplicity(atom_data.radical_num);
                atom->SetVector(atom_data.x, atom_data.y, atom_data.z);
            }

            for (const BondData &bond_data : data.bonds)
            {
                mol.AddBond(bond_data.begin_atom_idx, bond_data.end_atom_idx, bond_data.order);
            }

            mol.EndModify();
            return mol;
        }

        MoleculeData MoleculeDataFromOBMol(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol mol_copy(mol);
            MoleculeData data;
            data.total_charge = 0;
            data.total_radical_num = 0;

            data.atoms.reserve(mol_copy.NumAtoms());
            FOR_ATOMS_OF_MOL(atom, mol_copy)
            {
                AtomData atom_data;
                atom_data.atomic_num = atom->GetAtomicNum();
                atom_data.formal_charge = atom->GetFormalCharge();
                atom_data.radical_num = atom->GetSpinMultiplicity();
                atom_data.x = atom->GetX();
                atom_data.y = atom->GetY();
                atom_data.z = atom->GetZ();
                data.atoms.push_back(atom_data);

                data.total_charge += atom_data.formal_charge;
                data.total_radical_num += atom_data.radical_num;
            }

            data.bonds.reserve(mol_copy.NumBonds());
            FOR_BONDS_OF_MOL(bond, mol_copy)
            {
                BondData bond_data;
                bond_data.begin_atom_idx = bond->GetBeginAtom()->GetIdx();
                bond_data.end_atom_idx = bond->GetEndAtom()->GetIdx();
                bond_data.order = bond->GetBondOrder();
                data.bonds.push_back(bond_data);
            }

            return data;
        }
    }
}
