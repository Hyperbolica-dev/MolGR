#include "molgr/utils/conversions.h"
#include "molgr/utils/electrons.h"
#include "molgr/vendor/openbabel_threading.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"

namespace molgr
{
    namespace utils
    {
        // Electron bookkeeping: despite the historical name, clone topology plus
        // all three MolGR atom electron fields verbatim; do not infer one from another.
        OpenBabel::OBMol CloneMolTopologyOnly(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_source = const_cast<OpenBabel::OBMol &>(mol);
            const bool source_has_hybridization =
                mutable_source.HasHybridizationPerceived();
            OpenBabel::OBMol clone;
            clone.BeginModify();
            clone.ReserveAtoms(mol.NumAtoms());
            clone.SetDimension(mol.GetDimension());
            clone.SetTitle(mol.GetTitle());

            FOR_ATOMS_OF_MOL(atom, const_cast<OpenBabel::OBMol &>(mol))
            {
                OpenBabel::OBAtom *new_atom = clone.NewAtom();
                new_atom->SetAtomicNum(atom->GetAtomicNum());
                new_atom->SetFormalCharge(atom->GetFormalCharge());
                CopyElectronData(*atom, *new_atom);
                if (source_has_hybridization)
                {
                    new_atom->SetHyb(atom->GetHyb());
                }
                new_atom->SetAromatic(atom->IsAromatic());
                new_atom->SetInRing(atom->IsInRing());
                new_atom->SetVector(atom->GetX(), atom->GetY(), atom->GetZ());
            }

            FOR_BONDS_OF_MOL(bond, const_cast<OpenBabel::OBMol &>(mol))
            {
                const int begin_idx = bond->GetBeginAtom()->GetIdx();
                const int end_idx = bond->GetEndAtom()->GetIdx();
                clone.AddBond(
                    begin_idx,
                    end_idx,
                    bond->GetBondOrder());
                OpenBabel::OBBond *new_bond = clone.GetBond(
                    begin_idx,
                    end_idx);
                if (new_bond != nullptr)
                {
                    const bool source_bond_aromatic =
                        (bond->GetFlags() & OB_AROMATIC_BOND) != 0;
                    const bool source_bond_in_ring =
                        (bond->GetFlags() & OB_RING_BOND) != 0;
                    const bool source_bond_closure =
                        (bond->GetFlags() & OB_CLOSURE_BOND) != 0;
                    new_bond->SetAromatic(source_bond_aromatic);
                    new_bond->SetInRing(source_bond_in_ring);
                    new_bond->SetClosure(source_bond_closure);
                }
            }

            clone.EndModify();
            clone.SetRingAtomsAndBondsPerceived(
                mutable_source.HasRingAtomsAndBondsPerceived());
            clone.SetClosureBondsPerceived(mutable_source.HasClosureBondsPerceived());
            clone.SetSSSRPerceived(false);
            clone.SetLSSRPerceived(mutable_source.HasLSSRPerceived());
            clone.SetAromaticPerceived(mutable_source.HasAromaticPerceived());
            clone.SetHybridizationPerceived(source_has_hybridization);
            clone.SetAtomTypesPerceived(mutable_source.HasAtomTypesPerceived());
            return clone;
        }

        // Electron bookkeeping: restore radical_num as real unpaired electrons and
        // restore active lone-pair/unresolved fields independently.
        OpenBabel::OBMol MolFromMoleculeData(const MoleculeData &data)
        {
            OpenBabel::OBMol mol;
            mol.BeginModify();

            for (const AtomData &atom_data : data.atoms)
            {
                OpenBabel::OBAtom *atom = mol.NewAtom();
                atom->SetAtomicNum(atom_data.atomic_num);
                atom->SetFormalCharge(atom_data.formal_charge);
                SetUnpairedElectronCount(*atom, atom_data.radical_num);
                SetLonePairCount(*atom, atom_data.lone_pair_count);
                SetUnresolvedTwoElectronCenter(
                    *atom,
                    atom_data.unresolved_two_electron_center);
                atom->SetHyb(atom_data.hybridization);
                atom->SetVector(atom_data.x, atom_data.y, atom_data.z);
            }

            for (const BondData &bond_data : data.bonds)
            {
                mol.AddBond(bond_data.begin_atom_idx, bond_data.end_atom_idx, bond_data.order);
                if (bond_data.aromatic || bond_data.order == 5)
                {
                    OpenBabel::OBBond *bond = mol.GetBond(
                        bond_data.begin_atom_idx,
                        bond_data.end_atom_idx);
                    if (bond != nullptr)
                    {
                        bond->SetAromatic();
                    }
                }
            }

            mol.EndModify();
            return mol;
        }

        // Electron bookkeeping: serialize all three classifications independently;
        // only real unpaired electrons contribute to total_radical_num.
        MoleculeData MoleculeDataFromOBMol(const OpenBabel::OBMol &mol)
        {
            MoleculeData data;
            data.total_charge = 0;
            data.total_radical_num = 0;

            OpenBabel::OBMol &mutable_mol = const_cast<OpenBabel::OBMol &>(mol);
            molgr::vendor::openbabel_threading::EnsureHybridizationPerceived(mutable_mol);
            data.atoms.reserve(mol.NumAtoms());
            FOR_ATOMS_OF_MOL(atom, mutable_mol)
            {
                AtomData atom_data;
                atom_data.atomic_num = atom->GetAtomicNum();
                atom_data.formal_charge = atom->GetFormalCharge();
                atom_data.radical_num = GetUnpairedElectronCount(*atom);
                atom_data.lone_pair_count = GetLonePairCount(*atom);
                atom_data.unresolved_two_electron_center =
                    HasUnresolvedTwoElectronCenter(*atom);
                atom_data.hybridization = atom->GetHyb();
                atom_data.x = atom->GetX();
                atom_data.y = atom->GetY();
                atom_data.z = atom->GetZ();
                data.atoms.push_back(atom_data);

                data.total_charge += atom_data.formal_charge;
                data.total_radical_num += atom_data.radical_num;
            }

            data.bonds.reserve(mol.NumBonds());
            FOR_BONDS_OF_MOL(bond, mutable_mol)
            {
                BondData bond_data;
                bond_data.begin_atom_idx = bond->GetBeginAtom()->GetIdx();
                bond_data.end_atom_idx = bond->GetEndAtom()->GetIdx();
                bond_data.order = bond->GetBondOrder();
                bond_data.aromatic = (bond->GetFlags() & OB_AROMATIC_BOND) != 0;
                data.bonds.push_back(bond_data);
            }

            return data;
        }
    }
}
