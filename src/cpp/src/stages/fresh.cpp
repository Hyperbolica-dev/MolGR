#include "molgr/stages/fresh.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/consts.h"

#include <openbabel/atom.h>
#include <openbabel/obfunctions.h>
#include <openbabel/obiter.h>

#include <algorithm>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        const ElementInfo *GetElementInfo(int atomic_num)
        {
            auto it = kNonMetalDict.find(atomic_num);
            if (it != kNonMetalDict.end())
                return &(it->second);
            return nullptr;
        }

        int PythonModulo(int value, int modulus)
        {
            int remainder = value % modulus;
            if (remainder < 0)
                remainder += modulus;
            return remainder;
        }

        int AssignRadicalDots(OBAtom &atom)
        {
            const int total_valence = static_cast<int>(atom.GetTotalValence());
            const int typical_valence = GetTypicalValence(
                atom.GetAtomicNum(),
                total_valence,
                atom.GetFormalCharge());
            return std::max(0, typical_valence - total_valence);
        }

        void AssignChargeRadicalForAtom(OBAtom &atom)
        {
            const int rad = AssignRadicalDots(atom);
            if (rad > 0)
            {
                atom.SetSpinMultiplicity(rad);
                return;
            }

            const int atomic_num = atom.GetAtomicNum();
            const int current_val = atom.GetTotalValence();
            const int charge = atom.GetFormalCharge();
            const ElementInfo *info = GetElementInfo(atomic_num);
            if (!info)
            {
                return;
            }

            if (info->num_outer_electrons == 3 && current_val == 4)
            {
                atom.SetFormalCharge(-1);
                return;
            }

            const int spin = atom.GetSpinMultiplicity();
            const int total_elec = info->num_outer_electrons + current_val + spin - charge;

            const int low_valence_total = PythonModulo(total_elec, 8);
            const int high_valence_total = PythonModulo(info->num_outer_electrons - current_val + spin - charge, 2);

            if (low_valence_total == 0)
            {
                return;
            }

            if (low_valence_total <= high_valence_total)
            {
                atom.SetFormalCharge(low_valence_total);
            }
            else
            {
                const int new_spin = info->num_outer_electrons - current_val + spin - charge;
                atom.SetSpinMultiplicity(new_spin);
            }
        }

        void FreshOmolChargeRadical(OBMol &mol)
        {
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                AssignChargeRadicalForAtom(*atom_iter);
            }
        }
    }
}
