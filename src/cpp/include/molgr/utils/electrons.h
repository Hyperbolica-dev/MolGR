#pragma once

#include <openbabel/atom.h>
#include <openbabel/generic.h>

#include <stdexcept>
#include <string>

namespace molgr
{
    namespace utils
    {
        inline constexpr const char *kUnpairedElectronCountProp =
            "MOLGR_UNPAIRED_ELECTRON_COUNT";
        inline constexpr const char *kLonePairCountProp = "MOLGR_LONE_PAIR_COUNT";
        inline constexpr const char *kUnresolvedTwoElectronCenterProp =
            "MOLGR_UNRESOLVED_TWO_ELECTRON_CENTER";

        inline int GetAtomIntegerData(
            const OpenBabel::OBAtom &atom,
            const char *attribute)
        {
            auto &mutable_atom = const_cast<OpenBabel::OBAtom &>(atom);
            OpenBabel::OBGenericData *data = mutable_atom.GetData(attribute);
            if (data == nullptr)
            {
                return 0;
            }
            if (const auto *integer = dynamic_cast<const OpenBabel::OBPairInteger *>(data))
            {
                return integer->GetGenericValue();
            }
            if (const auto *pair = dynamic_cast<const OpenBabel::OBPairData *>(data))
            {
                try
                {
                    return std::stoi(pair->GetValue());
                }
                catch (const std::exception &)
                {
                    return 0;
                }
            }
            return 0;
        }

        inline void SetAtomIntegerData(
            OpenBabel::OBAtom &atom,
            const char *attribute,
            int value)
        {
            if (value < 0)
            {
                throw std::invalid_argument("atom electron count must be nonnegative");
            }
            atom.DeleteData(attribute);
            if (value == 0)
            {
                return;
            }
            auto *data = new OpenBabel::OBPairData();
            data->SetAttribute(attribute);
            data->SetValue(std::to_string(value));
            data->SetOrigin(OpenBabel::local);
            atom.SetData(data);
        }

        inline int GetUnpairedElectronCount(const OpenBabel::OBAtom &atom)
        {
            return GetAtomIntegerData(atom, kUnpairedElectronCountProp);
        }

        // Store only real unpaired electrons; active lone pairs must use their
        // independent field and never be folded into this count.
        inline void SetUnpairedElectronCount(OpenBabel::OBAtom &atom, int value)
        {
            SetAtomIntegerData(atom, kUnpairedElectronCountProp, value);
        }

        inline int GetLonePairCount(const OpenBabel::OBAtom &atom)
        {
            return GetAtomIntegerData(atom, kLonePairCountProp);
        }

        // Store reconstruction-active lone pairs, not the full Lewis lone-pair
        // inventory implied by ordinary valence and formal charge.
        inline void SetLonePairCount(OpenBabel::OBAtom &atom, int value)
        {
            SetAtomIntegerData(atom, kLonePairCountProp, value);
        }

        inline bool HasUnresolvedTwoElectronCenter(const OpenBabel::OBAtom &atom)
        {
            return GetAtomIntegerData(atom, kUnresolvedTwoElectronCenterProp) != 0;
        }

        // Mark a deferred two-electron occupancy without selecting singlet (0,1)
        // or triplet (2,0) classification.
        inline void SetUnresolvedTwoElectronCenter(OpenBabel::OBAtom &atom, bool value)
        {
            SetAtomIntegerData(atom, kUnresolvedTwoElectronCenterProp, value ? 1 : 0);
        }

        // Copy the three independent classifications verbatim; never reconstruct
        // lone-pair or unresolved state from unpaired-electron parity.
        inline void CopyElectronData(
            const OpenBabel::OBAtom &source,
            OpenBabel::OBAtom &target)
        {
            SetUnpairedElectronCount(target, GetUnpairedElectronCount(source));
            SetLonePairCount(target, GetLonePairCount(source));
            SetUnresolvedTwoElectronCenter(
                target,
                HasUnresolvedTwoElectronCenter(source));
        }
    }
}
