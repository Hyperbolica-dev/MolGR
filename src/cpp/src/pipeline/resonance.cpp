#include "molgr/pipeline/resonance.h"

#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>

#include <sstream>
#include <unordered_map>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        std::string SmilesFirstToken(const OBMol &mol)
        {
            thread_local OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetOutFormat("smi");
                initialized = true;
            }
            OBMol temp(mol);
            const std::string smi = conv.WriteString(&temp, true);
            std::istringstream iss(smi);
            std::string token;
            iss >> token;
            return token;
        }

        std::vector<OBMol> GetOneStepResonance(const OBMol &mol)
        {
            const std::string pattern = "[*]-,=,:[*]=,#,:[*]";
            OBMol query_mol(mol);
            const auto matches = molgr::utils::FindSmarts(query_mol, pattern);

            std::vector<OBMol> result;
            result.reserve(matches.size());

            for (const auto &idxs : matches)
            {
                if (idxs.size() < 3)
                {
                    continue;
                }

                OBAtom *atom1 = query_mol.GetAtom(idxs[0]);
                OBAtom *atom3 = query_mol.GetAtom(idxs[2]);
                OBBond *bond1 = query_mol.GetBond(idxs[0], idxs[1]);
                OBBond *bond2 = query_mol.GetBond(idxs[1], idxs[2]);

                if (atom1 == nullptr || atom3 == nullptr || bond1 == nullptr || bond2 == nullptr)
                {
                    continue;
                }

                if (atom1->GetSpinMultiplicity() == 1 &&
                    atom3->GetSpinMultiplicity() == 0 &&
                    bond1->GetBondOrder() <= 2 &&
                    bond2->GetBondOrder() >= 2)
                {
                    OBMol new_mol(mol);
                    OBAtom *atom1_clone = new_mol.GetAtom(idxs[0]);
                    OBAtom *atom3_clone = new_mol.GetAtom(idxs[2]);
                    OBBond *bond1_clone = new_mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2_clone = new_mol.GetBond(idxs[1], idxs[2]);
                    if (atom1_clone == nullptr || atom3_clone == nullptr ||
                        bond1_clone == nullptr || bond2_clone == nullptr)
                    {
                        continue;
                    }

                    bond1_clone->SetBondOrder(bond1_clone->GetBondOrder() + 1);
                    bond2_clone->SetBondOrder(bond2_clone->GetBondOrder() - 1);
                    atom1_clone->SetSpinMultiplicity(atom1_clone->GetSpinMultiplicity() - 1);
                    atom3_clone->SetSpinMultiplicity(atom3_clone->GetSpinMultiplicity() + 1);
                    result.push_back(new_mol);
                }
            }

            return result;
        }

        std::vector<OBMol> GetRadicalResonances(const OBMol &mol)
        {
            std::vector<OBMol> resonances;
            std::unordered_map<std::string, std::size_t> key_to_index;

            auto add_or_overwrite = [&resonances, &key_to_index](const OBMol &candidate)
            {
                const std::string key = SmilesFirstToken(candidate);
                const auto it = key_to_index.find(key);
                if (it == key_to_index.end())
                {
                    key_to_index.emplace(key, resonances.size());
                    resonances.emplace_back(candidate);
                    return;
                }

                resonances[it->second] = candidate;
            };

            add_or_overwrite(mol);

            const auto one_step = GetOneStepResonance(mol);
            for (const auto &candidate : one_step)
            {
                add_or_overwrite(candidate);
            }

            for (const auto &temp_mol : one_step)
            {
                const auto two_step = GetOneStepResonance(temp_mol);
                for (const auto &candidate : two_step)
                {
                    add_or_overwrite(candidate);
                }
            }

            return resonances;
        }

        std::pair<OBMol, int> ProcessResonance(const OBMol &mol, int charge)
        {
            OBMol processed(mol);
            int current_charge = charge;

            Eliminate13Dipole(processed, current_charge);
            EliminatePositiveCharges(processed, current_charge);
            EliminateNegativeCharges(processed, current_charge);
            CleanNeighborRadicals(processed);
            CleanResonances(processed);

            return {processed, current_charge};
        }

    }
}
