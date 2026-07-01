#include "molgr/stages/break_bond.h"

#include "molgr/utils/logger.h"
#include "molgr/utils/smarts.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        bool BreakDeformedEne(OBMol &mol, int allowed_charge, int allowed_radical, double tolerance)
        {
            bool hit = false;
            auto current_total_radical = [&]()
            {
                int sum = 0;
                FOR_ATOMS_OF_MOL(a, mol)
                sum += a->GetSpinMultiplicity();
                return sum;
            };

            struct EneCandidate
            {
                int idx1;
                int idx2;
                double deviation;
            };
            std::vector<EneCandidate> candidates;

            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::BREAK_DEFORMED_ENE_A);
                for (const auto &idxs : matches)
                {
                    OBBond *bond = mol.GetBond(idxs[1], idxs[2]);
                    if (!bond || bond->IsRotor() || bond->GetBondOrder() == 1)
                        continue;

                    double torsion = std::abs(mol.GetTorsion(idxs[0], idxs[1], idxs[2], idxs[3]));
                    double deviation = std::min(torsion, 180.0 - torsion);
                    if (deviation > tolerance)
                    {
                        candidates.push_back(EneCandidate{idxs[1], idxs[2], deviation});
                    }
                }
            }

            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::BREAK_DEFORMED_ENE_B);
                for (const auto &idxs : matches)
                {
                    OBBond *bond = mol.GetBond(idxs[1], idxs[2]);
                    if (!bond || bond->IsRotor() || bond->GetBondOrder() == 1)
                        continue;

                    double torsion = std::abs(mol.GetTorsion(idxs[0], idxs[1], idxs[2], idxs[3]));
                    double deviation = std::min(torsion, 180.0 - torsion);
                    if (deviation > tolerance)
                    {
                        candidates.push_back(EneCandidate{idxs[1], idxs[2], deviation});
                    }
                }
            }

            std::sort(candidates.begin(), candidates.end(), [](const EneCandidate &a, const EneCandidate &b)
                      { return a.deviation > b.deviation; });

            for (const auto &candidate : candidates)
            {
                if (current_total_radical() >= std::abs(allowed_charge) + allowed_radical)
                    return hit;
                OBBond *bond = mol.GetBond(candidate.idx1, candidate.idx2);
                if (!bond || bond->IsRotor() || bond->GetBondOrder() == 1)
                    continue;
                bond->SetBondOrder(bond->GetBondOrder() - 1);
                OBAtom *begin_atom = bond->GetBeginAtom();
                OBAtom *end_atom = bond->GetEndAtom();
                begin_atom->SetSpinMultiplicity(begin_atom->GetSpinMultiplicity() + 1);
                end_atom->SetSpinMultiplicity(end_atom->GetSpinMultiplicity() + 1);
                hit = true;
                LOG_DEBUG("[BreakDeformedEne] Broken sorted candidate");
            }
            return hit;
        }

        bool BreakOneBond(OBMol &mol, int &charge, int allowed_radical)
        {
            bool hit = false;
            auto check_cond = [&]()
            {
                int sum = 0;
                FOR_ATOMS_OF_MOL(a, mol)
                sum += a->GetSpinMultiplicity();
                return sum >= std::abs(charge) + allowed_radical;
            };

            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::BREAK_ONE_BOND_MULTIPLE);
                if (matches.empty())
                    break;
                if (check_cond())
                    return hit;
                const auto &idxs = matches.front();
                OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                if (!bond)
                    break;
                bond->SetBondOrder(bond->GetBondOrder() - 1);
                OBAtom *begin_atom = bond->GetBeginAtom();
                OBAtom *end_atom = bond->GetEndAtom();
                begin_atom->SetSpinMultiplicity(begin_atom->GetSpinMultiplicity() + 1);
                end_atom->SetSpinMultiplicity(end_atom->GetSpinMultiplicity() + 1);
                hit = true;
                LOG_DEBUG("[BreakOneBond] Broken triple/double candidate");
            }

            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::BREAK_ONE_BOND_CATION);
                if (matches.empty())
                    break;
                if (check_cond())
                    return hit;
                const auto &idxs = matches.front();
                OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                if (!bond)
                    break;
                bond->SetBondOrder(bond->GetBondOrder() - 1);
                OBAtom *idx0_atom = mol.GetAtom(idxs[0]);
                OBAtom *idx1_atom = mol.GetAtom(idxs[1]);
                if (!idx0_atom || !idx1_atom)
                    break;
                idx1_atom->SetSpinMultiplicity(idx1_atom->GetSpinMultiplicity() + 1);
                idx0_atom->SetFormalCharge(idx0_atom->GetFormalCharge() - 1);
                charge += 1;
                hit = true;
                LOG_DEBUG("[BreakOneBond] Broken charge-transfer candidate");
            }

            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::BREAK_ONE_BOND_AROMATIC);
                if (!matches.empty())
                {
                    if (check_cond())
                        return hit;
                    for (const auto &idxs : matches)
                    {
                        OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                        if (!bond)
                            continue;
                        if (bond->GetBondOrder() == 1)
                            continue;
                        bond->SetBondOrder(bond->GetBondOrder() - 1);
                        OBAtom *begin_atom = bond->GetBeginAtom();
                        OBAtom *end_atom = bond->GetEndAtom();
                        begin_atom->SetSpinMultiplicity(begin_atom->GetSpinMultiplicity() + 1);
                        end_atom->SetSpinMultiplicity(end_atom->GetSpinMultiplicity() + 1);
                        hit = true;
                        LOG_DEBUG("[BreakOneBond] Broken aromatic candidate");
                        break;
                    }
                }
            }

            if (check_cond())
                return hit;

            bool all_single = true;
            FOR_BONDS_OF_MOL(b, mol)
            if (b->GetBondOrder() != 1) all_single = false;

            if (all_single)
            {
                std::vector<OBBond *> single_bonds;
                FOR_BONDS_OF_MOL(b, mol)
                {
                    single_bonds.push_back(&*b);
                }
                if (!single_bonds.empty())
                {
                    if (check_cond())
                        return hit;
                    OBBond *single_bond = single_bonds.front();
                    OBAtom *b_at = single_bond->GetBeginAtom();
                    OBAtom *e_at = single_bond->GetEndAtom();
                    b_at->SetSpinMultiplicity(b_at->GetSpinMultiplicity() + 1);
                    e_at->SetSpinMultiplicity(e_at->GetSpinMultiplicity() + 1);
                    mol.DeleteBond(single_bond);
                    hit = true;
                    LOG_DEBUG("[BreakOneBond] Deleted single bond");
                }
            }
            return hit;
        }
    }
}
