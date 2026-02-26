#include "molgr/stages/preprocess.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/logger.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>

#include <algorithm>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        std::vector<int> GetFlatAtomList(OBMol &mol, const std::string &smarts)
        {
            std::vector<int> atom_indices;
            auto matches = molgr::utils::FindSmarts(mol, smarts);
            for (const auto &match : matches)
            {
                for (int idx : match)
                    atom_indices.push_back(idx);
            }
            return atom_indices;
        }

        bool ContainsAtomIdx(const std::vector<int> &atom_indices, int atom_idx)
        {
            return std::find(atom_indices.begin(), atom_indices.end(), atom_idx) != atom_indices.end();
        }

        void MakeConnections(OBMol &mol, double factor)
        {
            std::string donate_smarts = "[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]";
            std::string accept_smarts = "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]";

            std::vector<int> donate_atoms = GetFlatAtomList(mol, donate_smarts);
            std::vector<int> accept_atoms = GetFlatAtomList(mol, accept_smarts);

            while (!donate_atoms.empty() && !accept_atoms.empty())
            {
                int donor_idx = donate_atoms.front();
                donate_atoms.erase(donate_atoms.begin());

                OBAtom *donor = mol.GetAtom(donor_idx);
                if (!donor)
                    continue;

                std::vector<std::pair<int, int>> pairs;
                for (int acceptor_idx : accept_atoms)
                {
                    OBAtom *acceptor = mol.GetAtom(acceptor_idx);
                    if (acceptor && acceptor_idx != donor_idx)
                    {
                        pairs.push_back({donor_idx, acceptor_idx});
                    }
                }

                std::sort(pairs.begin(), pairs.end(), [&mol](const auto &a, const auto &b)
                          {
                    OBAtom *a1 = mol.GetAtom(a.first);
                    OBAtom *a2 = mol.GetAtom(a.second);
                    OBAtom *b1 = mol.GetAtom(b.first);
                    OBAtom *b2 = mol.GetAtom(b.second);
                    return a1->GetDistance(a2) < b1->GetDistance(b2); });

                if (pairs.empty())
                {
                    continue;
                }

                for (const auto &pair : pairs)
                {
                    int p1 = pair.first;
                    int p2 = pair.second;
                    OBAtom *a1 = mol.GetAtom(p1);
                    OBAtom *a2 = mol.GetAtom(p2);
                    if (!a1 || !a2)
                    {
                        continue;
                    }
                    double dist = a1->GetDistance(a2);

                    double r1 = OBElements::GetCovalentRad(a1->GetAtomicNum());
                    double r2 = OBElements::GetCovalentRad(a2->GetAtomicNum());

                    if (dist < (r1 + r2) * factor &&
                        ContainsAtomIdx(donate_atoms, p1) &&
                        ContainsAtomIdx(accept_atoms, p2))
                    {
                        OBBond *bond = mol.GetBond(a1, a2);
                        if (!bond)
                        {
                            mol.AddBond(p1, p2, 1);
                            LOG_DEBUG("[MakeConnections] Add Bond " << p1 << "-" << p2);
                            continue;
                        }
                        if (bond->GetBondOrder() == 0)
                        {
                            bond->SetBondOrder(1);
                            LOG_DEBUG("[MakeConnections] Set Bond Order 1 " << p1 << "-" << p2);
                            donate_atoms = GetFlatAtomList(mol, donate_smarts);
                            accept_atoms = GetFlatAtomList(mol, accept_smarts);
                        }
                    }
                }
            }
        }

        void PreClean(OBMol &mol)
        {
            while (true)
            {
                auto matches1 = molgr::utils::FindSmarts(mol, "[Cv5,Nv5,Pv5,Siv5]=,#[*]");
                if (matches1.empty())
                {
                    break;
                }
                const auto &match = matches1.front();
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                {
                    bond->SetBondOrder(bond->GetBondOrder() - 1);
                }
            }

            while (true)
            {
                auto matches2 = molgr::utils::FindSmarts(mol, "[#6]1([#6]2)([#6]3)[#7]23[#6]1");
                if (matches2.empty())
                {
                    break;
                }
                const auto &idxs = matches2.front();
                int n_idx = -1, c_idx = -1;
                for (int idx : idxs)
                {
                    OBAtom *atom = mol.GetAtom(idx);
                    int internal_degree = 0;
                    for (int other : idxs)
                        if (idx != other && mol.GetBond(idx, other))
                            internal_degree++;

                    if (internal_degree >= 3)
                    {
                        if (atom->GetAtomicNum() == 7)
                            n_idx = idx;
                        if (atom->GetAtomicNum() == 6)
                            c_idx = idx;
                    }
                }
                if (n_idx != -1 && c_idx != -1)
                {
                    OBBond *bond = mol.GetBond(n_idx, c_idx);
                    if (bond)
                        mol.DeleteBond(bond);
                }
            }

            while (true)
            {
                auto matches3 = molgr::utils::FindSmarts(mol, "[#6]1([#6]2)[#7]2[#6]1");
                if (matches3.empty())
                {
                    break;
                }
                const auto &idxs = matches3.front();
                int amine_n = -1, butyl_c = -1;
                for (int idx : idxs)
                {
                    OBAtom *atom = mol.GetAtom(idx);
                    int internal_degree = 0;
                    for (int other : idxs)
                        if (idx != other && mol.GetBond(idx, other))
                            internal_degree++;
                    if (internal_degree >= 2)
                    {
                        if (atom->GetAtomicNum() == 7)
                            amine_n = idx;
                        if (atom->GetAtomicNum() == 6)
                            butyl_c = idx;
                    }
                }
                if (amine_n != -1 && butyl_c != -1)
                {
                    OBBond *bond = mol.GetBond(amine_n, butyl_c);
                    if (bond)
                        mol.DeleteBond(bond);
                }
            }

            while (true)
            {
                auto matches4 = molgr::utils::FindSmarts(mol, "[Siv5]-[O,F]");
                if (matches4.empty())
                {
                    break;
                }
                const auto &match = matches4.front();
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                {
                    mol.DeleteBond(bond);
                }
            }
        }
    }
}
