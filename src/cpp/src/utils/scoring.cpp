/**
 * @file scoring.cpp
 * @brief Implementation of molecular structure scoring logic with logging.
 * @namespace molgr::scoring
 * @author TMJ
 * @date 2025-12-25
 */

#include "molgr/utils/scoring.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/utils.h"
#include "molgr/utils/logger.h"

#include <openbabel/graphsym.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>

#include <cmath>
#include <set>
#include <iostream>
#include <algorithm>

namespace molgr
{
    namespace scoring
    {

        using namespace OpenBabel;

        // 辅助：获取易读的原子标识，如 "C(1)"
        std::string AtomId(OBAtom *atom)
        {
            std::string s = OBElements::GetSymbol(atom->GetAtomicNum());
            s += "(" + std::to_string(atom->GetIdx()) + ")";
            return s;
        }

        // =============================================================================
        // 1. Geometry Deviation Score
        // =============================================================================

        double GetDeviationScore(OBMol &mol, OBAtom *atom)
        {
            std::vector<OBAtom *> neighbors;
            FOR_NB_OF_ATOM(nbr, atom)
            {
                neighbors.push_back(&*nbr);
            }

            // Case 1: 2 Neighbors -> Angle check (Target 108.0)
            if (neighbors.size() == 2)
            {
                double angle = mol.GetAngle(neighbors[0], atom, neighbors[1]);
                double score = std::abs(angle - 108.0) / 108.0;
                return score;
            }

            // Case 2: 3 Neighbors -> Tetrahedral shape quality
            if (neighbors.size() == 3)
            {
                vector3 p1 = neighbors[0]->GetVector();
                vector3 p2 = neighbors[1]->GetVector();
                vector3 p3 = neighbors[2]->GetVector();
                vector3 p_atom = atom->GetVector();

                double quality = molgr::utils::CalculateShapeQuality(p1, p2, p3, p_atom);
                double score = 1.0 - quality;
                return score;
            }

            return 0.0;
        }

        // =============================================================================
        // 2. Symmetry Penalty
        // =============================================================================

        double CalcSymmetryPenalty(OBMol &mol)
        {
            OBGraphSym gs(&mol);
            std::vector<unsigned int> symmetry_ids;
            gs.GetSymmetry(symmetry_ids);

            std::set<unsigned int> unique_ids(symmetry_ids.begin(), symmetry_ids.end());

            // logic: (num_unique_classes - num_atoms) * 2.0
            double penalty = (static_cast<double>(unique_ids.size()) - static_cast<double>(mol.NumAtoms())) * 2.0;

            LOG_DEBUG("[Symmetry] Classes: " << unique_ids.size()
                                             << ", Atoms: " << mol.NumAtoms()
                                             << " -> Penalty: " << penalty);

            return penalty;
        }

        // =============================================================================
        // 3. PhysChem Penalty (Internal Helpers)
        // =============================================================================

        double CalculateChargePenalty(OBAtom *atom)
        {
            double penalty = 0.0;
            int charge = atom->GetFormalCharge();
            if (charge == 0)
                return penalty;

            int atomic_num = atom->GetAtomicNum();
            if (kNonMetalDict.find(atomic_num) == kNonMetalDict.end())
                return penalty;

            const ElementInfo &info = kNonMetalDict.at(atomic_num);
            int total_electrons = info.num_outer_electrons + atom->GetTotalValence() - charge + atom->GetSpinMultiplicity();

            // 稳定八隅体或二隅体
            if (total_electrons == 8 || total_electrons == 2)
                return penalty;

            double en = OBElements::GetElectroNeg(atomic_num);

            if (charge > 0)
            {
                penalty += std::abs(charge) * std::max(0.0, en - 2.0) * 3.0;
            }
            if (charge < 0)
            {
                penalty += std::abs(charge) * std::max(0.0, 4.0 - en) * 3.0;
            }

            if (penalty > 0.0)
            {
                LOG_DEBUG("[PhysChem] Charge Penalty on " << AtomId(atom)
                                                          << " (Q=" << charge << ", EN=" << en << "): " << penalty);
            }

            return penalty;
        }

        double CalculateRadicalPenalty(OBAtom *atom)
        {
            double penalty = 0.0;
            int radical_num = atom->GetSpinMultiplicity();
            if (radical_num == 0)
                return penalty;

            int atomic_num = atom->GetAtomicNum();

            if (kHeteroatoms.find(atomic_num) != kHeteroatoms.end())
            {
                penalty = radical_num * 2.0;
            }
            else
            {
                penalty = (3.0 - static_cast<double>(atom->GetHvyDegree())) * 1.5;
            }

            LOG_DEBUG("[PhysChem] Radical Penalty on " << AtomId(atom)
                                                       << " (Rad=" << radical_num << "): " << penalty);

            return penalty;
        }

        double CalculateCoulombicPenalty(OBBond *bond)
        {
            OBAtom *atom1 = bond->GetBeginAtom();
            OBAtom *atom2 = bond->GetEndAtom();
            int q1 = atom1->GetFormalCharge();
            int q2 = atom2->GetFormalCharge();

            if (q1 == 0 || q2 == 0)
                return 0.0;

            if (q1 * q2 > 0)
            {
                // 同性相斥
                LOG_DEBUG("[PhysChem] Coulombic Repulsion: " << AtomId(atom1) << " <-> " << AtomId(atom2) << " (+15.0)");
                return 15.0;
            }
            // 异性相吸
            return -0.5;
        }

        double CalculatePhysChemPenalty(OBMol &mol)
        {
            double total_penalty = 0.0;

            FOR_ATOMS_OF_MOL(atom, mol)
            {
                if (atom->IsMetal())
                    continue;
                total_penalty += CalculateChargePenalty(&*atom);
                total_penalty += CalculateRadicalPenalty(&*atom);
            }

            FOR_BONDS_OF_MOL(bond, mol)
            {
                total_penalty += CalculateCoulombicPenalty(&*bond);
            }

            return total_penalty;
        }

        // =============================================================================
        // 4. Metal Penalty (Internal Helpers)
        // =============================================================================

        std::vector<std::pair<OBAtom *, double>> GetMetalCoordinationSphere(OBMol &mol, OBAtom *metal_atom, double cutoff = 2.8)
        {
            std::vector<std::pair<OBAtom *, double>> neighbors;
            FOR_ATOMS_OF_MOL(other, mol)
            {
                if (other->GetIdx() == metal_atom->GetIdx())
                    continue;
                double dist = metal_atom->GetDistance(&*other);
                if (dist <= cutoff)
                {
                    neighbors.push_back({&*other, dist});
                }
            }
            return neighbors;
        }

        double CalculateMetalPenalty(OBMol &mol)
        {
            double penalty = 0.0;

            FOR_ATOMS_OF_MOL(atom, mol)
            {
                if (!atom->IsMetal())
                    continue;

                std::string symbol = OBElements::GetSymbol(atom->GetAtomicNum());
                int valence = atom->GetFormalCharge();

                double atom_penalty = 0.0;

                // Rule 1: Valence Validity
                if (valence <= 0)
                    atom_penalty += 10.0;

                bool is_prior = false;
                bool is_minor = false;

                if (kMetalValencePrior.count(symbol))
                {
                    const auto &vec = kMetalValencePrior.at(symbol);
                    for (int v : vec)
                        if (v == valence)
                            is_prior = true;
                }

                if (kMetalValenceMinor.count(symbol))
                {
                    const auto &vec = kMetalValenceMinor.at(symbol);
                    for (int v : vec)
                        if (v == valence)
                            is_minor = true;
                }

                if (!is_prior)
                {
                    if (is_minor)
                    {
                        atom_penalty += 2.0;
                        LOG_DEBUG("[Metal] " << symbol << " Minor Valence (" << valence << ") +2.0");
                    }
                    else
                    {
                        atom_penalty += 10.0;
                        LOG_DEBUG("[Metal] " << symbol << " Invalid/Rare Valence (" << valence << ") +10.0");
                    }
                }

                // Rule 2: Coordination Sphere Electrostatics
                auto neighbors = GetMetalCoordinationSphere(mol, &*atom, 2.6);
                for (const auto &pair : neighbors)
                {
                    OBAtom *ligand = pair.first;
                    double dist = pair.second;
                    int ligand_charge = ligand->GetFormalCharge();

                    if (valence > 0)
                    {
                        if (ligand_charge > 0)
                        {
                            double p = 10.0 * (ligand_charge * valence) / (dist * dist);
                            atom_penalty += p;
                            LOG_DEBUG("[Metal] Repulsion with " << AtomId(ligand) << " Dist=" << dist << " -> +" << p);
                        }
                        else if (ligand_charge < 0)
                        {
                            double p = 2.0 * (std::abs(ligand_charge) * valence) / dist;
                            atom_penalty -= p; // Bonus
                            // LOG_DEBUG("[Metal] Attraction with " << AtomId(ligand) << " -> -" << p); // 奖励通常不打印，除非调试
                        }
                    }
                }

                penalty += atom_penalty;
            }
            return penalty;
        }

        // =============================================================================
        // Main Score Function
        // =============================================================================

        double OmolScore(OBMol &mol)
        {
            LOG_DEBUG("=== Scoring Start (Atoms: " << mol.NumAtoms() << ") ===");
            mol.SetAromaticPerceived(false);
            double score = 0.0;

            // 1. Symmetry
            double sym = CalcSymmetryPenalty(mol);
            score += sym;

            // 2. Geometry Deviation
            double geo_penalty = 0.0;
            FOR_ATOMS_OF_MOL(atom, mol)
            {
                if (atom->IsMetal())
                    continue;

                if (!atom->IsAromatic())
                    score += 5.0;

                double dev = GetDeviationScore(mol, &*atom);
                int charge = atom->GetFormalCharge();
                int radical = atom->GetSpinMultiplicity();

                double term = 0.0;
                if (radical > 0)
                    term = dev * 10.0;
                else if (charge > 0)
                    term = dev * 10.0;
                else if (charge < 0)
                    term = (1.0 - dev) * 10.0;

                if (term > 0.5)
                { // 只记录显著的几何偏差罚分
                    LOG_DEBUG("[Geometry] " << AtomId(&*atom) << " Dev=" << dev << " (Q=" << charge << ",R=" << radical << ") -> +" << term);
                }
                geo_penalty += term;
            }
            score += geo_penalty;
            LOG_DEBUG(">> Geometry Total: " << geo_penalty);

            // 3. PhysChem
            double phys = CalculatePhysChemPenalty(mol);
            score += phys;
            LOG_DEBUG(">> PhysChem Total: " << phys);

            // 4. Metal
            double metal = CalculateMetalPenalty(mol);
            score += metal;
            LOG_DEBUG(">> Metal Total:    " << metal);

            LOG_DEBUG("=== Final Score: " << score << " ===");
            return score;
        }

    } // namespace scoring
} // namespace molgr
