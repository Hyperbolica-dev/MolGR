#include "molgr/utils/scoring.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/smarts.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/graphsym.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <unordered_set>

namespace
{
    constexpr double kCoordinateScale = 1000000.0;

    OpenBabel::OBMol &MutableMol(const OpenBabel::OBMol &mol)
    {
        return const_cast<OpenBabel::OBMol &>(mol);
    }

    std::int64_t QuantizedCoordinate(double value)
    {
        return static_cast<std::int64_t>(std::llround(value * kCoordinateScale));
    }

    double CalculateChargeInteractionPenalty(int metal_valence, int ligand_charge, double dist_sq)
    {
        if (dist_sq <= 0.0)
        {
            return 0.0;
        }
        if (ligand_charge > 0)
        {
            return 100.0 * (ligand_charge * metal_valence) / dist_sq;
        }
        if (ligand_charge < 0)
        {
            return -5.0 * (std::abs(ligand_charge) * metal_valence) / dist_sq;
        }
        return 0.0;
    }

    double CalculateChargePenaltyFromData(
        int atomic_num,
        int charge,
        int total_valence,
        int spin_multiplicity)
    {
        if (charge == 0)
        {
            return 0.0;
        }

        const auto info_it = molgr::kNonMetalDict.find(atomic_num);
        if (info_it == molgr::kNonMetalDict.end())
        {
            return 0.0;
        }

        const int total_electrons =
            info_it->second.num_outer_electrons + total_valence - charge + spin_multiplicity;
        if (total_electrons == 8 || total_electrons == 2)
        {
            return 0.0;
        }

        const double en = OpenBabel::OBElements::GetElectroNeg(atomic_num);
        if (charge > 0)
        {
            return std::abs(charge) * std::max(0.0, en - 2.0) * 3.0;
        }
        return std::abs(charge) * std::max(0.0, 4.0 - en) * 3.0;
    }

    double CalculateRadicalPenaltyFromData(int atomic_num, int radical_num, int heavy_degree)
    {
        if (radical_num == 0)
        {
            return 0.0;
        }
        if (molgr::kHeteroatoms.find(atomic_num) != molgr::kHeteroatoms.end())
        {
            return radical_num * 10.0;
        }
        return (3.0 - static_cast<double>(heavy_degree)) * 1.5;
    }

    double CalculateCoulombicPenalty(const OpenBabel::OBBond &bond)
    {
        const int q1 = bond.GetBeginAtom()->GetFormalCharge();
        const int q2 = bond.GetEndAtom()->GetFormalCharge();
        if (q1 == 0 || q2 == 0)
        {
            return 0.0;
        }
        if (q1 * q2 > 0)
        {
            return 15.0;
        }
        return -0.5;
    }

    double CalcRemainBondOrderReward(const OpenBabel::OBMol &mol)
    {
        OpenBabel::OBMol &mutable_mol = MutableMol(mol);
        double remain_valence = 0.0;
        double bond_order_sum = 0.0;
        FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (atom.IsMetal())
            {
                continue;
            }
            const auto info_it = molgr::kNonMetalDict.find(atom.GetAtomicNum());
            if (info_it != molgr::kNonMetalDict.end())
            {
                remain_valence += info_it->second.default_valence;
            }
        }
        FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
        {
            bond_order_sum += bond_iter->GetBondOrder();
        }
        return (remain_valence - bond_order_sum * 2.0) * 5.0;
    }

    double MetalStateSymmetryPenalty(const std::vector<molgr::MetalAtomPosition> &metal_states)
    {
        std::unordered_set<std::string> unique_states;
        for (const auto &metal_state : metal_states)
        {
            unique_states.insert(
                std::to_string(metal_state.element_idx) + ":" +
                std::to_string(metal_state.valence) + ":" +
                std::to_string(metal_state.radical_num));
        }
        return (static_cast<double>(unique_states.size()) - static_cast<double>(metal_states.size())) * 2.0;
    }
}

namespace molgr
{
    namespace scoring
    {
        std::string BuildScoreKey(const OpenBabel::OBMol &mol)
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetOutFormat("molreport");
                initialized = true;
            }
            OpenBabel::OBMol mol_copy(mol);
            return conv.WriteString(&mol_copy, true);
        }

        std::string BuildMetalStateKey(const std::vector<MetalAtomPosition> &metal_states)
        {
            std::ostringstream oss;
            for (const auto &metal_state : metal_states)
            {
                oss << metal_state.idx << ','
                    << metal_state.element_idx << ','
                    << metal_state.valence << ','
                    << metal_state.radical_num << ','
                    << QuantizedCoordinate(metal_state.position_x) << ','
                    << QuantizedCoordinate(metal_state.position_y) << ','
                    << QuantizedCoordinate(metal_state.position_z) << ';';
            }
            return oss.str();
        }

        double GetDeviationScore(const OpenBabel::OBMol &mol, const OpenBabel::OBAtom *atom)
        {
            if (atom == nullptr)
            {
                return 0.0;
            }

            std::vector<const OpenBabel::OBAtom *> neighbors;
            FOR_NB_OF_ATOM(neighbor_iter, const_cast<OpenBabel::OBAtom *>(atom))
            {
                neighbors.push_back(&(*neighbor_iter));
            }

            if (neighbors.size() == 2)
            {
                const double angle = MutableMol(mol).GetAngle(
                    const_cast<OpenBabel::OBAtom *>(neighbors[0]),
                    const_cast<OpenBabel::OBAtom *>(atom),
                    const_cast<OpenBabel::OBAtom *>(neighbors[1]));
                return std::abs(angle - 108.0) / 108.0;
            }

            if (neighbors.size() == 3)
            {
                return 1.0 - molgr::utils::CalculateShapeQuality(
                                 neighbors[0]->GetVector(),
                                 neighbors[1]->GetVector(),
                                 neighbors[2]->GetVector(),
                                 atom->GetVector());
            }

            return 0.0;
        }

        double CalcSymmetryPenalty(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBGraphSym graph_sym(&MutableMol(mol));
            std::vector<unsigned int> symmetry_ids;
            graph_sym.GetSymmetry(symmetry_ids);
            std::unordered_set<unsigned int> unique_ids(symmetry_ids.begin(), symmetry_ids.end());
            return (static_cast<double>(unique_ids.size()) - static_cast<double>(mol.NumAtoms())) * 2.0;
        }

        double CalculatePhysChemPenalty(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            double total_penalty = 0.0;
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (atom.IsMetal())
                {
                    continue;
                }
                total_penalty += CalculateChargePenaltyFromData(
                    atom.GetAtomicNum(),
                    atom.GetFormalCharge(),
                    atom.GetTotalValence(),
                    atom.GetSpinMultiplicity());
                total_penalty += CalculateRadicalPenaltyFromData(
                    atom.GetAtomicNum(),
                    atom.GetSpinMultiplicity(),
                    atom.GetHvyDegree());
            }

            FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
            {
                total_penalty += CalculateCoulombicPenalty(*bond_iter);
            }
            return total_penalty;
        }

        double CalculateMetalPenalty(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            double penalty = 0.0;
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                OpenBabel::OBAtom &metal_atom = *atom_iter;
                if (!metal_atom.IsMetal())
                {
                    continue;
                }

                const std::string symbol = OpenBabel::OBElements::GetSymbol(metal_atom.GetAtomicNum());
                const int valence = metal_atom.GetFormalCharge();
                if (valence <= 0)
                {
                    penalty += 10.0 * std::max(std::abs(valence), 1);
                }

                const auto prior_it = molgr::kMetalValencePrior.find(symbol);
                const auto minor_it = molgr::kMetalValenceMinor.find(symbol);
                const bool in_prior = prior_it != molgr::kMetalValencePrior.end() &&
                                      std::find(prior_it->second.begin(), prior_it->second.end(), valence) !=
                                          prior_it->second.end();
                const bool in_minor = minor_it != molgr::kMetalValenceMinor.end() &&
                                      std::find(minor_it->second.begin(), minor_it->second.end(), valence) !=
                                          minor_it->second.end();
                if (!in_prior)
                {
                    if (in_minor)
                    {
                        penalty += 10.0;
                    }
                    else
                    {
                        penalty += 20.0;
                    }
                }

                if (valence <= 0)
                {
                    continue;
                }

                FOR_ATOMS_OF_MOL(other_iter, mutable_mol)
                {
                    OpenBabel::OBAtom &other = *other_iter;
                    if (other.GetIdx() == metal_atom.GetIdx())
                    {
                        continue;
                    }
                    const double dist = metal_atom.GetDistance(&other);
                    if (dist > 2.6)
                    {
                        continue;
                    }
                    penalty += CalculateChargeInteractionPenalty(
                        valence,
                        other.GetFormalCharge(),
                        dist * dist);
                }
            }
            return penalty;
        }

        double CalculateMetalPenaltyFromMetalStates(
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states,
            double cutoff)
        {
            double penalty = 0.0;
            const double cutoff_sq = cutoff * cutoff;

            for (std::size_t i = 0; i < metal_states.size(); ++i)
            {
                const auto &metal_state = metal_states[i];
                const int valence = metal_state.valence;
                if (valence <= 0)
                {
                    penalty += 10.0 * std::max(std::abs(valence), 1);
                }

                const auto prior_it = molgr::kMetalValencePrior.find(metal_state.symbol);
                const auto minor_it = molgr::kMetalValenceMinor.find(metal_state.symbol);
                const bool in_prior = prior_it != molgr::kMetalValencePrior.end() &&
                                      std::find(prior_it->second.begin(), prior_it->second.end(), valence) !=
                                          prior_it->second.end();
                const bool in_minor = minor_it != molgr::kMetalValenceMinor.end() &&
                                      std::find(minor_it->second.begin(), minor_it->second.end(), valence) !=
                                          minor_it->second.end();
                if (!in_prior)
                {
                    if (in_minor)
                    {
                        penalty += 10.0;
                    }
                    else
                    {
                        penalty += 20.0;
                    }
                }

                if (valence <= 0)
                {
                    continue;
                }

                for (const auto &[ligand_charge, x, y, z] : charged_atom_snapshots)
                {
                    const double dx = metal_state.position_x - x;
                    const double dy = metal_state.position_y - y;
                    const double dz = metal_state.position_z - z;
                    const double dist_sq = dx * dx + dy * dy + dz * dz;
                    if (dist_sq > cutoff_sq)
                    {
                        continue;
                    }
                    penalty += CalculateChargeInteractionPenalty(valence, ligand_charge, dist_sq);
                }

                for (std::size_t j = 0; j < metal_states.size(); ++j)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    const auto &other_metal = metal_states[j];
                    if (other_metal.valence == 0)
                    {
                        continue;
                    }
                    const double dx = metal_state.position_x - other_metal.position_x;
                    const double dy = metal_state.position_y - other_metal.position_y;
                    const double dz = metal_state.position_z - other_metal.position_z;
                    const double dist_sq = dx * dx + dy * dy + dz * dz;
                    if (dist_sq > cutoff_sq)
                    {
                        continue;
                    }
                    penalty += CalculateChargeInteractionPenalty(valence, other_metal.valence, dist_sq);
                }
            }

            return penalty;
        }

        static double CalculateHeteroatomPenalty(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            double penalty = 0.0;
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (molgr::kHeteroatoms.find(atom.GetAtomicNum()) == molgr::kHeteroatoms.end())
                {
                    continue;
                }
                penalty += 10.0 * (atom.GetFormalCharge() - atom.GetTotalValence());
            }
            return penalty;
        }

        static double CalculateConjugationReward(const OpenBabel::OBMol &mol)
        {
            auto matches = molgr::smarts::Match(
                MutableMol(mol),
                molgr::smarts::PatternId::SCORING_CONJUGATION);
            return static_cast<double>(matches.size()) * 2.0;
        }

        double OrganicCoreScore(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            double score = CalcRemainBondOrderReward(mol);

            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (atom.IsMetal())
                {
                    continue;
                }

                if (atom.IsAromatic())
                {
                    score -= 5.0 - std::abs(atom.GetFormalCharge()) * 3.0;
                }

                const bool needs_deviation = atom.GetSpinMultiplicity() > 0 || atom.GetFormalCharge() != 0;
                if (needs_deviation)
                {
                    const double deviation = GetDeviationScore(mol, &atom);
                    if (atom.GetSpinMultiplicity() > 0)
                    {
                        score += deviation * 10.0;
                    }
                    if (atom.GetFormalCharge() > 0)
                    {
                        score += deviation * 10.0;
                    }
                    if (atom.GetFormalCharge() < 0)
                    {
                        score += (1.0 - deviation) * 10.0;
                    }
                }

                if (atom.GetAtomicNum() == 6)
                {
                    bool all_double_bond = true;
                    FOR_BONDS_OF_ATOM(bond_iter, const_cast<OpenBabel::OBAtom *>(&atom))
                    {
                        if (bond_iter->GetBondOrder() != 2)
                        {
                            all_double_bond = false;
                            break;
                        }
                    }
                    if (all_double_bond)
                    {
                        score += 5.0;
                    }
                }
            }

            score += CalculatePhysChemPenalty(mol);
            score += CalculateHeteroatomPenalty(mol);
            score -= CalculateConjugationReward(mol);
            return score;
        }

        double PostReinsertionScore(const OpenBabel::OBMol &mol)
        {
            return CalcSymmetryPenalty(mol) + CalculateMetalPenalty(mol);
        }

        std::pair<double, ChargedAtomSnapshotList> BuildPostReinsertionBaseComponents(
            const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            ChargedAtomSnapshotList charged_atom_snapshots;
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (atom.GetFormalCharge() == 0)
                {
                    continue;
                }
                charged_atom_snapshots.emplace_back(
                    atom.GetFormalCharge(),
                    atom.GetX(),
                    atom.GetY(),
                    atom.GetZ());
            }
            return {CalcSymmetryPenalty(mol), charged_atom_snapshots};
        }

        double PostReinsertionScoreFromMetalStates(
            double base_symmetry_penalty,
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states)
        {
            return base_symmetry_penalty +
                   MetalStateSymmetryPenalty(metal_states) +
                   CalculateMetalPenaltyFromMetalStates(charged_atom_snapshots, metal_states);
        }

        double CombinedCandidateScoreFromMetalStates(
            double organic_score,
            const std::string &post_reinsertion_base_key,
            double base_symmetry_penalty,
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states)
        {
            (void)post_reinsertion_base_key;
            return organic_score +
                   PostReinsertionScoreFromMetalStates(
                       base_symmetry_penalty,
                       charged_atom_snapshots,
                       metal_states);
        }

        double OmolScore(const OpenBabel::OBMol &mol)
        {
            MutableMol(mol).SetAromaticPerceived(false);
            return OrganicCoreScore(mol) + PostReinsertionScore(mol);
        }
    }
}
