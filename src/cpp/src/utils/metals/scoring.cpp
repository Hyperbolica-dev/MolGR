#include "molgr/utils/metals/scoring.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/organic_topology.h"
#include "molgr/utils/scoring.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace
{
    constexpr std::size_t kWeightedSelectionFieldCount = 11;

    using OrganicStateKey = std::tuple<int, int, double, int, int, double>;
    using RawCandidateSelectionKey =
        std::tuple<int, int, double, int, int, double, double, double, double, double, double, int, double, int>;
    using WeightedSelectionKey =
        std::tuple<double, int, int, double, int, int, double, double, double, double, double, double, int, double, int>;

    struct Point3D
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    struct CoordinationBlocker
    {
        int atom_idx = 0;
        Point3D coordinates;
        double radius = 0.0;
    };

    struct MetalSiteEnvironmentProfile
    {
        double electrostatic_support = 0.0;
        double visible_anionic_donor_support = 0.0;
        double visible_neutral_donor_support = 0.0;
        double visible_effective_donor_support = 0.0;
        double obstructed_negative_effective_donor_support = 0.0;
    };

    struct MetalSiteHeuristicScore
    {
        double electrostatic_support = 0.0;
        double anionic_donor_support = 0.0;
        double neutral_donor_support = 0.0;
        double coordination_access_penalty = 0.0;
        double visible_coordination_reward = 0.0;
        double negative_metal_visible_coordination_penalty = 0.0;
        double obstructed_opposite_charge_penalty = 0.0;
        double electrostatic_penalty = 0.0;
        double donor_penalty = 0.0;
    };

    struct OrganicElectronicStateMetrics
    {
        int aromatic_atom_count = 0;
        int aromatic_ring_count = 0;
        int conjugated_atom_count = 0;
        int conjugated_bond_count = 0;
        int max_conjugated_component_size = 0;
        double radical_localization_penalty = 0.0;
        double charge_localization_penalty = 0.0;
    };

    struct OrganicElectronicStateSelectionContext
    {
        int max_aromatic_atom_count = 0;
        int max_aromatic_ring_count = 0;
        int max_conjugated_atom_count = 0;
        int max_conjugated_component_size = 0;
    };

    struct WeightedSelectionContext
    {
        std::array<double, kWeightedSelectionFieldCount> best_values{};
        std::array<double, kWeightedSelectionFieldCount> weights{};
        std::array<double, kWeightedSelectionFieldCount> scales{};
    };

    std::optional<double> LookupMapValue(
        const std::map<int, double> &values,
        int atomic_num)
    {
        const auto it = values.find(atomic_num);
        if (it == values.end())
        {
            return std::nullopt;
        }
        return it->second;
    }

    double NormalizeValue(
        const std::optional<double> &value,
        double lower,
        double upper,
        double fallback)
    {
        if (!value.has_value() || upper <= lower)
        {
            return fallback;
        }
        const double clamped = std::min(std::max(*value, lower), upper);
        return (clamped - lower) / (upper - lower);
    }

    bool AtomHasOddSpin(const OpenBabel::OBAtom &atom)
    {
        return static_cast<int>(atom.GetSpinMultiplicity()) % 2 == 1;
    }

    int BondOrder(const OpenBabel::OBBond &bond)
    {
        return bond.IsAromatic() ? 2 : static_cast<int>(bond.GetBondOrder());
    }

    bool HasMultipleBondToAtomicNum(
        OpenBabel::OBAtom &atom,
        const std::set<int> &atomic_nums)
    {
        FOR_BONDS_OF_ATOM(bond_iter, &atom)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (BondOrder(bond) < 2)
            {
                continue;
            }
            OpenBabel::OBAtom *neighbor = bond.GetNbrAtom(&atom);
            if (neighbor != nullptr &&
                atomic_nums.find(static_cast<int>(neighbor->GetAtomicNum())) != atomic_nums.end())
            {
                return true;
            }
        }
        return false;
    }

    bool IsCarbonylLikeOxygen(OpenBabel::OBAtom &atom)
    {
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        if (atomic_num != 8 && atomic_num != 16 && atomic_num != 34 && atomic_num != 52)
        {
            return false;
        }
        return HasMultipleBondToAtomicNum(atom, {6});
    }

    bool IsAmideLikeNitrogen(OpenBabel::OBAtom &atom)
    {
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        if (atomic_num != 7 && atomic_num != 15 && atomic_num != 33 && atomic_num != 51)
        {
            return false;
        }

        FOR_BONDS_OF_ATOM(bond_iter, &atom)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (BondOrder(bond) != 1)
            {
                continue;
            }
            OpenBabel::OBAtom *neighbor = bond.GetNbrAtom(&atom);
            if (neighbor == nullptr || static_cast<int>(neighbor->GetAtomicNum()) != 6)
            {
                continue;
            }
            if (HasMultipleBondToAtomicNum(*neighbor, {8, 16, 34, 52}))
            {
                return true;
            }
        }
        return false;
    }

    bool IsNitrileLikeNitrogen(OpenBabel::OBAtom &atom)
    {
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        if (atomic_num != 7 && atomic_num != 15 && atomic_num != 33 && atomic_num != 51)
        {
            return false;
        }

        FOR_BONDS_OF_ATOM(bond_iter, &atom)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (BondOrder(bond) >= 3)
            {
                return true;
            }
        }
        return false;
    }

    std::pair<double, double> AtomDonorSupport(OpenBabel::OBAtom &atom)
    {
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        const int formal_charge = static_cast<int>(atom.GetFormalCharge());
        if (atomic_num == 1 || atom.IsMetal() || formal_charge > 0)
        {
            return {0.0, 0.0};
        }

        if (formal_charge < 0)
        {
            const double magnitude = static_cast<double>(std::abs(formal_charge));
            if (atomic_num == 9 || atomic_num == 17 || atomic_num == 35 || atomic_num == 53)
            {
                return {2.2 * magnitude, 0.0};
            }
            if (atomic_num == 8 || atomic_num == 16 || atomic_num == 34 || atomic_num == 52)
            {
                return {2.0 * magnitude, 0.0};
            }
            if (atomic_num == 7 || atomic_num == 15 || atomic_num == 33 || atomic_num == 51)
            {
                return {1.6 * magnitude, 0.0};
            }
            if (atomic_num == 6)
            {
                return {1.0 * magnitude, 0.0};
            }
            return {1.2 * magnitude, 0.0};
        }

        if (atomic_num == 9 || atomic_num == 17 || atomic_num == 35 || atomic_num == 53)
        {
            return {0.0, 0.0};
        }
        if (atomic_num == 8 || atomic_num == 16 || atomic_num == 34 || atomic_num == 52)
        {
            return {0.0, IsCarbonylLikeOxygen(atom) ? 0.7 : 1.0};
        }
        if (atomic_num == 7 || atomic_num == 15 || atomic_num == 33 || atomic_num == 51)
        {
            if (IsAmideLikeNitrogen(atom))
            {
                return {0.0, 0.2};
            }
            if (IsNitrileLikeNitrogen(atom))
            {
                return {0.0, 0.4};
            }
            if (atom.IsAromatic())
            {
                return {0.0, 0.7};
            }
            return {0.0, 0.9};
        }
        return {0.0, 0.0};
    }

    double ChargeLocalizationPenaltyForAtom(
        OpenBabel::OBAtom &atom,
        bool is_conjugated)
    {
        const int formal_charge = static_cast<int>(atom.GetFormalCharge());
        if (formal_charge == 0)
        {
            return 0.0;
        }

        const double magnitude = static_cast<double>(std::abs(formal_charge));
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        const bool is_aromatic = atom.IsAromatic();
        const int radical_electrons = static_cast<int>(atom.GetSpinMultiplicity()) % 2;

        const auto generic_penalty = [&]()
        {
            double penalty = 0.0;
            if (atomic_num == 1)
            {
                penalty = 4.0 * magnitude;
            }
            else if (formal_charge < 0)
            {
                if (atomic_num == 9 || atomic_num == 17 || atomic_num == 35 || atomic_num == 53)
                {
                    penalty = 0.2 * magnitude;
                }
                else if (atomic_num == 8 || atomic_num == 16 || atomic_num == 34 || atomic_num == 52)
                {
                    penalty = 0.3 * magnitude;
                }
                else if (atomic_num == 7 || atomic_num == 15 || atomic_num == 33 || atomic_num == 51)
                {
                    penalty = 0.6 * magnitude;
                }
                else if (atomic_num == 6)
                {
                    penalty = (is_conjugated || is_aromatic ? 1.5 : 4.0) * magnitude;
                }
                else
                {
                    penalty = (is_conjugated || is_aromatic ? 1.0 : 2.0) * magnitude;
                }
            }
            else
            {
                if (atomic_num == 7 || atomic_num == 15 || atomic_num == 33 || atomic_num == 51)
                {
                    penalty = 0.4 * magnitude;
                }
                else if (atomic_num == 8 || atomic_num == 16 || atomic_num == 34 || atomic_num == 52)
                {
                    penalty = 1.0 * magnitude;
                }
                else if (atomic_num == 6)
                {
                    penalty = (is_conjugated || is_aromatic ? 1.2 : 3.0) * magnitude;
                }
                else
                {
                    penalty = (is_conjugated || is_aromatic ? 0.8 : 1.8) * magnitude;
                }
            }
            if (!is_conjugated && !is_aromatic)
            {
                penalty += 0.5 * magnitude;
            }
            return penalty;
        };

        const auto main_group_charge_penalty = [&]() -> std::optional<double>
        {
            const auto info_it = molgr::kNonMetalDict.find(atomic_num);
            if (info_it == molgr::kNonMetalDict.end())
            {
                return std::nullopt;
            }

            const molgr::ElementInfo &element_info = info_it->second;
            const int local_electron_count =
                static_cast<int>(element_info.num_outer_electrons) - formal_charge +
                static_cast<int>(atom.GetTotalValence());
            const int total_valence = static_cast<int>(atom.GetTotalValence());
            const int default_valence = static_cast<int>(element_info.default_valence);
            const int shell_target = atomic_num == 1 ? 2 : 8;
            const int coordination_excess = std::max(total_valence - default_valence, 0);
            const int shell_deficiency = std::max(shell_target - local_electron_count, 0);
            const int shell_surplus = std::max(local_electron_count - shell_target, 0);
            const double ionization_norm = NormalizeValue(
                LookupMapValue(molgr::kNonMetalFirstIonizationEnergyEv, atomic_num),
                7.5,
                18.0,
                0.5);
            const double electronegativity_norm = NormalizeValue(
                LookupMapValue(molgr::kNonMetalPaulingElectronegativity, atomic_num),
                1.9,
                4.0,
                0.5);
            const double local_environment_penalty =
                !is_conjugated && !is_aromatic ? 0.05 * magnitude : 0.0;

            if (formal_charge > 0 && shell_deficiency > 0)
            {
                return magnitude *
                       (0.45 + 0.55 * static_cast<double>(shell_deficiency) + 0.9 * ionization_norm +
                        local_environment_penalty);
            }

            if (local_electron_count >= shell_target && coordination_excess > 0)
            {
                return magnitude *
                       (0.18 + 0.12 * static_cast<double>(coordination_excess) +
                        0.04 * static_cast<double>(shell_surplus) + local_environment_penalty);
            }

            if (formal_charge < 0 && local_electron_count >= shell_target)
            {
                return magnitude *
                       (0.10 + 0.60 * (1.0 - electronegativity_norm) + local_environment_penalty);
            }
            return std::nullopt;
        };

        double penalty = main_group_charge_penalty().value_or(generic_penalty());
        if (radical_electrons > 0 && atomic_num == 6 && !is_conjugated && !is_aromatic)
        {
            penalty += 2.0 * static_cast<double>(radical_electrons);
        }
        return penalty;
    }

    double RadicalLocalizationPenaltyForAtom(
        const OpenBabel::OBAtom &atom,
        bool is_conjugated)
    {
        const int radical_electrons = static_cast<int>(atom.GetSpinMultiplicity()) % 2;
        if (radical_electrons <= 0)
        {
            return 0.0;
        }

        const double magnitude = static_cast<double>(radical_electrons);
        const int atomic_num = static_cast<int>(atom.GetAtomicNum());
        const bool is_aromatic = atom.IsAromatic();

        if (atomic_num == 1)
        {
            return 4.0 * magnitude;
        }
        if (atomic_num == 6)
        {
            return (is_conjugated || is_aromatic ? 0.6 : 2.5) * magnitude;
        }
        if (atomic_num == 7 || atomic_num == 15 || atomic_num == 33 || atomic_num == 51)
        {
            return (is_conjugated || is_aromatic ? 1.0 : 2.5) * magnitude;
        }
        if (atomic_num == 8 || atomic_num == 16 || atomic_num == 34 || atomic_num == 52)
        {
            return (is_conjugated || is_aromatic ? 1.5 : 3.0) * magnitude;
        }
        return (is_conjugated || is_aromatic ? 1.2 : 2.5) * magnitude;
    }

    OrganicElectronicStateMetrics ComputeOrganicElectronicStateMetrics(const OpenBabel::OBMol &mol)
    {
        try
        {
            const auto topology_metrics = molgr::organic_topology::ComputeOrganicTopologyMetrics(mol);
            std::set<int> conjugated_atom_indices(
                topology_metrics.conjugated_atom_indices.begin(),
                topology_metrics.conjugated_atom_indices.end());

            OrganicElectronicStateMetrics metrics;
            metrics.aromatic_atom_count = topology_metrics.aromatic_atom_count;
            metrics.aromatic_ring_count = topology_metrics.aromatic_ring_count;
            metrics.conjugated_atom_count = topology_metrics.conjugated_atom_count;
            metrics.conjugated_bond_count = topology_metrics.conjugated_bond_count;
            metrics.max_conjugated_component_size = topology_metrics.max_conjugated_component_size;

            FOR_ATOMS_OF_MOL(atom_iter, const_cast<OpenBabel::OBMol &>(mol))
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                const int atom_idx = static_cast<int>(atom.GetIdx()) - 1;
                const bool is_conjugated =
                    conjugated_atom_indices.find(atom_idx) != conjugated_atom_indices.end();
                metrics.radical_localization_penalty +=
                    RadicalLocalizationPenaltyForAtom(atom, is_conjugated);
                metrics.charge_localization_penalty +=
                    ChargeLocalizationPenaltyForAtom(atom, is_conjugated);
            }
            return metrics;
        }
        catch (...)
        {
            OrganicElectronicStateMetrics metrics;
            metrics.radical_localization_penalty = std::numeric_limits<double>::infinity();
            metrics.charge_localization_penalty = std::numeric_limits<double>::infinity();
            return metrics;
        }
    }

    Point3D AtomCoordinates(const OpenBabel::OBAtom &atom)
    {
        return {atom.GetX(), atom.GetY(), atom.GetZ()};
    }

    Point3D MetalCoordinates(const molgr::MetalAtomPosition &metal_state)
    {
        return {metal_state.position_x, metal_state.position_y, metal_state.position_z};
    }

    double DistanceToMetal(
        const OpenBabel::OBAtom &atom,
        const molgr::MetalAtomPosition &metal_state)
    {
        const double dx = static_cast<double>(atom.GetX()) - metal_state.position_x;
        const double dy = static_cast<double>(atom.GetY()) - metal_state.position_y;
        const double dz = static_cast<double>(atom.GetZ()) - metal_state.position_z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    double DistanceWeight(double distance, double cutoff, double min_distance_angstrom)
    {
        if (distance <= 0.0 || distance >= cutoff)
        {
            return 0.0;
        }
        const double scaled = distance / cutoff;
        const double attenuation = std::max(0.0, 1.0 - scaled * scaled);
        return (attenuation * attenuation) / std::max(distance, min_distance_angstrom);
    }

    double DistancePointToSegment(
        const Point3D &point,
        const Point3D &segment_start,
        const Point3D &segment_end)
    {
        const double vx = segment_end.x - segment_start.x;
        const double vy = segment_end.y - segment_start.y;
        const double vz = segment_end.z - segment_start.z;
        const double seg_len_sq = vx * vx + vy * vy + vz * vz;
        if (seg_len_sq <= 1.0e-12)
        {
            const double dx = point.x - segment_start.x;
            const double dy = point.y - segment_start.y;
            const double dz = point.z - segment_start.z;
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        const double wx = point.x - segment_start.x;
        const double wy = point.y - segment_start.y;
        const double wz = point.z - segment_start.z;
        const double projection = std::max(
            0.0,
            std::min(1.0, (wx * vx + wy * vy + wz * vz) / seg_len_sq));
        const Point3D closest{
            segment_start.x + projection * vx,
            segment_start.y + projection * vy,
            segment_start.z + projection * vz};
        const double dx = point.x - closest.x;
        const double dy = point.y - closest.y;
        const double dz = point.z - closest.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    std::vector<CoordinationBlocker> BuildCoordinationBlockers(
        const OpenBabel::OBMol &mol,
        const molgr::config::MetalScoringConfig &config)
    {
        std::vector<CoordinationBlocker> blockers;
        FOR_ATOMS_OF_MOL(atom_iter, const_cast<OpenBabel::OBMol &>(mol))
        {
            OpenBabel::OBAtom &blocker = *atom_iter;
            const double blocker_radius =
                config.metal_access_radius_scale *
                    OpenBabel::OBElements::GetCovalentRad(static_cast<int>(blocker.GetAtomicNum())) +
                config.metal_access_clearance_angstrom;
            if (blocker_radius <= 0.0)
            {
                continue;
            }
            blockers.push_back(
                CoordinationBlocker{
                    static_cast<int>(blocker.GetIdx()),
                    AtomCoordinates(blocker),
                    blocker_radius});
        }
        return blockers;
    }

    bool HasUnobstructedCoordinationPathFromBlockers(
        int atom_idx,
        const Point3D &atom_coordinates,
        const Point3D &segment_start,
        const std::vector<CoordinationBlocker> &blockers)
    {
        for (const auto &blocker : blockers)
        {
            if (blocker.atom_idx == atom_idx)
            {
                continue;
            }
            if (DistancePointToSegment(blocker.coordinates, segment_start, atom_coordinates) <
                blocker.radius)
            {
                return false;
            }
        }
        return true;
    }

    MetalSiteEnvironmentProfile BuildMetalSiteEnvironmentProfile(
        const OpenBabel::OBMol &mol,
        const molgr::MetalAtomPosition &metal_state,
        const molgr::config::MetalScoringConfig &config)
    {
        const auto blockers = BuildCoordinationBlockers(mol, config);
        const Point3D segment_start = MetalCoordinates(metal_state);
        MetalSiteEnvironmentProfile profile;

        FOR_ATOMS_OF_MOL(atom_iter, const_cast<OpenBabel::OBMol &>(mol))
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (atom.IsMetal())
            {
                continue;
            }

            const int atom_idx = static_cast<int>(atom.GetIdx());
            const Point3D atom_coordinates = AtomCoordinates(atom);
            const double distance = DistanceToMetal(atom, metal_state);
            const double electrostatic_weight = DistanceWeight(
                distance,
                config.metal_local_potential_cutoff_angstrom,
                config.min_distance_angstrom);
            const double formal_charge = static_cast<double>(static_cast<int>(atom.GetFormalCharge()));
            if (electrostatic_weight > 0.0)
            {
                profile.electrostatic_support += -formal_charge * electrostatic_weight;
            }

            const double coordination_weight = DistanceWeight(
                distance,
                config.metal_donor_cutoff_angstrom,
                config.min_distance_angstrom);
            if (coordination_weight <= 0.0)
            {
                continue;
            }

            auto [atom_anionic_support, atom_neutral_support] = AtomDonorSupport(atom);
            const double atom_effective_donor_support =
                atom_anionic_support + config.local_neutral_donor_weight * atom_neutral_support;
            if (atom_effective_donor_support <= 0.0)
            {
                continue;
            }

            const double weighted_effective_support = atom_effective_donor_support * coordination_weight;
            if (HasUnobstructedCoordinationPathFromBlockers(
                    atom_idx,
                    atom_coordinates,
                    segment_start,
                    blockers))
            {
                profile.visible_anionic_donor_support += atom_anionic_support * coordination_weight;
                profile.visible_neutral_donor_support += atom_neutral_support * coordination_weight;
                profile.visible_effective_donor_support += weighted_effective_support;
                continue;
            }

            if (formal_charge < 0.0)
            {
                profile.obstructed_negative_effective_donor_support += weighted_effective_support;
            }
        }

        return profile;
    }

    MetalSiteHeuristicScore ScoreMetalSiteEnvironmentFromProfile(
        const MetalSiteEnvironmentProfile &profile,
        const molgr::MetalAtomPosition &metal_state,
        const molgr::config::MetalScoringConfig &config)
    {
        MetalSiteHeuristicScore score;
        score.electrostatic_support = profile.electrostatic_support;
        score.anionic_donor_support = profile.visible_anionic_donor_support;
        score.neutral_donor_support = profile.visible_neutral_donor_support;

        const double metal_valence = static_cast<double>(metal_state.valence);
        if (metal_valence >= 0.0 && profile.visible_effective_donor_support > 0.0)
        {
            score.visible_coordination_reward =
                config.visible_coordination_reward_weight *
                profile.visible_effective_donor_support;
            score.coordination_access_penalty -= score.visible_coordination_reward;
        }

        if (metal_valence < 0.0 && profile.visible_effective_donor_support > 0.0)
        {
            score.negative_metal_visible_coordination_penalty =
                config.negative_metal_visible_coordination_penalty_weight *
                std::max(std::abs(metal_valence), 1.0) * profile.visible_effective_donor_support;
            score.coordination_access_penalty += score.negative_metal_visible_coordination_penalty;
        }

        if (metal_valence > 0.0 && profile.obstructed_negative_effective_donor_support > 0.0)
        {
            score.obstructed_opposite_charge_penalty =
                config.obstructed_opposite_charge_penalty_weight *
                std::max(std::abs(metal_valence), 1.0) *
                profile.obstructed_negative_effective_donor_support;
            score.coordination_access_penalty += score.obstructed_opposite_charge_penalty;
        }

        const double target_valence = static_cast<double>(std::max(metal_state.valence, 0));
        const double electrostatic_target =
            config.local_potential_target_per_valence * target_valence;
        const double electrostatic_under =
            std::max(electrostatic_target - profile.electrostatic_support, 0.0);
        const double electrostatic_over =
            std::max(profile.electrostatic_support - electrostatic_target, 0.0);
        score.electrostatic_penalty =
            electrostatic_under + config.local_potential_oversupport_weight * electrostatic_over;

        const double effective_donor_support =
            score.anionic_donor_support + config.local_neutral_donor_weight * score.neutral_donor_support;
        const double donor_target = config.local_donor_target_per_valence * target_valence;
        const double donor_under = std::max(donor_target - effective_donor_support, 0.0);
        const double donor_over = std::max(effective_donor_support - donor_target, 0.0);
        score.donor_penalty = donor_under + config.local_donor_oversupport_weight * donor_over;
        return score;
    }

    double SameElementValenceSpreadPenalty(
        const std::vector<molgr::MetalAtomPosition> &metal_states,
        const molgr::config::MetalScoringConfig &config)
    {
        std::map<std::string, std::vector<int>> grouped_valences;
        for (const auto &metal_state : metal_states)
        {
            grouped_valences[metal_state.symbol].push_back(metal_state.valence);
        }

        double penalty = 0.0;
        for (const auto &entry : grouped_valences)
        {
            const std::vector<int> &valences = entry.second;
            if (valences.size() < 2)
            {
                continue;
            }
            const auto [min_it, max_it] = std::minmax_element(valences.begin(), valences.end());
            penalty += config.same_element_valence_spread_weight *
                       static_cast<double>(*max_it - *min_it);
        }
        return penalty;
    }

    double MetalStateAssignmentPenaltyValue(const molgr::MetalAtomPosition &metal_state)
    {
        double penalty = 0.0;
        if (metal_state.valence <= 0)
        {
            penalty += 10.0 * std::max(std::abs(metal_state.valence), 1);
        }

        const auto prior_it = molgr::kMetalValencePrior.find(metal_state.symbol);
        const auto minor_it = molgr::kMetalValenceMinor.find(metal_state.symbol);
        const auto contains = [](const std::vector<int> &values, int target)
        {
            return std::find(values.begin(), values.end(), target) != values.end();
        };

        const bool in_prior =
            prior_it != molgr::kMetalValencePrior.end() &&
            contains(prior_it->second, metal_state.valence);
        const bool in_minor =
            minor_it != molgr::kMetalValenceMinor.end() &&
            contains(minor_it->second, metal_state.valence);

        if (!in_prior)
        {
            penalty += in_minor ? 10.0 : 20.0;
        }
        return penalty;
    }

    void AnnotateMetalEnvironmentConsistency(
        molgr::state::MetalCandidateState *candidate,
        const molgr::config::MolGRConfig &config)
    {
        if (candidate == nullptr || !candidate->no_metal_state)
        {
            throw std::runtime_error(
                "MetalCandidateState requires no_metal_state before metal scoring");
        }

        const auto &metal_scoring_config = config.metal_scoring;
        double total_prior_penalty = 0.0;
        double total_coordination_access_penalty = 0.0;
        double total_electrostatic_penalty = 0.0;
        double total_donor_penalty = 0.0;
        for (const auto &metal_state : candidate->metal_states)
        {
            const double prior_penalty = MetalStateAssignmentPenaltyValue(metal_state);
            const auto site_score = ScoreMetalSiteEnvironmentFromProfile(
                BuildMetalSiteEnvironmentProfile(
                    candidate->no_metal_state->Mol(),
                    metal_state,
                    metal_scoring_config),
                metal_state,
                metal_scoring_config);
            total_prior_penalty += prior_penalty;
            total_coordination_access_penalty += site_score.coordination_access_penalty;
            total_electrostatic_penalty += site_score.electrostatic_penalty;
            total_donor_penalty += site_score.donor_penalty;
        }

        candidate->metadata["metal_prior_penalty"] = total_prior_penalty;
        candidate->metadata["metal_coordination_access_penalty"] = total_coordination_access_penalty;
        candidate->metadata["metal_same_element_valence_spread_penalty"] =
            SameElementValenceSpreadPenalty(candidate->metal_states, metal_scoring_config);
        candidate->metadata["metal_electrostatic_penalty"] = total_electrostatic_penalty;
        candidate->metadata["metal_donor_penalty"] = total_donor_penalty;
    }

    void AnnotateOrganicElectronicStateConsistency(molgr::state::MetalCandidateState *candidate)
    {
        if (candidate == nullptr || !candidate->no_metal_state)
        {
            throw std::runtime_error(
                "MetalCandidateState requires no_metal_state before organic-state scoring");
        }

        const auto metrics = ComputeOrganicElectronicStateMetrics(candidate->no_metal_state->Mol());
        candidate->metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count;
        candidate->metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count;
        candidate->metadata["organic_conjugated_atom_count"] = metrics.conjugated_atom_count;
        candidate->metadata["organic_conjugated_bond_count"] = metrics.conjugated_bond_count;
        candidate->metadata["organic_max_conjugated_component_size"] =
            metrics.max_conjugated_component_size;
        candidate->metadata["organic_radical_localization_penalty"] =
            metrics.radical_localization_penalty;
        candidate->metadata["organic_charge_localization_penalty"] =
            metrics.charge_localization_penalty;
    }

    OrganicElectronicStateSelectionContext BuildOrganicElectronicStateSelectionContext(
        const std::vector<molgr::state::MetalCandidateState> &candidates,
        const std::vector<std::size_t> &candidate_indices)
    {
        OrganicElectronicStateSelectionContext context;
        for (const std::size_t candidate_index : candidate_indices)
        {
            const auto &candidate = candidates[candidate_index];
            context.max_aromatic_atom_count = std::max(
                context.max_aromatic_atom_count,
                molgr::metal::scoring::MetadataInt(candidate, "organic_aromatic_atom_count"));
            context.max_aromatic_ring_count = std::max(
                context.max_aromatic_ring_count,
                molgr::metal::scoring::MetadataInt(candidate, "organic_aromatic_ring_count"));
            context.max_conjugated_atom_count = std::max(
                context.max_conjugated_atom_count,
                molgr::metal::scoring::MetadataInt(candidate, "organic_conjugated_atom_count"));
            context.max_conjugated_component_size = std::max(
                context.max_conjugated_component_size,
                molgr::metal::scoring::MetadataInt(candidate, "organic_max_conjugated_component_size"));
        }
        return context;
    }

    OrganicStateKey OrganicElectronicStateKey(
        molgr::state::MetalCandidateState *candidate,
        const OrganicElectronicStateSelectionContext &selection_context)
    {
        const int aromatic_atom_loss = selection_context.max_aromatic_atom_count -
                                       molgr::metal::scoring::MetadataInt(
                                           *candidate,
                                           "organic_aromatic_atom_count");
        const int aromatic_ring_loss = selection_context.max_aromatic_ring_count -
                                       molgr::metal::scoring::MetadataInt(
                                           *candidate,
                                           "organic_aromatic_ring_count");
        const int conjugated_atom_loss = selection_context.max_conjugated_atom_count -
                                         molgr::metal::scoring::MetadataInt(
                                             *candidate,
                                             "organic_conjugated_atom_count");
        const int max_conjugated_component_loss =
            selection_context.max_conjugated_component_size -
            molgr::metal::scoring::MetadataInt(*candidate, "organic_max_conjugated_component_size");
        const double radical_localization_penalty =
            molgr::metal::scoring::MetadataDouble(
                *candidate,
                "organic_radical_localization_penalty",
                std::numeric_limits<double>::infinity());
        const double charge_localization_penalty =
            molgr::metal::scoring::MetadataDouble(
                *candidate,
                "organic_charge_localization_penalty",
                std::numeric_limits<double>::infinity());

        candidate->metadata["organic_aromatic_atom_loss"] = aromatic_atom_loss;
        candidate->metadata["organic_aromatic_ring_loss"] = aromatic_ring_loss;
        candidate->metadata["organic_conjugated_atom_loss"] = conjugated_atom_loss;
        candidate->metadata["organic_max_conjugated_component_loss"] =
            max_conjugated_component_loss;

        return {
            aromatic_ring_loss,
            max_conjugated_component_loss,
            charge_localization_penalty,
            aromatic_atom_loss,
            conjugated_atom_loss,
            radical_localization_penalty};
    }

    RawCandidateSelectionKey RawCandidateSelectionTuple(
        molgr::state::MetalCandidateState *candidate,
        const OrganicElectronicStateSelectionContext &organic_selection_context,
        double best_force_field_score,
        const molgr::config::MolGRConfig &config)
    {
        const double score_value =
            candidate->score.has_value() ? *candidate->score : candidate->CombinedScore(config);
        const int organic_bucket =
            molgr::metal::scoring::OrganicScoreBucketIndex(score_value, best_force_field_score, config);
        candidate->metadata["organic_score_bucket"] = organic_bucket;

        const OrganicStateKey organic_state_key =
            OrganicElectronicStateKey(candidate, organic_selection_context);
        return {
            std::get<0>(organic_state_key),
            std::get<1>(organic_state_key),
            std::get<2>(organic_state_key),
            std::get<3>(organic_state_key),
            std::get<4>(organic_state_key),
            std::get<5>(organic_state_key),
            molgr::metal::scoring::MetadataDouble(*candidate, "metal_coordination_access_penalty", 0.0),
            molgr::metal::scoring::MetadataDouble(
                *candidate,
                "metal_same_element_valence_spread_penalty",
                0.0),
            molgr::metal::scoring::MetadataDouble(*candidate, "metal_electrostatic_penalty", 0.0),
            molgr::metal::scoring::MetadataDouble(*candidate, "metal_donor_penalty", 0.0),
            molgr::metal::scoring::MetadataDouble(
                *candidate,
                "metal_prior_penalty",
                molgr::metal::scoring::MetadataDouble(*candidate, "metal_assignment_rank", 0.0)),
            organic_bucket,
            score_value,
            molgr::metal::scoring::CandidateCombinationIndex(*candidate)};
    }

    std::array<double, kWeightedSelectionFieldCount> WeightedMetricValues(
        const RawCandidateSelectionKey &raw_key)
    {
        return {
            static_cast<double>(std::get<0>(raw_key)),
            static_cast<double>(std::get<1>(raw_key)),
            std::get<2>(raw_key),
            static_cast<double>(std::get<3>(raw_key)),
            static_cast<double>(std::get<4>(raw_key)),
            std::get<5>(raw_key),
            std::get<6>(raw_key),
            std::get<7>(raw_key),
            std::get<8>(raw_key),
            std::get<9>(raw_key),
            std::get<10>(raw_key)};
    }

    std::array<double, kWeightedSelectionFieldCount> NormalizeWeightedSelectionValues(
        const std::vector<double> &raw_values,
        double default_value)
    {
        std::array<double, kWeightedSelectionFieldCount> normalized{};
        if (raw_values.empty())
        {
            normalized.fill(default_value);
            return normalized;
        }

        const double last_value = std::max(raw_values.back(), 0.0);
        normalized.fill(last_value);
        for (std::size_t idx = 0; idx < std::min(raw_values.size(), kWeightedSelectionFieldCount); ++idx)
        {
            normalized[idx] = std::max(raw_values[idx], 0.0);
        }
        return normalized;
    }

    WeightedSelectionContext BuildWeightedSelectionContext(
        const std::vector<RawCandidateSelectionKey> &raw_keys,
        const molgr::config::MolGRConfig &config)
    {
        WeightedSelectionContext context;
        context.weights = NormalizeWeightedSelectionValues(
            config.metal_scoring.selection_weight_values,
            1.0);
        context.scales = NormalizeWeightedSelectionValues(
            config.metal_scoring.selection_scale_values,
            1.0);
        context.best_values.fill(std::numeric_limits<double>::infinity());
        for (const auto &raw_key : raw_keys)
        {
            const auto metric_values = WeightedMetricValues(raw_key);
            for (std::size_t idx = 0; idx < kWeightedSelectionFieldCount; ++idx)
            {
                context.best_values[idx] =
                    std::min(context.best_values[idx], metric_values[idx]);
            }
        }
        return context;
    }

    WeightedSelectionKey BuildWeightedSelectionKey(
        molgr::state::MetalCandidateState *candidate,
        const RawCandidateSelectionKey &raw_key,
        const WeightedSelectionContext &context)
    {
        const auto metric_values = WeightedMetricValues(raw_key);
        double weighted_score = 0.0;
        for (std::size_t idx = 0; idx < kWeightedSelectionFieldCount; ++idx)
        {
            const double regret = std::max(
                0.0,
                (metric_values[idx] - context.best_values[idx]) /
                    std::max(context.scales[idx], 1e-12));
            weighted_score += regret * context.weights[idx];
        }
        candidate->metadata["weighted_selection_score"] = weighted_score;
        return {
            weighted_score,
            std::get<0>(raw_key),
            std::get<1>(raw_key),
            std::get<2>(raw_key),
            std::get<3>(raw_key),
            std::get<4>(raw_key),
            std::get<5>(raw_key),
            std::get<6>(raw_key),
            std::get<7>(raw_key),
            std::get<8>(raw_key),
            std::get<9>(raw_key),
            std::get<10>(raw_key),
            std::get<11>(raw_key),
            std::get<12>(raw_key),
            std::get<13>(raw_key)};
    }
}

namespace molgr
{
    namespace metal
    {
        namespace scoring
        {
            int CandidateCombinationIndex(const molgr::state::MetalCandidateState &candidate)
            {
                const auto metadata_it = candidate.metadata.find("combination_index");
                if (metadata_it == candidate.metadata.end())
                {
                    return 0;
                }
                if (const auto *value = std::get_if<int>(&metadata_it->second))
                {
                    return *value;
                }
                if (const auto *value = std::get_if<double>(&metadata_it->second))
                {
                    return static_cast<int>(*value);
                }
                return 0;
            }

            int MetadataInt(
                const molgr::state::MetalCandidateState &candidate,
                const std::string &key,
                int fallback)
            {
                const auto metadata_it = candidate.metadata.find(key);
                if (metadata_it == candidate.metadata.end())
                {
                    return fallback;
                }
                if (const auto *value = std::get_if<int>(&metadata_it->second))
                {
                    return *value;
                }
                if (const auto *value = std::get_if<double>(&metadata_it->second))
                {
                    return static_cast<int>(*value);
                }
                return fallback;
            }

            double MetadataDouble(
                const molgr::state::MetalCandidateState &candidate,
                const std::string &key,
                double fallback)
            {
                const auto metadata_it = candidate.metadata.find(key);
                if (metadata_it == candidate.metadata.end())
                {
                    return fallback;
                }
                if (const auto *value = std::get_if<double>(&metadata_it->second))
                {
                    return *value;
                }
                if (const auto *value = std::get_if<int>(&metadata_it->second))
                {
                    return static_cast<double>(*value);
                }
                return fallback;
            }

            int OrganicScoreBucketIndex(
                double score_value,
                double best_score,
                const molgr::config::MolGRConfig &config)
            {
                if (score_value <= best_score)
                {
                    return 0;
                }
                const double ratio =
                    std::max(config.metal_scoring.organic_score_bucket_relative_ratio, 1e-12);
                const double baseline_scale = std::max(std::abs(best_score), 1.0);
                const double relative_excess = (score_value - best_score) / baseline_scale;
                return static_cast<int>(std::floor(relative_excess / ratio));
            }

            bool PassesOrganicForceFieldGuard(
                double score_value,
                double best_score,
                const molgr::config::MolGRConfig &config)
            {
                const double hard_max_ratio = config.metal_scoring.organic_force_field_hard_max_ratio;
                if (hard_max_ratio <= 0.0 || best_score <= 0.0)
                {
                    return true;
                }
                return score_value <= best_score * hard_max_ratio;
            }

            std::optional<molgr::state::MetalCandidateState> SelectBestCandidate(
                const std::vector<molgr::state::MetalCandidateState> &candidates,
                const molgr::config::MolGRConfig &config)
            {
                if (candidates.empty())
                {
                    return std::nullopt;
                }

                double best_force_field_score = std::numeric_limits<double>::infinity();
                for (const auto &candidate : candidates)
                {
                    const double score_value =
                        candidate.score.has_value() ? *candidate.score : candidate.CombinedScore(config);
                    best_force_field_score = std::min(best_force_field_score, score_value);
                }

                std::vector<molgr::state::MetalCandidateState> eligible_candidates;
                eligible_candidates.reserve(candidates.size());
                for (auto candidate : candidates)
                {
                    const double score_value =
                        candidate.score.has_value() ? *candidate.score : candidate.CombinedScore(config);
                    const bool passes_force_field_guard =
                        PassesOrganicForceFieldGuard(score_value, best_force_field_score, config);
                    candidate.metadata["passes_organic_force_field_guard"] = passes_force_field_guard;
                    if (passes_force_field_guard)
                    {
                        eligible_candidates.push_back(std::move(candidate));
                    }
                }
                if (eligible_candidates.empty())
                {
                    return std::nullopt;
                }

                std::vector<std::size_t> eligible_indices(eligible_candidates.size());
                for (std::size_t idx = 0; idx < eligible_indices.size(); ++idx)
                {
                    eligible_indices[idx] = idx;
                }
                const OrganicElectronicStateSelectionContext organic_selection_context =
                    BuildOrganicElectronicStateSelectionContext(eligible_candidates, eligible_indices);

                std::vector<RawCandidateSelectionKey> raw_keys;
                raw_keys.reserve(eligible_candidates.size());
                for (auto &candidate : eligible_candidates)
                {
                    raw_keys.push_back(
                        RawCandidateSelectionTuple(
                            &candidate,
                            organic_selection_context,
                            best_force_field_score,
                            config));
                }
                const WeightedSelectionContext weighted_context =
                    BuildWeightedSelectionContext(raw_keys, config);

                std::optional<WeightedSelectionKey> best_selection_key;
                std::optional<molgr::state::MetalCandidateState> best_candidate;
                for (std::size_t idx = 0; idx < eligible_candidates.size(); ++idx)
                {
                    auto selection_key = BuildWeightedSelectionKey(
                        &eligible_candidates[idx],
                        raw_keys[idx],
                        weighted_context);
                    if (best_selection_key.has_value() && selection_key >= *best_selection_key)
                    {
                        continue;
                    }
                    best_selection_key = selection_key;
                    best_candidate = eligible_candidates[idx];
                }
                return best_candidate;
            }

            molgr::state::MetalCandidateState ScoreCandidateWithNoMetalState(
                const molgr::state::MetalCandidateState &candidate,
                const std::shared_ptr<molgr::state::ReconstructionState> &no_metal_state,
                const molgr::config::MolGRConfig &config)
            {
                auto machine = molgr::state::MetalCandidateStateMachine::FromCandidateState(candidate);
                machine.SetNoMetalState("reconstruct_no_metal", no_metal_state);
                machine.Annotate("score_candidate");
                auto scored_candidate = machine.Freeze();
                const double score = scored_candidate.CombinedScore(config);
                scored_candidate.score = score;
                scored_candidate.metadata["score"] = score;
                AnnotateOrganicElectronicStateConsistency(&scored_candidate);
                AnnotateMetalEnvironmentConsistency(&scored_candidate, config);
                return scored_candidate;
            }
        }
    }
}
