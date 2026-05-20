#include "molgr/utils/metals/scoring.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/organic_topology.h"
#include "molgr/utils/scoring.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/obfunctions.h>
#include <openbabel/obiter.h>

#include <algorithm>
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

    struct OrganicElectronicStateMetrics
    {
        int aromatic_atom_count = 0;
        int aromatic_ring_count = 0;
        double aromatic_stability_score = 0.0;
        int conjugated_atom_count = 0;
        int conjugated_bond_count = 0;
        int max_conjugated_component_size = 0;
        double radical_localization_penalty = 0.0;
        double charge_localization_penalty = 0.0;
    };

    struct NegativeMetalDiscordance
    {
        int count = 0;
        bool has_outer_sphere_cation_exception = false;
        bool has_positive_metal_counterion_exception = false;
    };

    constexpr double kNegativeMetalDiscordancePenalty = 0.5;

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

    bool IsInnerVisibleDiradicalDiscordanceAtom(const OpenBabel::OBAtom &atom)
    {
        switch (atom.GetAtomicNum())
        {
        case 15:
        case 16:
        case 17:
        case 35:
        case 53:
            return false;
        default:
            return static_cast<int>(atom.GetSpinMultiplicity()) >= 2;
        }
    }

    int BondOrder(const OpenBabel::OBBond &bond)
    {
        return bond.IsAromatic() ? 2 : static_cast<int>(bond.GetBondOrder());
    }

    OpenBabel::OBAtom *OtherBondAtom(OpenBabel::OBBond &bond, OpenBabel::OBAtom &atom)
    {
        return bond.GetNbrAtom(&atom);
    }

    OpenBabel::OBBond *BondBetweenAtoms(OpenBabel::OBAtom &lhs, OpenBabel::OBAtom &rhs)
    {
        const int rhs_idx = static_cast<int>(rhs.GetIdx());
        FOR_BONDS_OF_ATOM(bond_iter, &lhs)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            OpenBabel::OBAtom *other = OtherBondAtom(bond, lhs);
            if (other != nullptr && static_cast<int>(other->GetIdx()) == rhs_idx)
            {
                return &bond;
            }
        }
        return nullptr;
    }

    bool HasConjugatedBridgeBetweenChargedCarbons(
        OpenBabel::OBAtom &begin_atom,
        OpenBabel::OBAtom &end_atom)
    {
        FOR_BONDS_OF_ATOM(begin_bond_iter, &begin_atom)
        {
            OpenBabel::OBBond &begin_bond = *begin_bond_iter;
            if (BondOrder(begin_bond) != 1)
            {
                continue;
            }
            OpenBabel::OBAtom *begin_neighbor = OtherBondAtom(begin_bond, begin_atom);
            if (begin_neighbor == nullptr || begin_neighbor->GetAtomicNum() == 1)
            {
                continue;
            }
            FOR_BONDS_OF_ATOM(middle_bond_iter, begin_neighbor)
            {
                OpenBabel::OBBond &middle_bond = *middle_bond_iter;
                if (static_cast<int>(middle_bond.GetIdx()) == static_cast<int>(begin_bond.GetIdx()) ||
                    BondOrder(middle_bond) != 2)
                {
                    continue;
                }
                OpenBabel::OBAtom *middle_neighbor = OtherBondAtom(middle_bond, *begin_neighbor);
                if (middle_neighbor == nullptr ||
                    static_cast<int>(middle_neighbor->GetIdx()) == static_cast<int>(begin_atom.GetIdx()))
                {
                    continue;
                }
                OpenBabel::OBBond *bridge_end_bond = BondBetweenAtoms(*middle_neighbor, end_atom);
                if (bridge_end_bond != nullptr && BondOrder(*bridge_end_bond) == 1)
                {
                    return true;
                }
            }
        }
        return false;
    }

    int InnerVisibleConjugatedChargedCarbonPairCount(
        OpenBabel::OBMol &mol,
        const std::set<int> &visible_inner_atom_indices)
    {
        std::vector<OpenBabel::OBAtom *> charged_carbons;
        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (visible_inner_atom_indices.find(static_cast<int>(atom.GetIdx())) ==
                    visible_inner_atom_indices.end() ||
                atom.GetAtomicNum() != 6 ||
                atom.GetFormalCharge() == 0)
            {
                continue;
            }
            charged_carbons.push_back(&atom);
        }

        int count = 0;
        for (std::size_t i = 0; i < charged_carbons.size(); ++i)
        {
            for (std::size_t j = i + 1; j < charged_carbons.size(); ++j)
            {
                if ((charged_carbons[i]->GetFormalCharge() > 0) !=
                    (charged_carbons[j]->GetFormalCharge() > 0))
                {
                    continue;
                }
                if (BondBetweenAtoms(*charged_carbons[i], *charged_carbons[j]) != nullptr)
                {
                    continue;
                }
                if (HasConjugatedBridgeBetweenChargedCarbons(*charged_carbons[i], *charged_carbons[j]))
                {
                    ++count;
                }
            }
        }
        return count;
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

    OrganicElectronicStateMetrics ComputeOrganicElectronicStateMetrics(
        const OpenBabel::OBMol &mol,
        const molgr::config::MolGRConfig &config)
    {
        try
        {
            const auto topology_metrics = molgr::organic_topology::ComputeOrganicTopologyMetrics(
                mol,
                config.organic_topology);
            std::set<int> conjugated_atom_indices(
                topology_metrics.conjugated_atom_indices.begin(),
                topology_metrics.conjugated_atom_indices.end());

            OrganicElectronicStateMetrics metrics;
            metrics.aromatic_atom_count = topology_metrics.aromatic_atom_count;
            metrics.aromatic_ring_count = topology_metrics.aromatic_ring_count;
            metrics.aromatic_stability_score = topology_metrics.aromatic_stability_score;
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

    int ChargeSign(int charge)
    {
        if (charge > 0)
        {
            return 1;
        }
        if (charge < 0)
        {
            return -1;
        }
        return 0;
    }

    int NearestNonzeroMetalChargeSignToBond(
        const OpenBabel::OBAtom &begin_atom,
        const OpenBabel::OBAtom &end_atom,
        const std::vector<molgr::MetalAtomPosition> &metal_states)
    {
        const double midpoint_x =
            (static_cast<double>(begin_atom.GetX()) + static_cast<double>(end_atom.GetX())) * 0.5;
        const double midpoint_y =
            (static_cast<double>(begin_atom.GetY()) + static_cast<double>(end_atom.GetY())) * 0.5;
        const double midpoint_z =
            (static_cast<double>(begin_atom.GetZ()) + static_cast<double>(end_atom.GetZ())) * 0.5;
        double best_distance_sq = std::numeric_limits<double>::infinity();
        int best_charge_sign = 0;
        for (const auto &metal_state : metal_states)
        {
            const int metal_charge_sign = ChargeSign(metal_state.valence);
            if (metal_charge_sign == 0)
            {
                continue;
            }
            const double dx = midpoint_x - metal_state.position_x;
            const double dy = midpoint_y - metal_state.position_y;
            const double dz = midpoint_z - metal_state.position_z;
            const double distance_sq = dx * dx + dy * dy + dz * dz;
            if (distance_sq < best_distance_sq)
            {
                best_distance_sq = distance_sq;
                best_charge_sign = metal_charge_sign;
            }
        }
        return best_charge_sign;
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

    double CoordinationRadiusCutoff(
        const OpenBabel::OBAtom &atom,
        const molgr::MetalAtomPosition &metal_state,
        const molgr::config::MetalScoringConfig &config)
    {
        const double atom_radius =
            OpenBabel::OBElements::GetCovalentRad(static_cast<int>(atom.GetAtomicNum()));
        const double metal_radius =
            OpenBabel::OBElements::GetCovalentRad(static_cast<int>(metal_state.element_idx));
        return config.metal_access_radius_scale * (atom_radius + metal_radius) +
               config.metal_coordination_extra_tolerance_angstrom;
    }

    bool IsInnerSphereAtom(
        const OpenBabel::OBAtom &atom,
        const molgr::MetalAtomPosition &metal_state,
        const molgr::config::MetalScoringConfig &config)
    {
        const double distance = DistanceToMetal(atom, metal_state);
        return distance > 0.0 && distance <= CoordinationRadiusCutoff(atom, metal_state, config);
    }

    bool IsVisibleToMetalAtom(
        OpenBabel::OBAtom &atom,
        const molgr::MetalAtomPosition &metal_state,
        const std::vector<CoordinationBlocker> &blockers)
    {
        const double distance = DistanceToMetal(atom, metal_state);
        if (distance <= 0.0)
        {
            return false;
        }
        return HasUnobstructedCoordinationPathFromBlockers(
            static_cast<int>(atom.GetIdx()),
            AtomCoordinates(atom),
            MetalCoordinates(metal_state),
            blockers);
    }

    bool HasOuterSphereProton(
        OpenBabel::OBMol &mol,
        const std::vector<molgr::MetalAtomPosition> &metal_states,
        const molgr::config::MetalScoringConfig &config)
    {
        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (atom.IsMetal() || atom.GetAtomicNum() != 1 || atom.GetFormalCharge() <= 0)
            {
                continue;
            }

            bool is_inner_sphere_to_any_metal = false;
            for (const auto &metal_state : metal_states)
            {
                if (IsInnerSphereAtom(atom, metal_state, config))
                {
                    is_inner_sphere_to_any_metal = true;
                    break;
                }
            }
            if (!is_inner_sphere_to_any_metal)
            {
                return true;
            }
        }
        return false;
    }

    NegativeMetalDiscordance NegativeMetalDiscordanceCount(
        OpenBabel::OBMol &mol,
        const std::vector<molgr::MetalAtomPosition> &metal_states,
        const molgr::config::MetalScoringConfig &config)
    {
        int negative_metal_count = 0;
        bool has_positive_metal_counterion = false;
        for (const auto &metal_state : metal_states)
        {
            if (metal_state.valence < 0)
            {
                negative_metal_count += std::abs(metal_state.valence);
            }
            else if (metal_state.valence > 0)
            {
                has_positive_metal_counterion = true;
            }
        }
        if (negative_metal_count == 0)
        {
            return {};
        }

        const bool has_outer_sphere_proton = HasOuterSphereProton(mol, metal_states, config);
        if (has_outer_sphere_proton || has_positive_metal_counterion)
        {
            return NegativeMetalDiscordance{
                0,
                has_outer_sphere_proton,
                has_positive_metal_counterion};
        }
        return NegativeMetalDiscordance{negative_metal_count, false, false};
    }

    bool IsUnsaturatedOrganicCation(OpenBabel::OBAtom &atom)
    {
        if (atom.IsMetal() || static_cast<int>(atom.GetFormalCharge()) <= 0)
        {
            return false;
        }

        const int total_degree = static_cast<int>(atom.GetTotalDegree());
        const int total_valence = static_cast<int>(atom.GetTotalValence());
        const int typical_valence = OpenBabel::GetTypicalValence(
            static_cast<int>(atom.GetAtomicNum()),
            total_valence,
            static_cast<int>(atom.GetFormalCharge()));
        return total_degree < total_valence || total_valence < typical_valence;
    }

    bool HasAdjacentFormalChargeCancellation(OpenBabel::OBAtom &atom)
    {
        const int formal_charge = static_cast<int>(atom.GetFormalCharge());
        if (formal_charge == 0)
        {
            return false;
        }
        int adjacent_charge = 0;
        FOR_BONDS_OF_ATOM(bond_iter, &atom)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            OpenBabel::OBAtom *neighbor = OtherBondAtom(bond, atom);
            if (neighbor != nullptr && !neighbor->IsMetal())
            {
                adjacent_charge += static_cast<int>(neighbor->GetFormalCharge());
            }
        }
        return formal_charge + adjacent_charge == 0;
    }

    bool IsLocallyZwitterionicOrganicCation(OpenBabel::OBAtom &atom)
    {
        return static_cast<int>(atom.GetFormalCharge()) > 0 &&
               HasAdjacentFormalChargeCancellation(atom);
    }

    int ZeroValentMetalsWithOrganicCationCount(
        OpenBabel::OBMol &mol,
        const std::vector<molgr::MetalAtomPosition> &metal_states,
        int total_charge)
    {
        if (total_charge > 0)
        {
            return 0;
        }
        if (metal_states.empty())
        {
            return 0;
        }
        for (const auto &metal_state : metal_states)
        {
            if (metal_state.valence != 0)
            {
                return 0;
            }
        }
        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (!atom.IsMetal() &&
                atom.GetFormalCharge() > 0 &&
                !IsLocallyZwitterionicOrganicCation(atom))
            {
                return 1;
            }
        }
        return 0;
    }

    int NonnegativeMetalUnsaturatedOrganicCationCount(
        OpenBabel::OBMol &mol,
        const std::vector<molgr::MetalAtomPosition> &metal_states)
    {
        bool has_nonnegative_metal = false;
        for (const auto &metal_state : metal_states)
        {
            if (metal_state.valence >= 0)
            {
                has_nonnegative_metal = true;
                break;
            }
        }
        if (!has_nonnegative_metal)
        {
            return 0;
        }
        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            OpenBabel::OBAtom &atom = *atom_iter;
            if (IsUnsaturatedOrganicCation(atom) &&
                !IsLocallyZwitterionicOrganicCation(atom))
            {
                return 1;
            }
        }
        return 0;
    }

    double AnnotateCandidateDiscordanceFeatures(
        molgr::state::MetalCandidateState *candidate,
        const molgr::config::MolGRConfig &config)
    {
        if (candidate == nullptr || !candidate->no_metal_state)
        {
            throw std::runtime_error(
                "MetalCandidateState requires no_metal_state before discordance scoring");
        }

        auto &mol = const_cast<OpenBabel::OBMol &>(candidate->no_metal_state->Mol());
        const auto blockers = BuildCoordinationBlockers(mol, config.metal_scoring);

        int inner_visible_diradical_count = 0;
        int inner_visible_same_sign_charge_count = 0;
        int outer_or_invisible_adjacent_same_sign_double_charge_count = 0;
        int outer_or_invisible_adjacent_opposite_sign_double_charge_count = 0;
        int outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count = 0;
        int inner_visible_adjacent_carbanion_pair_count = 0;
        std::set<int> visible_inner_atom_indices;

        for (const auto &metal_state : candidate->metal_states)
        {
            const int metal_charge_sign = ChargeSign(metal_state.valence);
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (atom.IsMetal())
                {
                    continue;
                }
                if (!IsInnerSphereAtom(atom, metal_state, config.metal_scoring))
                {
                    continue;
                }
                if (!IsVisibleToMetalAtom(atom, metal_state, blockers))
                {
                    continue;
                }

                const int atom_idx = static_cast<int>(atom.GetIdx());
                visible_inner_atom_indices.insert(atom_idx);
                if (IsInnerVisibleDiradicalDiscordanceAtom(atom))
                {
                    ++inner_visible_diradical_count;
                }

                const int formal_charge = static_cast<int>(atom.GetFormalCharge());
                const int atom_charge_sign = ChargeSign(formal_charge);
                const bool has_inner_same_sign_charge =
                    (metal_charge_sign == 0 && formal_charge > 0) ||
                    (metal_charge_sign != 0 && atom_charge_sign == metal_charge_sign);
                if (has_inner_same_sign_charge &&
                    !HasAdjacentFormalChargeCancellation(atom))
                {
                    ++inner_visible_same_sign_charge_count;
                }
            }
        }

        int outer_or_invisible_adjacent_double_charge_count = 0;
        const int inner_visible_conjugated_carbanion_pair_count =
            InnerVisibleConjugatedChargedCarbonPairCount(mol, visible_inner_atom_indices);
        const auto negative_metal_discordance =
            NegativeMetalDiscordanceCount(mol, candidate->metal_states, config.metal_scoring);
        int total_charge = candidate->no_metal_state->total_charge;
        for (const auto &metal_state : candidate->metal_states)
        {
            total_charge += metal_state.valence;
        }
        const int zero_valent_metals_with_organic_cation_count =
            ZeroValentMetalsWithOrganicCationCount(mol, candidate->metal_states, total_charge);
        const int nonnegative_metal_unsaturated_organic_cation_count =
            NonnegativeMetalUnsaturatedOrganicCationCount(mol, candidate->metal_states);
        FOR_BONDS_OF_MOL(bond_iter, mol)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
            if (begin_atom == nullptr || end_atom == nullptr)
            {
                continue;
            }
            const int begin_charge = static_cast<int>(begin_atom->GetFormalCharge());
            const int end_charge = static_cast<int>(end_atom->GetFormalCharge());
            const int begin_charge_sign = ChargeSign(begin_charge);
            const int end_charge_sign = ChargeSign(end_charge);
            if (begin_charge_sign == 0 || end_charge_sign == 0)
            {
                continue;
            }
            if (begin_charge_sign != end_charge_sign)
            {
                continue;
            }

            const int begin_idx = static_cast<int>(begin_atom->GetIdx());
            const int end_idx = static_cast<int>(end_atom->GetIdx());
            const bool pair_both_inner_visible =
                visible_inner_atom_indices.find(begin_idx) != visible_inner_atom_indices.end() &&
                visible_inner_atom_indices.find(end_idx) != visible_inner_atom_indices.end();
            const bool is_inner_visible_adjacent_carbanion_pair =
                pair_both_inner_visible &&
                begin_atom->GetAtomicNum() == 6 &&
                end_atom->GetAtomicNum() == 6 &&
                begin_charge_sign == end_charge_sign;
            if (is_inner_visible_adjacent_carbanion_pair)
            {
                ++inner_visible_adjacent_carbanion_pair_count;
            }
            else if (!pair_both_inner_visible)
            {
                ++outer_or_invisible_adjacent_double_charge_count;
                const int metal_charge_sign = NearestNonzeroMetalChargeSignToBond(
                    *begin_atom,
                    *end_atom,
                    candidate->metal_states);
                if (metal_charge_sign == 0)
                {
                    ++outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count;
                }
                else if (begin_charge_sign == metal_charge_sign)
                {
                    ++outer_or_invisible_adjacent_same_sign_double_charge_count;
                }
                else
                {
                    ++outer_or_invisible_adjacent_opposite_sign_double_charge_count;
                }
            }
        }

        const double negative_metal_penalty =
            kNegativeMetalDiscordancePenalty * static_cast<double>(negative_metal_discordance.count);
        const double discordance_count =
            static_cast<double>(
                inner_visible_diradical_count +
                outer_or_invisible_adjacent_double_charge_count +
                inner_visible_adjacent_carbanion_pair_count +
                inner_visible_conjugated_carbanion_pair_count +
                inner_visible_same_sign_charge_count +
                zero_valent_metals_with_organic_cation_count +
                nonnegative_metal_unsaturated_organic_cation_count) +
            negative_metal_penalty;
        candidate->metadata["metal_discordance_structural_count"] = discordance_count;
        candidate->metadata["metal_discordance_aromatic_ring_deficit_count"] = 0;
        candidate->metadata["metal_discordance_count"] = discordance_count;
        candidate->metadata["metal_discordance_inner_visible_diradical_count"] =
            inner_visible_diradical_count;
        candidate->metadata["metal_discordance_outer_or_invisible_adjacent_double_charge_count"] =
            outer_or_invisible_adjacent_double_charge_count;
        candidate->metadata["metal_discordance_outer_or_invisible_adjacent_same_sign_double_charge_count"] =
            outer_or_invisible_adjacent_same_sign_double_charge_count;
        candidate->metadata["metal_discordance_outer_or_invisible_adjacent_opposite_sign_double_charge_count"] =
            outer_or_invisible_adjacent_opposite_sign_double_charge_count;
        candidate->metadata["metal_discordance_outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count"] =
            outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count;
        candidate->metadata["metal_discordance_inner_visible_adjacent_carbanion_pair_count"] =
            inner_visible_adjacent_carbanion_pair_count;
        candidate->metadata["metal_discordance_inner_visible_conjugated_carbanion_pair_count"] =
            inner_visible_conjugated_carbanion_pair_count;
        candidate->metadata["metal_discordance_inner_visible_same_sign_charge_count"] =
            inner_visible_same_sign_charge_count;
        candidate->metadata["metal_discordance_negative_metal_count"] =
            negative_metal_discordance.count;
        candidate->metadata["metal_discordance_negative_metal_penalty"] =
            negative_metal_penalty;
        candidate->metadata["metal_discordance_negative_metal_outer_sphere_cation_exception"] =
            negative_metal_discordance.has_outer_sphere_cation_exception;
        candidate->metadata["metal_discordance_negative_metal_positive_metal_counterion_exception"] =
            negative_metal_discordance.has_positive_metal_counterion_exception;
        candidate->metadata["metal_discordance_zero_valent_metals_with_organic_cation_count"] =
            zero_valent_metals_with_organic_cation_count;
        candidate->metadata["metal_discordance_nonnegative_metal_unsaturated_organic_cation_count"] =
            nonnegative_metal_unsaturated_organic_cation_count;
        return discordance_count;
    }

    void AnnotateOrganicElectronicStateConsistency(
        molgr::state::MetalCandidateState *candidate,
        const molgr::config::MolGRConfig &config)
    {
        if (candidate == nullptr || !candidate->no_metal_state)
        {
            throw std::runtime_error(
                "MetalCandidateState requires no_metal_state before organic-state scoring");
        }

        const auto metrics = ComputeOrganicElectronicStateMetrics(
            candidate->no_metal_state->Mol(),
            config);
        candidate->metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count;
        candidate->metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count;
        candidate->metadata["organic_aromatic_stability_score"] =
            metrics.aromatic_stability_score;
        candidate->metadata["organic_conjugated_atom_count"] = metrics.conjugated_atom_count;
        candidate->metadata["organic_conjugated_bond_count"] = metrics.conjugated_bond_count;
        candidate->metadata["organic_max_conjugated_component_size"] =
            metrics.max_conjugated_component_size;
        candidate->metadata["organic_radical_localization_penalty"] =
            metrics.radical_localization_penalty;
        candidate->metadata["organic_charge_localization_penalty"] =
            metrics.charge_localization_penalty;
    }

    void AnnotateSelectedCandidateMetrics(
        molgr::state::MetalCandidateState *candidate,
        const molgr::config::MolGRConfig &config)
    {
        if (candidate == nullptr)
        {
            throw std::runtime_error("candidate is null");
        }
        if (!candidate->score.has_value())
        {
            const double score = candidate->CombinedScore(config);
            candidate->score = score;
            candidate->metadata["score"] = score;
        }
        if (candidate->metadata.find("organic_aromatic_atom_count") == candidate->metadata.end())
        {
            AnnotateOrganicElectronicStateConsistency(candidate, config);
        }
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

            std::optional<molgr::state::MetalCandidateState> SelectBestCandidate(
                const std::vector<molgr::state::MetalCandidateState> &candidates,
                const molgr::config::MolGRConfig &config)
            {
                if (candidates.empty())
                {
                    return std::nullopt;
                }

                std::vector<molgr::state::MetalCandidateState> scored_candidates = candidates;
                int max_aromatic_ring_count = 0;
                for (auto &candidate : scored_candidates)
                {
                    if (candidate.metadata.find("organic_aromatic_ring_count") ==
                            candidate.metadata.end() ||
                        candidate.metadata.find("organic_aromatic_stability_score") ==
                            candidate.metadata.end())
                    {
                        AnnotateOrganicElectronicStateConsistency(&candidate, config);
                    }
                    if (candidate.metadata.find("organic_aromatic_stability_score") ==
                        candidate.metadata.end())
                    {
                        candidate.metadata["organic_aromatic_stability_score"] =
                            MetadataDouble(
                                candidate,
                                "organic_aromatic_ring_count",
                                0.0);
                    }
                    max_aromatic_ring_count = std::max(
                        max_aromatic_ring_count,
                        MetadataInt(candidate, "organic_aromatic_ring_count", 0));
                }
                double max_aromatic_stability_score = 0.0;
                for (const auto &candidate : scored_candidates)
                {
                    max_aromatic_stability_score = std::max(
                        max_aromatic_stability_score,
                        MetadataDouble(candidate, "organic_aromatic_stability_score", 0.0));
                }

                for (auto &candidate : scored_candidates)
                {
                    const double structural_discordance_count = MetadataDouble(
                        candidate,
                        "metal_discordance_structural_count",
                        MetadataDouble(candidate, "metal_discordance_count", 0.0));
                    const int aromatic_ring_deficit_count = std::max(
                        0,
                        max_aromatic_ring_count -
                            MetadataInt(candidate, "organic_aromatic_ring_count", 0));
                    const double aromatic_stability_deficit = std::max(
                        0.0,
                        max_aromatic_stability_score -
                            MetadataDouble(candidate, "organic_aromatic_stability_score", 0.0));
                    candidate.metadata["metal_discordance_structural_count"] =
                        structural_discordance_count;
                    candidate.metadata["metal_discordance_max_aromatic_ring_count"] =
                        max_aromatic_ring_count;
                    candidate.metadata["metal_discordance_max_aromatic_stability_score"] =
                        max_aromatic_stability_score;
                    candidate.metadata["metal_discordance_aromatic_ring_deficit_count"] =
                        aromatic_ring_deficit_count;
                    candidate.metadata["metal_discordance_aromatic_stability_deficit"] =
                        aromatic_stability_deficit;
                    candidate.metadata["metal_discordance_count"] =
                        structural_discordance_count + aromatic_stability_deficit;
                }

                double min_discordance_count = std::numeric_limits<double>::infinity();
                for (const auto &candidate : scored_candidates)
                {
                    min_discordance_count = std::min(
                        min_discordance_count,
                        MetadataDouble(candidate, "metal_discordance_count", 0.0));
                }

                std::vector<molgr::state::MetalCandidateState> discordance_filtered_candidates;
                discordance_filtered_candidates.reserve(scored_candidates.size());
                for (auto candidate : scored_candidates)
                {
                    const bool passes_discordance_filter =
                        MetadataDouble(candidate, "metal_discordance_count", 0.0) ==
                        min_discordance_count;
                    candidate.metadata["passes_metal_discordance_filter"] = passes_discordance_filter;
                    if (passes_discordance_filter)
                    {
                        discordance_filtered_candidates.push_back(std::move(candidate));
                    }
                }
                if (discordance_filtered_candidates.empty())
                {
                    return std::nullopt;
                }

                for (auto &candidate : discordance_filtered_candidates)
                {
                    AnnotateSelectedCandidateMetrics(&candidate, config);
                }

                std::optional<std::tuple<double, int>> best_selection_key;
                std::optional<molgr::state::MetalCandidateState> best_candidate;
                for (auto &candidate : discordance_filtered_candidates)
                {
                    const double score_value =
                        candidate.score.has_value() ? *candidate.score : candidate.CombinedScore(config);
                    const auto selection_key =
                        std::make_tuple(score_value, CandidateCombinationIndex(candidate));
                    candidate.metadata["selection_key"] =
                        std::to_string(min_discordance_count) + "," +
                        std::to_string(score_value) + "," +
                        std::to_string(CandidateCombinationIndex(candidate));
                    if (best_selection_key.has_value() && selection_key >= *best_selection_key)
                    {
                        continue;
                    }
                    best_selection_key = selection_key;
                    best_candidate = candidate;
                }
                return best_candidate;
            }

            molgr::state::MetalCandidateState PrepareCandidateWithNoMetalState(
                const molgr::state::MetalCandidateState &candidate,
                const std::shared_ptr<molgr::state::ReconstructionState> &no_metal_state,
                const molgr::config::MolGRConfig &config)
            {
                auto machine = molgr::state::MetalCandidateStateMachine::FromCandidateState(candidate);
                machine.SetNoMetalState("reconstruct_no_metal", no_metal_state);
                machine.Annotate("score_candidate");
                auto prepared_candidate = machine.Freeze();
                const double score = prepared_candidate.CombinedScore(config);
                prepared_candidate.score = score;
                prepared_candidate.metadata["score"] = score;
                AnnotateCandidateDiscordanceFeatures(&prepared_candidate, config);
                return prepared_candidate;
            }

        }
    }
}
