#include "molgr/utils/metal_radical_inference.h"

#include "molgr/utils/consts.h"

#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
    using Vector3 = std::tuple<double, double, double>;

    struct ShellOccupation
    {
        int remaining_f = 0;
        int effective_d = 0;
        int residual_sp = 0;
    };

    struct DonorSample
    {
        int atom_idx = 0;
        int atomic_num = 0;
        double distance_angstrom = 0.0;
        Vector3 vector;
    };

    const std::set<std::string> kDBlockRearrangementMetals = {
        "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y",  "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
        "Hg"};

    const std::map<int, double> kDonorFieldStrength = {
        {1, 1.20},
        {6, 1.30},  {7, 1.00},  {8, 0.85},  {9, 0.65}, {15, 1.15},
        {16, 0.75}, {17, 0.55}, {35, 0.50}, {53, 0.45},
    };

    const std::map<std::string, double> kGeometryFieldAdjustment = {
        {"free_ion", 0.0},
        {"terminal", 0.0},
        {"bent", 0.0},
        {"linear", 0.20},
        {"trigonal_planar", 0.15},
        {"trigonal_pyramidal", -0.05},
        {"tetrahedral", -0.25},
        {"square_planar", 0.35},
        {"octahedral_like", 0.05},
    };

    constexpr double kMetalHydrideCovalentToleranceAngstrom = 0.45;
    constexpr double kMinMetalHydrideCutoffAngstrom = 1.25;

    Vector3 VectorFromAtoms(OpenBabel::OBAtom &metal_atom, OpenBabel::OBAtom &donor_atom)
    {
        return {
            donor_atom.GetX() - metal_atom.GetX(),
            donor_atom.GetY() - metal_atom.GetY(),
            donor_atom.GetZ() - metal_atom.GetZ(),
        };
    }

    double VectorNorm(const Vector3 &vector)
    {
        return std::sqrt(
            std::get<0>(vector) * std::get<0>(vector) +
            std::get<1>(vector) * std::get<1>(vector) +
            std::get<2>(vector) * std::get<2>(vector));
    }

    double ClampedCovalentRadius(int atomic_num)
    {
        const double value = OpenBabel::OBElements::GetCovalentRad(atomic_num);
        if (value <= 0.0 || !std::isfinite(value))
        {
            return 0.0;
        }
        return value;
    }

    double MetalHydrideCutoffAngstrom(int metal_atomic_num)
    {
        return std::max(
            kMinMetalHydrideCutoffAngstrom,
            ClampedCovalentRadius(metal_atomic_num) + ClampedCovalentRadius(1) +
                kMetalHydrideCovalentToleranceAngstrom);
    }

    bool IsDirectMetalHydride(
        OpenBabel::OBAtom &metal_atom,
        OpenBabel::OBAtom &donor_atom,
        double distance)
    {
        return static_cast<int>(donor_atom.GetAtomicNum()) == 1 &&
               distance <= MetalHydrideCutoffAngstrom(static_cast<int>(metal_atom.GetAtomicNum()));
    }

    Vector3 Cross(const Vector3 &lhs, const Vector3 &rhs)
    {
        return {
            std::get<1>(lhs) * std::get<2>(rhs) - std::get<2>(lhs) * std::get<1>(rhs),
            std::get<2>(lhs) * std::get<0>(rhs) - std::get<0>(lhs) * std::get<2>(rhs),
            std::get<0>(lhs) * std::get<1>(rhs) - std::get<1>(lhs) * std::get<0>(rhs),
        };
    }

    double Dot(const Vector3 &lhs, const Vector3 &rhs)
    {
        return std::get<0>(lhs) * std::get<0>(rhs) +
               std::get<1>(lhs) * std::get<1>(rhs) +
               std::get<2>(lhs) * std::get<2>(rhs);
    }

    double AngleDegrees(const Vector3 &lhs, const Vector3 &rhs)
    {
        constexpr double kPi = 3.14159265358979323846;
        const double lhs_norm = VectorNorm(lhs);
        const double rhs_norm = VectorNorm(rhs);
        if (lhs_norm <= 1e-8 || rhs_norm <= 1e-8)
        {
            return 0.0;
        }
        double cos_value = Dot(lhs, rhs) / (lhs_norm * rhs_norm);
        cos_value = std::max(-1.0, std::min(1.0, cos_value));
        return std::acos(cos_value) * 180.0 / kPi;
    }

    double PlanarityDistance(const std::vector<Vector3> &vectors)
    {
        if (vectors.size() < 3)
        {
            return std::numeric_limits<double>::infinity();
        }

        std::optional<Vector3> best_normal;
        double best_norm = 0.0;
        for (std::size_t i = 0; i < vectors.size(); ++i)
        {
            for (std::size_t j = i + 1; j < vectors.size(); ++j)
            {
                const Vector3 normal = Cross(vectors[i], vectors[j]);
                const double normal_norm = VectorNorm(normal);
                if (normal_norm > best_norm)
                {
                    best_normal = normal;
                    best_norm = normal_norm;
                }
            }
        }
        if (!best_normal.has_value() || best_norm <= 1e-8)
        {
            return std::numeric_limits<double>::infinity();
        }

        const Vector3 normal = {
            std::get<0>(*best_normal) / best_norm,
            std::get<1>(*best_normal) / best_norm,
            std::get<2>(*best_normal) / best_norm,
        };
        double total = 0.0;
        for (const auto &vector : vectors)
        {
            total += std::abs(Dot(vector, normal));
        }
        return total / static_cast<double>(vectors.size());
    }

    std::vector<DonorSample> CollectCoordinationEnvironment(
        OpenBabel::OBAtom &metal_atom,
        const molgr::config::MetalRadicalInferenceConfig &config)
    {
        OpenBabel::OBMol *parent = metal_atom.GetParent();
        if (parent == nullptr)
        {
            return {};
        }

        std::vector<DonorSample> donors;
        FOR_ATOMS_OF_MOL(neighbor_iter, *parent)
        {
            OpenBabel::OBAtom &neighbor = *neighbor_iter;
            if (neighbor.GetIdx() == metal_atom.GetIdx() || neighbor.IsMetal())
            {
                continue;
            }
            const int atomic_num = static_cast<int>(neighbor.GetAtomicNum());

            const Vector3 vector = VectorFromAtoms(metal_atom, neighbor);
            const double distance = VectorNorm(vector);
            if (atomic_num <= 1)
            {
                if (!IsDirectMetalHydride(metal_atom, neighbor, distance))
                {
                    continue;
                }
            }
            else if (distance > config.coordination_cutoff_angstrom)
            {
                continue;
            }

            donors.push_back(DonorSample{
                static_cast<int>(neighbor.GetIdx()),
                atomic_num,
                distance,
                vector,
            });
        }

        std::sort(
            donors.begin(),
            donors.end(),
            [](const DonorSample &lhs, const DonorSample &rhs)
            {
                return lhs.distance_angstrom < rhs.distance_angstrom;
            });
        if (static_cast<int>(donors.size()) > config.max_considered_donors)
        {
            donors.resize(static_cast<std::size_t>(config.max_considered_donors));
        }
        return donors;
    }

    std::string ClassifyGeometry(
        const std::vector<DonorSample> &donors,
        const molgr::config::MetalRadicalInferenceConfig &config)
    {
        const int coordination_number = static_cast<int>(donors.size());
        std::vector<Vector3> vectors;
        vectors.reserve(donors.size());
        for (const auto &donor : donors)
        {
            vectors.push_back(donor.vector);
        }

        if (coordination_number == 0)
        {
            return "free_ion";
        }
        if (coordination_number == 1)
        {
            return "terminal";
        }
        if (coordination_number == 2)
        {
            const double angle = AngleDegrees(vectors[0], vectors[1]);
            return angle >= config.linear_angle_min_degrees ? "linear" : "bent";
        }
        if (coordination_number == 3)
        {
            return PlanarityDistance(vectors) <= config.trigonal_planar_planarity_tolerance_angstrom
                       ? "trigonal_planar"
                       : "trigonal_pyramidal";
        }
        if (coordination_number == 4)
        {
            return PlanarityDistance(vectors) <= config.square_planar_planarity_tolerance_angstrom
                       ? "square_planar"
                       : "tetrahedral";
        }
        return "octahedral_like";
    }

    double DonorFieldScore(
        const std::vector<DonorSample> &donors,
        const std::string &geometry)
    {
        if (donors.empty())
        {
            const auto geometry_it = kGeometryFieldAdjustment.find(geometry);
            return geometry_it == kGeometryFieldAdjustment.end() ? 0.0 : geometry_it->second;
        }

        double weighted_sum = 0.0;
        double total_weight = 0.0;
        for (const auto &donor : donors)
        {
            const auto base_strength_it = kDonorFieldStrength.find(donor.atomic_num);
            const double base_strength =
                base_strength_it == kDonorFieldStrength.end() ? 0.60 : base_strength_it->second;
            const double distance_weight =
                1.0 / std::max(donor.distance_angstrom * donor.distance_angstrom, 1.0);
            weighted_sum += base_strength * distance_weight;
            total_weight += distance_weight;
        }

        const auto geometry_it = kGeometryFieldAdjustment.find(geometry);
        const double geometry_adjustment =
            geometry_it == kGeometryFieldAdjustment.end() ? 0.0 : geometry_it->second;
        if (total_weight <= 0.0)
        {
            return geometry_adjustment;
        }
        return weighted_sum / total_weight + geometry_adjustment;
    }

    std::string ClassifyFieldStrength(
        double field_score,
        const molgr::config::MetalRadicalInferenceConfig &config)
    {
        if (field_score >= config.strong_field_threshold)
        {
            return "strong";
        }
        if (field_score <= config.weak_field_threshold)
        {
            return "weak";
        }
        return "intermediate";
    }

    std::optional<ShellOccupation> ShellOccupationAfterOxidation(
        const std::string &metal,
        int valence)
    {
        const auto fdsp_it = molgr::kMetalFDSP.find(metal);
        if (fdsp_it == molgr::kMetalFDSP.end())
        {
            return std::nullopt;
        }

        const molgr::FDSP &fdsp = fdsp_it->second;
        const int total_outer = fdsp.f + fdsp.d + fdsp.s + fdsp.p;
        if (valence > total_outer)
        {
            return std::nullopt;
        }

        const int remaining_sp = std::max(0, fdsp.s + fdsp.p - valence);
        const int removed_from_d = std::min(fdsp.d, std::max(0, valence - (fdsp.s + fdsp.p)));
        const int remaining_d = fdsp.d - removed_from_d;
        const int removed_from_f =
            std::min(fdsp.f, std::max(0, valence - (fdsp.s + fdsp.p + fdsp.d)));
        const int remaining_f = fdsp.f - removed_from_f;

        int promoted_to_d = 0;
        if (kDBlockRearrangementMetals.find(metal) != kDBlockRearrangementMetals.end())
        {
            promoted_to_d = std::min(remaining_sp, std::max(0, 10 - remaining_d));
        }

        return ShellOccupation{
            remaining_f,
            remaining_d + promoted_to_d,
            remaining_sp - promoted_to_d,
        };
    }

    std::vector<int> CandidateDUnpairedCounts(
        int effective_d,
        const std::string &geometry,
        const std::string &field_strength)
    {
        if (effective_d < 0 || effective_d >= static_cast<int>(molgr::kDElectronsSpin.size()))
        {
            return {0};
        }

        std::vector<int> free_ion_candidates = molgr::kDElectronsSpin[effective_d];
        std::sort(free_ion_candidates.begin(), free_ion_candidates.end());
        free_ion_candidates.erase(
            std::unique(free_ion_candidates.begin(), free_ion_candidates.end()),
            free_ion_candidates.end());

        if (geometry == "square_planar")
        {
            if (effective_d == 8)
            {
                return {0};
            }
            if (effective_d == 7 || effective_d == 9)
            {
                return {1};
            }
            return {free_ion_candidates.front(), free_ion_candidates.back()};
        }

        if (geometry == "tetrahedral")
        {
            return {free_ion_candidates.back()};
        }
        if (free_ion_candidates.size() == 1)
        {
            return {free_ion_candidates.front()};
        }
        return {free_ion_candidates.front(), free_ion_candidates.back()};
    }
}

namespace molgr
{
    namespace metal_radical_inference
    {
        MetalRadicalInferenceResult InferMetalRadicalState(
            OpenBabel::OBAtom &metal_atom,
            int valence,
            const molgr::config::MolGRConfig &config)
        {
            const auto &metal_radical_config = config.metal_radical_inference;
            const std::string symbol =
                OpenBabel::OBElements::GetSymbol(static_cast<int>(metal_atom.GetAtomicNum()));

            const auto occupation = ShellOccupationAfterOxidation(symbol, valence);
            if (!occupation.has_value())
            {
                return {};
            }

            const auto donors = CollectCoordinationEnvironment(metal_atom, metal_radical_config);
            const std::string geometry = ClassifyGeometry(donors, metal_radical_config);
            const double field_score = DonorFieldScore(donors, geometry);
            const std::string field_strength = ClassifyFieldStrength(field_score, metal_radical_config);
            const int base_unpaired = (occupation->remaining_f + occupation->residual_sp) % 2;

            const auto d_candidates =
                CandidateDUnpairedCounts(occupation->effective_d, geometry, field_strength);
            std::set<int> radical_counts_set;
            for (const int candidate : d_candidates)
            {
                radical_counts_set.insert(base_unpaired + candidate);
            }

            MetalRadicalInferenceResult result;
            result.radical_counts.assign(radical_counts_set.begin(), radical_counts_set.end());
            result.effective_d_electrons = occupation->effective_d;
            result.residual_sp_electrons = occupation->residual_sp;
            result.remaining_f_electrons = occupation->remaining_f;
            result.coordination_number = static_cast<int>(donors.size());
            result.geometry = geometry;
            result.field_score = field_score;
            result.field_strength = field_strength;
            return result;
        }

        std::vector<int> InferMetalRadicalCounts(
            OpenBabel::OBAtom &metal_atom,
            int valence,
            const molgr::config::MolGRConfig &config)
        {
            return InferMetalRadicalState(metal_atom, valence, config).radical_counts;
        }
    }
}
