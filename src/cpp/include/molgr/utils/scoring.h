#pragma once

#include "molgr/types.h"

#include <openbabel/mol.h>

#include <string>
#include <utility>

namespace molgr
{
    namespace scoring
    {
        using ChargedAtomSnapshot = molgr::ChargedAtomSnapshot;
        using ChargedAtomSnapshotList = molgr::ChargedAtomSnapshotList;
        using MetalAtomPosition = molgr::MetalAtomPosition;

        std::string BuildScoreKey(const OpenBabel::OBMol &mol);
        std::string BuildMetalStateKey(const std::vector<MetalAtomPosition> &metal_states);

        double OmolScore(const OpenBabel::OBMol &mol);
        double OrganicCoreScore(const OpenBabel::OBMol &mol);
        double PostReinsertionScore(const OpenBabel::OBMol &mol);
        std::pair<double, ChargedAtomSnapshotList> BuildPostReinsertionBaseComponents(
            const OpenBabel::OBMol &mol);
        double PostReinsertionScoreFromMetalStates(
            double base_symmetry_penalty,
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states);
        double CombinedCandidateScoreFromMetalStates(
            double organic_score,
            const std::string &post_reinsertion_base_key,
            double base_symmetry_penalty,
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states);

        double CalcSymmetryPenalty(const OpenBabel::OBMol &mol);
        double CalculatePhysChemPenalty(const OpenBabel::OBMol &mol);
        double CalculateMetalPenalty(const OpenBabel::OBMol &mol);
        double CalculateMetalPenaltyFromMetalStates(
            const ChargedAtomSnapshotList &charged_atom_snapshots,
            const std::vector<MetalAtomPosition> &metal_states,
            double cutoff = 2.6);
        double GetDeviationScore(const OpenBabel::OBMol &mol, const OpenBabel::OBAtom *atom);
    }
}
