#pragma once

#include "molgr/config.h"
#include "molgr/types.h"

#include <openbabel/mol.h>

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace scoring
    {
        struct ForceFieldEvaluation
        {
            double raw_energy = 0.0;
            std::string raw_unit;
            double energy_kj_mol = 0.0;
            int atom_count = 0;
            int heavy_atom_count = 0;
            bool contains_metals = false;
        };

        std::string BuildScoreKey(const OpenBabel::OBMol &mol);
        std::string BuildMetalStateKey(const std::vector<molgr::MetalAtomPosition> &metal_states);

        bool ContainsMetalAtoms(const OpenBabel::OBMol &mol);
        OpenBabel::OBMol StripMetalAtoms(const OpenBabel::OBMol &mol);

        ForceFieldEvaluation EvaluateForceField(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
        std::tuple<std::size_t, std::size_t, std::size_t> ForceFieldEvaluationCacheInfo();
        void ForceFieldEvaluationCacheClear();
        ForceFieldEvaluation OrganicForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
        ForceFieldEvaluation OrganicForceFieldEvaluationWithScoreKey(
            const OpenBabel::OBMol &mol,
            const std::string &score_key,
            const molgr::config::MolGRConfig &config);
        ForceFieldEvaluation CombinedForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
        ForceFieldEvaluation SelectionForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);

        double OrganicForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
        double CombinedForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
        double SelectionForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config);
    }
}
