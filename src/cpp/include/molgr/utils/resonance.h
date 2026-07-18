#pragma once

#include "molgr/config.h"

#include <openbabel/mol.h>

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace molgr
{
    namespace resonance
    {
        using ResonanceAtomKey = std::tuple<int, int, int, bool>;
        using ResonanceBondKey = std::tuple<int, int, int, bool>;
        using ResonanceBondIndexMap = std::map<std::pair<int, int>, std::size_t>;
        using UffLiteGainMetrics = std::array<double, 4>;

        struct ResonanceStateKey
        {
            std::vector<ResonanceAtomKey> atom_keys;
            std::vector<ResonanceBondKey> bond_keys;

            bool operator==(const ResonanceStateKey &other) const;
            bool operator<(const ResonanceStateKey &other) const;
        };

        using ProcessedResonanceKey = std::string;

        struct ResonanceTraversalMove
        {
            std::tuple<int, int, int> path;
            ResonanceStateKey next_state_key;
        };

        struct ResonanceSearchNode
        {
            const OpenBabel::OBMol &omol;
            const ResonanceStateKey &state_key;
            int depth;
            int discrepancy;
        };

        struct IndexedResonanceTraversalMove
        {
            std::tuple<int, int, int> idxs;
            ResonanceStateKey next_state_key;
        };

        using UffLiteGainBoundCacheKey = std::pair<ResonanceStateKey, int>;
        using UffLiteGainBoundCache = std::map<UffLiteGainBoundCacheKey, UffLiteGainMetrics>;

        struct LimitedDiscrepancyTraversalConfig
        {
            int max_discrepancy = 1;
        };

        std::string SmilesFirstToken(const OpenBabel::OBMol &mol);

        ResonanceStateKey BuildResonanceStateKey(const OpenBabel::OBMol &mol);
        ProcessedResonanceKey BuildProcessedResonanceKey(const OpenBabel::OBMol &mol);
        std::pair<ResonanceStateKey, ResonanceBondIndexMap> BuildResonanceSearchContext(
            const OpenBabel::OBMol &mol);
        ResonanceBondIndexMap BuildBondIndexMapFromStateKey(
            const ResonanceStateKey &state_key);
        ResonanceStateKey IncrementResonanceStateKey(
            const ResonanceStateKey &state_key,
            const ResonanceBondIndexMap &bond_index_map,
            const std::tuple<int, int, int> &idxs);

        std::vector<IndexedResonanceTraversalMove> EnumerateOneStepResonanceMoves(
            const OpenBabel::OBMol &mol,
            const ResonanceStateKey &state_key,
            const ResonanceBondIndexMap &bond_index_map);
        OpenBabel::OBMol MaterializeOneStepResonance(
            const OpenBabel::OBMol &mol,
            const std::tuple<int, int, int> &idxs);
        std::shared_ptr<OpenBabel::OBMol> MaterializeOneStepResonancePtr(
            const OpenBabel::OBMol &mol,
            const std::tuple<int, int, int> &idxs);
        std::vector<OpenBabel::OBMol> GetOneStepResonance(const OpenBabel::OBMol &mol);

        std::vector<IndexedResonanceTraversalMove> SelectLimitedDiscrepancyMoves(
            const OpenBabel::OBMol &mol,
            const std::vector<IndexedResonanceTraversalMove> &moves,
            const LimitedDiscrepancyTraversalConfig &traversal_config,
            const molgr::config::MolGRConfig &config);

        double EstimateRemainingResonanceScoreImprovementUpperBound(
            const OpenBabel::OBMol &mol,
            const ResonanceStateKey &state_key,
            int remaining_steps,
            UffLiteGainBoundCache *cache = nullptr);

        std::tuple<OpenBabel::OBMol, int, bool> ProcessResonanceDetailed(
            const OpenBabel::OBMol &mol,
            int charge);
        std::pair<OpenBabel::OBMol, int> ProcessResonance(
            const OpenBabel::OBMol &mol,
            int charge);
    }
}
