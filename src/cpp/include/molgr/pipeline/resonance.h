#pragma once

#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"

#include <openbabel/mol.h>

#include <array>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        using ResonanceAtomKey = std::tuple<int, int, int, bool>;
        using ResonanceBondKey = std::tuple<int, int, int, bool>;
        using ResonanceBondIndexMap = std::map<std::pair<int, int>, std::size_t>;
        using DirectGainMetrics = std::array<double, 4>;

        struct ResonanceStateKey
        {
            std::vector<ResonanceAtomKey> atom_keys;
            std::vector<ResonanceBondKey> bond_keys;

            bool operator==(const ResonanceStateKey &other) const;
            bool operator<(const ResonanceStateKey &other) const;
        };

        using ProcessedResonanceKey = ResonanceStateKey;

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
        };

        using ResonanceVisitCallback = std::function<bool(const ResonanceSearchNode &)>;
        using DirectGainBoundCacheKey = std::pair<ResonanceStateKey, int>;
        using DirectGainBoundCache = std::map<DirectGainBoundCacheKey, DirectGainMetrics>;

        struct LimitedDiscrepancyTraversalConfig
        {
            int max_discrepancy = 1;
            bool fallback_to_full_frontier = true;
        };

        ResonanceStateKey BuildResonanceStateKey(const OpenBabel::OBMol &mol);
        ProcessedResonanceKey BuildProcessedResonanceKey(const OpenBabel::OBMol &mol);

        std::vector<OpenBabel::OBMol> GetOneStepResonance(const OpenBabel::OBMol &mol);
        std::vector<OpenBabel::OBMol> GetRadicalResonances(
            const OpenBabel::OBMol &mol,
            int max_depth = 2);
        std::vector<OpenBabel::OBMol> GetRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            const LimitedDiscrepancyTraversalConfig &config = {});

        void WalkRadicalResonances(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            ResonanceVisitCallback visit = {});
        void WalkRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            ResonanceVisitCallback visit = {},
            const LimitedDiscrepancyTraversalConfig &config = {});

        double EstimateRemainingResonanceScoreImprovementUpperBound(
            const OpenBabel::OBMol &mol,
            const ResonanceStateKey &state_key,
            int remaining_steps,
            DirectGainBoundCache *cache = nullptr);

        std::tuple<OpenBabel::OBMol, int, bool> ProcessResonanceDetailed(
            const OpenBabel::OBMol &mol,
            int charge);
        std::pair<OpenBabel::OBMol, int> ProcessResonance(
            const OpenBabel::OBMol &mol,
            int charge);
        std::string SmilesFirstToken(const OpenBabel::OBMol &mol);
    }
}
