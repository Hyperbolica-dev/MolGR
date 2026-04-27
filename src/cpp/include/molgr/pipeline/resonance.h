#pragma once

#include "molgr/config.h"
#include "molgr/utils/resonance.h"

#include <openbabel/mol.h>

#include <functional>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        using molgr::resonance::BuildBondIndexMapFromStateKey;
        using molgr::resonance::BuildProcessedResonanceKey;
        using molgr::resonance::BuildResonanceSearchContext;
        using molgr::resonance::BuildResonanceStateKey;
        using molgr::resonance::DirectGainBoundCache;
        using molgr::resonance::DirectGainBoundCacheKey;
        using molgr::resonance::DirectGainMetrics;
        using molgr::resonance::EnumerateOneStepResonanceMoves;
        using molgr::resonance::EstimateRemainingResonanceScoreImprovementUpperBound;
        using molgr::resonance::GetOneStepResonance;
        using molgr::resonance::IncrementResonanceStateKey;
        using molgr::resonance::IndexedResonanceTraversalMove;
        using molgr::resonance::LimitedDiscrepancyTraversalConfig;
        using molgr::resonance::MaterializeOneStepResonance;
        using molgr::resonance::MaterializeOneStepResonancePtr;
        using molgr::resonance::OrderForceFieldMoves;
        using molgr::resonance::ProcessedResonanceKey;
        using molgr::resonance::ProcessResonance;
        using molgr::resonance::ProcessResonanceDetailed;
        using molgr::resonance::ResonanceAtomKey;
        using molgr::resonance::ResonanceBondIndexMap;
        using molgr::resonance::ResonanceBondKey;
        using molgr::resonance::ResonanceSearchNode;
        using molgr::resonance::ResonanceStateKey;
        using molgr::resonance::ResonanceTraversalMove;
        using molgr::resonance::ScoreOneStepResonanceWithForceField;
        using molgr::resonance::SelectLimitedDiscrepancyMoves;
        using molgr::resonance::SmilesFirstToken;

        using ResonanceVisitCallback = std::function<bool(const ResonanceSearchNode &)>;

        std::vector<OpenBabel::OBMol> GetRadicalResonances(
            const OpenBabel::OBMol &mol,
            int max_depth = 2);
        std::vector<OpenBabel::OBMol> GetRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            const LimitedDiscrepancyTraversalConfig &traversal_config = {},
            const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

        void WalkRadicalResonances(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            ResonanceVisitCallback visit = {});
        void WalkRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth = 2,
            ResonanceVisitCallback visit = {},
            const LimitedDiscrepancyTraversalConfig &traversal_config = {},
            const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
    }
}
