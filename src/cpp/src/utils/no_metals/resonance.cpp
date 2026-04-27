#include "molgr/utils/no_metals/resonance.h"

#include "molgr/pipeline/resonance.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/no_metals/selection.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace
{
    constexpr bool kDefaultResonanceFallbackToFullFrontier = true;

    struct RawResonanceCandidate
    {
        std::size_t resonance_index = 0;
        std::shared_ptr<OpenBabel::OBMol> omol;
    };

    struct PreparedResonanceCandidate
    {
        molgr::resonance::ProcessedResonanceKey processed_state_key;
        std::optional<molgr::state::ReconstructionState> candidate;
    };

    PreparedResonanceCandidate PrepareResonanceCandidate(
        RawResonanceCandidate raw_candidate,
        const molgr::state::ReconstructionState &state,
        const molgr::state::OmolStateMachine &base_machine)
    {
        auto branched_machine = base_machine.Branch(
            "branch_resonance_candidate",
            std::move(raw_candidate.omol));

        branched_machine.RunOmolChargeStage(
            std::nullopt,
            molgr::reconstruct::Eliminate13Dipole);
        branched_machine.RunOmolChargeStage(
            std::nullopt,
            molgr::reconstruct::EliminatePositiveCharges);
        branched_machine.RunOmolChargeStage(
            std::nullopt,
            molgr::reconstruct::EliminateNegativeCharges);
        branched_machine.RunOmolStage(
            std::nullopt,
            molgr::reconstruct::CleanNeighborRadicals);
        branched_machine.RunOmolStage(
            std::nullopt,
            molgr::reconstruct::CleanResonances);
        branched_machine.Annotate("process_resonance");

        PreparedResonanceCandidate prepared;
        prepared.processed_state_key =
            molgr::resonance::BuildProcessedResonanceKey(branched_machine.EnsureUniqueMol());

        if (molgr::reconstruct::ValidateOmol(
                branched_machine.EnsureUniqueMol(),
                state.total_charge,
                state.total_radical_electrons))
        {
            branched_machine.Annotate("validate_resonance_candidate");
            branched_machine.metadata["resonance_index"] =
                static_cast<int>(raw_candidate.resonance_index);
            prepared.candidate = branched_machine.FreezeLike(state);
        }

        return prepared;
    }

    void AnnotatePreparedCandidateTopology(
        molgr::state::ReconstructionState &candidate)
    {
        molgr::no_metals::selection::AnnotateNoMetalCandidateTopology(candidate);
    }
}

namespace molgr
{
    namespace no_metals
    {
        namespace resonance
        {
            std::vector<molgr::state::ReconstructionState> RecoverResonanceCandidates(
                const molgr::state::ReconstructionState &state,
                const molgr::config::MolGRConfig &config)
            {
                std::vector<RawResonanceCandidate> raw_candidates;
                auto base_machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                std::size_t resonance_index = 0;
                const int resonance_max_depth = std::max(0, config.resonance.max_depth);
                const int max_discrepancy = std::max(
                    0,
                    config.resonance.limited_discrepancy_max_discrepancy);

                reconstruct::WalkRadicalResonancesLimitedDiscrepancy(
                    state.Mol(),
                    resonance_max_depth,
                    [&](const reconstruct::ResonanceSearchNode &node) -> bool
                    {
                        raw_candidates.push_back(
                            RawResonanceCandidate{
                                resonance_index++,
                                std::make_shared<OpenBabel::OBMol>(node.omol),
                            });
                        return true;
                    },
                    reconstruct::LimitedDiscrepancyTraversalConfig{
                        max_discrepancy,
                        kDefaultResonanceFallbackToFullFrontier,
                    },
                    config);

                std::vector<std::optional<PreparedResonanceCandidate>> prepared_candidates(
                    raw_candidates.size());
                for (std::size_t candidate_index = 0;
                     candidate_index < raw_candidates.size();
                     ++candidate_index)
                {
                    prepared_candidates[candidate_index] = PrepareResonanceCandidate(
                        std::move(raw_candidates[candidate_index]),
                        state,
                        base_machine);
                }

                std::vector<molgr::state::ReconstructionState> candidates;
                candidates.reserve(prepared_candidates.size());
                std::set<reconstruct::ProcessedResonanceKey> seen_processed_states;
                for (auto &prepared_candidate : prepared_candidates)
                {
                    if (!prepared_candidate.has_value())
                    {
                        continue;
                    }
                    if (seen_processed_states.find(prepared_candidate->processed_state_key) !=
                        seen_processed_states.end())
                    {
                        continue;
                    }
                    seen_processed_states.insert(prepared_candidate->processed_state_key);
                    if (prepared_candidate->candidate.has_value())
                    {
                        candidates.push_back(std::move(*prepared_candidate->candidate));
                    }
                }

                for (auto &candidate : candidates)
                {
                    AnnotatePreparedCandidateTopology(candidate);
                }

                return candidates;
            }
        }
    }
}
