#include "molgr/utils/no_metals/resonance.h"

#include "molgr/pipeline/resonance.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/fresh.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/no_metals/selection.h"
#include "molgr/utils/resonance.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
    using RawStateKey = molgr::no_metals::resonance::RawStateKey;
    using ProcessedStateKey = molgr::no_metals::resonance::ProcessedStateKey;

    struct RawResonanceCandidate
    {
        std::size_t seed_index = 0;
        std::size_t resonance_index = 0;
        std::size_t raw_index = 0;
        std::shared_ptr<OpenBabel::OBMol> omol;
    };

    std::vector<molgr::state::ReconstructionState> DeduplicateStates(
        std::vector<molgr::state::ReconstructionState> states)
    {
        std::set<RawStateKey> seen;
        std::vector<molgr::state::ReconstructionState> unique;
        for (auto &state : states)
        {
            RawStateKey key{
                molgr::resonance::BuildResonanceStateKey(state.Mol()),
                state.given_charge,
            };
            if (seen.insert(std::move(key)).second)
            {
                unique.push_back(std::move(state));
            }
        }
        return unique;
    }

    template <typename Stage>
    void AppendOmolStageVariants(
        std::vector<molgr::state::ReconstructionState> &pool,
        const std::string &phase,
        Stage stage)
    {
        std::vector<molgr::state::ReconstructionState> additions;
        for (const auto &state : pool)
        {
            auto machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
            auto next_mol = std::make_shared<OpenBabel::OBMol>(
                molgr::utils::CloneMolTopologyOnly(state.Mol()));
            machine = machine.Branch(std::nullopt, std::move(next_mol));
            if (machine.RunOmolStage(phase, stage))
            {
                additions.push_back(machine.FreezeLike(state));
            }
        }
        for (auto &state : additions)
        {
            pool.push_back(std::move(state));
        }
        pool = DeduplicateStates(std::move(pool));
    }

    template <typename Stage>
    void AppendChargeStageVariants(
        std::vector<molgr::state::ReconstructionState> &pool,
        const std::string &phase,
        Stage stage)
    {
        std::vector<molgr::state::ReconstructionState> additions;
        for (const auto &state : pool)
        {
            auto machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
            auto next_mol = std::make_shared<OpenBabel::OBMol>(
                molgr::utils::CloneMolTopologyOnly(state.Mol()));
            machine = machine.Branch(std::nullopt, std::move(next_mol));
            if (machine.RunOmolChargeStage(phase, stage))
            {
                additions.push_back(machine.FreezeLike(state));
            }
        }
        for (auto &state : additions)
        {
            pool.push_back(std::move(state));
        }
        pool = DeduplicateStates(std::move(pool));
    }

    bool RegisterTraversalLabel(
        std::map<RawStateKey, std::vector<std::pair<int, int>>> &labels_by_state,
        const RawStateKey &state_key,
        const std::pair<int, int> &label)
    {
        auto &labels = labels_by_state[state_key];
        for (const auto &known : labels)
        {
            if (known.first <= label.first && known.second <= label.second)
            {
                return false;
            }
        }
        labels.erase(
            std::remove_if(
                labels.begin(),
                labels.end(),
                [&label](const auto &known)
                {
                    return label.first <= known.first && label.second <= known.second;
                }),
            labels.end());
        labels.push_back(label);
        return true;
    }

    std::optional<molgr::state::ReconstructionState> PrepareCandidateVariant(
        const RawResonanceCandidate &raw,
        const molgr::state::ReconstructionState &seed,
        const std::string &normalization,
        std::set<ProcessedStateKey> &seen_processed_states,
        const molgr::config::MolGRConfig &config)
    {
        auto machine = molgr::state::OmolStateMachine::FromReconstructionState(seed);
        auto candidate_mol = std::make_shared<OpenBabel::OBMol>(
            molgr::utils::CloneMolTopologyOnly(*raw.omol));
        machine = machine.Branch("branch_resonance_candidate", std::move(candidate_mol));
        if (normalization == "resonance_rule_normalization")
        {
            machine.RunOmolStage(normalization, molgr::reconstruct::CleanResonances);
        }
        else
        {
            machine.RunOmolChargeStage(std::nullopt, molgr::reconstruct::Eliminate13Dipole);
            machine.RunOmolChargeStage(
                std::nullopt,
                molgr::reconstruct::EliminateCPLikeRadicalAnion);
            machine.RunOmolChargeStage(
                std::nullopt,
                molgr::reconstruct::EliminatePositiveCharges);
            machine.RunOmolChargeStage(
                std::nullopt,
                molgr::reconstruct::EliminateNegativeCharges);
            machine.RunOmolChargeStage(
                std::nullopt,
                molgr::reconstruct::EliminatePositiveCharges);
            machine.RunOmolStage(std::nullopt, molgr::reconstruct::CleanNeighborRadicals);
            machine.RunOmolStage(std::nullopt, molgr::reconstruct::CleanResonances);
            machine.Annotate(normalization);
        }

        ProcessedStateKey processed_key{
            molgr::resonance::BuildProcessedResonanceKey(machine.EnsureUniqueMol()),
            machine.given_charge,
        };
        if (!seen_processed_states.insert(std::move(processed_key)).second)
        {
            return std::nullopt;
        }
        if (!molgr::reconstruct::ValidateOmol(
                machine.EnsureUniqueMol(),
                seed.total_charge,
                seed.total_radical_electrons))
        {
            return std::nullopt;
        }

        machine.Annotate("validate_no_metal_candidate");
        machine.metadata["resonance_seed_index"] = static_cast<int>(raw.seed_index);
        machine.metadata["resonance_index"] = static_cast<int>(raw.resonance_index);
        machine.metadata["resonance_raw_index"] = static_cast<int>(raw.raw_index);
        machine.metadata["resonance_normalization"] = normalization;
        auto candidate = machine.FreezeLike(seed);
        try
        {
            molgr::no_metals::selection::ScoreReconstructionCandidate(candidate, config);
        }
        catch (const std::exception &)
        {
            return std::nullopt;
        }
        molgr::no_metals::selection::AnnotateNoMetalCandidateTopology(candidate, config);
        return candidate;
    }
}

namespace molgr
{
    namespace no_metals
    {
        namespace resonance
        {
            std::vector<molgr::state::ReconstructionState> BuildResonanceSeedPool(
                std::vector<molgr::state::ReconstructionState> neighbor_seeds)
            {
                auto pool = DeduplicateStates(std::move(neighbor_seeds));
                AppendOmolStageVariants(
                    pool,
                    "relocate_carbene_radical_for_resonance",
                    reconstruct::CleanCarbeneNeighborUnsaturated);
                AppendChargeStageVariants(
                    pool,
                    "assign_negative_charges_from_radicals",
                    reconstruct::AssignNegativeChargesFromRadicals);
                AppendOmolStageVariants(
                    pool,
                    "refresh_electronic_labels_for_resonance",
                    reconstruct::FreshOmolChargeRadical);
                return pool;
            }

            std::vector<molgr::state::ReconstructionState> SearchResonanceCandidates(
                const std::vector<molgr::state::ReconstructionState> &states,
                const molgr::config::MolGRConfig &config,
                molgr::pipeline::perf::RunTimingReducer *timing_reducer,
                ResonanceSearchSession *session)
            {
                ResonanceSearchSession local_session;
                auto &search_session = session == nullptr ? local_session : *session;
                std::vector<RawResonanceCandidate> raw_candidates;
                const int resonance_max_depth = std::max(0, config.resonance.max_depth);
                const int max_discrepancy = std::max(
                    0,
                    config.resonance.limited_discrepancy_max_discrepancy);

                const auto walk_started = std::chrono::steady_clock::now();
                std::size_t pruned_expansion_count = 0;
                for (std::size_t seed_index = 0; seed_index < states.size(); ++seed_index)
                {
                    std::size_t resonance_index = 0;
                    const auto &state = states[seed_index];
                    reconstruct::WalkRadicalResonancesLimitedDiscrepancy(
                        state.Mol(),
                        resonance_max_depth,
                        [&](const reconstruct::ResonanceSearchNode &node) -> bool
                        {
                            RawStateKey raw_key{node.state_key, state.given_charge};
                            bool should_expand = true;
                            if (node.depth < resonance_max_depth)
                            {
                                should_expand = RegisterTraversalLabel(
                                    search_session.labels_by_state,
                                    raw_key,
                                    {node.depth, node.discrepancy});
                                if (!should_expand)
                                {
                                    ++pruned_expansion_count;
                                }
                            }
                            if (search_session.seen_raw_states.insert(raw_key).second)
                            {
                                raw_candidates.push_back(
                                    RawResonanceCandidate{
                                        seed_index,
                                        resonance_index,
                                        search_session.next_raw_index++,
                                        std::make_shared<OpenBabel::OBMol>(
                                            molgr::utils::CloneMolTopologyOnly(node.omol)),
                                    });
                            }
                            ++resonance_index;
                            return should_expand;
                        },
                        reconstruct::LimitedDiscrepancyTraversalConfig{max_discrepancy},
                        config);
                }
                if (timing_reducer != nullptr)
                {
                    const auto now = std::chrono::steady_clock::now();
                    timing_reducer->AddResonanceWalkMs(
                        std::chrono::duration<double, std::milli>(now - walk_started).count());
                    timing_reducer->AddResonanceRawCandidates(
                        static_cast<double>(raw_candidates.size()));
                    timing_reducer->AddResonancePrunedExpansions(
                        static_cast<double>(pruned_expansion_count));
                }

                std::vector<molgr::state::ReconstructionState> candidates;
                const auto prepare_started = std::chrono::steady_clock::now();
                for (const auto &raw : raw_candidates)
                {
                    const auto &seed = states[raw.seed_index];
                    for (const std::string normalization : {
                             "resonance_rule_normalization",
                             "full_resonance_normalization",
                         })
                    {
                        auto candidate = PrepareCandidateVariant(
                            raw,
                            seed,
                            normalization,
                            search_session.seen_processed_states,
                            config);
                        if (candidate.has_value())
                        {
                            candidates.push_back(std::move(*candidate));
                        }
                    }
                }
                if (timing_reducer != nullptr)
                {
                    const auto now = std::chrono::steady_clock::now();
                    timing_reducer->AddResonancePrepareMs(
                        std::chrono::duration<double, std::milli>(now - prepare_started).count());
                    timing_reducer->AddResonancePreparedCandidates(
                        static_cast<double>(raw_candidates.size() * 2));
                    timing_reducer->AddResonanceValidCandidates(
                        static_cast<double>(candidates.size()));
                    timing_reducer->AddResonanceDedupCandidates(
                        static_cast<double>(candidates.size()));
                }
                return candidates;
            }
        }
    }
}
