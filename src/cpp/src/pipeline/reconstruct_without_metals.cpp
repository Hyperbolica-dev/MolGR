#include "molgr/pipeline/reconstruct_without_metals.h"

#include "molgr/pipeline/resonance.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/no_metals/preparation.h"
#include "molgr/utils/no_metals/resonance.h"
#include "molgr/utils/no_metals/selection.h"
#include "molgr/utils/perf.h"
#include "molgr/utils/utils.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <variant>

namespace molgr
{
    namespace pipeline
    {
        namespace reconstruct_without_metals
        {
            namespace
            {
                void RecordNoMetalElapsed(
                    const std::chrono::steady_clock::time_point &started,
                    perf::RunTimingReducer *timing_reducer)
                {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_ms =
                        std::chrono::duration<double, std::milli>(now - started).count();
                    if (timing_reducer != nullptr)
                    {
                        timing_reducer->AddNoMetalPipelineMs(elapsed_ms);
                    }
                }

                void RecordResonanceElapsed(
                    const std::chrono::steady_clock::time_point &started,
                    perf::RunTimingReducer *timing_reducer)
                {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_ms =
                        std::chrono::duration<double, std::milli>(now - started).count();
                    if (timing_reducer != nullptr)
                    {
                        timing_reducer->AddResonanceHandlingEnumerationMs(elapsed_ms);
                    }
                }
            }

            std::optional<molgr::state::ReconstructionState> RunNoMetalPipelineFromSeedState(
                molgr::state::ReconstructionState state,
                const molgr::config::MolGRConfig &config,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle,
                const std::chrono::steady_clock::time_point &no_metal_started)
            {
                if (state.total_radical_electrons < 0)
                {
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                if (!state.omol)
                {
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }
                const auto linear_started = std::chrono::steady_clock::now();
                auto candidate_states =
                    molgr::no_metals::preparation::EnumerateNoMetalCandidateStates(state);
                if (timing_reducer != nullptr)
                {
                    const auto linear_now = std::chrono::steady_clock::now();
                    timing_reducer->AddNoMetalLinearPipelineMs(
                        std::chrono::duration<double, std::milli>(linear_now - linear_started).count());
                }

                std::vector<molgr::state::ReconstructionState> direct_candidates;
                std::vector<molgr::state::ReconstructionState> resonance_candidates;
                const auto resonance_started = std::chrono::steady_clock::now();
                bool recovered_any_resonance = false;
                for (auto &linear_state : candidate_states)
                {
                    const auto validate_started = std::chrono::steady_clock::now();
                    const bool direct_candidate_valid =
                        reconstruct::ValidateOmol(
                            linear_state.MutableMol(),
                            linear_state.total_charge,
                            linear_state.total_radical_electrons);
                    if (timing_reducer != nullptr)
                    {
                        const auto validate_now = std::chrono::steady_clock::now();
                        timing_reducer->AddNoMetalValidateMs(
                            std::chrono::duration<double, std::milli>(
                                validate_now - validate_started)
                                .count());
                    }
                    if (direct_candidate_valid)
                    {
                        auto result_machine =
                            molgr::state::OmolStateMachine::FromReconstructionState(linear_state);
                        auto result_mol = std::make_shared<OpenBabel::OBMol>(
                            molgr::utils::CloneMolTopologyOnly(linear_state.Mol()));
                        result_machine = result_machine.Branch(
                            std::nullopt,
                            std::move(result_mol));
                        result_machine.Annotate("validate_direct_candidate");
                        result_machine.RunOmolStage("clean_resonances", reconstruct::CleanResonances);
                        auto result_state = result_machine.FreezeLike(linear_state);
                        bool direct_candidate_scored = false;
                        try
                        {
                            molgr::no_metals::selection::ScoreReconstructionCandidate(
                                result_state,
                                config);
                            direct_candidate_scored = true;
                        }
                        catch (const std::exception &)
                        {
                            // Drop only the direct fallback; this branch still enters resonance recovery.
                        }
                        if (direct_candidate_scored)
                        {
                            molgr::no_metals::selection::AnnotateNoMetalCandidateTopology(
                                result_state,
                                config);
                            direct_candidates.push_back(std::move(result_state));
                        }
                    }

                    auto recovered_resonances =
                        molgr::no_metals::resonance::RecoverResonanceCandidates(
                            linear_state,
                            config,
                            timing_reducer);
                    recovered_any_resonance = true;
                    for (auto &candidate : recovered_resonances)
                    {
                        resonance_candidates.push_back(std::move(candidate));
                    }
                }
                std::vector<molgr::state::ReconstructionState> *candidate_pool = nullptr;
                if (!resonance_candidates.empty())
                {
                    candidate_pool = &resonance_candidates;
                }
                else if (!direct_candidates.empty())
                {
                    candidate_pool = &direct_candidates;
                }

                if (candidate_pool == nullptr)
                {
                    if (recovered_any_resonance)
                    {
                        RecordResonanceElapsed(resonance_started, timing_reducer);
                    }
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                std::size_t best_idx = 0;
                std::optional<std::tuple<double, int, double, double, double, double>>
                    best_selection_key;
                for (std::size_t i = 0; i < candidate_pool->size(); ++i)
                {
                    const auto selection_key =
                        molgr::no_metals::selection::NoMetalCandidateSelectionKey(
                            (*candidate_pool)[i],
                            config);
                    if (!best_selection_key.has_value() || selection_key < *best_selection_key)
                    {
                        best_selection_key = selection_key;
                        best_idx = i;
                    }
                }

                if (!best_selection_key.has_value())
                {
                    if (recovered_any_resonance)
                    {
                        RecordResonanceElapsed(resonance_started, timing_reducer);
                    }
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                auto result_state = std::move((*candidate_pool)[best_idx]);
                if (std::find(
                        result_state.phase_history.begin(),
                        result_state.phase_history.end(),
                        "validate_resonance_candidate") != result_state.phase_history.end())
                {
                    auto result_machine =
                        molgr::state::OmolStateMachine::FromReconstructionState(result_state);
                    result_machine.Annotate("select_best_resonance_candidate");
                    result_state = result_machine.FreezeLike(result_state);
                }

                if (recovered_any_resonance)
                {
                    RecordResonanceElapsed(resonance_started, timing_reducer);
                }
                if (preheat_score_bundle)
                {
                    result_state.PreheatScoreBundle(config);
                }
                RecordNoMetalElapsed(no_metal_started, timing_reducer);
                return result_state;
            }

            std::optional<molgr::state::ReconstructionState> XyzToOmolNoMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle)
            {
                const auto no_metal_started = std::chrono::steady_clock::now();
                auto state = molgr::no_metals::preparation::SeedState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                return RunNoMetalPipelineFromSeedState(
                    std::move(state),
                    config,
                    timing_reducer,
                    preheat_score_bundle,
                    no_metal_started);
            }

            std::optional<molgr::state::ReconstructionState> SeedOmolToOmolNoMetalState(
                const OpenBabel::OBMol &seed_omol,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle)
            {
                const auto no_metal_started = std::chrono::steady_clock::now();
                auto state = molgr::no_metals::preparation::SeedStateFromOmol(
                    seed_omol,
                    total_charge,
                    total_radical_electrons);
                return RunNoMetalPipelineFromSeedState(
                    std::move(state),
                    config,
                    timing_reducer,
                    preheat_score_bundle,
                    no_metal_started);
            }

            std::optional<molgr::state::ReconstructionState> SeedOmolCopyToOmolNoMetalState(
                std::shared_ptr<OpenBabel::OBMol> seed_omol,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle)
            {
                if (!seed_omol)
                {
                    return std::nullopt;
                }
                const auto no_metal_started = std::chrono::steady_clock::now();
                auto state = molgr::no_metals::preparation::BuildSeedState(
                    std::move(seed_omol),
                    total_charge,
                    total_radical_electrons);
                return RunNoMetalPipelineFromSeedState(
                    std::move(state),
                    config,
                    timing_reducer,
                    preheat_score_bundle,
                    no_metal_started);
            }

            std::unique_ptr<molgr::utils::MoleculeData> XyzToMolDataNoMetal(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config)
            {
                molgr::pipeline::perf::RunTimingScope timing_scope;
                auto state = XyzToOmolNoMetalState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons,
                    config,
                    &timing_scope.Reducer(),
                    false);
                if (!state.has_value())
                {
                    return nullptr;
                }
                return std::make_unique<molgr::utils::MoleculeData>(
                    molgr::utils::MoleculeDataFromOBMol(state->Mol()));
            }

            std::vector<DebugNoMetalCandidateSummary> DebugNoMetalResonanceCandidateSummaries(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config)
            {
                std::vector<DebugNoMetalCandidateSummary> summaries;
                if (total_radical_electrons < 0)
                {
                    return summaries;
                }

                auto state = molgr::no_metals::preparation::SeedState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                if (!state.omol)
                {
                    return summaries;
                }
                auto linear_states =
                    molgr::no_metals::preparation::EnumerateNoMetalCandidateStates(state);
                std::vector<molgr::state::ReconstructionState> candidates;
                for (auto &linear_state : linear_states)
                {
                    auto branch_candidates =
                        molgr::no_metals::resonance::RecoverResonanceCandidates(
                            linear_state,
                            config);
                    for (auto &candidate : branch_candidates)
                    {
                        candidates.push_back(std::move(candidate));
                    }
                }
                summaries.reserve(candidates.size());
                for (auto &candidate : candidates)
                {
                    molgr::no_metals::selection::NoMetalCandidateSelectionKey(candidate, config);
                    const auto get_int = [&](const std::string &key, int fallback = 0)
                    {
                        const auto it = candidate.metadata.find(key);
                        if (it == candidate.metadata.end())
                        {
                            return fallback;
                        }
                        if (const auto *value = std::get_if<int>(&it->second))
                        {
                            return *value;
                        }
                        if (const auto *value = std::get_if<double>(&it->second))
                        {
                            return static_cast<int>(*value);
                        }
                        return fallback;
                    };
                    const auto get_double = [&](const std::string &key, double fallback = 0.0)
                    {
                        const auto it = candidate.metadata.find(key);
                        if (it == candidate.metadata.end())
                        {
                            return fallback;
                        }
                        if (const auto *value = std::get_if<double>(&it->second))
                        {
                            return *value;
                        }
                        if (const auto *value = std::get_if<int>(&it->second))
                        {
                            return static_cast<double>(*value);
                        }
                        return fallback;
                    };

                    summaries.push_back(DebugNoMetalCandidateSummary{
                        reconstruct::SmilesFirstToken(candidate.Mol()),
                        get_int("resonance_index", -1),
                        get_double("score"),
                        get_double("organic_aromatic_stability_score"),
                        get_int("organic_aromatic_atom_count"),
                        get_int("organic_max_conjugated_component_size"),
                        get_int("organic_conjugated_atom_count"),
                        get_int("organic_conjugated_bond_count"),
                        get_int("organic_formal_charge_absolute_sum"),
                        get_double("organic_conjugation_charge_penalty"),
                        get_double("organic_adjusted_max_conjugated_component_size"),
                        get_double("organic_adjusted_conjugated_atom_count"),
                        get_double("organic_adjusted_conjugated_bond_count"),
                    });
                }
                return summaries;
            }
        }
    }
}
