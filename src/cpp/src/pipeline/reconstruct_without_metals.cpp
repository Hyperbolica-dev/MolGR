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

            std::optional<molgr::state::ReconstructionState> XyzToOmolNoMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle)
            {
                const auto no_metal_started = std::chrono::steady_clock::now();
                if (total_radical_electrons < 0)
                {
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                auto state = molgr::no_metals::preparation::SeedState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                if (!state.omol)
                {
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }
                state = molgr::no_metals::preparation::RunLinearPipeline(state);

                if (reconstruct::ValidateOmol(
                        state.MutableMol(),
                        total_charge,
                        total_radical_electrons))
                {
                    auto result_machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                    result_machine.Annotate("validate_direct_candidate");
                    result_machine.RunOmolStage("clean_resonances", reconstruct::CleanResonances);
                    auto result_state = result_machine.FreezeLike(state);
                    molgr::no_metals::selection::AnnotateNoMetalCandidateTopology(
                        result_state,
                        config);
                    if (preheat_score_bundle)
                    {
                        result_state.PreheatScoreBundle(config);
                    }
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return result_state;
                }

                const auto resonance_started = std::chrono::steady_clock::now();
                auto recovered_resonances =
                    molgr::no_metals::resonance::RecoverResonanceCandidates(state, config);

                if (recovered_resonances.empty())
                {
                    RecordResonanceElapsed(resonance_started, timing_reducer);
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                std::size_t best_idx = 0;
                std::optional<molgr::no_metals::selection::NoMetalTopologySelectionKey>
                    best_topology_key;
                std::vector<std::size_t> best_topology_indices;
                for (std::size_t i = 0; i < recovered_resonances.size(); ++i)
                {
                    const auto topology_key =
                        molgr::no_metals::selection::NoMetalCandidateTopologySelectionKey(
                            recovered_resonances[i],
                            config);
                    if (!best_topology_key.has_value() || topology_key < *best_topology_key)
                    {
                        best_topology_key = topology_key;
                        best_topology_indices.clear();
                        best_topology_indices.push_back(i);
                    }
                    else if (topology_key == *best_topology_key)
                    {
                        best_topology_indices.push_back(i);
                    }
                }

                if (best_topology_indices.empty())
                {
                    RecordResonanceElapsed(resonance_started, timing_reducer);
                    RecordNoMetalElapsed(no_metal_started, timing_reducer);
                    return std::nullopt;
                }

                double best_score = std::numeric_limits<double>::infinity();
                for (const std::size_t candidate_index : best_topology_indices)
                {
                    const double score = molgr::no_metals::selection::ScoreReconstructionCandidate(
                        recovered_resonances[candidate_index],
                        config);
                    if (score < best_score)
                    {
                        best_score = score;
                        best_idx = candidate_index;
                    }
                }

                auto result_machine =
                    molgr::state::OmolStateMachine::FromReconstructionState(recovered_resonances[best_idx]);
                result_machine.Annotate("select_best_resonance_candidate");
                auto result_state = result_machine.FreezeLike(recovered_resonances[best_idx]);

                RecordResonanceElapsed(resonance_started, timing_reducer);
                if (preheat_score_bundle)
                {
                    result_state.PreheatScoreBundle(config);
                }
                RecordNoMetalElapsed(no_metal_started, timing_reducer);
                return result_state;
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
                state = molgr::no_metals::preparation::RunLinearPipeline(state);
                auto candidates = molgr::no_metals::resonance::RecoverResonanceCandidates(state, config);
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
