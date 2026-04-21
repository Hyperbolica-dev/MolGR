/**
 * @file initial_reconstructor.cpp
 * @brief Implementation of initial reconstruction logic.
 * @details STRICTLY aligned with Python 'GraphReconstruction.py'.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/pipeline/reconstruct_without_metals.h"

#include "molgr/state.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <optional>

namespace molgr
{
    namespace pipeline
    {
        namespace perf
        {
            namespace
            {
                std::mutex t_last_run_timing_breakdown_mutex;
                molgr::pipeline::perf::RunTimingBreakdown t_last_run_timing_breakdown;
            }

            void RunTimingReducer::AddNoMetalPipelineMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.no_metal_pipeline_ms += delta_ms;
            }

            void RunTimingReducer::AddResonanceHandlingEnumerationMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.resonance_handling_enumeration_ms += delta_ms;
            }

            void RunTimingReducer::AddMetalEnumerationCombinationMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.metal_enumeration_combination_ms += delta_ms;
            }

            molgr::pipeline::perf::RunTimingBreakdown RunTimingReducer::Snapshot() const
            {
                std::lock_guard<std::mutex> lock(mutex_);
                return timing_;
            }

            RunTimingReducer &RunTimingScope::Reducer()
            {
                return reducer_;
            }

            const RunTimingReducer &RunTimingScope::Reducer() const
            {
                return reducer_;
            }

            RunTimingScope::~RunTimingScope()
            {
                SetRunTimingBreakdown(reducer_.Snapshot());
            }

            molgr::pipeline::perf::RunTimingBreakdown GetRunTimingBreakdown()
            {
                std::lock_guard<std::mutex> lock(t_last_run_timing_breakdown_mutex);
                return t_last_run_timing_breakdown;
            }

            void SetRunTimingBreakdown(const RunTimingBreakdown &timing)
            {
                std::lock_guard<std::mutex> lock(t_last_run_timing_breakdown_mutex);
                t_last_run_timing_breakdown = timing;
            }
        }
    }

    namespace reconstruct
    {
        using namespace OpenBabel;

        OBConversion &ThreadLocalXyzInConversion()
        {
            thread_local OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetInFormat("xyz");
                initialized = true;
            }
            return conv;
        }

        bool ValidateOmol(OBMol &mol, int total_charge, int total_radical, bool emit_warnings)
        {
            int charge_sum = 0;
            int radical_sum = 0;
            int radical_sum_singlet = 0;

            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom *atom = &(*atom_iter);
                charge_sum += atom->GetFormalCharge();
                int spin = atom->GetSpinMultiplicity();
                radical_sum += spin;
                radical_sum_singlet += (spin % 2);
            }

            if (charge_sum != total_charge)
            {
                if (emit_warnings)
                {
                    LOG_WARN("[Validate] Charge mismatch. Target: " << total_charge << ", Actual: " << charge_sum);
                }
                return false;
            }

            if (radical_sum_singlet == total_radical)
            {
                radical_sum = radical_sum_singlet;
            }

            if (radical_sum != total_radical)
            {
                if (emit_warnings)
                {
                    LOG_WARN("[Validate] Radical mismatch. Target: " << total_radical << ", Actual: " << radical_sum);
                }
                return false;
            }
            return true;
        }

    }

    namespace pipeline
    {
        namespace reconstruct_without_metals
        {
            namespace
            {
                constexpr int kDefaultResonanceSearchMaxDepth = 2;
                constexpr int kDefaultResonanceMaxDiscrepancy = 1;
                constexpr bool kDefaultResonanceFallbackToFullFrontier = true;
                constexpr double kResonanceIncumbentPruneMargin = 5.0;

                molgr::state::ReconstructionState SeedState(
                    const std::string &xyz_block,
                    int total_charge,
                    int total_radical_electrons)
                {
                    auto omol = std::make_shared<OpenBabel::OBMol>();
                    OpenBabel::OBConversion &conv = reconstruct::ThreadLocalXyzInConversion();
                    if (!conv.ReadString(omol.get(), xyz_block))
                    {
                        return {};
                    }
                    return molgr::state::ReconstructionState(
                        omol,
                        0,
                        total_charge,
                        total_radical_electrons,
                        {"read_xyz"},
                        {{"source", std::string("xyz_to_omol_no_metal_state")}},
                        0);
                }

                molgr::state::ReconstructionState RunLinearPipeline(
                    const molgr::state::ReconstructionState &state)
                {
                    auto machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                    machine.RunOmolStage("make_connections", reconstruct::MakeConnections, 1.4);
                    machine.RunOmolStage("pre_clean", reconstruct::PreClean);
                    machine.RunOmolStage(
                        "fresh_omol_charge_radical_initial",
                        reconstruct::FreshOmolChargeRadical);

                    int formal_charge_sum = 0;
                    FOR_ATOMS_OF_MOL(atom_iter, machine.EnsureUniqueMol())
                    {
                        formal_charge_sum += atom_iter->GetFormalCharge();
                    }
                    machine.SetGivenCharge(
                        "initialize_charge_budget",
                        state.total_charge - formal_charge_sum);

                    machine.RunOmolChargeStage("eliminate_NNN_negative", reconstruct::EliminateNNN, false);
                    machine.RunOmolChargeStage(
                        "eliminate_high_positive_charge_atoms",
                        reconstruct::EliminateHighPositiveChargeAtoms);
                    machine.RunOmolChargeStage(
                        "eliminate_CN_in_doubt",
                        reconstruct::EliminateCNInDoubt);
                    machine.RunOmolChargeStage("eliminate_NNN_positive", reconstruct::EliminateNNN, true);
                    machine.RunOmolChargeStage("eliminate_carboxyl", reconstruct::EliminateCarboxyl);
                    machine.RunOmolStage(
                        "clean_carbene_neighbor_unsaturated_first",
                        reconstruct::CleanCarbeneNeighborUnsaturated);
                    machine.RunOmolChargeStage(
                        "eliminate_carbene_neighbor_heteroatom",
                        reconstruct::EliminateCarbeneNeighborHeteroatom);
                    machine.RunOmolStage("clean_neighbor_radicals", reconstruct::CleanNeighborRadicals);
                    machine.RunOmolStage(
                        "clean_carbene_neighbor_unsaturated_second",
                        reconstruct::CleanCarbeneNeighborUnsaturated);
                    machine.RunOmolChargeStage(
                        "eliminate_charge_spliting",
                        reconstruct::EliminateChargeSpliting);
                    machine.RunOmolStage(
                        "break_deformed_ene",
                        reconstruct::BreakDeformedEne,
                        machine.given_charge,
                        state.total_radical_electrons,
                        5.0);
                    machine.RunOmolChargeStage(
                        "break_one_bond",
                        reconstruct::BreakOneBond,
                        state.total_radical_electrons);
                    machine.RunOmolStage(
                        "fresh_omol_charge_radical_final",
                        reconstruct::FreshOmolChargeRadical);
                    return machine.FreezeLike(state);
                }

                std::vector<molgr::state::ReconstructionState> RecoverResonanceCandidates(
                    const molgr::state::ReconstructionState &state)
                {
                    std::vector<molgr::state::ReconstructionState> candidates;
                    auto base_machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                    std::map<reconstruct::ProcessedResonanceKey, std::optional<double>> processed_state_scores;
                    reconstruct::DirectGainBoundCache bound_cache;
                    double best_score = std::numeric_limits<double>::infinity();
                    std::size_t resonance_index = 0;

                    reconstruct::WalkRadicalResonancesLimitedDiscrepancy(
                        state.Mol(),
                        kDefaultResonanceSearchMaxDepth,
                        [&](const reconstruct::ResonanceSearchNode &node) -> bool
                        {
                            const std::size_t current_resonance_index = resonance_index++;
                            auto branched_machine = base_machine.Branch(
                                "branch_resonance_candidate",
                                std::make_shared<OpenBabel::OBMol>(node.omol));

                            bool process_hit = false;
                            process_hit = branched_machine.RunOmolChargeStage(
                                std::nullopt,
                                reconstruct::Eliminate13Dipole) || process_hit;
                            process_hit = branched_machine.RunOmolChargeStage(
                                std::nullopt,
                                reconstruct::EliminatePositiveCharges) || process_hit;
                            process_hit = branched_machine.RunOmolChargeStage(
                                std::nullopt,
                                reconstruct::EliminateNegativeCharges) || process_hit;
                            process_hit = branched_machine.RunOmolStage(
                                std::nullopt,
                                reconstruct::CleanNeighborRadicals) || process_hit;
                            process_hit = branched_machine.RunOmolStage(
                                std::nullopt,
                                reconstruct::CleanResonances) || process_hit;
                            branched_machine.Annotate("process_resonance");

                            const reconstruct::ProcessedResonanceKey processed_state_key =
                                reconstruct::BuildProcessedResonanceKey(branched_machine.EnsureUniqueMol());

                            auto score_it = processed_state_scores.find(processed_state_key);
                            std::optional<double> cached_score;
                            if (score_it != processed_state_scores.end())
                            {
                                cached_score = score_it->second;
                            }
                            else
                            {
                                processed_state_scores.emplace(processed_state_key, std::nullopt);
                                if (reconstruct::ValidateOmol(
                                        branched_machine.EnsureUniqueMol(),
                                        state.total_charge,
                                        state.total_radical_electrons))
                                {
                                    branched_machine.Annotate("validate_resonance_candidate");
                                    branched_machine.metadata["resonance_index"] =
                                        static_cast<int>(current_resonance_index);
                                    auto candidate = branched_machine.FreezeLike(state);
                                    const double score = candidate.FullScore();
                                    processed_state_scores[processed_state_key] = score;
                                    cached_score = score;
                                    candidates.push_back(std::move(candidate));
                                    if (score < best_score)
                                    {
                                        best_score = score;
                                    }
                                }
                            }

                            const int remaining_steps =
                                kDefaultResonanceSearchMaxDepth - node.depth;
                            if (remaining_steps <= 0 || !cached_score.has_value())
                            {
                                return node.depth < kDefaultResonanceSearchMaxDepth;
                            }
                            if (*cached_score < best_score + kResonanceIncumbentPruneMargin)
                            {
                                return true;
                            }

                            const double optimistic_improvement =
                                reconstruct::EstimateRemainingResonanceScoreImprovementUpperBound(
                                    node.omol,
                                    node.state_key,
                                    remaining_steps,
                                    &bound_cache);
                            return *cached_score - optimistic_improvement < best_score;
                        },
                        reconstruct::LimitedDiscrepancyTraversalConfig{
                            kDefaultResonanceMaxDiscrepancy,
                            kDefaultResonanceFallbackToFullFrontier,
                        });

                    return candidates;
                }
            }

            std::optional<molgr::state::ReconstructionState> XyzToOmolNoMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                perf::RunTimingReducer *timing_reducer,
                bool preheat_score_bundle)
            {
                const auto no_metal_started = std::chrono::steady_clock::now();
                const auto record_no_metal_elapsed = [&]()
                {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_ms = std::chrono::duration<double, std::milli>(now - no_metal_started).count();
                    if (timing_reducer != nullptr)
                    {
                        timing_reducer->AddNoMetalPipelineMs(elapsed_ms);
                    }
                };

                if (total_radical_electrons < 0)
                {
                    record_no_metal_elapsed();
                    return std::nullopt;
                }

                auto state = SeedState(xyz_block, total_charge, total_radical_electrons);
                if (!state.omol)
                {
                    record_no_metal_elapsed();
                    return std::nullopt;
                }
                state = RunLinearPipeline(state);

                if (reconstruct::ValidateOmol(
                        state.MutableMol(),
                        total_charge,
                        total_radical_electrons))
                {
                    auto result_machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                    result_machine.Annotate("validate_direct_candidate");
                    result_machine.RunOmolStage("clean_resonances", reconstruct::CleanResonances);
                    auto result_state = result_machine.FreezeLike(state);
                    if (preheat_score_bundle)
                    {
                        result_state.PreheatScoreBundle();
                    }
                    record_no_metal_elapsed();
                    return result_state;
                }

                const auto resonance_started = std::chrono::steady_clock::now();
                auto recovered_resonances = RecoverResonanceCandidates(state);

                if (recovered_resonances.empty())
                {
                    const auto resonance_now = std::chrono::steady_clock::now();
                    const double resonance_ms = std::chrono::duration<double, std::milli>(resonance_now - resonance_started).count();
                    if (timing_reducer != nullptr)
                    {
                        timing_reducer->AddResonanceHandlingEnumerationMs(resonance_ms);
                    }
                    record_no_metal_elapsed();
                    return std::nullopt;
                }

                std::size_t best_idx = 0;
                double best_score = std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < recovered_resonances.size(); ++i)
                {
                    const double score = recovered_resonances[i].FullScore();
                    if (score < best_score)
                    {
                        best_idx = i;
                        best_score = score;
                    }
                }

                auto result_machine =
                    molgr::state::OmolStateMachine::FromReconstructionState(recovered_resonances[best_idx]);
                result_machine.Annotate("select_best_resonance_candidate");
                auto result_state = result_machine.FreezeLike(recovered_resonances[best_idx]);

                const auto resonance_now = std::chrono::steady_clock::now();
                const double resonance_ms =
                    std::chrono::duration<double, std::milli>(resonance_now - resonance_started).count();
                if (timing_reducer != nullptr)
                {
                    timing_reducer->AddResonanceHandlingEnumerationMs(resonance_ms);
                }
                if (preheat_score_bundle)
                {
                    result_state.PreheatScoreBundle();
                }
                record_no_metal_elapsed();
                return result_state;
            }

            std::unique_ptr<molgr::utils::MoleculeData> XyzToMolDataNoMetal(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons)
            {
                molgr::pipeline::perf::RunTimingScope timing_scope;
                auto state = XyzToOmolNoMetalState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons,
                    &timing_scope.Reducer(),
                    false);
                if (!state.has_value())
                {
                    return nullptr;
                }
                return std::make_unique<molgr::utils::MoleculeData>(
                    molgr::utils::MoleculeDataFromOBMol(state->Mol()));
            }
        }
    }
}
