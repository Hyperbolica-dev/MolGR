/**
 * @file initial_reconstructor.cpp
 * @brief Implementation of initial reconstruction logic.
 * @details STRICTLY aligned with Python 'GraphReconstruction.py'.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/pipeline/reconstruct_without_metals.h"

#include "molgr/utils/logger.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <chrono>
#include <memory>

namespace molgr
{
    namespace pipeline
    {
        namespace perf
        {
            namespace
            {
                thread_local molgr::pipeline::perf::RunTimingBreakdown t_run_timing_breakdown;
            }

            void ResetRunTimingBreakdown()
            {
                t_run_timing_breakdown = molgr::pipeline::perf::RunTimingBreakdown{};
            }

            molgr::pipeline::perf::RunTimingBreakdown GetRunTimingBreakdown()
            {
                return t_run_timing_breakdown;
            }

            void AddNoMetalPipelineMs(double delta_ms)
            {
                t_run_timing_breakdown.no_metal_pipeline_ms += delta_ms;
            }

            void AddResonanceHandlingEnumerationMs(double delta_ms)
            {
                t_run_timing_breakdown.resonance_handling_enumeration_ms += delta_ms;
            }

            void AddMetalEnumerationCombinationMs(double delta_ms)
            {
                t_run_timing_breakdown.metal_enumeration_combination_ms += delta_ms;
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

        std::unique_ptr<OBMol> ReconstructFromXYZNoMetal(const std::string &xyz_block, int total_charge, int total_radical)
        {
            LOG_DEBUG("[ReconstructNoMetal] Start. Target Charge=" << total_charge << " Radical=" << total_radical);

            if (total_radical < 0)
                return nullptr;

            auto mol = std::make_unique<OBMol>();
            OBConversion &conv = ThreadLocalXyzInConversion();
            if (!conv.ReadString(mol.get(), xyz_block))
            {
                return nullptr;
            }

            MakeConnections(*mol);
            PreClean(*mol);

            FreshOmolChargeRadical(*mol);

            int current_charge_sum = 0;
            FOR_ATOMS_OF_MOL(a, *mol)
            current_charge_sum += a->GetFormalCharge();
            int given_charge = total_charge - current_charge_sum;

            EliminateNNN(*mol, given_charge, false);
            EliminateHighPositiveChargeAtoms(*mol, given_charge);
            EliminateCNInDoubt(*mol, given_charge);
            EliminateNNN(*mol, given_charge, true);
            EliminateCarboxyl(*mol, given_charge);
            CleanCarbeneNeighborUnsaturated(*mol);
            EliminateCarbeneNeighborHeteroatom(*mol, given_charge);
            CleanNeighborRadicals(*mol);
            CleanCarbeneNeighborUnsaturated(*mol);

            EliminateChargeSpliting(*mol, given_charge);

            BreakDeformedEne(*mol, given_charge, total_radical);
            BreakOneBond(*mol, given_charge, total_radical);

            FreshOmolChargeRadical(*mol);

            return mol;
        }
    }

    namespace pipeline
    {
        namespace reconstruct_without_metals
        {
            std::unique_ptr<molgr::utils::MoleculeData> XyzToMolDataNoMetal(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons)
            {
                molgr::pipeline::perf::ResetRunTimingBreakdown();
                const auto no_metal_started = std::chrono::steady_clock::now();
                const auto record_no_metal_elapsed = [&]()
                {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_ms = std::chrono::duration<double, std::milli>(now - no_metal_started).count();
                    molgr::pipeline::perf::AddNoMetalPipelineMs(elapsed_ms);
                };

                if (total_radical_electrons < 0)
                {
                    record_no_metal_elapsed();
                    return nullptr;
                }

                OpenBabel::OBMol mol;
                OpenBabel::OBConversion &conv = reconstruct::ThreadLocalXyzInConversion();
                if (!conv.ReadString(&mol, xyz_block))
                {
                    record_no_metal_elapsed();
                    return nullptr;
                }

                reconstruct::MakeConnections(mol);
                reconstruct::PreClean(mol);
                reconstruct::FreshOmolChargeRadical(mol);

                int formal_charge_sum = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    formal_charge_sum += atom_iter->GetFormalCharge();
                }
                int given_charge = total_charge - formal_charge_sum;

                reconstruct::EliminateNNN(mol, given_charge, false);
                reconstruct::EliminateHighPositiveChargeAtoms(mol, given_charge);
                reconstruct::EliminateCNInDoubt(mol, given_charge);
                reconstruct::EliminateNNN(mol, given_charge, true);
                reconstruct::EliminateCarboxyl(mol, given_charge);
                reconstruct::CleanCarbeneNeighborUnsaturated(mol);
                reconstruct::EliminateCarbeneNeighborHeteroatom(mol, given_charge);
                reconstruct::CleanNeighborRadicals(mol);
                reconstruct::CleanCarbeneNeighborUnsaturated(mol);
                reconstruct::EliminateChargeSpliting(mol, given_charge);
                reconstruct::BreakDeformedEne(mol, given_charge, total_radical_electrons);
                reconstruct::BreakOneBond(mol, given_charge, total_radical_electrons);
                reconstruct::FreshOmolChargeRadical(mol);

                if (reconstruct::ValidateOmol(mol, total_charge, total_radical_electrons))
                {
                    reconstruct::CleanResonances(mol);
                    record_no_metal_elapsed();
                    return std::make_unique<molgr::utils::MoleculeData>(molgr::utils::MoleculeDataFromOBMol(mol));
                }

                const auto resonance_started = std::chrono::steady_clock::now();
                const auto possible_resonances = reconstruct::GetRadicalResonances(mol);
                std::vector<OpenBabel::OBMol> recovered_resonances;
                recovered_resonances.reserve(possible_resonances.size());

                for (const auto &resonance : possible_resonances)
                {
                    int charge = given_charge;
                    auto processed = reconstruct::ProcessResonance(resonance, charge);
                    if (reconstruct::ValidateOmol(processed.first, total_charge, total_radical_electrons, false))
                    {
                        recovered_resonances.push_back(std::move(processed.first));
                    }
                }

                if (recovered_resonances.empty())
                {
                    const auto resonance_now = std::chrono::steady_clock::now();
                    const double resonance_ms = std::chrono::duration<double, std::milli>(resonance_now - resonance_started).count();
                    molgr::pipeline::perf::AddResonanceHandlingEnumerationMs(resonance_ms);
                    record_no_metal_elapsed();
                    return nullptr;
                }

                size_t best_idx = 0;
                double best_score = molgr::scoring::OmolScore(recovered_resonances[0]);
                for (size_t i = 1; i < recovered_resonances.size(); ++i)
                {
                    const double score = molgr::scoring::OmolScore(recovered_resonances[i]);
                    if (score < best_score)
                    {
                        best_idx = i;
                        best_score = score;
                    }
                }

                const auto resonance_now = std::chrono::steady_clock::now();
                const double resonance_ms = std::chrono::duration<double, std::milli>(resonance_now - resonance_started).count();
                molgr::pipeline::perf::AddResonanceHandlingEnumerationMs(resonance_ms);
                record_no_metal_elapsed();
                return std::make_unique<molgr::utils::MoleculeData>(
                    molgr::utils::MoleculeDataFromOBMol(recovered_resonances[best_idx]));
            }
        }
    }
}
