#include "molgr/utils/no_metals/preparation.h"

#include "molgr/stages/break_bond.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/fresh.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/xyz.h"

#include "molgr/compat/openbabel_iter.h"

#include <openbabel/bond.h>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace preparation
        {
            namespace
            {
                constexpr const char *kNeighborRadicalResolutionStrategyKey =
                    "neighbor_radical_resolution_strategy";

                void NormalizeSeedElectronicLabels(OpenBabel::OBMol &mol)
                {
                    FOR_ATOMS_OF_MOL(atom_iter, mol)
                    {
                        atom_iter->SetFormalCharge(0);
                        atom_iter->SetSpinMultiplicity(0);
                    }
                }

                std::vector<std::pair<int, int>> NeighborRadicalBondPairs(OpenBabel::OBMol &mol)
                {
                    std::vector<std::pair<int, int>> pairs;
                    FOR_BONDS_OF_MOL(bond_iter, mol)
                    {
                        OpenBabel::OBBond *bond = &(*bond_iter);
                        OpenBabel::OBAtom *begin_atom = bond->GetBeginAtom();
                        OpenBabel::OBAtom *end_atom = bond->GetEndAtom();
                        if (begin_atom != nullptr && end_atom != nullptr &&
                            begin_atom->GetSpinMultiplicity() > 0 &&
                            end_atom->GetSpinMultiplicity() > 0)
                        {
                            pairs.emplace_back(begin_atom->GetIdx(), end_atom->GetIdx());
                        }
                    }
                    return pairs;
                }

                bool CleanNeighborRadicalsChargeSplit(
                    OpenBabel::OBMol &mol,
                    int begin_charge_sign)
                {
                    bool hit = false;
                    const auto pairs = NeighborRadicalBondPairs(mol);
                    for (const auto &[begin_idx, end_idx] : pairs)
                    {
                        OpenBabel::OBAtom *begin_atom = mol.GetAtom(begin_idx);
                        OpenBabel::OBAtom *end_atom = mol.GetAtom(end_idx);
                        if (begin_atom == nullptr || end_atom == nullptr)
                        {
                            continue;
                        }
                        const int spin_to_consume = std::min(
                            begin_atom->GetSpinMultiplicity(),
                            end_atom->GetSpinMultiplicity());
                        if (spin_to_consume <= 0)
                        {
                            continue;
                        }
                        begin_atom->SetSpinMultiplicity(
                            begin_atom->GetSpinMultiplicity() - spin_to_consume);
                        end_atom->SetSpinMultiplicity(
                            end_atom->GetSpinMultiplicity() - spin_to_consume);
                        begin_atom->SetFormalCharge(
                            begin_atom->GetFormalCharge() + begin_charge_sign * spin_to_consume);
                        end_atom->SetFormalCharge(
                            end_atom->GetFormalCharge() - begin_charge_sign * spin_to_consume);
                        hit = true;
                    }
                    return hit;
                }

                molgr::state::OmolStateMachine RunDeterministicPreResolutionStages(
                    const molgr::state::ReconstructionState &state)
                {
                    auto machine = molgr::state::OmolStateMachine::FromReconstructionState(state);
                    machine.RunOmolStage("make_connections", reconstruct::MakeConnections, 0.15);
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
                    return machine;
                }

                molgr::state::OmolStateMachine RunDirectNeighborRadicalResolution(
                    molgr::state::OmolStateMachine machine)
                {
                    machine.metadata[kNeighborRadicalResolutionStrategyKey] =
                        std::string("direct");
                    machine.RunOmolStage(
                        "clean_neighbor_radicals",
                        reconstruct::CleanNeighborRadicals);
                    return machine;
                }

                molgr::state::OmolStateMachine RunChargeSplitNeighborRadicalResolution(
                    molgr::state::OmolStateMachine machine,
                    const std::string &strategy,
                    const std::string &phase,
                    int begin_charge_sign)
                {
                    machine.metadata[kNeighborRadicalResolutionStrategyKey] = strategy;
                    machine.RunOmolStage(
                        phase,
                        CleanNeighborRadicalsChargeSplit,
                        begin_charge_sign);
                    return machine;
                }

                molgr::state::ReconstructionState RunDeterministicPostResolutionStages(
                    molgr::state::OmolStateMachine machine,
                    const molgr::state::ReconstructionState &state)
                {
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

                std::vector<molgr::state::OmolStateMachine> EnumerateNeighborRadicalResolutionMachines(
                    molgr::state::OmolStateMachine machine)
                {
                    if (NeighborRadicalBondPairs(machine.EnsureUniqueMol()).empty())
                    {
                        return {RunDirectNeighborRadicalResolution(std::move(machine))};
                    }

                    auto direct_mol = std::make_shared<OpenBabel::OBMol>(
                        molgr::utils::CloneMolTopologyOnly(machine.EnsureUniqueMol()));
                    auto direct_machine = machine.Branch(
                        std::nullopt,
                        std::move(direct_mol),
                        std::nullopt);
                    direct_machine = RunDirectNeighborRadicalResolution(std::move(direct_machine));

                    auto begin_positive_mol = std::make_shared<OpenBabel::OBMol>(
                        molgr::utils::CloneMolTopologyOnly(machine.EnsureUniqueMol()));
                    auto begin_positive_machine = machine.Branch(
                        std::nullopt,
                        std::move(begin_positive_mol),
                        std::nullopt);
                    begin_positive_machine = RunChargeSplitNeighborRadicalResolution(
                        std::move(begin_positive_machine),
                        "charge_begin_positive",
                        "clean_neighbor_radicals_charge_begin_positive",
                        1);

                    auto begin_negative_mol = std::make_shared<OpenBabel::OBMol>(
                        molgr::utils::CloneMolTopologyOnly(machine.EnsureUniqueMol()));
                    auto begin_negative_machine = machine.Branch(
                        std::nullopt,
                        std::move(begin_negative_mol),
                        std::nullopt);
                    begin_negative_machine = RunChargeSplitNeighborRadicalResolution(
                        std::move(begin_negative_machine),
                        "charge_begin_negative",
                        "clean_neighbor_radicals_charge_begin_negative",
                        -1);

                    std::vector<molgr::state::OmolStateMachine> branches;
                    branches.reserve(3);
                    branches.push_back(std::move(direct_machine));
                    branches.push_back(std::move(begin_positive_machine));
                    branches.push_back(std::move(begin_negative_machine));
                    return branches;
                }
            }

            molgr::state::ReconstructionState BuildSeedState(
                std::shared_ptr<OpenBabel::OBMol> omol,
                int total_charge,
                int total_radical_electrons)
            {
                return molgr::state::ReconstructionState(
                    std::move(omol),
                    0,
                    total_charge,
                    total_radical_electrons,
                    {"read_xyz", "normalize_seed_electronic_labels"},
                    {{"source", std::string("xyz_to_omol_no_metal_state")}},
                    0);
            }

            std::shared_ptr<OpenBabel::OBMol> SeedOmolFromXyzBlock(
                const std::string &xyz_block)
            {
                auto omol = std::make_shared<OpenBabel::OBMol>();
                if (!molgr::utils::ReadXyzBlockToMol(xyz_block, omol.get()))
                {
                    return nullptr;
                }
                NormalizeSeedElectronicLabels(*omol);
                return omol;
            }

            std::shared_ptr<OpenBabel::OBMol> NormalizeSeedOmolCopy(
                const OpenBabel::OBMol &seed_omol)
            {
                auto omol = std::make_shared<OpenBabel::OBMol>(
                    molgr::utils::CloneMolTopologyOnly(seed_omol));
                NormalizeSeedElectronicLabels(*omol);
                return omol;
            }

            molgr::state::ReconstructionState SeedStateFromOmol(
                const OpenBabel::OBMol &seed_omol,
                int total_charge,
                int total_radical_electrons)
            {
                auto omol = NormalizeSeedOmolCopy(seed_omol);
                return BuildSeedState(
                    std::move(omol),
                    total_charge,
                    total_radical_electrons);
            }

            molgr::state::ReconstructionState SeedState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons)
            {
                auto omol = SeedOmolFromXyzBlock(xyz_block);
                if (!omol)
                {
                    return {};
                }
                return BuildSeedState(
                    std::move(omol),
                    total_charge,
                    total_radical_electrons);
            }

            molgr::state::ReconstructionState RunLinearPipeline(
                const molgr::state::ReconstructionState &state)
            {
                auto machine = RunDeterministicPreResolutionStages(state);
                machine = RunDirectNeighborRadicalResolution(std::move(machine));
                return RunDeterministicPostResolutionStages(std::move(machine), state);
            }

            std::vector<molgr::state::ReconstructionState> EnumerateNoMetalCandidateStates(
                const molgr::state::ReconstructionState &state)
            {
                auto machine = RunDeterministicPreResolutionStages(state);
                auto branches = EnumerateNeighborRadicalResolutionMachines(std::move(machine));
                std::vector<molgr::state::ReconstructionState> candidates;
                candidates.reserve(branches.size());
                for (auto &branch : branches)
                {
                    candidates.push_back(RunDeterministicPostResolutionStages(std::move(branch), state));
                }
                return candidates;
            }
        }
    }
}
