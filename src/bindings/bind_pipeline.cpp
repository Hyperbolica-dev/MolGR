#include "bindings.h"

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/pipeline/reconstruct_without_metals.h"
#include "molgr/pipeline/resonance.h"
#include "molgr/python_config.h"
#include "molgr/stages/break_bond.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/fresh.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/force_field.h"
#include "molgr/utils/metals/scoring.h"
#include "molgr/utils/metals/search.h"
#include "molgr/utils/metals/preparation.h"
#include "molgr/utils/no_metals/preparation.h"
#include "molgr/utils/utils.h"
#include "molgr/vendor/forcefielduff.h"

#include <openbabel/obconversion.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace bind
    {
        namespace
        {
            OpenBabel::OBMol *require_obmol_ptr(intptr_t mol_ptr)
            {
                if (mol_ptr == 0)
                {
                    throw std::runtime_error("null OBMol pointer");
                }
                return reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
            }

            py::tuple metal_state_signature(const molgr::metal::MetalAtomPosition &metal_state)
            {
                py::tuple item(8);
                item[0] = metal_state.idx;
                item[1] = metal_state.symbol;
                item[2] = metal_state.element_idx;
                item[3] = metal_state.valence;
                item[4] = metal_state.radical_num;
                item[5] = static_cast<long long>(std::llround(metal_state.position_x * 1000000.0));
                item[6] = static_cast<long long>(std::llround(metal_state.position_y * 1000000.0));
                item[7] = static_cast<long long>(std::llround(metal_state.position_z * 1000000.0));
                return item;
            }

            py::list metal_state_choice_signature(
                const std::vector<molgr::metal::MetalAtomPosition> &choice)
            {
                py::list out;
                for (const auto &metal_state : choice)
                {
                    out.append(metal_state_signature(metal_state));
                }
                return out;
            }

            py::list metal_state_search_group_signature(
                const molgr::metal::search::MetalStateChoiceGroup &group)
            {
                py::list out;
                for (const auto &choice : group)
                {
                    out.append(metal_state_choice_signature(choice));
                }
                return out;
            }

            py::list metal_state_search_layer_signature(
                const molgr::metal::search::MetalStateSearchLayer &layer)
            {
                py::list out;
                for (const auto &group : layer)
                {
                    out.append(metal_state_search_group_signature(group));
                }
                return out;
            }

            py::list phase_history_signature(const std::vector<std::string> &phase_history)
            {
                py::list out;
                for (const auto &phase : phase_history)
                {
                    out.append(phase);
                }
                return out;
            }

            std::unique_ptr<OpenBabel::OBMol> mol_from_smiles(const std::string &smiles)
            {
                auto mol = std::make_unique<OpenBabel::OBMol>();
                OpenBabel::OBConversion conv;
                conv.SetInFormat("smi");
                if (!conv.ReadString(mol.get(), smiles))
                {
                    return nullptr;
                }
                return mol;
            }

            void bind_reconstruct_with_metals_public_ns(py::module_ &ns)
            {
                ns.def(
                    "xyz2omol",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons,
                       py::object config)
                        -> std::unique_ptr<molgr::utils::MoleculeData>
                    {
                        auto runtime_config = molgr::config::FromPython(config);
                        py::gil_scoped_release release;
                        auto mol_data = molgr::pipeline::reconstruct_with_metals::Xyz2OmolMolData(
                            xyz_block,
                            total_charge,
                            total_radical_electrons,
                            runtime_config);
                        if (!mol_data)
                        {
                            return nullptr;
                        }
                        return std::move(mol_data);
                    },
                    "Reconstruct molecule data from XYZ with metal-aware pipeline.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0,
                    py::kw_only(),
                    py::arg("config") = py::none());
            }

            void bind_reconstruct_with_metals_dev_ns(py::module_ &ns)
            {
                py::class_<molgr::metal::MetalAtomPosition>(ns, "MetalAtomPosition")
                    .def(py::init<>())
                    .def_readwrite("idx", &molgr::metal::MetalAtomPosition::idx)
                    .def_readwrite("symbol", &molgr::metal::MetalAtomPosition::symbol)
                    .def_readwrite("element_idx", &molgr::metal::MetalAtomPosition::element_idx)
                    .def_readwrite("valence", &molgr::metal::MetalAtomPosition::valence)
                    .def_readwrite("radical_num", &molgr::metal::MetalAtomPosition::radical_num)
                    .def_readwrite("position_x", &molgr::metal::MetalAtomPosition::position_x)
                    .def_readwrite("position_y", &molgr::metal::MetalAtomPosition::position_y)
                    .def_readwrite("position_z", &molgr::metal::MetalAtomPosition::position_z)
                    .def("__repr__", [](const molgr::metal::MetalAtomPosition &metal)
                         { return "<MetalAtomPosition " + metal.symbol +
                                   " val=" + std::to_string(metal.valence) +
                                   " rad=" + std::to_string(metal.radical_num) + ">"; },
                         "Return a concise debug representation of MetalAtomPosition.");

                ns.def(
                    "build_metal_states_ptr",
                    [](intptr_t mol_ptr, int atom_idx) -> std::vector<molgr::metal::MetalAtomPosition>
                    {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        auto *atom = mol->GetAtom(atom_idx);
                        if (!atom)
                        {
                            throw std::runtime_error("invalid atom index");
                        }
                        return molgr::metal::preparation::BuildMetalStates(*atom);
                    },
                    "Build candidate metal states for a specific atom index on an OBMol pointer.",
                    py::arg("mol_ptr"),
                    py::arg("atom_idx"),
                    py::call_guard<py::gil_scoped_release>());

                ns.def(
                    "combine_metal_with_omol_ptr",
                    [](intptr_t mol_ptr, const std::vector<molgr::metal::MetalAtomPosition> &metals)
                    {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        molgr::metal::preparation::CombineMetalWithOmol(*mol, metals);
                    },
                    "Combine metal states with an existing OBMol pointer in place.",
                    py::arg("mol_ptr"),
                    py::arg("metals"),
                    py::call_guard<py::gil_scoped_release>());

                ns.def(
                    "debug_metal_search_summaries",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons,
                       py::object config)
                    {
                        auto runtime_config = molgr::config::FromPython(config);
                        py::dict out;
                        const auto base_state = molgr::metal::preparation::PrepareMetalState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons,
                            runtime_config);

                        py::list available_states;
                        for (const auto &state_options : base_state.available_valence_radical_states)
                        {
                            py::list state_options_out;
                            for (const auto &metal_state : state_options)
                            {
                                state_options_out.append(metal_state_signature(metal_state));
                            }
                            available_states.append(std::move(state_options_out));
                        }
                        out["available_valence_radical_states"] = std::move(available_states);
                        out["base_phase_history"] = phase_history_signature(base_state.phase_history);

                        const auto state_search_groups =
                            molgr::metal::search::BuildMetalStateSearchGroups(
                                base_state.available_valence_radical_states,
                                runtime_config);
                        py::list state_search_groups_out;
                        for (const auto &group : state_search_groups)
                        {
                            state_search_groups_out.append(metal_state_search_group_signature(group));
                        }
                        out["state_search_groups"] = std::move(state_search_groups_out);

                        const auto layered_state_search_groups =
                            molgr::metal::search::BuildLayeredMetalStateSearchGroups(
                                state_search_groups,
                                total_radical_electrons,
                                runtime_config);
                        py::list layers_out;
                        py::list target_buckets_by_layer;
                        for (const auto &layer : layered_state_search_groups)
                        {
                            layers_out.append(metal_state_search_layer_signature(layer));

                            auto grouped_candidates =
                                molgr::metal::search::GroupCandidatesByTargetDp(
                                    base_state.phase_history,
                                    layer,
                                    total_charge,
                                    total_radical_electrons,
                                    runtime_config);
                            py::list layer_buckets;
                            for (const auto &target_entry : grouped_candidates)
                            {
                                py::dict bucket;
                                bucket["target"] = py::make_tuple(
                                    target_entry.first.no_metal_charge,
                                    target_entry.first.no_metal_radicals);
                                py::list candidates_out;
                                for (const auto &candidate : target_entry.second)
                                {
                                    py::dict candidate_out;
                                    candidate_out["combination_index"] =
                                        molgr::metal::scoring::CandidateCombinationIndex(candidate);
                                    candidate_out["no_metal_charge_target"] =
                                        candidate.no_metal_charge_target;
                                    candidate_out["no_metal_radical_target"] =
                                        candidate.no_metal_radical_target;
                                    candidate_out["metal_assignment_rank"] =
                                        molgr::metal::scoring::MetadataDouble(
                                            candidate,
                                            "metal_assignment_rank");
                                    candidate_out["metal_states"] =
                                        metal_state_choice_signature(candidate.metal_states);
                                    candidate_out["phase_history"] =
                                        phase_history_signature(candidate.phase_history);
                                    candidates_out.append(std::move(candidate_out));
                                }
                                bucket["candidates"] = std::move(candidates_out);
                                layer_buckets.append(std::move(bucket));
                            }
                            target_buckets_by_layer.append(std::move(layer_buckets));
                        }
                        out["layered_state_search_groups"] = std::move(layers_out);
                        out["target_buckets_by_layer"] = std::move(target_buckets_by_layer);
                        return out;
                    },
                    "Return C++ metal-search groups, layers, and target buckets for parity debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0,
                    py::kw_only(),
                    py::arg("config") = py::none());

                ns.def(
                    "debug_scored_candidate_summaries",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons,
                       py::object config)
                    {
                        auto runtime_config = molgr::config::FromPython(config);
                        py::dict out;
                        const auto base_state = molgr::metal::preparation::PrepareMetalState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons,
                            runtime_config);
                        if (base_state.phase_history.empty())
                        {
                            out["layer_index"] = py::none();
                            out["candidates"] = py::list();
                            return out;
                        }

                        const auto state_search_groups =
                            molgr::metal::search::BuildMetalStateSearchGroups(
                                base_state.available_valence_radical_states,
                                runtime_config);
                        const auto layered_state_search_groups =
                            molgr::metal::search::BuildLayeredMetalStateSearchGroups(
                                state_search_groups,
                                total_radical_electrons,
                                runtime_config);

                        std::vector<molgr::state::MetalCandidateState> selected_layer_candidates;
                        std::optional<std::size_t> selected_layer_index;
                        for (std::size_t layer_index = 0; layer_index < layered_state_search_groups.size();
                             ++layer_index)
                        {
                            auto grouped_candidates = molgr::metal::search::GroupCandidatesByTargetDp(
                                base_state.phase_history,
                                layered_state_search_groups[layer_index],
                                total_charge,
                                total_radical_electrons,
                                runtime_config);
                            if (grouped_candidates.empty())
                            {
                                continue;
                            }

                            std::vector<molgr::state::MetalCandidateState> current_layer_candidates;
                            for (auto &target_entry : grouped_candidates)
                            {
                                if (target_entry.second.empty())
                                {
                                    continue;
                                }
                                const auto &prototype = target_entry.second.front();
                                auto no_metal_state =
                                    molgr::pipeline::reconstruct_without_metals::XyzToOmolNoMetalState(
                                        base_state.no_metal_xyz_block,
                                        prototype.no_metal_charge_target,
                                        prototype.no_metal_radical_target,
                                        runtime_config);
                                if (!no_metal_state.has_value())
                                {
                                    continue;
                                }
                                auto shared_no_metal_state =
                                    std::make_shared<molgr::state::ReconstructionState>(
                                        std::move(*no_metal_state));
                                for (const auto &candidate : target_entry.second)
                                {
                                    current_layer_candidates.push_back(
                                        molgr::metal::scoring::PrepareCandidateWithNoMetalState(
                                            candidate,
                                            shared_no_metal_state,
                                            runtime_config));
                                }
                            }
                            if (current_layer_candidates.empty())
                            {
                                continue;
                            }
                            selected_layer_candidates = std::move(current_layer_candidates);
                            selected_layer_index = layer_index;
                            break;
                        }

                        out["layer_index"] = selected_layer_index.has_value()
                                                 ? py::object(py::int_(static_cast<int>(*selected_layer_index)))
                                                 : py::object(py::none());
                        py::list candidates_out;
                        auto selected_candidate =
                            molgr::metal::scoring::SelectBestCandidateInPlace(
                                &selected_layer_candidates,
                                runtime_config);
                        const int selected_combination_index =
                            selected_candidate.has_value()
                                ? molgr::metal::scoring::CandidateCombinationIndex(*selected_candidate)
                                : -1;
                        for (const auto &candidate : selected_layer_candidates)
                        {
                            py::dict item;
                            item["combination_index"] =
                                molgr::metal::scoring::CandidateCombinationIndex(candidate);
                            item["selected"] =
                                molgr::metal::scoring::CandidateCombinationIndex(candidate) ==
                                selected_combination_index;
                            item["no_metal_charge_target"] = candidate.no_metal_charge_target;
                            item["no_metal_radical_target"] = candidate.no_metal_radical_target;
                            if (candidate.no_metal_state)
                            {
                                item["no_metal_smiles"] = molgr::reconstruct::SmilesFirstToken(
                                    candidate.no_metal_state->Mol());
                                item["no_metal_score_key"] = molgr::scoring::BuildScoreKey(
                                    candidate.no_metal_state->Mol());
                                item["no_metal_total_charge"] = candidate.no_metal_state->total_charge;
                                item["no_metal_total_radical_electrons"] =
                                    candidate.no_metal_state->total_radical_electrons;
                                py::list atom_signature;
                                OpenBabel::OBMol &no_metal_mol =
                                    const_cast<OpenBabel::OBMol &>(candidate.no_metal_state->Mol());
                                FOR_ATOMS_OF_MOL(atom_iter, no_metal_mol)
                                {
                                    py::tuple atom_item(8);
                                    atom_item[0] = static_cast<int>(atom_iter->GetAtomicNum());
                                    atom_item[1] = atom_iter->GetFormalCharge();
                                    atom_item[2] = atom_iter->GetSpinMultiplicity();
                                    atom_item[3] = atom_iter->GetHyb();
                                    atom_item[4] = static_cast<long long>(
                                        std::llround(atom_iter->GetX() * 1000000.0));
                                    atom_item[5] = static_cast<long long>(
                                        std::llround(atom_iter->GetY() * 1000000.0));
                                    atom_item[6] = static_cast<long long>(
                                        std::llround(atom_iter->GetZ() * 1000000.0));
                                    atom_item[7] = atom_iter->IsAromatic();
                                    atom_signature.append(std::move(atom_item));
                                }
                                item["no_metal_atom_signature"] = std::move(atom_signature);

                                std::vector<std::tuple<int, int, int, bool>> bond_signature;
                                FOR_BONDS_OF_MOL(bond_iter, no_metal_mol)
                                {
                                    int begin_idx = bond_iter->GetBeginAtom()->GetIdx();
                                    int end_idx = bond_iter->GetEndAtom()->GetIdx();
                                    if (begin_idx > end_idx)
                                    {
                                        std::swap(begin_idx, end_idx);
                                    }
                                    bond_signature.emplace_back(
                                        begin_idx,
                                        end_idx,
                                        bond_iter->GetBondOrder(),
                                        bond_iter->IsAromatic());
                                }
                                std::sort(bond_signature.begin(), bond_signature.end());
                                py::list bond_signature_out;
                                for (const auto &bond_item_source : bond_signature)
                                {
                                    py::tuple bond_item(4);
                                    bond_item[0] = std::get<0>(bond_item_source);
                                    bond_item[1] = std::get<1>(bond_item_source);
                                    bond_item[2] = std::get<2>(bond_item_source);
                                    bond_item[3] = std::get<3>(bond_item_source);
                                    bond_signature_out.append(std::move(bond_item));
                                }
                                item["no_metal_bond_signature"] = std::move(bond_signature_out);
                            }
                            else
                            {
                                item["no_metal_smiles"] = py::none();
                                item["no_metal_score_key"] = py::none();
                                item["no_metal_total_charge"] = py::none();
                                item["no_metal_total_radical_electrons"] = py::none();
                                item["no_metal_atom_signature"] = py::none();
                                item["no_metal_bond_signature"] = py::none();
                            }
                            item["score"] = candidate.score.has_value()
                                                ? py::object(py::float_(*candidate.score))
                                                : py::object(py::none());
                            const auto selection_key_it = candidate.metadata.find("selection_key");
                            if (selection_key_it == candidate.metadata.end())
                            {
                                item["selection_key"] = py::none();
                            }
                            else if (const auto *selection_key =
                                         std::get_if<std::string>(&selection_key_it->second))
                            {
                                py::tuple parsed_selection_key(3);
                                std::size_t first_comma = selection_key->find(',');
                                std::size_t second_comma = first_comma == std::string::npos
                                                               ? std::string::npos
                                                               : selection_key->find(',', first_comma + 1);
                                if (first_comma != std::string::npos && second_comma != std::string::npos)
                                {
                                    parsed_selection_key[0] =
                                        std::strtod(selection_key->substr(0, first_comma).c_str(), nullptr);
                                    parsed_selection_key[1] =
                                        std::strtod(
                                            selection_key
                                                ->substr(first_comma + 1, second_comma - first_comma - 1)
                                                .c_str(),
                                            nullptr);
                                    parsed_selection_key[2] =
                                        std::atoi(selection_key->substr(second_comma + 1).c_str());
                                    item["selection_key"] = std::move(parsed_selection_key);
                                }
                                else
                                {
                                    item["selection_key"] = *selection_key;
                                }
                            }
                            else
                            {
                                item["selection_key"] = py::none();
                            }
                            item["metal_assignment_rank"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "metal_assignment_rank");
                            item["organic_aromatic_atom_count"] =
                                molgr::metal::scoring::MetadataInt(candidate, "organic_aromatic_atom_count");
                            item["organic_aromatic_ring_count"] =
                                molgr::metal::scoring::MetadataInt(candidate, "organic_aromatic_ring_count");
                            item["organic_aromatic_stability_score"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "organic_aromatic_stability_score");
                            item["organic_conjugated_atom_count"] =
                                molgr::metal::scoring::MetadataInt(candidate, "organic_conjugated_atom_count");
                            item["organic_max_conjugated_component_size"] =
                                molgr::metal::scoring::MetadataInt(candidate, "organic_max_conjugated_component_size");
                            item["organic_charge_localization_penalty"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "organic_charge_localization_penalty");
                            item["organic_radical_localization_penalty"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "organic_radical_localization_penalty");
                            item["metal_discordance_structural_count"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "metal_discordance_structural_count");
                            item["metal_discordance_count"] =
                                molgr::metal::scoring::MetadataDouble(candidate, "metal_discordance_count");
                            item["metal_discordance_aromatic_ring_deficit_count"] =
                                molgr::metal::scoring::MetadataInt(
                                    candidate,
                                    "metal_discordance_aromatic_ring_deficit_count");
                            item["metal_discordance_aromatic_stability_deficit"] =
                                molgr::metal::scoring::MetadataDouble(
                                    candidate,
                                    "metal_discordance_aromatic_stability_deficit");
                            const auto passes_filter_it =
                                candidate.metadata.find("passes_metal_discordance_filter");
                            if (passes_filter_it == candidate.metadata.end())
                            {
                                item["passes_metal_discordance_filter"] = py::none();
                            }
                            else if (const auto *passes_filter =
                                         std::get_if<bool>(&passes_filter_it->second))
                            {
                                item["passes_metal_discordance_filter"] = *passes_filter;
                            }
                            else
                            {
                                item["passes_metal_discordance_filter"] = py::none();
                            }
                            py::list metal_states;
                            for (const auto &metal_state : candidate.metal_states)
                            {
                                py::dict metal_item;
                                metal_item["idx"] = metal_state.idx;
                                metal_item["symbol"] = metal_state.symbol;
                                metal_item["valence"] = metal_state.valence;
                                metal_item["radical_num"] = metal_state.radical_num;
                                metal_states.append(std::move(metal_item));
                            }
                            item["metal_states"] = std::move(metal_states);
                            candidates_out.append(std::move(item));
                        }
                        out["candidates"] = std::move(candidates_out);
                        return out;
                    },
                    "Return scored metal candidates for the first successful search layer.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0,
                    py::kw_only(),
                    py::arg("config") = py::none());
            }

            void bind_reconstruct_without_metals_dev_ns(py::module_ &ns)
            {
                ns.def(
                    "debug_linear_pipeline_state",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons)
                        -> py::object
                    {
                        auto seed_state = molgr::no_metals::preparation::SeedState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons);
                        if (!seed_state.omol)
                        {
                            return py::none();
                        }
                        auto state = molgr::no_metals::preparation::RunLinearPipeline(seed_state);
                        py::dict item;
                        item["smiles"] = molgr::reconstruct::SmilesFirstToken(state.Mol());
                        item["given_charge"] = state.given_charge;
                        item["total_charge"] = state.total_charge;
                        item["total_radical_electrons"] = state.total_radical_electrons;
                        item["valid"] = molgr::reconstruct::ValidateOmol(
                            state.MutableMol(),
                            state.total_charge,
                            state.total_radical_electrons);
                        py::list phases;
                        for (const auto &phase : state.phase_history)
                        {
                            phases.append(phase);
                        }
                        item["phase_history"] = std::move(phases);
                        return item;
                    },
                    "Return the C++ no-metal linear-pipeline state for parity debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0);

                ns.def(
                    "debug_no_metal_candidate_states",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons)
                        -> py::object
                    {
                        auto seed_state = molgr::no_metals::preparation::SeedState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons);
                        if (!seed_state.omol)
                        {
                            return py::none();
                        }
                        auto states = molgr::no_metals::preparation::EnumerateNoMetalCandidateStates(
                            seed_state);
                        py::list out;
                        for (auto &state : states)
                        {
                            py::dict item;
                            item["smiles"] = molgr::reconstruct::SmilesFirstToken(state.Mol());
                            item["given_charge"] = state.given_charge;
                            item["total_charge"] = state.total_charge;
                            item["total_radical_electrons"] = state.total_radical_electrons;
                            item["valid"] = molgr::reconstruct::ValidateOmol(
                                state.MutableMol(),
                                state.total_charge,
                                state.total_radical_electrons);
                            py::list phases;
                            for (const auto &phase : state.phase_history)
                            {
                                phases.append(phase);
                            }
                            item["phase_history"] = std::move(phases);
                            const auto strategy_it =
                                state.metadata.find("neighbor_radical_resolution_strategy");
                            if (strategy_it != state.metadata.end())
                            {
                                if (const auto *strategy =
                                        std::get_if<std::string>(&strategy_it->second))
                                {
                                    item["neighbor_radical_resolution_strategy"] = *strategy;
                                }
                            }
                            out.append(std::move(item));
                        }
                        return out;
                    },
                    "Return C++ no-metal candidate states for parity debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0);

                ns.def(
                    "debug_linear_pipeline_trace",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons) -> py::object
                    {
                        auto seed_state = molgr::no_metals::preparation::SeedState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons);
                        if (!seed_state.omol)
                        {
                            return py::none();
                        }

                        auto machine = molgr::state::OmolStateMachine::FromReconstructionState(seed_state);
                        py::list trace;
                        const auto append_snapshot = [&](const std::string &phase)
                        {
                            py::dict item;
                            item["phase"] = phase;
                            item["smiles"] = molgr::reconstruct::SmilesFirstToken(machine.EnsureUniqueMol());
                            item["score_key"] = molgr::scoring::BuildScoreKey(machine.EnsureUniqueMol());
                            item["given_charge"] = machine.given_charge;
                            item["valid"] = molgr::reconstruct::ValidateOmol(
                                machine.EnsureUniqueMol(),
                                total_charge,
                                total_radical_electrons);
                            trace.append(std::move(item));
                        };

                        append_snapshot("read_xyz");
                        machine.RunOmolStage("make_connections", reconstruct::MakeConnections, 0.15);
                        append_snapshot("make_connections");
                        machine.RunOmolStage("pre_clean", reconstruct::PreClean);
                        append_snapshot("pre_clean");
                        machine.RunOmolStage(
                            "fresh_omol_charge_radical_initial",
                            reconstruct::FreshOmolChargeRadical);
                        append_snapshot("fresh_omol_charge_radical_initial");

                        int formal_charge_sum = 0;
                        FOR_ATOMS_OF_MOL(atom_iter, machine.EnsureUniqueMol())
                        {
                            formal_charge_sum += atom_iter->GetFormalCharge();
                        }
                        machine.SetGivenCharge(
                            "initialize_charge_budget",
                            total_charge - formal_charge_sum);
                        append_snapshot("initialize_charge_budget");

                        machine.RunOmolChargeStage("eliminate_NNN_negative", reconstruct::EliminateNNN, false);
                        append_snapshot("eliminate_NNN_negative");
                        machine.RunOmolChargeStage(
                            "eliminate_high_positive_charge_atoms",
                            reconstruct::EliminateHighPositiveChargeAtoms);
                        append_snapshot("eliminate_high_positive_charge_atoms");
                        machine.RunOmolChargeStage(
                            "eliminate_CN_in_doubt",
                            reconstruct::EliminateCNInDoubt);
                        append_snapshot("eliminate_CN_in_doubt");
                        machine.RunOmolChargeStage("eliminate_NNN_positive", reconstruct::EliminateNNN, true);
                        append_snapshot("eliminate_NNN_positive");
                        machine.RunOmolChargeStage("eliminate_carboxyl", reconstruct::EliminateCarboxyl);
                        append_snapshot("eliminate_carboxyl");
                        machine.RunOmolStage(
                            "clean_carbene_neighbor_unsaturated_first",
                            reconstruct::CleanCarbeneNeighborUnsaturated);
                        append_snapshot("clean_carbene_neighbor_unsaturated_first");
                        machine.RunOmolChargeStage(
                            "eliminate_carbene_neighbor_heteroatom",
                            reconstruct::EliminateCarbeneNeighborHeteroatom);
                        append_snapshot("eliminate_carbene_neighbor_heteroatom");
                        machine.RunOmolStage("clean_neighbor_radicals", reconstruct::CleanNeighborRadicals);
                        append_snapshot("clean_neighbor_radicals");
                        machine.RunOmolStage(
                            "clean_carbene_neighbor_unsaturated_second",
                            reconstruct::CleanCarbeneNeighborUnsaturated);
                        append_snapshot("clean_carbene_neighbor_unsaturated_second");
                        machine.RunOmolChargeStage(
                            "eliminate_charge_spliting",
                            reconstruct::EliminateChargeSpliting);
                        append_snapshot("eliminate_charge_spliting");
                        machine.RunOmolStage(
                            "break_deformed_ene",
                            reconstruct::BreakDeformedEne,
                            machine.given_charge,
                            total_radical_electrons,
                            5.0);
                        append_snapshot("break_deformed_ene");
                        machine.RunOmolChargeStage(
                            "break_one_bond",
                            reconstruct::BreakOneBond,
                            total_radical_electrons);
                        append_snapshot("break_one_bond");
                        machine.RunOmolStage(
                            "fresh_omol_charge_radical_final",
                            reconstruct::FreshOmolChargeRadical);
                        append_snapshot("fresh_omol_charge_radical_final");
                        return std::move(trace);
                    },
                    "Return per-stage C++ no-metal linear-pipeline snapshots for parity debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0);

                ns.def(
                    "debug_resonance_candidate_summaries",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons,
                       py::object config)
                    {
                        auto runtime_config = molgr::config::FromPython(config);
                        std::vector<molgr::pipeline::reconstruct_without_metals::DebugNoMetalCandidateSummary>
                            summaries;
                        {
                            py::gil_scoped_release release;
                            summaries = molgr::pipeline::reconstruct_without_metals::
                                DebugNoMetalResonanceCandidateSummaries(
                                    xyz_block,
                                    total_charge,
                                    total_radical_electrons,
                                    runtime_config);
                        }
                        py::list out;
                        for (const auto &summary : summaries)
                        {
                            py::dict item;
                            item["smiles"] = summary.smiles;
                            item["resonance_index"] = summary.resonance_index;
                            item["score"] = summary.score;
                            item["aromatic_stability_score"] = summary.aromatic_stability_score;
                            item["aromatic_atom_count"] = summary.aromatic_atom_count;
                            item["max_conjugated_component_size"] =
                                summary.max_conjugated_component_size;
                            item["conjugated_atom_count"] = summary.conjugated_atom_count;
                            item["conjugated_bond_count"] = summary.conjugated_bond_count;
                            item["formal_charge_absolute_sum"] = summary.formal_charge_absolute_sum;
                            item["conjugation_charge_penalty"] = summary.conjugation_charge_penalty;
                            item["adjusted_max_conjugated_component_size"] =
                                summary.adjusted_max_conjugated_component_size;
                            item["adjusted_conjugated_atom_count"] =
                                summary.adjusted_conjugated_atom_count;
                            item["adjusted_conjugated_bond_count"] =
                                summary.adjusted_conjugated_bond_count;
                            out.append(std::move(item));
                        }
                        return out;
                    },
                    "Return C++ no-metal resonance candidates for parity debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0,
                    py::kw_only(),
                    py::arg("config") = py::none());

                ns.def(
                    "debug_processed_root_resonance",
                    [](const std::string &xyz_block,
                       int total_charge,
                       int total_radical_electrons)
                        -> py::object
                    {
                        auto seed_state = molgr::no_metals::preparation::SeedState(
                            xyz_block,
                            total_charge,
                            total_radical_electrons);
                        if (!seed_state.omol)
                        {
                            return py::none();
                        }
                        auto state = molgr::no_metals::preparation::RunLinearPipeline(seed_state);
                        auto processed = molgr::resonance::ProcessResonanceDetailed(
                            state.Mol(),
                            state.given_charge);
                        py::dict item;
                        item["smiles"] = molgr::reconstruct::SmilesFirstToken(std::get<0>(processed));
                        item["given_charge"] = std::get<1>(processed);
                        item["hit"] = std::get<2>(processed);
                        auto &processed_mol = std::get<0>(processed);
                        item["valid"] = molgr::reconstruct::ValidateOmol(
                            processed_mol,
                            state.total_charge,
                            state.total_radical_electrons);
                        item["processed_key"] =
                            molgr::resonance::BuildProcessedResonanceKey(processed_mol);
                        return item;
                    },
                    "Process the linear no-metal state as a root resonance candidate for debugging.",
                    py::arg("xyz_block"),
                    py::arg("total_charge") = 0,
                    py::arg("total_radical_electrons") = 0);
            }

            void bind_resonance_dev_ns(py::module_ &ns)
            {
                ns.def(
                    "get_radical_resonances_smi",
                    [](const std::string &smiles) -> std::vector<std::string>
                    {
                        auto mol = mol_from_smiles(smiles);
                        if (!mol)
                        {
                            throw std::runtime_error("failed to parse SMILES");
                        }
                        const auto resonances = molgr::reconstruct::GetRadicalResonances(*mol);
                        std::vector<std::string> out;
                        out.reserve(resonances.size());
                        for (const auto &resonance : resonances)
                        {
                            out.push_back(molgr::reconstruct::SmilesFirstToken(resonance));
                        }
                        return out;
                    },
                    "Enumerate radical resonance structures from a SMILES string.",
                    py::arg("smiles"),
                    py::call_guard<py::gil_scoped_release>());

                ns.def(
                    "process_resonance_smi",
                    [](const std::string &smiles, int charge) -> std::tuple<std::string, int>
                    {
                        auto mol = mol_from_smiles(smiles);
                        if (!mol)
                        {
                            throw std::runtime_error("failed to parse SMILES");
                        }
                        auto processed = molgr::reconstruct::ProcessResonance(*mol, charge);
                        return std::make_tuple(
                            molgr::reconstruct::SmilesFirstToken(processed.first),
                            processed.second);
                    },
                    "Process one resonance step on SMILES and return updated token plus charge.",
                    py::arg("smiles"),
                    py::arg("charge"),
                    py::call_guard<py::gil_scoped_release>());

                ns.def(
                    "get_radical_resonances_ptr",
                    [](intptr_t mol_ptr) -> std::vector<intptr_t>
                    {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        py::gil_scoped_release release;
                        const auto resonances = molgr::reconstruct::GetRadicalResonances(*mol);
                        std::vector<intptr_t> out;
                        out.reserve(resonances.size());
                        for (const auto &resonance : resonances)
                        {
                            auto *res_ptr = new OpenBabel::OBMol(resonance);
                            out.push_back(reinterpret_cast<intptr_t>(res_ptr));
                        }
                        return out;
                    },
                    "Get radical resonance OBMol pointers for an input OBMol pointer.",
                    py::arg("mol_ptr"));

                ns.def(
                    "process_resonance_ptr",
                    [](intptr_t mol_ptr, int charge) -> py::tuple
                    {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        std::pair<OpenBabel::OBMol, int> processed;
                        {
                            py::gil_scoped_release release;
                            processed = molgr::reconstruct::ProcessResonance(*mol, charge);
                        }
                        auto *res_ptr = new OpenBabel::OBMol(processed.first);
                        return py::make_tuple(reinterpret_cast<intptr_t>(res_ptr), processed.second);
                    },
                    "Process resonance on an OBMol pointer and return (new_ptr, updated_charge).",
                    py::arg("mol_ptr"),
                    py::arg("charge"));

                ns.def(
                    "smiles_token_ptr",
                    [](intptr_t mol_ptr) -> std::string
                    {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        py::gil_scoped_release release;
                        return molgr::reconstruct::SmilesFirstToken(*mol);
                    },
                    "Return canonical first-token SMILES string for an OBMol pointer.",
                    py::arg("mol_ptr"));
            }
        }

        void bind_pipeline(py::module_ &m)
        {
            auto reconstruct_with_metals = m.def_submodule(
                "reconstruct_with_metals",
                "Fallback-aligned reconstruction helpers with metals");

            auto reconstruct_without_metals = m.def_submodule(
                "reconstruct_without_metals",
                "Fallback-aligned no-metal reconstruction helpers");

            bind_reconstruct_with_metals_public_ns(reconstruct_with_metals);

            m.def(
                "xyz2omol",
                [](const std::string &xyz_block,
                   int total_charge,
                   int total_radical_electrons,
                   py::object config)
                    -> std::unique_ptr<molgr::utils::MoleculeData>
                {
                    auto runtime_config = molgr::config::FromPython(config);
                    py::gil_scoped_release release;
                    auto mol_data = molgr::pipeline::reconstruct_with_metals::Xyz2OmolMolData(
                        xyz_block,
                        total_charge,
                        total_radical_electrons,
                        runtime_config);
                    if (!mol_data)
                    {
                        return nullptr;
                    }
                    return std::move(mol_data);
                },
                "Reconstruct molecule data from XYZ with metal-aware pipeline.",
                py::arg("xyz_block"),
                py::arg("total_charge") = 0,
                py::arg("total_radical_electrons") = 0,
                py::kw_only(),
                py::arg("config") = py::none());

            reconstruct_without_metals.def(
                "xyz_to_omol_no_metal",
                [](const std::string &xyz_block,
                   int total_charge,
                   int total_radical_electrons,
                   py::object config)
                    -> std::unique_ptr<molgr::utils::MoleculeData>
                {
                    auto runtime_config = molgr::config::FromPython(config);
                    py::gil_scoped_release release;
                    auto mol_data = molgr::pipeline::reconstruct_without_metals::XyzToMolDataNoMetal(
                        xyz_block,
                        total_charge,
                        total_radical_electrons,
                        runtime_config);
                    if (!mol_data)
                    {
                        return nullptr;
                    }
                    return std::move(mol_data);
                },
                "Reconstruct molecule data from XYZ without metal handling.",
                py::arg("xyz_block"),
                py::arg("total_charge") = 0,
                py::arg("total_radical_electrons") = 0,
                py::kw_only(),
                py::arg("config") = py::none());

            m.def(
                "get_last_run_timing_breakdown_ms",
                []() -> py::dict
                {
                    const auto timing = molgr::pipeline::perf::GetRunTimingBreakdown();
                    py::dict out;
                    out["no_metal_pipeline_ms"] = timing.no_metal_pipeline_ms;
                    out["no_metal_linear_pipeline_ms"] = timing.no_metal_linear_pipeline_ms;
                    out["no_metal_validate_ms"] = timing.no_metal_validate_ms;
                    out["resonance_handling_enumeration_ms"] = timing.resonance_handling_enumeration_ms;
                    out["resonance_walk_ms"] = timing.resonance_walk_ms;
                    out["resonance_prepare_ms"] = timing.resonance_prepare_ms;
                    out["resonance_dedup_score_ms"] = timing.resonance_dedup_score_ms;
                    out["resonance_score_ms"] = timing.resonance_score_ms;
                    out["resonance_topology_ms"] = timing.resonance_topology_ms;
                    out["resonance_raw_candidates"] = timing.resonance_raw_candidates;
                    out["resonance_prepared_candidates"] = timing.resonance_prepared_candidates;
                    out["resonance_valid_candidates"] = timing.resonance_valid_candidates;
                    out["resonance_dedup_candidates"] = timing.resonance_dedup_candidates;
                    out["metal_enumeration_combination_ms"] = timing.metal_enumeration_combination_ms;
                    out["force_field_total_ms"] = timing.force_field_total_ms;
                    out["force_field_cache_key_ms"] = timing.force_field_cache_key_ms;
                    out["force_field_prepare_ms"] = timing.force_field_prepare_ms;
                    out["force_field_setup_key_ms"] = timing.force_field_setup_key_ms;
                    out["force_field_setup_ms"] = timing.force_field_setup_ms;
                    out["force_field_energy_ms"] = timing.force_field_energy_ms;
                    out["force_field_calls"] = timing.force_field_calls;
                    return out;
                },
                "Return timing breakdown (milliseconds) for the most recent reconstruction run.");
        }

        void bind_dev_pipeline(py::module_ &m)
        {
            auto reconstruct_with_metals = m.def_submodule(
                "reconstruct_with_metals",
                "Development-only helpers for metal reconstruction internals");
            auto reconstruct_without_metals = m.def_submodule(
                "reconstruct_without_metals",
                "Development-only helpers for no-metal reconstruction internals");
            auto resonance = m.def_submodule(
                "resonance",
                "Development-only helpers for resonance internals");

            bind_reconstruct_with_metals_dev_ns(reconstruct_with_metals);
            bind_reconstruct_without_metals_dev_ns(reconstruct_without_metals);
            bind_resonance_dev_ns(resonance);

            m.def(
                "clear_force_field_evaluation_cache",
                []()
                {
                    molgr::scoring::ForceFieldEvaluationCacheClear();
                },
                "Clear the C++ force-field evaluation cache.");
            m.def(
                "clear_uff_atom_typing_cache",
                []()
                {
                    OpenBabel::ClearMolgrUffAtomTypeAssignmentCache();
                },
                "Clear the C++ UFF atom-typing assignment cache.");
            m.def(
                "get_uff_atom_typing_cache_info",
                []()
                {
                    const auto [hits, misses, size] =
                        OpenBabel::MolgrUffAtomTypeAssignmentCacheInfo();
                    py::dict out;
                    out["hits"] = hits;
                    out["misses"] = misses;
                    out["size"] = size;
                    return out;
                },
                "Return C++ UFF atom-typing cache hit/miss/size counters.");
        }
    }
}
