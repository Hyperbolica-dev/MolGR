#include "bindings.h"

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/pipeline/reconstruct_without_metals.h"
#include "molgr/pipeline/resonance.h"
#include "molgr/utils/utils.h"

#include <openbabel/obconversion.h>

#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace bind
    {

        static OpenBabel::OBMol *require_obmol_ptr(intptr_t mol_ptr)
        {
            if (mol_ptr == 0)
            {
                throw std::runtime_error("null OBMol pointer");
            }
            return reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
        }

        static std::unique_ptr<OpenBabel::OBMol> mol_from_smiles(const std::string &smiles)
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

        struct ReconstructWithMetalsNS
        {
        };

        struct ReconstructWithoutMetalsNS
        {
        };

        struct ResonanceNS
        {
        };

        static void bind_reconstruct_with_metals_ns(py::class_<ReconstructWithMetalsNS> &ns)
        {
            py::class_<molgr::metal::MetalAtomPosition>(ns, "MetalAtomPosition")
                .def_readwrite("idx", &molgr::metal::MetalAtomPosition::idx)
                .def_readwrite("symbol", &molgr::metal::MetalAtomPosition::symbol)
                .def_readwrite("element_idx", &molgr::metal::MetalAtomPosition::element_idx)
                .def_readwrite("valence", &molgr::metal::MetalAtomPosition::valence)
                .def_readwrite("radical_num", &molgr::metal::MetalAtomPosition::radical_num)
                .def_readwrite("position_x", &molgr::metal::MetalAtomPosition::x)
                .def_readwrite("position_y", &molgr::metal::MetalAtomPosition::y)
                .def_readwrite("position_z", &molgr::metal::MetalAtomPosition::z)
                .def("__repr__", [](const molgr::metal::MetalAtomPosition &metal)
                     { return "<MetalAtomPosition " + metal.symbol +
                              " val=" + std::to_string(metal.valence) +
                              " rad=" + std::to_string(metal.radical_num) + ">"; });

            py::class_<molgr::metal::MetalHandler>(ns, "MetalHandler")
                .def(py::init([](intptr_t mol_ptr)
                              {
                          auto *mol = require_obmol_ptr(mol_ptr);
                          return new molgr::metal::MetalHandler(*mol); }),
                     py::arg("mol_ptr"))
                .def("strip_metals", [](molgr::metal::MetalHandler &self, intptr_t mol_ptr)
                     {
                 auto *mol = require_obmol_ptr(mol_ptr);
                 return self.StripMetals(*mol); }, py::arg("mol_ptr"))
                .def("generate_combinations", &molgr::metal::MetalHandler::GenerateCombinations, py::arg("total_radical_electrons"))
                .def_static("combine_metal_with_mol", [](intptr_t mol_ptr, const std::vector<molgr::metal::MetalAtomPosition> &metals)
                            {
                        auto *mol = require_obmol_ptr(mol_ptr);
                        molgr::metal::MetalHandler::CombineMetalWithMol(*mol, metals); }, py::arg("mol_ptr"), py::arg("metals"));

            ns.def_static(
                "get_possible_metal_radicals",
                &molgr::pipeline::reconstruct_with_metals::get_possible_metal_radicals,
                py::arg("metal"),
                py::arg("valence"),
                py::call_guard<py::gil_scoped_release>());

            ns.def_static(
                "build_metal_states_ptr",
                [](intptr_t mol_ptr, int atom_idx) -> std::vector<molgr::metal::MetalAtomPosition>
                {
                    auto *mol = require_obmol_ptr(mol_ptr);
                    auto *atom = mol->GetAtom(atom_idx);
                    if (!atom)
                    {
                        throw std::runtime_error("invalid atom index");
                    }
                    return molgr::pipeline::reconstruct_with_metals::build_metal_states(*atom);
                },
                py::arg("mol_ptr"),
                py::arg("atom_idx"),
                py::call_guard<py::gil_scoped_release>());

            ns.def_static(
                "combine_metal_with_omol_ptr",
                [](intptr_t mol_ptr, const std::vector<molgr::metal::MetalAtomPosition> &metals)
                {
                    auto *mol = require_obmol_ptr(mol_ptr);
                    molgr::pipeline::reconstruct_with_metals::combine_metal_with_omol(*mol, metals);
                },
                py::arg("mol_ptr"),
                py::arg("metals"),
                py::call_guard<py::gil_scoped_release>());
        }

        void bind_pipeline(py::module_ &m)
        {
            auto reconstruct_with_metals = py::class_<ReconstructWithMetalsNS>(
                m,
                "reconstruct_with_metals",
                "Fallback-aligned reconstruction helpers with metals");

            auto reconstruct_without_metals = py::class_<ReconstructWithoutMetalsNS>(
                m,
                "reconstruct_without_metals",
                "Fallback-aligned no-metal reconstruction helpers");

            auto resonance = py::class_<ResonanceNS>(
                m,
                "resonance",
                "Fallback-aligned resonance helpers");

            bind_reconstruct_with_metals_ns(reconstruct_with_metals);

            reconstruct_with_metals.def_static(
                "xyz2omol",
                [](const std::string &xyz_block, int total_charge, int total_radical_electrons)
                    -> std::unique_ptr<molgr::utils::MoleculeData>
                {
                    auto mol_data = molgr::pipeline::reconstruct_with_metals::Xyz2OmolMolData(
                        xyz_block,
                        total_charge,
                        total_radical_electrons);
                    if (!mol_data)
                    {
                        return nullptr;
                    }
                    return std::move(mol_data);
                },
                py::arg("xyz_block"),
                py::arg("total_charge") = 0,
                py::arg("total_radical_electrons") = 0,
                py::call_guard<py::gil_scoped_release>());

            reconstruct_without_metals.def_static(
                "xyz_to_omol_no_metal",
                [](const std::string &xyz_block, int total_charge, int total_radical_electrons)
                    -> std::unique_ptr<molgr::utils::MoleculeData>
                {
                    auto mol_data = molgr::pipeline::reconstruct_without_metals::XyzToMolDataNoMetal(
                        xyz_block,
                        total_charge,
                        total_radical_electrons);
                    if (!mol_data)
                    {
                        return nullptr;
                    }
                    return std::move(mol_data);
                },
                py::arg("xyz_block"),
                py::arg("total_charge") = 0,
                py::arg("total_radical_electrons") = 0,
                py::call_guard<py::gil_scoped_release>());

            resonance.def_static(
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
                py::arg("smiles"),
                py::call_guard<py::gil_scoped_release>());

            resonance.def_static(
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
                py::arg("smiles"),
                py::arg("charge"),
                py::call_guard<py::gil_scoped_release>());

            m.def(
                "get_last_run_timing_breakdown_ms",
                []() -> py::dict
                {
                    const auto timing = molgr::pipeline::perf::GetRunTimingBreakdown();
                    py::dict out;
                    out["no_metal_pipeline_ms"] = timing.no_metal_pipeline_ms;
                    out["resonance_handling_enumeration_ms"] = timing.resonance_handling_enumeration_ms;
                    out["metal_enumeration_combination_ms"] = timing.metal_enumeration_combination_ms;
                    return out;
                });

            m.def(
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
                py::arg("mol_ptr"));

            m.def(
                "process_resonance_ptr",
                [](intptr_t mol_ptr, int charge) -> py::tuple
                {
                    auto *mol = require_obmol_ptr(mol_ptr);
                    py::gil_scoped_release release;
                    auto processed = molgr::reconstruct::ProcessResonance(*mol, charge);
                    auto *res_ptr = new OpenBabel::OBMol(processed.first);
                    return py::make_tuple(reinterpret_cast<intptr_t>(res_ptr), processed.second);
                },
                py::arg("mol_ptr"),
                py::arg("charge"));

            m.def(
                "smiles_token_ptr",
                [](intptr_t mol_ptr) -> std::string
                {
                    auto *mol = require_obmol_ptr(mol_ptr);
                    py::gil_scoped_release release;
                    return molgr::reconstruct::SmilesFirstToken(*mol);
                },
                py::arg("mol_ptr"));
        }

    }
}
