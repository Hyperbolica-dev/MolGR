#include "bindings.h"

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/pipeline/reconstruct_without_metals.h"
#include "molgr/pipeline/resonance.h"

#include <openbabel/obconversion.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
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

        static void bind_reconstruct_with_metals(py::module_ &m)
        {
            py::class_<molgr::metal::MetalAtomPosition>(m, "MetalAtomPosition")
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

            py::class_<molgr::metal::MetalHandler>(m, "MetalHandler")
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

            m.def(
                "get_possible_metal_radicals",
                &molgr::pipeline::reconstruct_with_metals::get_possible_metal_radicals,
                py::arg("metal"),
                py::arg("valence"));

            m.def(
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
                py::arg("atom_idx"));

            m.def(
                "combine_metal_with_omol_ptr",
                [](intptr_t mol_ptr, const std::vector<molgr::metal::MetalAtomPosition> &metals)
                {
                    auto *mol = require_obmol_ptr(mol_ptr);
                    molgr::pipeline::reconstruct_with_metals::combine_metal_with_omol(*mol, metals);
                },
                py::arg("mol_ptr"),
                py::arg("metals"));
        }

        void bind_pipeline(py::module_ &m)
        {
            const auto xyz2omol = [](const std::string &xyz_block, int total_charge, int total_radical_electrons) -> py::object
            {
                py::gil_scoped_release release;
                auto mol_data = molgr::pipeline::reconstruct_with_metals::Xyz2OmolMolData(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                py::gil_scoped_acquire acquire;
                if (!mol_data)
                {
                    return py::none();
                }
                return py::cast(*mol_data);
            };

            const auto xyz_to_omol_no_metal = [](const std::string &xyz_block, int total_charge, int total_radical_electrons) -> py::object
            {
                py::gil_scoped_release release;
                auto mol_data = molgr::pipeline::reconstruct_without_metals::XyzToMolDataNoMetal(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                py::gil_scoped_acquire acquire;
                if (!mol_data)
                {
                    return py::none();
                }
                return py::cast(*mol_data);
            };

            const auto get_radical_resonances_smi = [](const std::string &smiles) -> std::vector<std::string>
            {
                auto mol = mol_from_smiles(smiles);
                if (!mol)
                {
                    throw std::runtime_error("failed to parse SMILES");
                }
                py::gil_scoped_release release;
                const auto resonances = molgr::reconstruct::GetRadicalResonances(*mol);
                std::vector<std::string> out;
                out.reserve(resonances.size());
                for (const auto &resonance : resonances)
                {
                    out.push_back(molgr::reconstruct::SmilesFirstToken(resonance));
                }
                return out;
            };

            const auto process_resonance_smi = [](const std::string &smiles, int charge) -> py::tuple
            {
                auto mol = mol_from_smiles(smiles);
                if (!mol)
                {
                    throw std::runtime_error("failed to parse SMILES");
                }
                py::gil_scoped_release release;
                auto processed = molgr::reconstruct::ProcessResonance(*mol, charge);
                return py::make_tuple(
                    molgr::reconstruct::SmilesFirstToken(processed.first),
                    processed.second);
            };

            auto m_reconstruct_without_metals = m.def_submodule(
                "reconstruct_without_metals",
                "Fallback-aligned no-metal reconstruction helpers");
            auto m_reconstruct_with_metals = m.def_submodule(
                "reconstruct_with_metals",
                "Fallback-aligned reconstruction helpers with metals");
            auto m_resonance = m.def_submodule(
                "resonance",
                "Fallback-aligned resonance helpers");

            bind_reconstruct_with_metals(m_reconstruct_with_metals);

            m_reconstruct_with_metals.def(
                "xyz2omol",
                xyz2omol,
                py::arg("xyz_block"),
                py::arg("total_charge"),
                py::arg("total_radical_electrons"));

            m_reconstruct_without_metals.def(
                "xyz_to_omol_no_metal",
                xyz_to_omol_no_metal,
                py::arg("xyz_block"),
                py::arg("total_charge"),
                py::arg("total_radical_electrons"));

            m_resonance.def(
                "get_radical_resonances_smi",
                get_radical_resonances_smi,
                py::arg("smiles"));

            m_resonance.def(
                "process_resonance_smi",
                process_resonance_smi,
                py::arg("smiles"),
                py::arg("charge"));

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
