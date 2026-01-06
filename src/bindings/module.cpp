/**
 * @file module.cpp
 * @brief Pybind11 bindings for the molgr C++ library.
 * @author TMJ
 * @date 2025-12-25
 */

#pragma once
#include <pybind11/pybind11.h>
#include <cstdint> // 用于 intptr_t

#include "molgr/logger.h" // 日志模块
#include "bindings.h"
#include <openbabel/mol.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        MolGR Core C++ Implementation
        -----------------------------
        Exposes optimized C++ algorithms for MolOP Graph Reconstruction.
    )pbdoc";

    py::enum_<molgr::LogLevel>(m, "LogLevel")
        .value("DEBUG", molgr::LogLevel::DEBUG)
        .value("INFO", molgr::LogLevel::INFO)
        .value("WARN", molgr::LogLevel::WARN)
        .value("ERROR", molgr::LogLevel::ERROR)
        .value("OFF", molgr::LogLevel::OFF)
        .export_values();

    m.def("set_log_level", &molgr::SetLogLevel,
          "Set the logging level for the C++ core (DEBUG=0, INFO=1, WARN=2, ERROR=3, OFF=4)",
          py::arg("level"));

    m.def("free_obmol_ptr", [](intptr_t mol_ptr)
          {
        if (mol_ptr != 0) {
            delete reinterpret_cast<OpenBabel::OBMol*>(mol_ptr);
        } }, "Manually delete the OBMol pointer");

    // 1. Consts 子模块 -> _core.consts
    auto m_consts = m.def_submodule("consts", "Chemical constants");
    molgr::bind::bind_consts(m_consts);

    // 2. Utils 子模块 -> _core.utils
    // 包含 AtomData, MoleculeData 等结构体
    auto m_utils = m.def_submodule("utils", "Utilities and Data Structures");
    molgr::bind::bind_utils(m_utils);

    // 3. Metal 子模块 -> _core.metal
    auto m_metal = m.def_submodule("metal", "Metal handling logic");
    molgr::bind::bind_metal(m_metal);

    // 4. Scoring 子模块 -> _core.scoring
    auto m_scoring = m.def_submodule("scoring", "Scoring functions");
    molgr::bind::bind_scoring(m_scoring);

    // 5. Reconstruct 子模块 -> _core.reconstruct
    auto m_recon = m.def_submodule("reconstruct", "Graph reconstruction algorithms");
    molgr::bind::bind_reconstruct(m_recon);
}