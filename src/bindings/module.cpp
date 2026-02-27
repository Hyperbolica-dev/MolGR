/*
 * @Author: TMJ
 * @Date: 2026-02-23 00:02:13
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-27 18:02:24
 * @Description: 请填写简介
 */
/**
 * @file module.cpp
 * @brief Pybind11 bindings for the molgr C++ library.
 * @author TMJ
 * @date 2026-02-27
 */

#include <pybind11/pybind11.h>
#include <cstdint> // 用于 intptr_t

#include "molgr/utils/logger.h" // 日志模块
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

    auto m_utils = m.def_submodule("utils", "Utilities and Data Structures");
    molgr::bind::bind_utils(m_utils);

    auto m_pipeline = m.def_submodule("pipeline", "Pipeline-level helpers");
    molgr::bind::bind_pipeline(m_pipeline);

    auto m_stages = m.def_submodule("stages", "Stage-level helpers for parity testing");
    molgr::bind::bind_stages(m_stages);
}
