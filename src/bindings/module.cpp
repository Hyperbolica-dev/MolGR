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

#include "molgr/python_config.h"
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

    m.def("set_default_config", &molgr::config::SetDefaultConfigFromPython,
          "Set the default runtime configuration used by the C++ backend.",
          py::arg("config"));

    m.def("get_default_config", &molgr::config::DefaultConfigSummary,
          "Return a summary of the default runtime configuration used by the C++ backend.");

    m.def("free_obmol_ptr", [](intptr_t mol_ptr)
          {
        if (mol_ptr != 0) {
            delete reinterpret_cast<OpenBabel::OBMol*>(mol_ptr);
        } }, "Manually delete the OBMol pointer");

    auto m_utils = m.def_submodule("utils", "Utilities and Data Structures");
    molgr::bind::bind_utils(m_utils);

    auto m_pipeline = m.def_submodule("pipeline", "Pipeline-level helpers");
    molgr::bind::bind_pipeline(m_pipeline);

    auto m_dev = m.def_submodule("dev", "Development and parity-testing helpers");
    auto m_dev_utils = m_dev.def_submodule("utils", "Development-only utility helpers");
    molgr::bind::bind_dev_utils(m_dev_utils);
    auto m_dev_pipeline = m_dev.def_submodule("pipeline", "Development-only pipeline helpers");
    molgr::bind::bind_dev_pipeline(m_dev_pipeline);
    auto m_dev_stages = m_dev.def_submodule("stages", "Development-only stage helpers");
    molgr::bind::bind_stages(m_dev_stages);
}
