/**
 * @file module.cpp
 * @brief Pybind11 bindings for the molgr C++ library.
 * @author TMJ
 * @date 2025-12-25
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // 核心：自动转换 STL 容器
#include <cstdint>        // 用于 intptr_t

#include "molgr/consts.h"  // 常数模块
#include "molgr/utils.h"   // 工具模块
#include "molgr/scoring.h" // 评分模块
#include "molgr/logger.h"  // 日志模块

// OpenBabel 头文件 (用于指针转换和内部测试构建)
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

namespace py = pybind11;

// 引用工具函数
using molgr::utils::ToVector3;

// =============================================================================
// 内部辅助函数 (Internal Helpers)
// =============================================================================

// 宏：将整数地址强制转换为 OBMol 指针
// 这是连接 Python SWIG 对象和 C++ Pybind11 的桥梁
#define CAST_TO_OBMOL(ptr) reinterpret_cast<OpenBabel::OBMol *>(ptr)

// 辅助：从 SMILES 创建 OBMol (用于 C++ 独立单元测试)
std::unique_ptr<OpenBabel::OBMol> MolFromSmiles(const std::string &smiles)
{
  auto mol = std::make_unique<OpenBabel::OBMol>();
  OpenBabel::OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(mol.get(), smiles);
  return mol;
}

// 辅助：从 XYZ Block 创建 OBMol (用于 C++ 独立单元测试)
std::unique_ptr<OpenBabel::OBMol> MolFromXYZ(const std::string &xyz)
{
  auto mol = std::make_unique<OpenBabel::OBMol>();
  OpenBabel::OBConversion conv;
  conv.SetInFormat("xyz");
  conv.ReadString(mol.get(), xyz);
  return mol;
}

// =============================================================================
// 模块定义
// =============================================================================

PYBIND11_MODULE(_core, m)
{
  m.doc() = R"pbdoc(
        MolGR Core C++ Implementation
        -----------------------------
        Exposes optimized C++ algorithms for molecular graph reconstruction.
    )pbdoc";

  // =========================================================================
  // 0. Logging Module (新增)
  // =========================================================================

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

  // =========================================================================
  // 1. Consts 模块绑定
  // =========================================================================

  m.def("get_possible_metal_radicals", &molgr::GetPossibleMetalRadicals,
        R"pbdoc(
            Get possible radical electron counts for a metal given its valence.
            
            Args:
                metal (str): The chemical symbol (e.g., "Fe").
                valence (int): The oxidation state.
            
            Returns:
                set[int]: A set of possible unpaired electron counts.
        )pbdoc",
        py::arg("metal"), py::arg("valence"));

  // =========================================================================
  // 2. Utils 模块绑定
  // =========================================================================

  m.def("calculate_tetrahedron_volume", [](const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &p3, const std::vector<double> &p4)
        { return molgr::utils::CalculateTetrahedronVolume(
              ToVector3(p1), ToVector3(p2), ToVector3(p3), ToVector3(p4)); },
        R"pbdoc(
            Calculate the volume of a tetrahedron defined by 4 points.
            
            Args:
                p1, p2, p3, p4 (list[float]): Coordinates [x, y, z].
            
            Returns:
                float: The volume.
        )pbdoc",
        py::arg("p1"), py::arg("p2"), py::arg("p3"), py::arg("p4"));

  m.def("calculate_shape_quality", [](const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &p3, const std::vector<double> &p4)
        { return molgr::utils::CalculateShapeQuality(
              ToVector3(p1), ToVector3(p2), ToVector3(p3), ToVector3(p4)); },
        R"pbdoc(
            Calculate the shape quality score of a tetrahedron.
            
            Args:
                p1, p2, p3 (list[float]): Coordinates of neighbor atoms.
                p4 (list[float]): Coordinates of the central atom.
            
            Returns:
                float: Quality score between 0.0 (coplanar/bad) and 1.0 (ideal).
        )pbdoc",
        py::arg("p1"), py::arg("p2"), py::arg("p3"), py::arg("p4"));

  // =========================================================================
  // 3. Scoring 模块 API (生产环境用 - 指针传递)
  // =========================================================================

  m.def("omol_score_from_ptr", [](intptr_t mol_ptr)
        {
        OpenBabel::OBMol* mol = CAST_TO_OBMOL(mol_ptr);
        if (!mol) throw std::runtime_error("Received null pointer for OBMol");
        return molgr::scoring::OmolScore(*mol); }, R"pbdoc(
        Calculate total OMolScore using a memory pointer to an OpenBabel::OBMol.
        This allows compatibility with SWIG-wrapped OpenBabel objects.
        
        Args:
            mol_ptr (int): The memory address of the OBMol object (use `int(mol.this)` in Python).
            
        Returns:
            float: Total score.
    )pbdoc",
        py::arg("mol_ptr"));

  // =========================================================================
  // 4. Scoring 模块 API (测试环境用 - 字符串构建)
  // =========================================================================
  // 这些函数用于 pytest 单元测试，独立于 Python 的 openbabel 包

  m.def("test_symmetry_penalty", [](const std::string &smiles)
        {
        auto mol = MolFromSmiles(smiles);
        return molgr::scoring::CalcSymmetryPenalty(*mol); }, "Calculate symmetry penalty from SMILES (For Testing)", py::arg("smiles"));

  m.def("test_physchem_penalty", [](const std::string &smiles)
        {
        auto mol = MolFromSmiles(smiles);
        return molgr::scoring::CalculatePhysChemPenalty(*mol); }, "Calculate PhysChem penalty from SMILES (For Testing)", py::arg("smiles"));

  m.def("test_deviation_score", [](const std::string &xyz_block, int atom_idx)
        {
        auto mol = MolFromXYZ(xyz_block);
        OpenBabel::OBAtom* atom = mol->GetAtom(atom_idx);
        if (!atom) return -1.0; 
        return molgr::scoring::GetDeviationScore(*mol, atom); }, "Calculate geometry deviation for atom (1-based index) from XYZ (For Testing)", py::arg("xyz_block"), py::arg("atom_idx"));

  m.def("test_total_score", [](const std::string &xyz_block)
        {
        auto mol = MolFromXYZ(xyz_block);
        return molgr::scoring::OmolScore(*mol); }, "Calculate total OMolScore from XYZ block (For Testing)", py::arg("xyz_block"));
}