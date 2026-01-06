/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:26:22
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:27:43
 * @Description: 请填写简介
 */
#include "bindings.h"
#include "molgr/scoring.h"
#include <openbabel/obconversion.h>

#define CAST_TO_OBMOL(ptr) reinterpret_cast<OpenBabel::OBMol *>(ptr)

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
void molgr::bind::bind_scoring(py::module_ &m)
{
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
        
        Parameters:
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