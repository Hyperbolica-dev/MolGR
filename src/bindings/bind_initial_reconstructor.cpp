/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:21:59
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:24:07
 * @Description: 请填写简介
 */
#include "bindings.h"
#include "molgr/initial_reconstructor.h"

void molgr::bind::bind_reconstruct(py::module_ &m)
{
    m.def("reconstruct_from_xyz_no_metal", [](const std::string &xyz_block, int total_charge, int total_radical) -> intptr_t
          {
        // 调用 C++ 重建逻辑
        auto mol_ptr = molgr::reconstruct::ReconstructFromXYZNoMetal(xyz_block, total_charge, total_radical);
        
        if (!mol_ptr) {
            return 0; // 返回空指针表示失败
        }

        // 关键：释放 unique_ptr 的所有权，防止 C++ 自动删除对象
        // 返回裸指针地址给 Python，Python 端需负责将其包装为 OBMol 并管理生命周期
        return reinterpret_cast<intptr_t>(mol_ptr.release()); },
          R"pbdoc(
        Reconstruct molecule topology and state from XYZ block (No Metals).
        
        Parameters:
            xyz_block (str): The XYZ content.
            total_charge (int): Target charge.
            total_radical (int): Target radical electrons.
            
        Returns:
            int: Memory address (pointer) of the created OBMol, or 0 if failed.
    )pbdoc",
          py::arg("xyz_block"), py::arg("total_charge"), py::arg("total_radical"));
}