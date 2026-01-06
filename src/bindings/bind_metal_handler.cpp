/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:19:28
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:26:56
 * @Description: 请填写简介
 */
#include "bindings.h"
#include "molgr/metal_handler.h"
#include <pybind11/stl.h>

void molgr::bind::bind_metal(py::module_ &m)
{
    // 绑定 MetalAtomPosition 结构体
    py::class_<molgr::metal::MetalAtomPosition>(m, "MetalAtomPosition")
        .def_readwrite("idx", &molgr::metal::MetalAtomPosition::idx)
        .def_readwrite("symbol", &molgr::metal::MetalAtomPosition::symbol)
        .def_readwrite("element_idx", &molgr::metal::MetalAtomPosition::element_idx)
        .def_readwrite("valence", &molgr::metal::MetalAtomPosition::valence)
        .def_readwrite("radical_num", &molgr::metal::MetalAtomPosition::radical_num)
        .def_readwrite("x", &molgr::metal::MetalAtomPosition::x)
        .def_readwrite("y", &molgr::metal::MetalAtomPosition::y)
        .def_readwrite("z", &molgr::metal::MetalAtomPosition::z)
        .def("__repr__", [](const molgr::metal::MetalAtomPosition &m)
             { return "<MetalAtomPosition " + m.symbol +
                      " val=" + std::to_string(m.valence) +
                      " rad=" + std::to_string(m.radical_num) + ">"; });

    // 绑定 MetalHandler 类
    py::class_<molgr::metal::MetalHandler>(m, "MetalHandler")
        .def(py::init([](intptr_t mol_ptr)
                      {
            // 构造函数接收指针，转为引用传给 C++
            OpenBabel::OBMol* mol = reinterpret_cast<OpenBabel::OBMol*>(mol_ptr);
            if (!mol) throw std::runtime_error("Null pointer passed to MetalHandler");
            return new molgr::metal::MetalHandler(*mol); }),
             py::arg("mol_ptr"))

        .def("strip_metals", [](molgr::metal::MetalHandler &self, intptr_t mol_ptr)
             {
            // strip_metals 也需要操作原始 OBMol
            OpenBabel::OBMol* mol = reinterpret_cast<OpenBabel::OBMol*>(mol_ptr);
            return self.StripMetals(*mol); }, py::arg("mol_ptr"))

        .def("generate_combinations", &molgr::metal::MetalHandler::GenerateCombinations, py::arg("total_radical_electrons"))

        .def_static("combine_metal_with_mol", [](intptr_t mol_ptr, const std::vector<molgr::metal::MetalAtomPosition> &metals)
                    {
            OpenBabel::OBMol* mol = reinterpret_cast<OpenBabel::OBMol*>(mol_ptr);
            molgr::metal::MetalHandler::CombineMetalWithMol(*mol, metals); }, py::arg("mol_ptr"), py::arg("metals"));
}