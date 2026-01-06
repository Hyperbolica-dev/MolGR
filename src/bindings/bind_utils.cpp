/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:15:31
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:28:31
 * @Description: 请填写简介
 */
#include "bindings.h"
#include "molgr/utils.h"
#include <pybind11/stl.h> // 必须包含

void molgr::bind::bind_utils(py::module_ &m)
{
    m.def("calculate_tetrahedron_volume", [](const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &p3, const std::vector<double> &p4)
          { return molgr::utils::CalculateTetrahedronVolume(
                molgr::utils::ToVector3(p1),
                molgr::utils::ToVector3(p2),
                molgr::utils::ToVector3(p3),
                molgr::utils::ToVector3(p4)); },
          R"pbdoc(
            Calculate the volume of a tetrahedron defined by 4 points.
            
            Parameters:
                p1 (list[float]): Coordinates of the first atom.
                p2 (list[float]): Coordinates of the second atom.
                p3 (list[float]): Coordinates of the third atom.
                p4 (list[float]): Coordinates of the fourth atom.
            
            Returns:
                float: The volume.
        )pbdoc",
          py::arg("p1"), py::arg("p2"), py::arg("p3"), py::arg("p4"));

    m.def("calculate_shape_quality", [](const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &p3, const std::vector<double> &p4)
          { return molgr::utils::CalculateShapeQuality(
                molgr::utils::ToVector3(p1),
                molgr::utils::ToVector3(p2),
                molgr::utils::ToVector3(p3),
                molgr::utils::ToVector3(p4)); },
          R"pbdoc(
            Calculate the shape quality score of a tetrahedron.
            
            Parameters:
                p1 (list[float]): Coordinates of the first atom.
                p2 (list[float]): Coordinates of the second atom.
                p3 (list[float]): Coordinates of the third atom.
                p4 (list[float]): Coordinates of the fourth atom.
            
            Returns:
                float: Quality score between 0.0 (coplanar/bad) and 1.0 (ideal).
        )pbdoc",
          py::arg("p1"), py::arg("p2"), py::arg("p3"), py::arg("p4"));

    py::class_<molgr::utils::AtomData>(m, "AtomData")
        .def_readwrite("atomic_num", &molgr::utils::AtomData::atomic_num)
        .def_readwrite("formal_charge", &molgr::utils::AtomData::formal_charge)
        .def_readwrite("radical_num", &molgr::utils::AtomData::radical_num)
        .def_readwrite("x", &molgr::utils::AtomData::x)
        .def_readwrite("y", &molgr::utils::AtomData::y)
        .def_readwrite("z", &molgr::utils::AtomData::z)
        .def("__repr__", [](const molgr::utils::AtomData &a)
             { return "<AtomData Z=" + std::to_string(a.atomic_num) +
                      " charge=" + std::to_string(a.formal_charge) +
                      " radical_num=" + std::to_string(a.radical_num) +
                      " pos=(" + std::to_string(a.x) + "," +
                      std::to_string(a.y) + "," + std::to_string(a.z) + ")>"; });

    py::class_<molgr::utils::BondData>(m, "BondData")
        .def_readwrite("begin_atom_idx", &molgr::utils::BondData::begin_atom_idx)
        .def_readwrite("end_atom_idx", &molgr::utils::BondData::end_atom_idx)
        .def_readwrite("order", &molgr::utils::BondData::order)
        .def("__repr__", [](const molgr::utils::BondData &b)
             { return "<BondData " + std::to_string(b.begin_atom_idx) + "-" +
                      std::to_string(b.end_atom_idx) + " order=" + std::to_string(b.order) + ">"; });

    py::class_<molgr::utils::MoleculeData>(m, "MoleculeData")
        .def_readwrite("atoms", &molgr::utils::MoleculeData::atoms)
        .def_readwrite("bonds", &molgr::utils::MoleculeData::bonds)
        .def_readwrite("total_charge", &molgr::utils::MoleculeData::total_charge)
        .def_readwrite("total_radical_num", &molgr::utils::MoleculeData::total_radical_num);

    // 绑定提取函数
    m.def("extract_molecule_data", &molgr::utils::ExtractMoleculeData,
          "Extracts OBMol content into a structured object.",
          py::arg("mol_ptr"));
}