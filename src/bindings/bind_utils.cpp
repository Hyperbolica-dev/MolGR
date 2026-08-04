/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:15:31
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:28:31
 * @Description: 请填写简介
 */
#include "bindings.h"
#include "molgr/python_config.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/force_field.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/organic_topology.h"
#include "molgr/utils/utils.h"
#include "molgr/utils/xyz.h"
#include "molgr/vendor/forcefielduff.h"
#include "molgr/vendor/openbabel_threading.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <pybind11/stl.h> // 必须包含

#include <memory>

static std::unique_ptr<OpenBabel::OBMol> mol_from_smiles(const std::string &smiles)
{
    auto mol = std::make_unique<OpenBabel::OBMol>();
    OpenBabel::OBConversion conv;
    conv.SetInFormat("smi");
    conv.ReadString(mol.get(), smiles);
    return mol;
}

static std::unique_ptr<OpenBabel::OBMol> mol_from_xyz(const std::string &xyz)
{
    auto mol = std::make_unique<OpenBabel::OBMol>();
    return molgr::utils::ReadXyzBlockToMol(xyz, mol.get()) ? std::move(mol) : nullptr;
}

void molgr::bind::bind_utils(py::module_ &m)
{
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
        .def_readwrite("lone_pair_count", &molgr::utils::AtomData::lone_pair_count)
        .def_readwrite("unresolved_two_electron_center", &molgr::utils::AtomData::unresolved_two_electron_center)
        .def_readwrite("hybridization", &molgr::utils::AtomData::hybridization)
        .def_readwrite("x", &molgr::utils::AtomData::x)
        .def_readwrite("y", &molgr::utils::AtomData::y)
        .def_readwrite("z", &molgr::utils::AtomData::z)
        .def("__repr__", [](const molgr::utils::AtomData &a)
             { return "<AtomData Z=" + std::to_string(a.atomic_num) +
                      " charge=" + std::to_string(a.formal_charge) +
                      " radical_num=" + std::to_string(a.radical_num) +
                      " lone_pair_count=" + std::to_string(a.lone_pair_count) +
                      " unresolved_two_electron_center=" +
                      (a.unresolved_two_electron_center ? "true" : "false") +
                      " hyb=" + std::to_string(a.hybridization) +
                      " pos=(" + std::to_string(a.x) + "," +
                      std::to_string(a.y) + "," + std::to_string(a.z) + ")>"; },
             "Return a concise debug representation of AtomData.");

    py::class_<molgr::utils::BondData>(m, "BondData")
        .def_readwrite("begin_atom_idx", &molgr::utils::BondData::begin_atom_idx)
        .def_readwrite("end_atom_idx", &molgr::utils::BondData::end_atom_idx)
        .def_readwrite("order", &molgr::utils::BondData::order)
        .def_readwrite("aromatic", &molgr::utils::BondData::aromatic)
        .def("__repr__", [](const molgr::utils::BondData &b)
             { return "<BondData " + std::to_string(b.begin_atom_idx) + "-" +
                      std::to_string(b.end_atom_idx) + " order=" + std::to_string(b.order) +
                      " aromatic=" + (b.aromatic ? "true" : "false") + ">"; },
             "Return a concise debug representation of BondData.");

    py::class_<molgr::utils::MoleculeData>(m, "MoleculeData")
        .def_readwrite("atoms", &molgr::utils::MoleculeData::atoms)
        .def_readwrite("bonds", &molgr::utils::MoleculeData::bonds)
        .def_readwrite("total_charge", &molgr::utils::MoleculeData::total_charge)
        .def_readwrite("total_radical_num", &molgr::utils::MoleculeData::total_radical_num);

    // 绑定提取函数
    m.def("extract_molecule_data", &molgr::utils::ExtractMoleculeData,
          "Extracts OBMol content into a structured object.",
          py::arg("mol_ptr"));

    m.def("molecule_data_to_obmol_ptr",
          [](const molgr::utils::MoleculeData &molecule_data)
          {
              OpenBabel::OBMol *mol = new OpenBabel::OBMol(molgr::utils::MolFromMoleculeData(molecule_data));
              return reinterpret_cast<intptr_t>(mol);
          },
          "Converts MoleculeData to a newly allocated OBMol pointer. Free it with _core.free_obmol_ptr.",
          py::arg("molecule_data"));
}

void molgr::bind::bind_dev_utils(py::module_ &m)
{
    m.def("debug_xyz_seed_molecule_data",
          [](const std::string &xyz_block)
          {
              OpenBabel::OBMol mol;
              molgr::utils::MoleculeData molecule_data;
              {
                  py::gil_scoped_release release;
                  if (!molgr::utils::ReadXyzBlockToMol(xyz_block, &mol))
                  {
                      throw std::runtime_error("failed to parse XYZ block");
                  }
                  molecule_data = molgr::utils::MoleculeDataFromOBMol(mol);
              }
              return molecule_data;
          },
          "Return the C++ vendor-perceived seed molecule data for an XYZ block.",
          py::arg("xyz_block"));

    m.def("compute_organic_topology_metrics_ptr",
          [](intptr_t mol_ptr, py::object config)
          {
              if (mol_ptr == 0)
              {
                  throw std::runtime_error("null OBMol pointer");
              }
              const auto runtime_config = molgr::config::FromPython(config);
              auto *mol = reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
              const auto metrics = molgr::organic_topology::ComputeOrganicTopologyMetrics(
                  *mol,
                  runtime_config.organic_topology);
              py::dict out;
              out["aromatic_atom_count"] = metrics.aromatic_atom_count;
              out["aromatic_ring_count"] = metrics.aromatic_ring_count;
              out["aromatic_stability_score"] = metrics.aromatic_stability_score;
              out["conjugated_atom_count"] = metrics.conjugated_atom_count;
              out["conjugated_bond_count"] = metrics.conjugated_bond_count;
              out["max_conjugated_component_size"] = metrics.max_conjugated_component_size;
              out["conjugated_atom_indices"] = metrics.conjugated_atom_indices;
              out["hyperconjugative_donor_count"] = metrics.hyperconjugative_donor_count;
              out["hyperconjugation_score"] = metrics.hyperconjugation_score;
              return out;
          },
          "Compute C++ organic topology metrics for an OBMol pointer.",
          py::arg("mol_ptr"),
          py::kw_only(),
          py::arg("config") = py::none());

    m.def("organic_force_field_energy_ptr",
          [](intptr_t mol_ptr, py::object config)
          {
              if (mol_ptr == 0)
              {
                  throw std::runtime_error("null OBMol pointer");
              }
              const auto runtime_config = molgr::config::FromPython(config);
              auto *mol = reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
              return molgr::scoring::OrganicForceFieldEvaluation(*mol, runtime_config).energy_kj_mol;
          },
          "Compute C++ organic force-field energy for an OBMol pointer.",
          py::arg("mol_ptr"),
          py::kw_only(),
          py::arg("config") = py::none());

    m.def("debug_vendor_uff_ptr",
          [](intptr_t mol_ptr)
          {
              if (mol_ptr == 0)
              {
                  throw std::runtime_error("null OBMol pointer");
              }
              auto *mol = reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
              OpenBabel::OBMol working = molgr::utils::CloneMolTopologyOnly(*mol);
              molgr::vendor::openbabel_threading::SetAromaticPerceived(working, false);
              OpenBabel::MolgrForceFieldUFF force_field("MolGR-UFF-debug", false);
              force_field.SetLogLevel(OBFF_LOGLVL_NONE);
              py::dict out;
              const bool setup_ok = force_field.Setup(working);
              out["setup_ok"] = setup_ok;
              py::list atom_types;
              if (setup_ok)
              {
                  for (const std::string &atom_type : force_field.DebugAtomTypes())
                  {
                      atom_types.append(atom_type);
                  }
                  out["energy"] = force_field.Energy(false);
                  out["bond"] = force_field.E_Bond(false);
                  out["angle"] = force_field.E_Angle(false);
                  out["torsion"] = force_field.E_Torsion(false);
                  out["oop"] = force_field.E_OOP(false);
                  out["vdw"] = force_field.E_VDW(false);
              }
              out["atom_types"] = atom_types;
              return out;
          },
          "Return MolGR vendor UFF atom types and energy terms for an OBMol pointer.",
          py::arg("mol_ptr"));

    m.def("test_symmetry_penalty", [](const std::string &smiles)
          {
              auto mol = mol_from_smiles(smiles);
              return molgr::scoring::CalcSymmetryPenalty(*mol);
          },
          "Calculate symmetry penalty from SMILES (For Testing)",
          py::arg("smiles"));

    m.def("test_physchem_penalty", [](const std::string &smiles)
          {
              auto mol = mol_from_smiles(smiles);
              return molgr::scoring::CalculatePhysChemPenalty(*mol);
          },
          "Calculate PhysChem penalty from SMILES (For Testing)",
          py::arg("smiles"));

    m.def("test_deviation_score", [](const std::string &xyz_block, int atom_idx)
          {
              auto mol = mol_from_xyz(xyz_block);
              if (!mol)
              {
                  return -1.0;
              }
              OpenBabel::OBAtom *atom = mol->GetAtom(atom_idx);
              if (!atom)
              {
                  return -1.0;
              }
              return molgr::scoring::GetDeviationScore(*mol, atom);
          },
          "Calculate geometry deviation for atom (1-based index) from XYZ (For Testing)",
          py::arg("xyz_block"),
          py::arg("atom_idx"));

    m.def("test_total_score", [](const std::string &xyz_block)
          {
              auto mol = mol_from_xyz(xyz_block);
              if (!mol)
              {
                  throw std::runtime_error("failed to parse XYZ");
              }
              double score = molgr::scoring::OmolScore(*mol);
              const int n = mol->NumAtoms();
              for (int i = 1; i <= n; ++i)
              {
                  OpenBabel::OBAtom *a = mol->GetAtom(i);
                  if (!a)
                  {
                      continue;
                  }
                  for (int j = i + 1; j <= n; ++j)
                  {
                      OpenBabel::OBAtom *b = mol->GetAtom(j);
                      if (!b)
                      {
                          continue;
                      }
                      const double dist = a->GetDistance(b);
                      if (dist < 0.5)
                      {
                          score += (0.5 - dist) * 1000.0;
                      }
                  }
              }
              return score;
          },
          "Calculate total OMolScore from XYZ block (For Testing)",
          py::arg("xyz_block"));
}
