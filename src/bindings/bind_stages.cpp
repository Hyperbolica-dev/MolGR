#include "bindings.h"

#include "molgr/stages/clean.h"
#include "molgr/stages/break_bond.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/fresh.h"
#include "molgr/stages/preprocess.h"

#include <cstdint>
#include <stdexcept>

namespace molgr {
namespace bind {

static OpenBabel::OBMol *require_obmol_ptr(intptr_t mol_ptr)
{
    if (mol_ptr == 0)
    {
        throw std::runtime_error("null OBMol pointer");
    }
    return reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
}

void bind_stages(py::module_ &m)
{
    // These dev helpers operate on OBMol pointers owned by Python OpenBabel wrappers.
    // Keep the GIL held while using those pointers and while constructing py::object results.
    const auto make_connections_ptr = [](intptr_t mol_ptr, double factor)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::MakeConnections(*mol, factor);
    };

    const auto pre_clean_ptr = [](intptr_t mol_ptr)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::PreClean(*mol);
    };

    const auto fresh_omol_charge_radical_ptr = [](intptr_t mol_ptr)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::FreshOmolChargeRadical(*mol);
    };

    const auto assign_radical_dots_ptr = [](intptr_t mol_ptr, int atom_idx) -> int
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        auto *atom = mol->GetAtom(atom_idx);
        if (!atom)
        {
            throw std::runtime_error("invalid atom index");
        }
        return molgr::reconstruct::AssignRadicalDots(*atom);
    };

    const auto assign_charge_radical_for_atom_ptr = [](intptr_t mol_ptr, int atom_idx)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        auto *atom = mol->GetAtom(atom_idx);
        if (!atom)
        {
            throw std::runtime_error("invalid atom index");
        }
        return molgr::reconstruct::AssignChargeRadicalForAtom(*atom);
    };

    const auto validate_omol_ptr = [](intptr_t mol_ptr, int total_charge, int total_radical) -> bool
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::ValidateOmol(*mol, total_charge, total_radical);
    };

    const auto eliminate_1_3_dipole_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::Eliminate13Dipole(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_positive_charges_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminatePositiveCharges(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_negative_charges_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateNegativeCharges(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_nnn_ptr = [](intptr_t mol_ptr, int given_charge, bool positive) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateNNN(*mol, given_charge, positive);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_high_positive_charge_atoms_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateHighPositiveChargeAtoms(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_cn_in_doubt_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateCNInDoubt(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_carboxyl_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateCarboxyl(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_carbene_neighbor_heteroatom_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateCarbeneNeighborHeteroatom(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto eliminate_charge_spliting_ptr = [](intptr_t mol_ptr, int given_charge) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::EliminateChargeSpliting(*mol, given_charge);
        return py::make_tuple(given_charge, hit);
    };

    const auto clean_resonances_ptr = [](intptr_t mol_ptr)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::CleanResonances(*mol);
    };

    const auto clean_neighbor_radicals_ptr = [](intptr_t mol_ptr)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::CleanNeighborRadicals(*mol);
    };

    const auto clean_carbene_neighbor_unsaturated_ptr = [](intptr_t mol_ptr)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::CleanCarbeneNeighborUnsaturated(*mol);
    };

    const auto break_deformed_ene_ptr = [](intptr_t mol_ptr, int given_charge, int given_radical, double tolerance)
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        return molgr::reconstruct::BreakDeformedEne(*mol, given_charge, given_radical, tolerance);
    };

    const auto break_one_bond_ptr = [](intptr_t mol_ptr, int given_charge, int given_radical) -> py::tuple
    {
        auto *mol = require_obmol_ptr(mol_ptr);
        const bool hit = molgr::reconstruct::BreakOneBond(*mol, given_charge, given_radical);
        return py::make_tuple(given_charge, hit);
    };

    auto m_preprocess = m.def_submodule("preprocess", "Fallback-aligned preprocess stage helpers");
    auto m_fresh = m.def_submodule("fresh", "Fallback-aligned fresh stage helpers");
    auto m_eliminate = m.def_submodule("eliminate", "Fallback-aligned eliminate stage helpers");
    auto m_clean = m.def_submodule("clean", "Fallback-aligned clean stage helpers");
    auto m_break_bond = m.def_submodule("break_bond", "Fallback-aligned break_bond stage helpers");

    m_preprocess.def(
        "make_connections_ptr",
        make_connections_ptr,
        R"pbdoc(
Apply preprocess.make_connections to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    factor: distance factor (default matches python fallback)
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("factor") = 1.4);

    m_preprocess.def(
        "pre_clean_ptr",
        pre_clean_ptr,
        R"pbdoc(
Apply preprocess.pre_clean to an existing OBMol.
)pbdoc",
        py::arg("mol_ptr"));

    m_preprocess.def(
        "validate_omol_ptr",
        validate_omol_ptr,
        R"pbdoc(
Validate conservation of total charge and radical electrons.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    total_charge: expected total formal charge
    total_radical: expected total radical electrons
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("total_charge"),
        py::arg("total_radical"));

    m_fresh.def(
        "fresh_omol_charge_radical_ptr",
        fresh_omol_charge_radical_ptr,
        R"pbdoc(
Apply fresh.fresh_omol_charge_radical to an existing OBMol.
)pbdoc",
        py::arg("mol_ptr"));

    m_fresh.def(
        "assign_radical_dots_ptr",
        assign_radical_dots_ptr,
        R"pbdoc(
Assign radical dots for a specific atom on an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    atom_idx: 1-based atom index
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("atom_idx"));

    m_fresh.def(
        "assign_charge_radical_for_atom_ptr",
        assign_charge_radical_for_atom_ptr,
        R"pbdoc(
Assign charge/radical state for a specific atom on an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    atom_idx: 1-based atom index
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("atom_idx"));

    m_eliminate.def(
        "eliminate_1_3_dipole_ptr",
        eliminate_1_3_dipole_ptr,
        R"pbdoc(
Apply eliminate.eliminate_1_3_dipole to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_positive_charges_ptr",
        eliminate_positive_charges_ptr,
        R"pbdoc(
Apply eliminate.eliminate_positive_charges to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_negative_charges_ptr",
        eliminate_negative_charges_ptr,
        R"pbdoc(
Apply eliminate.eliminate_negative_charges to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_nnn_ptr",
        eliminate_nnn_ptr,
        R"pbdoc(
Apply eliminate.eliminate_nnn to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
    positive: whether to run positive-direction elimination
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"),
        py::arg("positive") = false);

    m_eliminate.def(
        "eliminate_high_positive_charge_atoms_ptr",
        eliminate_high_positive_charge_atoms_ptr,
        R"pbdoc(
Apply eliminate.eliminate_high_positive_charge_atoms to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_cn_in_doubt_ptr",
        eliminate_cn_in_doubt_ptr,
        R"pbdoc(
Apply eliminate.eliminate_cn_in_doubt to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_carboxyl_ptr",
        eliminate_carboxyl_ptr,
        R"pbdoc(
Apply eliminate.eliminate_carboxyl to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_carbene_neighbor_heteroatom_ptr",
        eliminate_carbene_neighbor_heteroatom_ptr,
        R"pbdoc(
Apply eliminate.eliminate_carbene_neighbor_heteroatom to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_eliminate.def(
        "eliminate_charge_spliting_ptr",
        eliminate_charge_spliting_ptr,
        R"pbdoc(
Apply eliminate.eliminate_charge_spliting to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge deficit to be updated in place and returned
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge"));

    m_clean.def(
        "clean_resonances_ptr",
        clean_resonances_ptr,
        R"pbdoc(
Apply clean.clean_resonances to an existing OBMol.
)pbdoc",
        py::arg("mol_ptr"));

    m_clean.def(
        "clean_neighbor_radicals_ptr",
        clean_neighbor_radicals_ptr,
        R"pbdoc(
Apply clean.clean_neighbor_radicals to an existing OBMol.
)pbdoc",
        py::arg("mol_ptr"));

    m_clean.def(
        "clean_carbene_neighbor_unsaturated_ptr",
        clean_carbene_neighbor_unsaturated_ptr,
        R"pbdoc(
Apply clean.clean_carbene_neighbor_unsaturated to an existing OBMol.
)pbdoc",
        py::arg("mol_ptr"));

    m_break_bond.def(
        "break_deformed_ene_ptr",
        break_deformed_ene_ptr,
        R"pbdoc(
Apply break_bond.break_deformed_ene to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: target charge budget used by the stage
    given_radical: target radical budget used by the stage
    tolerance: torsion-angle tolerance threshold
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge") = 0,
        py::arg("given_radical") = 0,
        py::arg("tolerance") = 5.0);

    m_break_bond.def(
        "break_one_bond_ptr",
        break_one_bond_ptr,
        R"pbdoc(
Apply break_bond.break_one_bond to an existing OBMol.

Args:
    mol_ptr: int address of OpenBabel::OBMol
    given_charge: charge budget to be updated in place and returned
    given_radical: target radical budget used by the stage
)pbdoc",
        py::arg("mol_ptr"),
        py::arg("given_charge") = 0,
        py::arg("given_radical") = 0);

}

}
}
