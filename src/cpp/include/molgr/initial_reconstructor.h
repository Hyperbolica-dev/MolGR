/**
 * @file initial_reconstructor.h
 * @brief Module for initial graph building and rule-based cleaning (pre-resonance).
 * @details Corresponds to the first half of 'xyz_to_omol_no_metal' in Python.
 * @author TMJ
 * @date 2025-12-28
 */

#pragma once

#include <openbabel/mol.h>
#include <string>
#include <memory>

namespace molgr
{
    namespace reconstruct
    {

        /**
         * @brief Main entry point for initial reconstruction (No Metal).
         * @details
         * Python Correspondence: xyz_to_omol_no_metal (Part 1 - before resonance search)
         * * Flow:
         * 1. Read XYZ
         * 2. MakeConnections
         * 3. PreClean
         * 4. FreshOmolChargeRadical
         * 5. Run sequence of Eliminate* rules
         * 6. Validate (Check charge/radical conservation)
         * * @param xyz_block The XYZ string without metals.
         * @param total_charge Target total charge.
         * @param total_radical Target total radical electrons.
         * @return unique_ptr to OBMol if successful, nullptr if validation fails.
         */
        std::unique_ptr<OpenBabel::OBMol> ReconstructFromXYZNoMetal(
            const std::string &xyz_block,
            int total_charge,
            int total_radical);

        // =============================================================================
        // Individual Rules (Exposed for testing/fine-grained control)
        // =============================================================================

        // Python: make_connections
        void MakeConnections(OpenBabel::OBMol &mol, double factor = 1.4);

        // Python: pre_clean
        void PreClean(OpenBabel::OBMol &mol);

        // Python: fresh_omol_charge_radical
        void FreshOmolChargeRadical(OpenBabel::OBMol &mol);

        // --- Elimination Rules ---

        // Python: eliminate_NNN
        void EliminateNNN(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: eliminate_high_positive_charge_atoms
        void EliminateHighPositiveChargeAtoms(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: eliminate_CN_in_doubt
        void EliminateCNInDoubt(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: eliminate_carboxyl
        void EliminateCarboxyl(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: clean_carbene_neighbor_unsaturated
        void CleanCarbeneNeighborUnsaturated(OpenBabel::OBMol &mol);

        // Python: eliminate_carbene_neighbor_heteroatom
        void EliminateCarbeneNeighborHeteroatom(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: clean_neighbor_radicals
        void CleanNeighborRadicals(OpenBabel::OBMol &mol);

        // Python: eliminate_charge_spliting
        void EliminateChargeSpliting(OpenBabel::OBMol &mol, int &current_charge_deficit);

        // Python: break_deformed_ene
        void BreakDeformedEne(OpenBabel::OBMol &mol, int allowed_charge, int allowed_radical, double tolerance = 5.0);

        // Python: break_one_bond
        void BreakOneBond(OpenBabel::OBMol &mol, int &current_charge_deficit, int allowed_radical);

    } // namespace reconstruct
} // namespace molgr