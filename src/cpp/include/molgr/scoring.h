/**
 * @file scoring.h
 * @brief Scoring functions for molecular graph reconstruction.
 * @details Evaluates the plausibility of a molecular structure based on geometry,
 * physical chemistry rules (electronegativity, charge distribution), and metal coordination.
 * The lower the score, the more plausible the structure.
 * * @namespace molgr::scoring
 * @author TMJ
 * @date 2025-12-25
 */

#pragma once

#include <openbabel/mol.h>
#include <vector>

namespace molgr
{
    namespace scoring
    {

        /**
         * @brief Calculate the total resonance score of a molecule.
         * * @details This is the main entry point for the scoring system. It aggregates penalties from:
         * 1. Symmetry (Is the graph symmetry consistent with atom count?)
         * 2. Geometry Deviation (Are bond angles/shapes consistent with hybridization?)
         * 3. Physical Chemistry (Are charges and radicals on appropriate atoms?)
         * 4. Metal Coordination (Are metal valences and ligand environments reasonable?)
         * * @param mol The OpenBabel molecule object to evaluate.
         * @return double Total penalty score (Lower is better).
         */
        double OmolScore(OpenBabel::OBMol &mol);

        // --- Sub-functions exposed for unit testing ---

        /**
         * @brief Calculate penalty based on graph symmetry.
         * @param mol Target molecule.
         * @return double Penalty score.
         */
        double CalcSymmetryPenalty(OpenBabel::OBMol &mol);

        /**
         * @brief Calculate total physical chemistry penalty (Charge + Radical + Coulombic).
         * @param mol Target molecule.
         * @return double Penalty score.
         */
        double CalculatePhysChemPenalty(OpenBabel::OBMol &mol);

        /**
         * @brief Calculate metal-specific penalties (Valence state + Ligand interaction).
         * @param mol Target molecule.
         * @return double Penalty score.
         */
        double CalculateMetalPenalty(OpenBabel::OBMol &mol);

        /**
         * @brief Calculate geometric deviation score for a single atom.
         * @param mol Target molecule.
         * @param atom The atom to check.
         * @return double Deviation score (0.0 to 1.0).
         */
        double GetDeviationScore(OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom);

    } // namespace scoring
} // namespace molgr