/*
 * @Author: TMJ
 * @Date: 2025-12-25 20:13:30
 * @LastEditors: TMJ
 * @LastEditTime: 2025-12-26 16:20:45
 * @Description: 请填写简介
 */
/**
 * @file utils.h
 * @brief Utility functions for geometric calculations and vector math.
 * @namespace molgr::utils
 * @author TMJ
 * @date 2025-12-25
 */

#pragma once

#include <vector>
#include <openbabel/mol.h>
#include <openbabel/math/vector3.h>
#include <string>

namespace molgr
{
    namespace utils
    {

        /**
         * @brief Calculate the volume of a tetrahedron defined by 4 points.
         * * Uses the scalar triple product formula: V = |(a-d) . ((b-d) x (c-d))| / 6
         * * @param p1 Coordinates of point 1 (x, y, z).
         * @param p2 Coordinates of point 2 (x, y, z).
         * @param p3 Coordinates of point 3 (x, y, z).
         * @param p4 Coordinates of point 4 (x, y, z).
         * @return double The volume of the tetrahedron.
         */
        double CalculateTetrahedronVolume(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                          const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4);

        /**
         * @brief Calculate the shape quality of a tetrahedron.
         * * Used to score the deviation of an atom's geometry from an ideal tetrahedron.
         * Quality = (Volume / L_rms^3) * Normalization_Constant.
         * The result is clamped between 0.0 and 1.0.
         * * @param p1 Coordinates of neighbor 1.
         * @param p2 Coordinates of neighbor 2.
         * @param p3 Coordinates of neighbor 3.
         * @param p4 Coordinates of the central atom.
         * @return double Shape quality score (0.0 to 1.0).
         */
        double CalculateShapeQuality(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                     const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4);

        /**
         * @brief Helper to get an OpenBabel::vector3 from a standard vector.
         * Useful for interfacing with Python lists [x, y, z].
         * * @param coords A std::vector containing at least 3 doubles.
         * @return OpenBabel::vector3
         */
        OpenBabel::vector3 ToVector3(const std::vector<double> &coords);

        /**
         * @brief Find substructures matching a SMARTS pattern.
         * @param mol The molecule to search.
         * @param smarts The SMARTS pattern string.
         * @return A vector of matches, where each match is a vector of atom indices (1-based).
         */
        std::vector<std::vector<int>> FindSmarts(OpenBabel::OBMol &mol, const std::string &smarts);

    } // namespace utils
} // namespace molgr