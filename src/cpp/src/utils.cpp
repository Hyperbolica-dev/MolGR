/**
 * @file utils.cpp
 * @brief Implementation of geometric utility functions.
 * @namespace molgr::utils
 * @author TMJ
 * @date 2025-12-25
 */

#include "molgr/utils.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <openbabel/parsmart.h>

namespace molgr
{
    namespace utils
    {

        // -----------------------------------------------------------------------------
        // Internal Helper: Calculate Squared Length manually
        // OpenBabel::vector3 doesn't have LengthSq(), so we use dot(v, v)
        // -----------------------------------------------------------------------------
        inline double LengthSq(const OpenBabel::vector3 &v)
        {
            return OpenBabel::dot(v, v);
        }

        // -----------------------------------------------------------------------------
        // Implementation of CalculateTetrahedronVolume
        // -----------------------------------------------------------------------------
        double CalculateTetrahedronVolume(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                          const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4)
        {
            // 算法原理：
            // V = (1/6) * | (p1-p4) · ((p2-p4) x (p3-p4)) |

            OpenBabel::vector3 v1 = p1 - p4;
            OpenBabel::vector3 v2 = p2 - p4;
            OpenBabel::vector3 v3 = p3 - p4;

            // 注意：显式使用 OpenBabel::cross 和 OpenBabel::dot
            OpenBabel::vector3 cross_product = OpenBabel::cross(v2, v3);
            double scalar_triple_product = OpenBabel::dot(v1, cross_product);

            return std::abs(scalar_triple_product) / 6.0;
        }

        // -----------------------------------------------------------------------------
        // Implementation of CalculateShapeQuality
        // -----------------------------------------------------------------------------
        double CalculateShapeQuality(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                     const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4)
        {

            // 1. 计算体积
            double volume = CalculateTetrahedronVolume(p1, p2, p3, p4);

            // 如果体积接近 0（共面），直接返回 0
            if (volume < 1e-9)
            {
                return 0.0;
            }

            // 2. 计算 6 条边的长度平方和
            // 使用我们定义的 helper 函数 LengthSq
            double edges_sq_sum = 0.0;
            edges_sq_sum += LengthSq(p1 - p2);
            edges_sq_sum += LengthSq(p1 - p3);
            edges_sq_sum += LengthSq(p1 - p4);
            edges_sq_sum += LengthSq(p2 - p3);
            edges_sq_sum += LengthSq(p2 - p4);
            edges_sq_sum += LengthSq(p3 - p4);

            // 3. 计算均方根边长的立方 (L_rms^3)
            // L_rms = sqrt( sum(edges^2) / 6 )
            // L_rms^3 = ( sum(edges^2) / 6 ) ^ 1.5
            double l_rms_squared = edges_sq_sum / 6.0;
            double l_rms_cubed = std::pow(l_rms_squared, 1.5);

            // 4. 计算 Quality
            // Normalization Constant = 6 * sqrt(2) ≈ 8.48528
            const double NORMALIZATION_CONST = 6.0 * std::sqrt(2.0);

            // 避免除以零
            if (l_rms_cubed < 1e-9)
                return 0.0;

            double quality = NORMALIZATION_CONST * (volume / l_rms_cubed);

            // 5. Clamp result to [0.0, 1.0]
            return std::max(0.0, std::min(1.0, quality));
        }

        // -----------------------------------------------------------------------------
        // Implementation of Helper ToVector3
        // -----------------------------------------------------------------------------
        OpenBabel::vector3 ToVector3(const std::vector<double> &coords)
        {
            if (coords.size() < 3)
                return OpenBabel::vector3(0, 0, 0);
            return OpenBabel::vector3(coords[0], coords[1], coords[2]);
        }

        std::vector<std::vector<int>> FindSmarts(OpenBabel::OBMol &mol, const std::string &smarts)
        {
            OpenBabel::OBSmartsPattern sp;
            if (!sp.Init(smarts))
            {
                std::cerr << "Error: Invalid SMARTS pattern: " << smarts << std::endl;
                return {};
            }

            std::vector<std::vector<int>> results;
            if (sp.Match(mol))
            {
                results = sp.GetMapList();
            }
            return results;
        }

    } // namespace utils
} // namespace molgr