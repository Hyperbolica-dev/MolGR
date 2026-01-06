/**
 * @file metal_handler.h
 * @brief Handles metal atom removal, state combination generation, and restoration.
 * @author TMJ
 * @date 2025-12-27
 */

#pragma once

#include <openbabel/mol.h>
#include <string>
#include <vector>
#include <memory>

namespace molgr
{
    namespace metal
    {

        /**
         * @brief Represents the state and original position of a metal atom.
         * Corresponds to Python: MetalAtomPosition dataclass
         */
        struct MetalAtomPosition
        {
            int idx;            // 原始原子索引 (1-based)
            std::string symbol; // 元素符号
            int element_idx;    // 原子序数
            int valence;        // 化合价 (Formal Charge)
            int radical_num;    // 自旋多重度/单电子数
            double x, y, z;     // 坐标
        };

        /**
         * @brief Helper class to manage metal stripping and recombination.
         */
        class MetalHandler
        {
        public:
            /**
             * @brief Initialize with the raw molecule from XYZ.
             * Identifying and storing metal information.
             */
            explicit MetalHandler(OpenBabel::OBMol &mol);

            /**
             * @brief Remove all metal atoms from the molecule.
             * @return A string containing the XYZ block of the organic part (no metals).
             */
            std::string StripMetals(OpenBabel::OBMol &mol);

            /**
             * @brief Generate all valid combinations of metal states.
             * Corresponds to Python's itertools.product logic inside xyz2omol.
             * @param total_radical_electrons Max allowed total radical electrons for filtering.
             * @return A vector of combinations, where each combination is a vector of MetalAtomPosition.
             */
            std::vector<std::vector<MetalAtomPosition>> GenerateCombinations(int total_radical_electrons);

            /**
             * @brief Add metals back to an organic molecule.
             * Corresponds to Python: combine_metal_with_omol
             */
            static void CombineMetalWithMol(OpenBabel::OBMol &mol, const std::vector<MetalAtomPosition> &metals);

        private:
            // 存储从原始分子中提取的金属信息（只包含坐标和基本属性，不包含具体的价态组合）
            struct RawMetalInfo
            {
                int idx;
                int atomic_num;
                double x, y, z;
            };
            std::vector<RawMetalInfo> raw_metals_;
        };

    } // namespace metal
} // namespace molgr