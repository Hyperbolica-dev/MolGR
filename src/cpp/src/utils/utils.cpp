/**
 * @file utils.cpp
 * @brief Implementation of geometric utility functions and Thread-Safe SMARTS caching.
 * @details Uses thread_local storage to allow lock-free parallel execution of SMARTS matching.
 * @namespace molgr::utils
 * @author TMJ
 * @date 2025-12-27
 */

#include "molgr/utils/utils.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/logger.h"

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/parsmart.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <vector>

namespace molgr
{
    namespace utils
    {

        // =============================================================================
        // Thread-Safe SMARTS Caching System
        // =============================================================================

        // 使用 thread_local 保证每个线程有独立的缓存，无需加锁，实现真正的并行加速
        thread_local std::unordered_map<std::string, std::unique_ptr<OpenBabel::OBSmartsPattern>> t_smarts_cache;

        OpenBabel::OBSmartsPattern *GetCompiledSmarts(const std::string &smarts)
        {
            // 1. 查找当前线程的缓存
            auto it = t_smarts_cache.find(smarts);
            if (it != t_smarts_cache.end())
            {
                return it->second.get();
            }

            // 2. 未找到，编译新模式
            auto sp = std::make_unique<OpenBabel::OBSmartsPattern>();
            if (!sp->Init(smarts))
            {
                LOG_ERROR("Invalid SMARTS pattern: " << smarts);
                return nullptr;
            }

            // 3. 存入当前线程的缓存
            OpenBabel::OBSmartsPattern *ptr = sp.get();
            t_smarts_cache[smarts] = std::move(sp);
            return ptr;
        }

        std::vector<std::vector<int>> FindSmarts(OpenBabel::OBMol &mol, const std::string &smarts)
        {
            OpenBabel::OBSmartsPattern *sp = GetCompiledSmarts(smarts);
            if (!sp)
                return {};

            // Match() 修改内部状态，但在 thread_local 下是安全的
            if (sp->Match(mol))
            {
                // 返回结果是 1-based index 的 vector<vector<int>>
                return sp->GetMapList();
            }

            return {};
        }

        // =============================================================================
        // Geometric Functions
        // =============================================================================

        inline double LengthSq(const OpenBabel::vector3 &v)
        {
            return OpenBabel::dot(v, v);
        }

        double CalculateTetrahedronVolume(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                          const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4)
        {
            OpenBabel::vector3 v1 = p1 - p4;
            OpenBabel::vector3 v2 = p2 - p4;
            OpenBabel::vector3 v3 = p3 - p4;
            OpenBabel::vector3 cross_product = OpenBabel::cross(v2, v3);
            double scalar_triple_product = OpenBabel::dot(v1, cross_product);
            return std::abs(scalar_triple_product) / 6.0;
        }

        double CalculateShapeQuality(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                     const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4)
        {
            double volume = CalculateTetrahedronVolume(p1, p2, p3, p4);
            if (volume < 1e-9)
                return 0.0;

            double edges_sq_sum = 0.0;
            edges_sq_sum += LengthSq(p1 - p2);
            edges_sq_sum += LengthSq(p1 - p3);
            edges_sq_sum += LengthSq(p1 - p4);
            edges_sq_sum += LengthSq(p2 - p3);
            edges_sq_sum += LengthSq(p2 - p4);
            edges_sq_sum += LengthSq(p3 - p4);

            double l_rms_squared = edges_sq_sum / 6.0;
            double l_rms_cubed = std::pow(l_rms_squared, 1.5);

            const double NORMALIZATION_CONST = 6.0 * std::sqrt(2.0);
            if (l_rms_cubed < 1e-9)
                return 0.0;
            double quality = NORMALIZATION_CONST * (volume / l_rms_cubed);
            return std::max(0.0, std::min(1.0, quality));
        }

        OpenBabel::vector3 ToVector3(const std::vector<double> &coords)
        {
            if (coords.size() < 3)
                return OpenBabel::vector3(0, 0, 0);
            return OpenBabel::vector3(coords[0], coords[1], coords[2]);
        }

        MoleculeData ExtractMoleculeData(intptr_t mol_ptr)
        {
            OpenBabel::OBMol *mol = reinterpret_cast<OpenBabel::OBMol *>(mol_ptr);
            MoleculeData data;
            data.total_charge = 0;
            data.total_radical_num = 0;

            if (!mol)
            {
                return data;
            }
            return MoleculeDataFromOBMol(*mol);
        }

    } // namespace utils
} // namespace molgr
