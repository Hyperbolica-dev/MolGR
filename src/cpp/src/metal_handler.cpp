/**
 * @file metal_handler.cpp
 * @brief Implementation of metal handling logic.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/metal_handler.h"
#include "molgr/consts.h"
#include "molgr/logger.h"
#include "molgr/utils.h"

#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

namespace molgr
{
    namespace metal
    {

        using namespace OpenBabel;

        void CartesianProduct(
            const std::vector<std::vector<MetalAtomPosition>> &input_pools,
            std::vector<MetalAtomPosition> &current_combo,
            std::vector<std::vector<MetalAtomPosition>> &results,
            size_t depth,
            int total_radical_limit)
        {
            if (depth == input_pools.size())
            {
                int current_total_radical = 0;
                for (const auto &m : current_combo)
                    current_total_radical += m.radical_num;

                if (current_total_radical <= total_radical_limit)
                {
                    results.push_back(current_combo);
                }
                return;
            }

            for (const auto &item : input_pools[depth])
            {
                current_combo.push_back(item);
                CartesianProduct(input_pools, current_combo, results, depth + 1, total_radical_limit);
                current_combo.pop_back();
            }
        }

        MetalHandler::MetalHandler(OBMol &mol)
        {
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom *atom = &(*atom_iter);
                if (atom->IsMetal())
                {
                    raw_metals_.push_back({static_cast<int>(atom->GetIdx()),
                                           static_cast<int>(atom->GetAtomicNum()),
                                           atom->GetX(), atom->GetY(), atom->GetZ()});
                }
            }
        }

        std::string MetalHandler::StripMetals(OBMol &mol)
        {
            std::vector<OBAtom *> to_delete;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom *atom = &(*atom_iter);
                if (atom->IsMetal())
                {
                    to_delete.push_back(atom);
                }
            }

            // OpenBabel 的 DeleteAtom 是安全的，只要我们不重复删除
            for (OBAtom *atom : to_delete)
            {
                mol.DeleteAtom(atom);
            }

            OBConversion conv;
            conv.SetOutFormat("xyz");
            return conv.WriteString(&mol);
        }

        std::vector<std::vector<MetalAtomPosition>> MetalHandler::GenerateCombinations(int total_radical_electrons)
        {
            std::vector<std::vector<MetalAtomPosition>> all_pools;

            for (const auto &raw : raw_metals_)
            {
                std::vector<MetalAtomPosition> possible_states;
                std::string symbol = OBElements::GetSymbol(raw.atomic_num);

                // 使用 vector 保持插入顺序 (Priority First)
                std::vector<int> ordered_valences;
                // 使用 set 辅助去重
                std::set<int> seen_valences;

                // 定义一个 lambda 辅助函数来按顺序添加价态
                auto add_valences = [&](const std::vector<int> &source)
                {
                    for (int val : source)
                    {
                        if (seen_valences.find(val) == seen_valences.end())
                        {
                            ordered_valences.push_back(val);
                            seen_valences.insert(val);
                        }
                    }
                };

                // 1. 先添加 Prior (主要价态)
                if (kMetalValencePrior.count(symbol))
                {
                    add_valences(kMetalValencePrior.at(symbol));
                }

                // 2. 后添加 Minor (次要价态)
                if (kMetalValenceMinor.count(symbol))
                {
                    add_valences(kMetalValenceMinor.at(symbol));
                }

                // 3. 按照保留下来的顺序生成状态
                // 这里的 val 顺序就是 Prior + Minor，未经过数值排序
                for (int val : ordered_valences)
                {
                    std::set<int> radicals = GetPossibleMetalRadicals(symbol, val);
                    // radicals 本身是 set，这就已经按数值排序了，但电子数的顺序通常影响较小
                    // 如果 Python 中的 get_possible_metal_radicals 返回的顺序也很重要，
                    // 那么 Consts.cpp 中的 GetPossibleMetalRadicals 也需要改为返回 vector。
                    // 但根据物理意义，同一个价态下的不同自旋态优先级差异不大，这里先维持 set。
                    for (int rad : radicals)
                    {
                        possible_states.push_back({raw.idx, symbol, raw.atomic_num, val, rad, raw.x, raw.y, raw.z});
                    }
                }

                // 如果没有找到任何可能的价态，回退到默认
                if (possible_states.empty())
                {
                    possible_states.push_back({raw.idx, symbol, raw.atomic_num, 0, 0, raw.x, raw.y, raw.z});
                }

                all_pools.push_back(possible_states);
            }

            std::vector<std::vector<MetalAtomPosition>> results;
            std::vector<MetalAtomPosition> current;

            if (!all_pools.empty())
                CartesianProduct(all_pools, current, results, 0, total_radical_electrons);
            else
                results.push_back({}); // Handle edge case of no metals? Though unlikely here.

            return results;
        }

        void MetalHandler::CombineMetalWithMol(OBMol &mol, const std::vector<MetalAtomPosition> &metals)
        {
            mol.BeginModify(); // 这是一个好习惯，通知 OB 我们要大改结构了

            int num_organic = mol.NumAtoms();
            int num_metals = static_cast<int>(metals.size());
            int total_atoms = num_organic + num_metals;

            // 1. Append Metals
            for (const auto &m : metals)
            {
                OBAtom *atom = mol.NewAtom();
                atom->SetAtomicNum(m.element_idx);
                atom->SetFormalCharge(m.valence);
                atom->SetSpinMultiplicity(m.radical_num);
                atom->SetVector(m.x, m.y, m.z);
            }

            // 2. Build Renumber Map
            std::vector<int> new_order(total_atoms, 0);
            bool error_flag = false;

            // A. Place Metals
            for (int i = 0; i < num_metals; ++i)
            {
                int current_idx = num_organic + 1 + i; // 1-based index
                int target_slot = metals[i].idx - 1;   // 0-based slot

                if (target_slot >= 0 && target_slot < total_atoms)
                {
                    if (new_order[target_slot] != 0)
                    {
                        LOG_ERROR("[MetalHandler] Index collision at slot " << target_slot);
                        error_flag = true;
                    }
                    new_order[target_slot] = current_idx;
                }
                else
                {
                    LOG_ERROR("[MetalHandler] Original index out of bounds: " << metals[i].idx);
                    error_flag = true;
                }
            }

            // B. Place Organic
            int current_organic_idx = 1;
            for (int i = 0; i < total_atoms; ++i)
            {
                if (new_order[i] == 0)
                {
                    if (current_organic_idx <= num_organic)
                    {
                        new_order[i] = current_organic_idx;
                        current_organic_idx++;
                    }
                    else
                    {
                        LOG_ERROR("[MetalHandler] Not enough organic atoms to fill slots.");
                        error_flag = true;
                    }
                }
            }

            // C. Final Safety Check
            for (int idx : new_order)
            {
                if (idx == 0)
                {
                    LOG_ERROR("[MetalHandler] Invalid 0 index in renumber map. Aborting renumber.");
                    error_flag = true;
                    break;
                }
            }

            // 3. Execute Renumber
            if (!error_flag)
            {
                // DEBUG: Print the vector to see what we are passing
                // std::stringstream ss;
                // for(int k : new_order) ss << k << " ";
                // LOG_DEBUG("Renumbering with: " << ss.str());

                mol.RenumberAtoms(new_order);
            }

            mol.EndModify(); // 完成修改
        }

    } // namespace metal
} // namespace molgr