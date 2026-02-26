/**
 * @file metal_handler.cpp
 * @brief Implementation of metal handling logic.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/pipeline/reconstruct_without_metals.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/utils.h"

#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
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
            int current_radical_sum,
            int total_radical_limit)
        {
            if (depth == input_pools.size())
            {
                results.push_back(current_combo);
                return;
            }

            for (const auto &item : input_pools[depth])
            {
                const int next_radical_sum = current_radical_sum + item.radical_num;
                if (next_radical_sum > total_radical_limit)
                {
                    continue;
                }
                current_combo.push_back(item);
                CartesianProduct(input_pools, current_combo, results, depth + 1, next_radical_sum, total_radical_limit);
                current_combo.pop_back();
            }
        }

        OpenBabel::OBConversion &ThreadLocalXyzOutConversion()
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetOutFormat("xyz");
                initialized = true;
            }
            return conv;
        }

        OpenBabel::OBConversion &ThreadLocalXyzInConversion()
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetInFormat("xyz");
                initialized = true;
            }
            return conv;
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

            OpenBabel::OBConversion &conv = ThreadLocalXyzOutConversion();
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
                CartesianProduct(all_pools, current, results, 0, 0, total_radical_electrons);
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

namespace molgr
{
    namespace pipeline
    {
        namespace reconstruct_with_metals
        {
            std::set<int> get_possible_metal_radicals(const std::string &metal_symbol, int valence)
            {
                return molgr::GetPossibleMetalRadicals(metal_symbol, valence);
            }

            std::vector<molgr::metal::MetalAtomPosition> build_metal_states(const OpenBabel::OBAtom &obatom)
            {
                const int atomic_num = static_cast<int>(obatom.GetAtomicNum());
                const std::string symbol = OpenBabel::OBElements::GetSymbol(atomic_num);

                const auto default_state = [&]()
                {
                    return molgr::metal::MetalAtomPosition{
                        static_cast<int>(obatom.GetIdx()),
                        symbol,
                        atomic_num,
                        0,
                        0,
                        obatom.GetX(),
                        obatom.GetY(),
                        obatom.GetZ()};
                };

                std::vector<int> valences;
                std::set<int> seen_valences;

                const auto add_valences = [&](const std::vector<int> &source)
                {
                    for (const int valence : source)
                    {
                        if (seen_valences.insert(valence).second)
                        {
                            valences.push_back(valence);
                        }
                    }
                };

                if (kMetalValencePrior.count(symbol))
                {
                    add_valences(kMetalValencePrior.at(symbol));
                }
                if (kMetalValenceMinor.count(symbol))
                {
                    add_valences(kMetalValenceMinor.at(symbol));
                }
                if (valences.empty())
                {
                    valences.push_back(0);
                }

                if (!kMetalFDSP.count(symbol))
                {
                    return {default_state()};
                }

                std::vector<molgr::metal::MetalAtomPosition> states;
                for (const int valence : valences)
                {
                    const auto radicals = get_possible_metal_radicals(symbol, valence);
                    for (const int radical_num : radicals)
                    {
                        states.push_back(
                            molgr::metal::MetalAtomPosition{
                                static_cast<int>(obatom.GetIdx()),
                                symbol,
                                atomic_num,
                                valence,
                                radical_num,
                                obatom.GetX(),
                                obatom.GetY(),
                                obatom.GetZ()});
                    }
                }

                if (states.empty())
                {
                    return {default_state()};
                }
                return states;
            }

            void combine_metal_with_omol(OpenBabel::OBMol &mol, const std::vector<molgr::metal::MetalAtomPosition> &metals)
            {
                molgr::metal::MetalHandler::CombineMetalWithMol(mol, metals);
            }

            std::unique_ptr<molgr::utils::MoleculeData> Xyz2OmolMolData(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons)
            {
                molgr::pipeline::perf::ResetRunTimingBreakdown();

                if (total_radical_electrons < 0)
                {
                    return nullptr;
                }

                OpenBabel::OBMol mol;
                OpenBabel::OBConversion &conv = molgr::metal::ThreadLocalXyzInConversion();
                if (!conv.ReadString(&mol, xyz_block))
                {
                    return nullptr;
                }

                molgr::metal::MetalHandler handler(mol);
                const std::string no_metal_xyz = handler.StripMetals(mol);

                const auto metal_enum_started = std::chrono::steady_clock::now();
                const auto metal_combinations = handler.GenerateCombinations(total_radical_electrons);
                const auto metal_enum_now = std::chrono::steady_clock::now();
                const double metal_enum_ms = std::chrono::duration<double, std::milli>(metal_enum_now - metal_enum_started).count();
                molgr::pipeline::perf::AddMetalEnumerationCombinationMs(metal_enum_ms);

                bool has_best = false;
                double best_score = 0.0;
                std::unique_ptr<molgr::utils::MoleculeData> best_data;

                for (const auto &metal_states : metal_combinations)
                {
                    int metal_charge = 0;
                    int metal_radical = 0;
                    for (const auto &state : metal_states)
                    {
                        metal_charge += state.valence;
                        metal_radical += state.radical_num;
                    }

                    if (metal_radical > total_radical_electrons)
                    {
                        continue;
                    }

                    const auto organic_data = molgr::pipeline::reconstruct_without_metals::XyzToMolDataNoMetal(
                        no_metal_xyz,
                        total_charge - metal_charge,
                        total_radical_electrons - metal_radical);
                    if (!organic_data)
                    {
                        continue;
                    }

                    OpenBabel::OBMol candidate = molgr::utils::MolFromMoleculeData(*organic_data);
                    molgr::metal::MetalHandler::CombineMetalWithMol(candidate, metal_states);

                    const double score = molgr::scoring::OmolScore(candidate);
                    if (!has_best || score < best_score)
                    {
                        has_best = true;
                        best_score = score;
                        best_data = std::make_unique<molgr::utils::MoleculeData>(
                            molgr::utils::MoleculeDataFromOBMol(candidate));
                    }
                }

                if (!has_best)
                {
                    return nullptr;
                }
                return best_data;
            }
        }
    }
}
