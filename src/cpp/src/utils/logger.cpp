/**
 * @file logger.cpp
 * @brief Implementation of logging control.
 */

#include "molgr/utils/logger.h"

// OpenBabel Includes (仅在 .cpp 中引入)
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>

#include <iomanip>

namespace molgr
{

    // 默认日志级别：WARN (生产环境通常只看警告和错误)
    LogLevel g_current_log_level = LogLevel::WARN;

    void SetLogLevel(LogLevel level)
    {
        g_current_log_level = level;
    }

    void LogOmolInfos(OpenBabel::OBMol &mol, const std::string &comment)
    {
        using namespace OpenBabel;

        // 1. Log Comment (if exists)
        if (!comment.empty())
        {
            LOG_DEBUG("[Info] " << comment);
        }

        // 2. Generate SMILES
        // 注意：OBConversion 开销较大，如果频繁调用且有性能要求，考虑将其设为 thread_local static
        OBConversion conv;
        conv.SetOutFormat("smi");
        
        // WriteString 通常会包含换行符，我们需要去掉它以便日志整洁
        std::string smiles = conv.WriteString(&mol);
        if (!smiles.empty() && smiles.back() == '\n') {
            smiles.pop_back();
        }
        // 去除可能的首尾空白
        smiles.erase(0, smiles.find_first_not_of(" \t\r\n"));
        smiles.erase(smiles.find_last_not_of(" \t\r\n") + 1);

        // 3. Log General Info
        LOG_DEBUG("[Info] The omol contains " << mol.NumAtoms() << " atoms "
                  << "and " << mol.NumBonds() << " bonds. smiles: " << smiles);

        // 4. Log Atom Details
        // 对应 Python: f"{idx:03d} ... {symbol:<2s} ... {charge:+d} ... {spin:+d}"

        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            OBAtom *atom = &(*atom_iter);
            std::string symbol = OBElements::GetSymbol(atom->GetAtomicNum());
            
            std::stringstream ss;
            
            // Format: "Atom 001 is C  with formal charge +0 and spin multiplicity +0."
            ss << "[Info] | Atom " 
               << std::setw(3) << std::setfill('0') << atom->GetIdx() 
               << " is "
               << std::left << std::setw(2) << std::setfill(' ') << symbol 
               << " with formal charge " 
               << std::showpos << atom->GetFormalCharge() // showpos 显示 + 号
               << " and radical number " 
               << atom->GetSpinMultiplicity()
               << std::noshowpos; // 恢复不显示 + 号，以免影响后续
               
            LOG_DEBUG(ss.str());
        }
    }

} // namespace molgr
