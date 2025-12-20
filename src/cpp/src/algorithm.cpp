/*
 * @Author: TMJ
 * @Date: 2025-12-19 21:53:47
 * @LastEditors: TMJ
 * @LastEditTime: 2025-12-19 21:53:56
 * @Description: 请填写简介
 */
#include "mylib/algorithm.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

namespace mylib {

    int calculate_atom_count(const std::string& smiles) {
        OpenBabel::OBConversion conv;
        OpenBabel::OBMol mol;
        
        conv.SetInFormat("SMILES");
        
        // 这里的逻辑只关乎 C++ 和 OpenBabel，完全不知道 Python 的存在
        if(conv.ReadString(&mol, smiles)) {
            return mol.NumAtoms();
        }
        return -1;
    }

}