#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        bool EliminateNNN(OpenBabel::OBMol &mol, int &current_charge_deficit, bool positive = false);
        bool EliminateHighPositiveChargeAtoms(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminateCNInDoubt(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminateCarboxyl(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminateCarbeneNeighborHeteroatom(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool Eliminate13DipolePostive(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminatePossibleCPLikeRadicalAnion(
            OpenBabel::OBMol &mol,
            int &current_charge_deficit,
            int total_radical_electrons);
        bool EliminatePositiveCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminateNegativeCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
    }
}
