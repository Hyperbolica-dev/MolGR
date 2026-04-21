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
        bool EliminateChargeSpliting(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool Eliminate13Dipole(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminatePositiveCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
        bool EliminateNegativeCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
    }
}
