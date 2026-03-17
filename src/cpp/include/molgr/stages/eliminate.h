#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        void EliminateNNN(OpenBabel::OBMol &mol, int &current_charge_deficit, bool positive = false);
        void EliminateHighPositiveChargeAtoms(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminateCNInDoubt(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminateCarboxyl(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminateCarbeneNeighborHeteroatom(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminateChargeSpliting(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void Eliminate13Dipole(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminatePositiveCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
        void EliminateNegativeCharges(OpenBabel::OBMol &mol, int &current_charge_deficit);
    }
}
