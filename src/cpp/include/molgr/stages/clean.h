#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        bool CleanCarbeneNeighborUnsaturated(OpenBabel::OBMol &mol);
        bool CleanNeighborRadicals(
            OpenBabel::OBMol &mol,
            int given_charge,
            int total_radical_electrons);
        bool Clean14Radicals(
            OpenBabel::OBMol &mol,
            int given_charge,
            int total_radical_electrons);
        bool Clean16Radicals(
            OpenBabel::OBMol &mol,
            int given_charge,
            int total_radical_electrons);
        bool CleanPossible13Dipole(
            OpenBabel::OBMol &mol,
            int given_charge,
            int total_radical_electrons);
        bool CleanResonances14(OpenBabel::OBMol &mol);
        bool CleanResonances16(OpenBabel::OBMol &mol);
        bool CleanResonances17(OpenBabel::OBMol &mol);
        bool CleanResonances18(OpenBabel::OBMol &mol);
        bool CleanResonances(OpenBabel::OBMol &mol);
    }
}
