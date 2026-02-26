#pragma once

#include <openbabel/math/vector3.h>
#include <openbabel/mol.h>

#include <string>
#include <vector>

#define FOR_NB_OF_ATOM(a, p) for (OpenBabel::OBAtomAtomIter a(p); a; ++a)

namespace molgr
{
    namespace utils
    {
        double CalculateTetrahedronVolume(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                          const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4);

        double CalculateShapeQuality(const OpenBabel::vector3 &p1, const OpenBabel::vector3 &p2,
                                     const OpenBabel::vector3 &p3, const OpenBabel::vector3 &p4);

        OpenBabel::vector3 ToVector3(const std::vector<double> &coords);

        std::vector<std::vector<int>> FindSmarts(OpenBabel::OBMol &mol, const std::string &smarts);

        struct AtomData
        {
            int atomic_num;
            int formal_charge;
            int radical_num;
            double x, y, z;
        };

        struct BondData
        {
            int begin_atom_idx;
            int end_atom_idx;
            int order;
        };

        struct MoleculeData
        {
            std::vector<AtomData> atoms;
            std::vector<BondData> bonds;
            int total_charge;
            int total_radical_num;
        };

        MoleculeData ExtractMoleculeData(intptr_t mol_ptr);
        OpenBabel::OBMol MolFromMoleculeData(const MoleculeData &data);
        MoleculeData MoleculeDataFromOBMol(const OpenBabel::OBMol &mol);

    }
}
