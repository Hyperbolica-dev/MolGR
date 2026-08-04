#pragma once

#include <cstdint>

#include <openbabel/math/vector3.h>
#include <openbabel/mol.h>

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

        struct AtomData
        {
            int atomic_num;
            int formal_charge;
            int radical_num;
            int lone_pair_count = 0;
            bool unresolved_two_electron_center = false;
            int hybridization = 0;
            double x, y, z;
        };

        struct BondData
        {
            int begin_atom_idx;
            int end_atom_idx;
            int order;
            bool aromatic = false;
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
