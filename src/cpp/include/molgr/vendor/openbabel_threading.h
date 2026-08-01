#pragma once

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <openbabel/ring.h>

namespace molgr
{
    namespace vendor
    {
        namespace openbabel_threading
        {
            void ConnectTheDotsAndPerceiveBondOrders(OpenBabel::OBMol &mol);
            void ConnectTheDotsLocal(OpenBabel::OBMol &mol);
            void PerceiveBondOrdersLocal(OpenBabel::OBMol &mol);
            void FindRingAtomsAndBonds(OpenBabel::OBMol &mol);
            void PrepareForSmartsMatching(OpenBabel::OBMol &mol);
            void EnsureHybridizationPerceived(OpenBabel::OBMol &mol);
            void SetAromaticPerceived(OpenBabel::OBMol &mol, bool perceived);
            void ResetAndAssignAromaticFlags(OpenBabel::OBMol &mol);
            bool AtomIsAromatic(OpenBabel::OBAtom &atom);
            bool BondIsAromatic(OpenBabel::OBBond &bond);
            bool RingIsAromatic(OpenBabel::OBRing &ring);
        }
    }
}
