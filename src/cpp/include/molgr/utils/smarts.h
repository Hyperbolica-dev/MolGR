#pragma once

#include <openbabel/mol.h>

#include <cstddef>
#include <vector>

namespace molgr
{
    namespace smarts
    {
        enum class PatternId : std::size_t
        {
            PREPROCESS_DONATE = 0,
            PREPROCESS_ACCEPT,
            PRE_CLEAN_HYPERVALENT,
            PRE_CLEAN_BCP_RING_5,
            PRE_CLEAN_BCP_RING_4,
            PRE_CLEAN_SI_O_F,
            ELIM_HIGH_POSITIVE,
            ELIM_CN_IN_DOUBT,
            ELIM_CARBOXYL,
            ELIM_NNN_NEGATIVE,
            ELIM_NNN_POSITIVE,
            ELIM_1_3_DIPOLE,
            ELIM_POSITIVE_N,
            ELIM_POSITIVE_C_H,
            ELIM_NEGATIVE_C_V3,
            ELIM_NEGATIVE_H,
            ELIM_NEGATIVE_C_LOW,
            ELIM_NEGATIVE_CP,
            CLEAN_CARBENE_NEIGHBOR_UNSAT,
            CLEAN_RESONANCE_0,
            CLEAN_RESONANCE_1,
            CLEAN_RESONANCE_2,
            CLEAN_RESONANCE_3,
            CLEAN_RESONANCE_4,
            CLEAN_RESONANCE_5,
            CLEAN_RESONANCE_6,
            CLEAN_RESONANCE_7,
            CLEAN_RESONANCE_8,
            CLEAN_RESONANCE_9,
            CLEAN_RESONANCE_10,
            CLEAN_RESONANCE_11,
            CLEAN_RESONANCE_12,
            CLEAN_RESONANCE_13,
            CLEAN_RESONANCE_14,
            CLEAN_RESONANCE_15,
            CLEAN_RESONANCE_16,
            BREAK_DEFORMED_ENE_A,
            BREAK_DEFORMED_ENE_B,
            BREAK_ONE_BOND_MULTIPLE,
            BREAK_ONE_BOND_CATION,
            BREAK_ONE_BOND_AROMATIC,
            RESONANCE_ONE_STEP,
            SCORING_CONJUGATION,
            COUNT,
        };

        std::vector<std::vector<int>> Match(OpenBabel::OBMol &mol, PatternId pattern_id);
    }
}
