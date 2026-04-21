from __future__ import annotations

from openbabel import pybel


def _smarts(pattern: str) -> pybel.Smarts:
    return pybel.Smarts(pattern)


PREPROCESS_DONATE = _smarts("[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]")
PREPROCESS_ACCEPT = _smarts(
    "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]"
)
PRE_CLEAN_HYPERVALENT = _smarts("[Cv5,Nv5,Pv5,Siv5]=,#[*]")
PRE_CLEAN_BCP_RING_5 = _smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
PRE_CLEAN_BCP_RING_4 = _smarts("[#6]1([#6]2)[#7]2[#6]1")
PRE_CLEAN_SI_O_F = _smarts("[Siv5]-[O,F]")

ELIM_HIGH_POSITIVE = _smarts("[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]")
ELIM_CN_IN_DOUBT = _smarts("[#6v4+0]=,#[#7v4+1,#15v4+1]")
ELIM_CARBOXYL = _smarts("[Ov1+0]-C=O")
ELIM_NNN_NEGATIVE = _smarts("[#7v1+0]-[#7v2+0]-[#7v1+0]")
ELIM_NNN_POSITIVE = _smarts("[#7v3+0]-[#7v2+0]-[#7v3+0]")
ELIM_1_3_DIPOLE = _smarts("[*-1]-,=[N+0,O+0]-,=[*]")
ELIM_POSITIVE_N = _smarts("[Nv3+0]=[Nv2+0]")
ELIM_POSITIVE_C_H = _smarts("[#6v3+0,#6v2+0,#1v0+0]")
ELIM_NEGATIVE_C_V3 = _smarts("[#6v3+0]")
ELIM_NEGATIVE_H = _smarts("[#1v0+0]")
ELIM_NEGATIVE_C_LOW = _smarts("[#6v2+0,#6v1+0,#6v0+0]")

CLEAN_CARBENE_NEIGHBOR_UNSAT = _smarts("[*]-[*]=[*]")
CLEAN_RESONANCE_0 = _smarts("[*-]-[*]=[*]~[*+]")
CLEAN_RESONANCE_1 = _smarts("[*-]=[*+]=[*+0]")
CLEAN_RESONANCE_2 = _smarts("[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]")
CLEAN_RESONANCE_3 = _smarts("[#7v2+]=[*]-[*]=[*]-[#8-]")
CLEAN_RESONANCE_4 = _smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
CLEAN_RESONANCE_5 = _smarts("[#7+0,#8+0,#16+0]=[*+0]-[#6-,#7-]")
CLEAN_RESONANCE_6 = _smarts("[#6]=[#6]=[#6-,#7-]")
CLEAN_RESONANCE_7 = _smarts("[*-]1-,:[*](=,:[*])-,:[*]=,:[*]-,:[*]=,:[*]1")
CLEAN_RESONANCE_8 = _smarts("[*-]1-,:[*]=,:[*]-,:[*](=,:[*])-,:[*]=,:[*]1")
CLEAN_RESONANCE_9 = _smarts("[*+,*+2,*+3]-,=[*-,*-2,*-3]")
CLEAN_RESONANCE_10 = _smarts("[*]-[*]=,#[*]-[*]")
CLEAN_RESONANCE_11 = _smarts("[#7v3+0,#8v2+0,#16v2+0]-,=,:[*+1]")
CLEAN_RESONANCE_12 = _smarts("[#7v3+0,#8v2+0,#16v2+0]-,:[*]=,:[*]-,:[*+1]")
CLEAN_RESONANCE_13 = _smarts("[*-]:[*]=[#7+0,#8+0]")

BREAK_DEFORMED_ENE_A = _smarts("[*]~[*+0]=,:[*+0]~[*]")
BREAK_DEFORMED_ENE_B = _smarts("[*]~[*+0](=,:[*+0])~[*]")
BREAK_ONE_BOND_MULTIPLE = _smarts("[*+0]#,=[*+0]")
BREAK_ONE_BOND_CATION = _smarts("[#7+1,#15+1]=[*+0]")
BREAK_ONE_BOND_AROMATIC = _smarts("[*+0]:[*+0]")

RESONANCE_ONE_STEP = _smarts("[*]-,=,:[*]=,#,:[*]")
SCORING_CONJUGATION = _smarts("[*]=,#,:[*]-,:[*]=,#,:[*]")


__all__ = [
    "BREAK_DEFORMED_ENE_A",
    "BREAK_DEFORMED_ENE_B",
    "BREAK_ONE_BOND_AROMATIC",
    "BREAK_ONE_BOND_CATION",
    "BREAK_ONE_BOND_MULTIPLE",
    "CLEAN_CARBENE_NEIGHBOR_UNSAT",
    "CLEAN_RESONANCE_0",
    "CLEAN_RESONANCE_1",
    "CLEAN_RESONANCE_10",
    "CLEAN_RESONANCE_11",
    "CLEAN_RESONANCE_12",
    "CLEAN_RESONANCE_13",
    "CLEAN_RESONANCE_2",
    "CLEAN_RESONANCE_3",
    "CLEAN_RESONANCE_4",
    "CLEAN_RESONANCE_5",
    "CLEAN_RESONANCE_6",
    "CLEAN_RESONANCE_7",
    "CLEAN_RESONANCE_8",
    "CLEAN_RESONANCE_9",
    "ELIM_1_3_DIPOLE",
    "ELIM_CARBOXYL",
    "ELIM_CN_IN_DOUBT",
    "ELIM_HIGH_POSITIVE",
    "ELIM_NEGATIVE_C_LOW",
    "ELIM_NEGATIVE_C_V3",
    "ELIM_NEGATIVE_H",
    "ELIM_NNN_NEGATIVE",
    "ELIM_NNN_POSITIVE",
    "ELIM_POSITIVE_C_H",
    "ELIM_POSITIVE_N",
    "PREPROCESS_ACCEPT",
    "PREPROCESS_DONATE",
    "PRE_CLEAN_BCP_RING_4",
    "PRE_CLEAN_BCP_RING_5",
    "PRE_CLEAN_HYPERVALENT",
    "PRE_CLEAN_SI_O_F",
    "RESONANCE_ONE_STEP",
    "SCORING_CONJUGATION",
]
