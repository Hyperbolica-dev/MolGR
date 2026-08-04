from __future__ import annotations

from typing import Dict, List, Optional, cast

from rdkit import Chem
from rdkit.Chem import rdchem


# -----------------------------
# 基础工具
# -----------------------------
def _is_pi_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() in (rdchem.BondType.DOUBLE, rdchem.BondType.TRIPLE)


def _bond_order_int(bt: rdchem.BondType) -> int:
    if bt == rdchem.BondType.SINGLE:
        return 1
    if bt == rdchem.BondType.DOUBLE:
        return 2
    if bt == rdchem.BondType.TRIPLE:
        return 3
    return 0


def _order_to_bondtype(order: int) -> rdchem.BondType:
    if order == 1:
        return rdchem.BondType.SINGLE
    if order == 2:
        return rdchem.BondType.DOUBLE
    if order == 3:
        return rdchem.BondType.TRIPLE
    raise ValueError(f"unsupported bond order: {order}")


def _kekulize_copy(mol: Chem.Mol) -> Chem.Mol:
    m = Chem.Mol(mol)
    # 对芳香体系：先 Kekulize 成显式单双键再做迁移
    try:
        Chem.Kekulize(m, clearAromaticFlags=True)
    except Chem.KekulizeException:
        # Reconstruction may leave aromatic flags inconsistent with the
        # explicit bond orders.  Keep the supplied graph usable and let the
        # resonance walker skip transformations that require explicit pi
        # bonds instead of failing the whole reconstruction.
        m.ClearComputedProps()
        sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        try:
            Chem.SanitizeMol(m, sanitizeOps=sanitize_ops)
        except Exception:  # noqa: BLE001
            m.UpdatePropertyCache(strict=False)
        Chem.KekulizeIfPossible(m, clearAromaticFlags=False)
    return m


def _key(mol: Chem.Mol) -> str:
    # canonical SMILES 会包含自由基与形式电荷信息，用于去重足够好
    return Chem.MolToSmiles(mol, canonical=True)


def _sanitize_reject_valence(m: Chem.Mol) -> bool:
    try:
        Chem.SanitizeMol(m)
        return True
    except Chem.rdchem.AtomValenceException:
        # Explicit valence for atom # ... greater than permitted
        return False
    except Exception:
        return False


# -----------------------------
# 规则 1：一般自由基 1,3 迁移
#   r•-a=b  ->  r=a-b•
# 支持 a=b 为 double / triple（triple->double）
# -----------------------------
def _try_radical_shift_kekule(
    mol_k: Chem.Mol, r_idx: int, a_idx: int, b_idx: int
) -> Optional[Chem.Mol]:
    rw = Chem.RWMol(mol_k)
    r = rw.GetAtomWithIdx(r_idx)
    b = rw.GetAtomWithIdx(b_idx)

    bond_ra = rw.GetBondBetweenAtoms(r_idx, a_idx)
    bond_ab = rw.GetBondBetweenAtoms(a_idx, b_idx)
    if bond_ra is None or bond_ab is None:
        return None

    if r.GetNumRadicalElectrons() < 1:
        return None
    if bond_ra.GetBondType() != rdchem.BondType.SINGLE:
        return None
    if not _is_pi_bond(bond_ab):
        return None

    ra_order = _bond_order_int(bond_ra.GetBondType())
    ab_order = _bond_order_int(bond_ab.GetBondType())
    if ra_order != 1 or ab_order not in (2, 3):
        return None

    # 键级迁移：r-a +1；a-b -1
    bond_ra.SetBondType(_order_to_bondtype(ra_order + 1))  # single->double
    bond_ab.SetBondType(_order_to_bondtype(ab_order - 1))  # double->single / triple->double

    # 自由基迁移：r -> b（只移动 1 个自由基电子）
    r.SetNumRadicalElectrons(r.GetNumRadicalElectrons() - 1)
    b.SetNumRadicalElectrons(b.GetNumRadicalElectrons() + 1)

    new_m = rw.GetMol()
    if not _sanitize_reject_valence(new_m):
        return None
    return new_m


# -----------------------------
# 规则 2：同原子多自由基（>=2） -> 分离电荷贡献式
#   x(••) - a=b  ->  x=a-b
#   并生成：x(+1), b(-1)，同时消耗 2 个自由基电子（闭壳层化）
# -----------------------------
def _try_diradical_charge_sep_shift_kekule(
    mol_k: Chem.Mol, x_idx: int, a_idx: int, b_idx: int, positive: bool = True
) -> Optional[Chem.Mol]:
    rw = Chem.RWMol(mol_k)
    x = rw.GetAtomWithIdx(x_idx)
    b = rw.GetAtomWithIdx(b_idx)

    bond_xa = rw.GetBondBetweenAtoms(x_idx, a_idx)
    bond_ab = rw.GetBondBetweenAtoms(a_idx, b_idx)
    if bond_xa is None or bond_ab is None:
        return None

    nrad = x.GetNumRadicalElectrons()
    if nrad < 2:  # 改为 nrad != 2 可只处理“恰好二自由基”
        return None

    if bond_xa.GetBondType() != rdchem.BondType.SINGLE:
        return None
    if not _is_pi_bond(bond_ab):
        return None

    xa_order = _bond_order_int(bond_xa.GetBondType())
    ab_order = _bond_order_int(bond_ab.GetBondType())
    if xa_order != 1 or ab_order not in (2, 3):
        return None

    # 键级迁移：x-a +1；a-b -1
    bond_xa.SetBondType(_order_to_bondtype(xa_order + 1))
    bond_ab.SetBondType(_order_to_bondtype(ab_order - 1))

    # “同原子 +/- 对”分离成：x(+1), b(-1) 或 x(-1), b(+1)
    # 同时消耗 2 个自由基电子：让中心更接近闭壳层贡献式（避免生成大量不合理超价式）
    x.SetNumRadicalElectrons(nrad - 2)
    if positive:
        x.SetFormalCharge(x.GetFormalCharge() + 1)
        b.SetFormalCharge(b.GetFormalCharge() - 1)
    else:
        x.SetFormalCharge(x.GetFormalCharge() - 1)
        b.SetFormalCharge(b.GetFormalCharge() + 1)

    new_m = rw.GetMol()
    if not _sanitize_reject_valence(new_m):
        return None
    return new_m


# -----------------------------
# 一步扩展：把两类规则合并
# -----------------------------
def _one_step_unified(
    mol: Chem.Mol,
    enable_radical_shift: bool = True,
    enable_multi_radical_charge_sep: bool = True,
) -> List[Chem.Mol]:
    mol_k = _kekulize_copy(mol)
    out: List[Chem.Mol] = []

    for c_idx in range(mol_k.GetNumAtoms()):
        c = mol_k.GetAtomWithIdx(c_idx)
        nrad = c.GetNumRadicalElectrons()

        # 统一在 "c - a = b" 的三原子路径上做变换
        for a in c.GetNeighbors():
            a_idx = cast(rdchem.Atom, a).GetIdx()
            bond_ca = mol_k.GetBondBetweenAtoms(c_idx, a_idx)
            if bond_ca is None or bond_ca.GetBondType() != rdchem.BondType.SINGLE:
                continue

            for b in cast(rdchem.Atom, a).GetNeighbors():
                b_idx = cast(rdchem.Atom, b).GetIdx()
                if b_idx == c_idx:
                    continue
                bond_ab = mol_k.GetBondBetweenAtoms(a_idx, b_idx)
                if bond_ab is None or not _is_pi_bond(bond_ab):
                    continue

                # 规则 1：一般自由基（>=1）迁移
                if enable_radical_shift and nrad >= 1:
                    m2 = _try_radical_shift_kekule(mol_k, c_idx, a_idx, b_idx)
                    if m2 is not None:
                        out.append(m2)

                # 规则 2：同原子多自由基（>=2）分离电荷
                if enable_multi_radical_charge_sep and nrad >= 2:
                    m3 = _try_diradical_charge_sep_shift_kekule(
                        mol_k, c_idx, a_idx, b_idx, positive=True
                    )
                    if m3 is not None:
                        out.append(m3)
                    m3_negative = _try_diradical_charge_sep_shift_kekule(
                        mol_k, c_idx, a_idx, b_idx, positive=False
                    )
                    if m3_negative is not None:
                        out.append(m3_negative)

    return out


# -----------------------------
# 统一入口：BFS 多步枚举 + 去重
# -----------------------------
def enumerate_resonance_radical(
    mol: Chem.Mol,
    depth: int = 3,
    enable_radical_shift: bool = True,
    enable_multi_radical_charge_sep: bool = True,
) -> List[Chem.Mol]:
    if depth < 0:
        raise ValueError("depth must be >= 0")

    start = Chem.Mol(mol)
    if not _sanitize_reject_valence(start):
        return []

    seen: Dict[str, Chem.Mol] = {_key(start): start}
    frontier = [start]

    for _ in range(depth):
        nxt: List[Chem.Mol] = []
        for m in frontier:
            for m2 in _one_step_unified(
                m,
                enable_radical_shift=enable_radical_shift,
                enable_multi_radical_charge_sep=enable_multi_radical_charge_sep,
            ):
                k = _key(m2)
                if k not in seen:
                    seen[k] = m2
                    nxt.append(m2)
        frontier = nxt
        if not frontier:
            break

    return list(seen.values())


def enumerate_resonance_radical_from_smiles(
    smiles: str,
    depth: int = 3,
    enable_radical_shift: bool = True,
    enable_multi_radical_charge_sep: bool = True,
) -> List[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    res = enumerate_resonance_radical(
        mol,
        depth=depth,
        enable_radical_shift=enable_radical_shift,
        enable_multi_radical_charge_sep=enable_multi_radical_charge_sep,
    )
    return sorted({_key(m) for m in res})
