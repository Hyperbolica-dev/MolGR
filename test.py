import molgr._core as core
from openbabel import pybel


def test_pybel_to_cpp_extraction():
    # 1. 在 Python 中创建一个 Pybel 分子
    # 这个对象是由 Python (SWIG) 管理的
    pybel_mol = pybel.readstring("smi", "c1ccccc1[C+]C(C=C)C#C")
    pybel_mol.addh()
    pybel_mol.make3D()
    # 2. 获取它的底层 C++ 指针地址
    # .OBMol 是 openbabel.OBMol 对象
    # .this 通常是一个字符串 (如 "_0x123456_p_OpenBabel__OBMol") 或直接是 int
    # 使用 int() 强转可以兼容这两种情况拿到纯整数地址
    ptr_addr = int(pybel_mol.OBMol.this)

    # 3. 传给你的 C++ 扩展进行数据提取
    # C++ 并不在乎这个指针是谁创建的，只要它是有效的 OBMol 地址
    mol_data = core.extract_molecule_data(ptr_addr)

    # 4. 验证结果
    print(f"Atoms count from C++: {len(mol_data.atoms)}")
    print(f"Total charge from C++: {mol_data.total_charge}")

    # 验证原子信息
    for atom in mol_data.atoms:
        print(f"Atom: {atom.atomic_num} at ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})")
    for bond in mol_data.bonds:
        print(f"Bond: {bond.begin_atom_idx} - {bond.end_atom_idx} type {bond.order}")

    # 注意：千万不要调用 core.free_obmol_ptr(ptr_addr)！
    # 因为这个对象是 Python 创建的，Python 的垃圾回收器(GC)会负责释放它。
    # 只有 core.reconstruct_... 返回的指针才需要手动 free。


test_pybel_to_cpp_extraction()