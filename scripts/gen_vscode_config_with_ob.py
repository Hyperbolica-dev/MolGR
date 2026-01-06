import json
import os
import platform
import shutil
import sysconfig


def get_python_include_path():
    """获取 Python 头文件路径"""
    paths = sysconfig.get_paths()
    includes = [paths["include"]]
    if "platinclude" in paths and paths["platinclude"] not in includes:
        includes.append(paths["platinclude"])
    return includes


def get_pybind11_include_path():
    """获取 pybind11 头文件路径"""
    try:
        import pybind11

        return [pybind11.get_include()]
    except ImportError:
        return []


def get_openbabel_include_path():
    """
    尝试获取 openbabel 的 include 路径。
    策略：
    1. 检查 Python site-packages/openbabel/include (pip 安装常见)
    2. 检查 sys.prefix/include/openbabel3 (conda 或系统安装常见)
    """
    paths = []

    # 策略 1: 通过 import 也就是 site-packages 内部查找
    try:
        import openbabel

        if hasattr(openbabel, "__file__"):
            # 定位到 .../site-packages/openbabel
            ob_module_dir = os.path.dirname(openbabel.__file__)

            # 检查 .../site-packages/openbabel/include
            ob_include_dir = os.path.join(ob_module_dir, "include")

            if os.path.exists(ob_include_dir):
                # 很多时候里面还有一层 openbabel3
                ob3_subdir = os.path.join(ob_include_dir, "openbabel3")
                if os.path.exists(ob3_subdir):
                    paths.append(ob3_subdir)
                else:
                    paths.append(ob_include_dir)
    except ImportError:
        pass

    # 策略 2: 检查环境根目录 (常见于 Conda 环境或 Linux 系统级安装)
    # 目标通常是 {prefix}/include/openbabel3
    env_include_root = sysconfig.get_paths()["include"]  # 通常是 .../include/python3.8
    # 我们需要往上一级找通用的 include 目录
    base_include = os.path.dirname(env_include_root)  # .../include

    possible_ob_path = os.path.join(base_include, "openbabel3")
    if os.path.exists(possible_ob_path) and possible_ob_path not in paths:
        paths.append(possible_ob_path)

    return paths


def detect_compiler_and_mode():
    """根据系统检测编译器和 IntelliSense 模式"""
    system = platform.system()
    machine = platform.machine().lower()

    arch = "arm64" if "arm" in machine or "aarch64" in machine else "x64"

    compiler_path = ""
    intelli_sense_mode = ""
    define_os = ""

    if system == "Linux":
        define_os = "Linux"
        compiler_path = shutil.which("gcc") or "/usr/bin/gcc"
        intelli_sense_mode = f"linux-gcc-{arch}"

    elif system == "Darwin":
        define_os = "Mac"
        compiler_path = shutil.which("clang") or "/usr/bin/clang"
        intelli_sense_mode = f"macos-clang-{arch}"

    elif system == "Windows":
        define_os = "Win32"
        compiler_path = shutil.which("gcc") or shutil.which("cl") or ""
        if compiler_path and "cl.exe" in compiler_path:
            intelli_sense_mode = f"windows-msvc-{arch}"
        else:
            intelli_sense_mode = f"windows-gcc-{arch}"

    return define_os, compiler_path, intelli_sense_mode


def generate_config():
    # 1. 基础项目路径
    workspace_path = "${workspaceFolder}/**"
    project_include_path = "${workspaceFolder}/src/cpp/include"

    # 2. 自动探测路径
    python_includes = get_python_include_path()
    pybind_includes = get_pybind11_include_path()
    openbabel_includes = get_openbabel_include_path()  # 新增

    # 3. 合并去重
    # 顺序很重要：项目 > 第三方库 > 系统/Python
    raw_paths = (
        [workspace_path, project_include_path]
        + pybind_includes
        + openbabel_includes
        + python_includes
    )

    # 过滤无效路径
    include_paths = []
    seen = set()
    for p in raw_paths:
        if p and p not in seen:
            include_paths.append(p)
            seen.add(p)

    # 4. 获取环境配置
    config_name, compiler_path, mode = detect_compiler_and_mode()

    # 5. 生成 JSON
    config_data = {
        "configurations": [
            {
                "name": config_name,
                "includePath": include_paths,
                "defines": [],
                "compilerPath": compiler_path,
                "cStandard": "c17",
                "cppStandard": "c++14",
                "intelliSenseMode": mode,
            }
        ],
        "version": 4,
    }

    # 6. 写入
    vscode_dir = ".vscode"
    if not os.path.exists(vscode_dir):
        os.makedirs(vscode_dir)

    config_file = os.path.join(vscode_dir, "c_cpp_properties.json")

    with open(config_file, "w", encoding="utf-8") as f:
        json.dump(config_data, f, indent=4)

    # 7. 打印结果摘要
    print("-" * 40)
    print(f"Generated: {config_file}")
    print(f"Platform : {config_name} ({mode})")
    print("-" * 40)
    print("Included Paths:")
    for p in include_paths:
        label = "  "
        if "site-packages" in p:
            label = "[Lib]"
        elif "openbabel" in p:
            label = "[OB ]"
        elif "python" in p:
            label = "[Py ]"
        elif "workspace" in p:
            label = "[Prj]"
        print(f"{label} {p}")
    print("-" * 40)

    if not openbabel_includes:
        print("⚠️ Warning: OpenBabel headers not found automatically.")
        print("   Ensure openbabel is installed in this environment.")


if __name__ == "__main__":
    generate_config()
