import argparse
import json
import os
import platform
import shutil
import sysconfig
from pathlib import Path


DEFAULT_BUILD_DIR = os.environ.get("MOLGR_CMAKE_BUILD_DIR", "build")


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

    try:
        import openbabel

        if hasattr(openbabel, "__file__"):
            ob_module_dir = os.path.dirname(openbabel.__file__)
            ob_include_dir = os.path.join(ob_module_dir, "include")

            if os.path.exists(ob_include_dir):
                ob3_subdir = os.path.join(ob_include_dir, "openbabel3")
                if os.path.exists(ob3_subdir):
                    paths.append(ob3_subdir)
                else:
                    paths.append(ob_include_dir)
    except ImportError:
        pass

    env_include_root = sysconfig.get_paths()["include"]
    base_include = os.path.dirname(env_include_root)

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
        compiler_path = shutil.which("g++") or shutil.which("gcc") or "/usr/bin/g++"
        intelli_sense_mode = f"linux-gcc-{arch}"

    elif system == "Darwin":
        define_os = "Mac"
        compiler_path = shutil.which("clang++") or shutil.which("clang") or "/usr/bin/clang++"
        intelli_sense_mode = f"macos-clang-{arch}"

    elif system == "Windows":
        define_os = "Win32"
        compiler_path = shutil.which("cl") or shutil.which("g++") or shutil.which("gcc") or ""
        if compiler_path and "cl.exe" in compiler_path.lower():
            intelli_sense_mode = f"windows-msvc-{arch}"
        else:
            intelli_sense_mode = f"windows-gcc-{arch}"

    return define_os, compiler_path, intelli_sense_mode


def to_workspace_reference(workspace_dir: Path, path: Path) -> str:
    """将路径尽量转换为 VSCode 可移植的 ${workspaceFolder} 引用。"""
    try:
        relative = path.relative_to(workspace_dir)
    except ValueError:
        return str(path)

    if str(relative) == ".":
        return "${workspaceFolder}"
    return "${workspaceFolder}/" + relative.as_posix()


def resolve_paths(workspace_dir: Path, build_dir_arg: str):
    """解析 compilation database 目录和 compile_commands.json 路径。"""
    raw_path = Path(build_dir_arg)
    resolved_path = raw_path if raw_path.is_absolute() else workspace_dir / raw_path

    if resolved_path.name == "compile_commands.json":
        compile_commands_path = resolved_path
        database_dir = resolved_path.parent
    else:
        database_dir = resolved_path
        compile_commands_path = resolved_path / "compile_commands.json"

    return database_dir.resolve(), compile_commands_path.resolve()


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=4) + "\n", encoding="utf-8")


def write_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def ensure_compile_commands_link(workspace_dir: Path, compile_commands_path: Path):
    """
    在仓库根创建 compile_commands.json 链接。
    clangd 会优先自动发现这里，避免仅靠 IDE 私有配置。
    """
    workspace_link = workspace_dir / "compile_commands.json"
    if not compile_commands_path.exists():
        return False, f"compile_commands not found: {compile_commands_path}"

    try:
        if workspace_link.exists() or workspace_link.is_symlink():
            if workspace_link.resolve() == compile_commands_path.resolve():
                return True, f"compile_commands link already correct: {workspace_link}"
            workspace_link.unlink()

        relative_target = os.path.relpath(compile_commands_path, workspace_link.parent)
        workspace_link.symlink_to(relative_target)
        return True, f"linked compile_commands.json -> {relative_target}"
    except OSError:
        shutil.copy2(compile_commands_path, workspace_link)
        return True, f"copied compile_commands.json -> {workspace_link}"


def generate_config(build_dir: str):
    workspace_dir = Path(__file__).resolve().parent.parent
    database_dir, compile_commands_path = resolve_paths(workspace_dir, build_dir)

    workspace_path = "${workspaceFolder}/**"
    project_include_path = "${workspaceFolder}/src/cpp/include"

    python_includes = get_python_include_path()
    pybind_includes = get_pybind11_include_path()
    openbabel_includes = get_openbabel_include_path()

    raw_paths = (
        [workspace_path, project_include_path]
        + pybind_includes
        + openbabel_includes
        + python_includes
    )

    include_paths = []
    seen = set()
    for p in raw_paths:
        if p and p not in seen:
            include_paths.append(p)
            seen.add(p)

    config_name, compiler_path, mode = detect_compiler_and_mode()
    compile_commands_ref = to_workspace_reference(workspace_dir, compile_commands_path)

    config_data = {
        "configurations": [
            {
                "name": config_name,
                "includePath": include_paths,
                "compileCommands": compile_commands_ref,
                "defines": [],
                "compilerPath": compiler_path,
                "cStandard": "c17",
                "cppStandard": "c++17",
                "intelliSenseMode": mode,
            }
        ],
        "version": 4,
    }

    config_file = workspace_dir / ".vscode" / "c_cpp_properties.json"
    write_json(config_file, config_data)

    clangd_file = workspace_dir / ".clangd"
    clangd_database_path = to_workspace_reference(workspace_dir, database_dir)
    if clangd_database_path.startswith("${workspaceFolder}/"):
        clangd_database_path = clangd_database_path.replace("${workspaceFolder}/", "", 1)
    elif clangd_database_path == "${workspaceFolder}":
        clangd_database_path = "."
    write_text(
        clangd_file,
        "CompileFlags:\n"
        f"  CompilationDatabase: {clangd_database_path}\n",
    )

    link_ok, link_message = ensure_compile_commands_link(workspace_dir, compile_commands_path)

    print("-" * 40)
    print(f"Generated: {config_file}")
    print(f"Updated  : {clangd_file}")
    print(f"Platform : {config_name} ({mode})")
    print(f"Compiler : {compiler_path}")
    print(f"Build dir: {database_dir}")
    print(f"Compile DB: {compile_commands_path}")
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
    print(link_message)

    if not openbabel_includes:
        print("Warning: OpenBabel headers not found automatically.")
        print("Ensure openbabel is installed in this environment.")
    if not link_ok:
        print("Warning: compile_commands.json link was not created.")
        print("Run `cmake -S . -B build` first, or pass `--build-dir <dir>`.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate VSCode/clangd C++ IDE config for MolGR."
    )
    parser.add_argument(
        "--build-dir",
        default=DEFAULT_BUILD_DIR,
        help=(
            "CMake build directory or direct path to compile_commands.json. "
            f"Default: {DEFAULT_BUILD_DIR}"
        ),
    )
    args = parser.parse_args()
    generate_config(args.build_dir)


if __name__ == "__main__":
    main()
