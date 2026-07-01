# MolGR

[English](README.md) | [中文](README.zh-CN.md)

MolGR 名称来源于 Moleculer Graph Reconstructer。

MolGR 是一个从 XYZ 坐标重建分子图的 Python 包。它接收 XYZ 文本、总电荷和自旋多重度，
输出带有键级、三维构象、可选配位键和可选立体化学信息的 `rdkit.Chem.Mol`。

当前项目有两套可切换后端：

- `cpp`：默认后端，通过 `molgr._core` 暴露 C++ 加速实现，适合常规使用和性能敏感场景。
- `python`：位于 [`src/molgr/fallback/`](src/molgr/fallback/) 的 Python 参考实现，用于语义对齐、调试和回归测试。

## 快速使用

主要入口是 `molgr.interface.xyz_to_rdmol`：

```python
from molgr.interface import xyz_to_rdmol

xyz = """3
water
O 0.000000 0.000000 0.000000
H 0.758602 0.000000 0.504284
H -0.758602 0.000000 0.504284
"""

mol = xyz_to_rdmol(
    xyz,
    total_charge=0,
    spin_multiplicity=1,
    backend="cpp",  # 也可以使用 "python"
)
```

`spin_multiplicity` 会在入口处转换为 `total_radical_electrons = spin_multiplicity - 1`。
默认会补充可能的配位键和立体化学信息；需要关闭时可传入
`make_dative_bonds=False` 或 `make_stereochemistry=False`。

## 配置

运行时配置基于 dataclass。默认使用包级全局对象
[`molgr.config.CONFIG`](src/molgr/config.py)；调用方可以直接修改这个全局配置对象，也可以向
`xyz_to_rdmol(..., config=...)` 传入独立的 [`MolGRConfig`](src/molgr/config.py)。

```python
from molgr.config import CONFIG

CONFIG.cpp_backend.enable_target_bucket_parallelism = True
CONFIG.cpp_backend.target_bucket_parallel_threshold = 1
CONFIG.cpp_backend.target_bucket_parallel_max_threads = None
```

C++ 后端是 Python fallback 语义的默认加速实现。C++ 专属开关可以改变调度、缓存或线程安全
vendor 实现，但同一个 `MolGRConfig` 下不能改变最终入选分子。

## 安装

MolGR 需要 Python `>=3.8`。运行时依赖和构建依赖以 [`pyproject.toml`](pyproject.toml) 为准；
默认开发工作流使用 `uv`。

### C++ 编译环境要求

从源码构建会编译 `molgr._core` 扩展。本地编译环境需要提供：

- CMake `>=3.15`。
- 支持 C++17 的编译器：
  - Linux：GCC 9+ 或 Clang 10+。
  - macOS：当前 Xcode Command Line Tools 提供的 Apple Clang。
  - Windows：MSVC Build Tools 2019+，并安装 C++ workload。
- 当前 Python 解释器对应的 Python development headers。使用 `uv` 创建虚拟环境时，
  通常由所选 Python 安装提供。
- [`pyproject.toml`](pyproject.toml) 中声明的 Python 构建包：`scikit-build-core`、
  `pybind11`、`setuptools-scm` 和 `openbabel-wheel`。
- 目标平台可用的 OpenBabel wheel。CMake 会从已安装的 `openbabel-wheel` 包中定位
  OpenBabel 头文件和库文件。

从源码安装并构建 C++ 扩展：

```bash
cd MolGR
uv sync --dev
uv pip install -e . -v --no-build-isolation
```

可编辑安装会通过 `scikit-build-core`、`pybind11` 和当前 Python 环境中的 OpenBabel 构建
C++ 扩展。如果修改了 C++ 代码或 pybind11 绑定，运行相关测试前需要重新构建：

```bash
make cpp-build
```

如果需要模拟 CI 中直接使用 PyPI 的依赖解析方式，而不是使用 [`pyproject.toml`](pyproject.toml) 中配置的本地镜像，
可以显式指定 index：

```bash
uv sync --dev --index https://pypi.org/simple --default-index https://pypi.org/simple
uv pip install -e . -v --no-build-isolation \
  --index https://pypi.org/simple \
  --default-index https://pypi.org/simple
```

## 项目结构

- [`src/molgr/interface.py`](src/molgr/interface.py)：公开 API，负责配置解析、后端选择和 RDKit 后处理。
- [`src/molgr/fallback/`](src/molgr/fallback/)：Python 参考实现。
- [`src/cpp/`](src/cpp/)：C++ 重建实现。
- [`src/bindings/`](src/bindings/)：`molgr._core` 的 pybind11 绑定。
- [`src/molgr/utils/`](src/molgr/utils/)：转换、等价性检查和 RDKit 后处理工具。
- [`tests/`](tests/)：公开行为测试、fallback 测试、转换器测试和 C++ 一致性测试。
- [`benchmarks/`](benchmarks/)：benchmark 入口、专用依赖和运行说明。
- [`docs/`](docs/)：架构说明和发布流程文档。

## 算法概览

MolGR 把 Python fallback 作为语义参考，把 C++ 后端作为加速路径。两套后端遵循同一套
高层重建逻辑：

1. 解析 XYZ 输入，并根据总电荷、自旋多重度和运行时配置确定重建目标。
2. 根据 `backend="cpp"` 或 `backend="python"` 分发到对应后端。
3. 对不含金属的结构，执行有机骨架重建流程；如果直接候选无效，则进入共振候选恢复。
4. 对含金属的结构，先把金属从有机骨架中临时分离，枚举可能的金属状态，并把这些状态按照
   诱导出的有机部分电荷和自由基目标分桶。
5. 每个有机目标桶只重建一次，然后复用于对应的金属候选组合。
6. 对有机候选或“金属 + 有机骨架”组合候选评分，选择最优状态。
7. 将后端结果转换为 RDKit 分子，并按需补充配位键和立体化学信息。

更完整的调用图、金属搜索、共振恢复和 C++ 优化说明见
[`docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md`](docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md)。

## 开发

创建或刷新开发环境：

```bash
uv sync --dev
```

以可编辑模式构建扩展：

```bash
uv pip install -e . -v --no-build-isolation
```

常用检查命令：

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy src
uv run pytest
```

构建源码包和 wheel：

```bash
uv build
```

根目录 [`makefile`](makefile) 提供常用开发快捷命令：

- `make format`：运行 Ruff formatter，会修改文件。
- `make lint`：运行 `ruff check . --fix`，会修改文件。
- `make test`：运行 pytest。
- `make type-check`：对 `src` 运行 MyPy。
- `make cpp-build`：重新构建可编辑 C++ 扩展，并刷新 C++ IDE 配置。
- `make stubs`：重新生成编译扩展的类型存根。

`make stubs` 会对 `molgr._core` 运行 `pybind11-stubgen`，保留公开 `pipeline`
stubs 以及 [`src/molgr/_core/dev/`](src/molgr/_core/dev/) 下的开发/一致性测试 helper
stubs，然后格式化并 lint 生成的 `.pyi` 文件。

可选：根据当前环境生成 VSCode/clangd C++ 配置：

```bash
uv run python scripts/gen_vscode_config_with_ob.py --build-dir build
```

## 测试

运行完整测试：

```bash
uv run pytest
```

开发时可以运行更小范围的测试：

```bash
uv run pytest tests/test_main.py
uv run pytest tests/test_backend_reconstruction_regression.py
uv run pytest tests/test_fallback_resonance.py
uv run pytest tests/test_fallback_get_possible_metal_radicals.py
```

修改 C++、pybind11 绑定或 `_core` 公开接口时，推荐流程是：

1. 使用 `make cpp-build` 重新构建扩展。
2. 运行相关的 C++/Python 后端一致性测试，例如
   [`tests/test_backend_reconstruction_regression.py`](tests/test_backend_reconstruction_regression.py)
   [`tests/test_cpp_python_metal_candidate_parity.py`](tests/test_cpp_python_metal_candidate_parity.py)、
   [`tests/test_cpp_uff_atom_typing_cache.py`](tests/test_cpp_uff_atom_typing_cache.py)
   和 [`tests/test_force_field_scoring_policy.py`](tests/test_force_field_scoring_policy.py)。
3. 如果绑定接口变化，运行 `make stubs`，并检查生成的 `.pyi` 文件。

C++ 后端是 Python fallback 语义的加速实现。修改 SMARTS 匹配、UFF setup 状态、
共振候选选择或 C++ 专属加速开关时，需要同步检查
[`docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md`](docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md#cpppython-后端一致性护栏)
中的后端一致性护栏。

部分测试依赖 OpenBabel；如果当前环境缺少对应依赖，这些测试可能会被跳过。

## Benchmark

Benchmark 依赖与主项目依赖分离在 [`benchmarks/pyproject.toml`](benchmarks/pyproject.toml) 中；它们不参与发布包的依赖解析，
也不参与主 [`uv.lock`](uv.lock)。benchmark 使用独立环境，默认路径是 `.venv-benchmark`，默认 Python 是
`python3.10`。

创建或刷新 benchmark 环境：

```bash
bash scripts/benchmark_env.sh create
```

运行 SMILES/XYZ benchmark：

```bash
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/smiles-demo
```

运行 Molfile/SDF benchmark：

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

运行一个 tmQMg 后端对齐子集：

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --limit 1000 \
  --out benchmarks/_runs/tmqmg-parity \
  --case-timeout-seconds 1.0 \
  --cpp-accelerations all \
  --methods molgr_cpp,molgr_fallback
```

连续运行 benchmark 命令时，也可以切换当前 shell 环境：

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

恢复默认项目环境：

```bash
unset UV_PROJECT UV_PROJECT_ENVIRONMENT UV_PYTHON
```

更多说明见 [`benchmarks/README.md`](benchmarks/README.md) 或
[`benchmarks/README.zh-CN.md`](benchmarks/README.zh-CN.md)。

## 文档导航

- [`docs/README.md`](docs/README.md)：英文文档索引。
- [`docs/README.zh-CN.md`](docs/README.zh-CN.md)：中文文档索引。
- [`docs/architecture/ALGORITHM_ARCHITECTURE.md`](docs/architecture/ALGORITHM_ARCHITECTURE.md)：英文算法架构说明。
- [`docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md`](docs/architecture/ALGORITHM_ARCHITECTURE.zh-CN.md)：中文算法架构说明。
- [`docs/release/DEVELOPMENT_RELEASE_GUIDE.md`](docs/release/DEVELOPMENT_RELEASE_GUIDE.md)：英文开发与发布指南。
- [`docs/release/DEVELOPMENT_RELEASE_GUIDE.zh-CN.md`](docs/release/DEVELOPMENT_RELEASE_GUIDE.zh-CN.md)：中文开发与发布指南。
- [`benchmarks/README.md`](benchmarks/README.md)：英文 benchmark 环境、运行命令和输出说明。
- [`benchmarks/README.zh-CN.md`](benchmarks/README.zh-CN.md)：中文 benchmark 环境、运行命令和输出说明。

## 协议

MolGR 使用 GNU General Public License v2 发布。完整协议内容见 [`LICENSE`](LICENSE)。
