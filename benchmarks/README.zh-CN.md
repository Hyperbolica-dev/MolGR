# Benchmarks

[English](README.md) | [中文](README.zh-CN.md)

本目录包含 benchmark 入口、benchmark 专用依赖和 benchmark 说明。Benchmark 依赖与根项目依赖图隔离，
不会影响正式发布包的依赖解析。

## 可用 Benchmark

- `smiles_xyz_benchmark`：从 SMILES 输入生成 XYZ case，并比较不同重建方法。
- `molfile_xyz_benchmark`：读取 `.mol`、`.molfile` 和 `.sdf` 输入，转换为 XYZ case，
  并运行同一套方法注册表。
- `tmqmg_xyz_benchmark`：对 tmQMg CSV/XYZ 对运行同一套方法注册表，并支持按行号或 id
  选择子集。

所有 benchmark 入口都使用同一套共享方法注册表：

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`
- `molgr_cpp`
- `xyzgraph_cheminf_full`

`cell2mol_v2` 和 `xyzgraph_cheminf_full` 是竞品基线。`molgr_fallback` 使用 Python 参考后端，
`molgr_cpp` 使用默认 C++ 后端。

## 环境

- 项目运行时：Python `>=3.8`
- Benchmark 运行时：Python `>=3.10,<3.12`
- 推荐包管理器：`uv`
- Benchmark 依赖文件：[`pyproject.toml`](pyproject.toml)
- Benchmark lock 文件：[`uv.lock`](uv.lock)
- 默认 benchmark 虚拟环境：`.venv-benchmark`
- 默认 benchmark Python 可执行文件：`python3.10`

Benchmark 运行时范围比项目运行时更窄，因为 `xyzgraph`、`cell2mol` 和 `cosymlib`
等可选对比栈有较旧的依赖约束。Benchmark 依赖集也保留 `numpy<2`，以兼容这些栈。

## 快速开始

在仓库根目录创建或刷新专用 benchmark 环境：

```bash
bash scripts/benchmark_env.sh create
```

运行一个小规模 SMILES/XYZ benchmark：

```bash
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/smiles-demo
```

运行一个小规模 Molfile/SDF benchmark：

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

运行一个带子集筛选的小规模 tmQMg benchmark：

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --start-row 1 \
  --end-row 100 \
  --ids ABC123,DEF456 \
  --limit 10 \
  --out benchmarks/_runs/tmqmg-demo
```

验证后端行为对齐时，tmQMg 只需要比较 MolGR 的两个后端：

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

`tmqmg_xyz_benchmark` 默认使用 `--process-workers 1` 串行运行；C++ 单分子调用仍保留
内部 target-bucket 并行。对 `molgr_cpp` 指定更大的值时，单个子进程会把 worker 数交给
native batch 池，避免外部进程和 C++ 线程叠加。批量吞吐仍应按目标机器实测 worker 数。

如果需要连续运行 benchmark 命令，可以把当前 shell 切到 benchmark project：

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

然后使用普通 `uv run` 命令：

```bash
uv run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/smiles-demo
```

恢复默认项目环境：

```bash
unset UV_PROJECT UV_PROJECT_ENVIRONMENT UV_PYTHON
```

## 输入

SMILES/XYZ benchmark：

- `--input` 指向包含 SMILES 的文本或 CSV 风格文件。
- 接受名为 `smiles`、`canonicalsmiles` 或 `canonicalsmi` 的表头。
- 每个 SMILES 会由 RDKit 嵌入、UFF 优化、转换为 XYZ，并与原始 RDKit 分子比较。

Molfile/SDF benchmark：

- `--input` 接受单个 `.mol`、`.molfile` 或 `.sdf` 文件，也接受递归扫描的目录。
- 推荐 fixture 根目录：[`../tests/data/sdf/`](../tests/data/sdf/)。
- 支持 [`../tests/data/sdf/`](../tests/data/sdf/) 下的嵌套 fixture 分类。
- 每个输入结构会由 RDKit 读取、转换为 XYZ case，并与源分子比较。

共享参数：

- `--limit`：限制 case 数量，便于快速检查。
- `--out`：输出运行目录。
- `--methods`：tmQMg 入口可用，用 method id 限制共享方法注册表。
- `--process-workers`：tmQMg 串行/native batch worker 数；`1` 使用单分子串行路径，
  C++ 更大取值使用单个 native batch 池。
- `--cpp-accelerations`：tmQMg 入口可用，选择 C++ 加速 preset。

## 输出

每次运行会在 `--out` 目录下写入两个主文件：

- `results.csv`：每个 `(case, method)` 尝试一行，包含状态、预测、等价性结果和计时字段。
- `summary.csv`：按方法聚合的指标，包括数量和延迟统计。

`results.csv` 的主要列包括：

- `case_idx`
- `method_id`
- `input_smiles`
- `ground_truth_smiles`
- `status`
- `error`
- `predicted_smiles`
- `equivalent`
- `equivalence_method`
- `timing_ms_total`
- `method_ms`
- `equivalence_ms`
- `timing_ms_breakdown_json`

出现的计时列会被展开。`timing_ms_breakdown_json` 保留完整的方法专用计时字典，便于 profiling。

## 结果解读

- 使用 `summary.csv` 做方法级比较。
- 使用 `results.csv` 做逐 case 调试和失败分析。
- `equivalent` 为 `True` 的行计为成功。
- 方法报错或 `equivalent` 为 `False` 的行计为失败。
- 输入 case 准备失败的行计为 skipped。

## 可选 `cell2mol_v2` 说明

`cell2mol_v2` 是可选 benchmark 基线。Benchmark 依赖文件指向上游 `cell2mol` 的 `v2`
Git revision。如果启用或重新分发这个基线，请确认用法符合其许可证条款，包括可能的 GPL 义务。

## 子 Benchmark 文档

- [`smiles_xyz_benchmark/README.zh-CN.md`](smiles_xyz_benchmark/README.zh-CN.md)
- [`smiles_xyz_benchmark/README.md`](smiles_xyz_benchmark/README.md)
- [`molfile_xyz_benchmark/README.zh-CN.md`](molfile_xyz_benchmark/README.zh-CN.md)
- [`molfile_xyz_benchmark/README.md`](molfile_xyz_benchmark/README.md)
- [`tmqmg_xyz_benchmark/README.zh-CN.md`](tmqmg_xyz_benchmark/README.zh-CN.md)
- [`tmqmg_xyz_benchmark/README.md`](tmqmg_xyz_benchmark/README.md)
