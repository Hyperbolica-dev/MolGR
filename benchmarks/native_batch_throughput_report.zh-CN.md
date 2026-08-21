# MolGR 原生 Batch 吞吐报告

- 测量日期：2026-08-18
- 代码版本：`ab8db1e` + 工作区本轮 native batch 稳定性修复（未提交）

## 测量条件

- 平台：WSL2 Linux，28 个逻辑 CPU
- Python：3.8.20
- 输入：`tests/data/tmqmg/reconstruction/` 中 6 个 XYZ fixture
- fixture 原子数范围：27–71；1000 条任务按这 6 个 fixture 循环构成，共 51316 个原子
- 每个配置独立 Python 进程运行
- 预热 1 批，正式运行 3 次，表中使用 wall-time 中位数
- C++ batch 使用原生 worker queue，未使用 joblib；`ordered=False`
- 每次均验证 1000/1000 成功，三次原子数 checksum 均一致

可复现实验入口：

```bash
uv run --frozen python benchmarks/native_batch_throughput.py \
  --mode cpp_batch --workers 0 --items 1000 --repeats 3 --warmup 1
```

其中 `--workers 0` 表示由 MolGR 按硬件并行度自动选择 worker；本机实际为 28。

唯一 tmQMg 输入的命令：

```bash
uv run --frozen python benchmarks/native_batch_throughput.py \
  --mode cpp_batch --workers 0 --items 1000 --repeats 3 --warmup 1 \
  --tmqmg-xyz-dir /mnt/e/download/tmQMg_xyz/xyz \
  --tmqmg-csv /mnt/e/download/tmQMg_properties_and_targets.csv
```

## 1000 条固定回归任务结果

| 路径 | worker | 中位耗时 | 吞吐 | 成功率 |
| --- | ---: | ---: | ---: | ---: |
| C++ batch（高层，含 RDKit 转换） | 1 | 5.457 s | 183.2 items/s | 100% |
| C++ batch（高层，含 RDKit 转换） | 4 | 3.151 s | 317.3 items/s | 100% |
| C++ batch（高层，自动 worker） | 自动（28） | 1.539 s | **649.7 items/s** | 100% |
| C++ batch（native 数据层，不含 RDKit 转换） | 4 | 3.093 s | 323.3 items/s | 100% |
| C++ batch（native 数据层，不含 RDKit 转换） | 自动（28） | 0.883 s | **1132.5 items/s** | 100% |

作为 Python 参考路径，60 条任务顺序执行耗时 22.086 s，吞吐 2.72 items/s；该后端按设计不创建 Python 线程池。

## 1000 个唯一 tmQMg 输入

另外使用本机完整 tmQMg 数据目录中按文件名排序的前 1000 个唯一 XYZ 输入；电荷来自
`tmQMg_properties_and_targets.csv`，自旋多重度按总电子数奇偶性设置为 1/2。该组输入共
54907 个原子，包含更广的结构复杂度，因而吞吐低于固定六样本回归组。

| 路径 | worker | 中位耗时 | 处理吞吐 | 成功率 |
| --- | ---: | ---: | ---: | ---: |
| C++ batch（native 数据层） | 自动（28） | 3.744 s | **267.1 items/s** | 998/1000 (99.8%) |
| C++ batch（高层，含 RDKit 转换） | 自动（28） | 3.777 s | **264.8 items/s** | 998/1000 (99.8%) |

两次失败均为输入数据本身无法得到有机图候选，诊断为
`NO_VALID_ORGANIC_CANDIDATE`：`AHIKAD.xyz`（index 884）和 `AHIXEV.xyz`（index 891）。
吞吐按全部已处理任务数计算；如果只按成功结果计算，高层约为 264.3 molecules/s。

## 结论

1. 推荐批量调用使用 C++ backend，并让 `max_workers=None`（自动 worker）。在固定回归样本上，高层端到端吞吐约 **650 molecules/s**；在 1000 个唯一 tmQMg 输入上约 **265 molecules/s**，后者更接近实际复杂度。
2. 固定样本上 native 数据层约 1132 molecules/s，高层约 650 molecules/s；唯一 tmQMg 样本上两者分别约 267 和 265 molecules/s，说明复杂输入的主要成本在 native 重建，RDKit 后处理不是主要瓶颈。
3. 当前 batch worker 使用进程级固定线程池，单分子内部桶并行和 batch 外层并行互斥；不会为每个 batch 重复创建/销毁 Open Babel worker。生产环境仍应按目标机器实测 worker 数，不能直接照搬 28。
4. 无序流不会丢失输入对应关系：高层 `ReconstructionBatchResult` 现在直接返回可解包的 `(input, result, status)` 三元组，且携带原始输入；也可用 `as_pair()` 或 `as_dict()`。因此可以安全使用 `ordered=False`，避免为排序保留额外队列压力。

## 生命周期回归

独立 Python 子进程中执行了 100 轮、每轮 60 条任务，交替覆盖完整消费和首条结果后提前 `close()`；同一回归独立运行 5 次，均正常退出且无 `SIGABRT`、`SIGSEGV` 或 allocator corruption。另验证了两个同时存在的 batch iterator 可以共享 worker 并完整消费。

## 两个 C++ batch 接口的区别

- `core.pipeline.reconstruct_with_metals.batch_xyz2omol` 是底层 native 接口：返回 `MoleculeData`，不包含 RDKit 转换和高层后处理，适合 C++/native 数据链路或性能剖析。
- `molgr.iter_xyz_to_rdmol_batch(..., backend="cpp")` 是实际业务接口：在 native 重建后完成 RDKit 对象转换、后处理和结构化诊断，返回 `ReconstructionBatchResult`。
- 高层 iterator 的每项运行时类型是长度为 3 的 `NamedTuple`，可直接解包为
  `(input, result, status)`；`status=None` 表示成功，失败时为 `ReconstructionDiagnostics`。外层输入接受
  `ReconstructionBatchRequest`、三元 sequence，以及 list/generator 等任意 one-shot iterable。

因此 native 行吞吐更高，但不是端到端用户指标。业务代码最终需要 `rdkit.Chem.Mol` 时应使用高层接口；只有下游明确消费 `MoleculeData` 时才直接调用 native 接口。

## 限制

- 这是一组小型、固定 tmQMg fixture 的吞吐结果，不代表所有分子尺寸、金属数或失败比例。
- 本报告没有把外部 joblib 作为基线复现；此前已知外部进程/线程叠加可能造成过量并发，原生 batch 的目标正是把 worker 预算统一收敛到 C++ 层。
- wall-time 包含 batch 迭代、MoleculeData 到 RDKit 的转换和高层后处理；native 行仅用于定位转换开销，不是用户端到端指标。
