# tmQMg XYZ Benchmark

这个 benchmark 用 tmQMg 的 CSV/XYZ 对比较不同重建方法，复用同一套方法注册表。

## 输入

`--csv` 指向 tmQMg 元数据 CSV，至少需要 `id`、`smiles`、`charge` 列。
`--xyz-dir` 指向包含 `<id>.xyz` 文件的目录。

子集筛选参数：

- `--start-row` 和 `--end-row` 按 1-based CSV 行号筛选。
- `--ids` 按 tmQMg id 筛选，可以重复传入或使用逗号分隔。
- `--limit` 在筛选后限制 case 数量。

执行参数：

- `--methods` 限制共享方法注册表。验证 C++/Python 后端对齐时使用
  `molgr_cpp,molgr_fallback`。
- `--cpp-accelerations default|all` 选择 C++ 后端加速 preset。两个 preset 都保持
  target-bucket 并行开启；`all` 会额外开启 vendor UFF atom-typing cache 等可选 C++ 加速项。
- `--enable-uff-atom-typing-cache` 在 default preset 下开启 vendor UFF atom-typing cache。
- `--process-workers N` 把每个方法拆到 `N` 个 worker 子进程中运行。它可以和 C++ 内部
  target-bucket 线程叠加，但过高取值可能竞争 CPU。
- `--case-timeout-seconds` 设置每个方法、每个 case 的 wall-time 限制；传 `0` 表示关闭。

后端对齐检查示例：

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --limit 1000 \
  --out benchmarks/_runs/tmqmg-parity \
  --progress-every 50 \
  --case-timeout-seconds 1.0 \
  --cpp-accelerations all \
  --methods molgr_cpp,molgr_fallback
```

## 输出

benchmark 会在 `--out` 下写出 `results.csv` 和 `summary.csv`。
