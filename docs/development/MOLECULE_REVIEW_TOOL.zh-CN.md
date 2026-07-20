# 分子图重建审核与 Trace 工具

`tools/molgr_review/` 用于审核 MolGR 的分子图重建结果。服务读取 XYZ、电子态、
参考 SMILES 和候选结果，提供二维/三维结构对照与完整 Trace，并将人工确认的答案
固化为回归测试 fixture。

审核服务有意使用当前工作区中的 MolGR，而不是独立镜像或已发布版本。这样可以保证
页面展示、Trace 和 fixture 生成均对应正在开发的 Python/C++ 实现。数据集相关逻辑
则限制在队列生成脚本中；服务本身只读取统一格式的审核队列。

## 数据与状态

| 路径 | 内容 | 是否提交 |
| --- | --- | --- |
| `.local/molgr_review/review.sqlite` | 当前审核状态 | 否 |
| `.local/molgr_review/tmqmg/` | tmQMg benchmark 结果与审核队列 | 否 |
| `.local/datasets/tmqmg/` | tmQMg 原始数据 | 否 |
| `tests/data/reviewed/` | 已确认的标准答案 | 是 |

SQLite、队列 CSV、Trace、benchmark 输出和 JSONL 备份均属于本地开发状态。
仓库只保存能够作为测试标准答案的 fixture。

## 准备项目环境

首次使用时安装开发依赖并构建当前 C++ 扩展：

```bash
uv sync --dev
uv pip install -e . -v --no-build-isolation
```

Python 代码变化后重启审核服务即可。修改 `src/cpp` 或 `src/bindings` 后，需要重新
执行可编辑安装以更新扩展模块。

服务启动时会验证 `molgr` 的导入路径，并比较 C++ 扩展与源码的修改时间。若当前
解释器加载的不是本工作区，或扩展早于源码，服务会拒绝启动。

## 获取 tmQMg 数据

tmQMg 数据来自 [uiocompcat/tmQMg](https://github.com/uiocompcat/tmQMg)，对应论文
为 [D2DD00129B](https://doi.org/10.1039/D2DD00129B)。为保证 benchmark 和 fixture
可复现，当前固定使用提交 `e1dc9887b8f20a217a1db6ca972d726bcbaab45b`。

下载属性表与 XYZ 数据，并验证文件校验和：

```bash
TMQMG_REV=e1dc9887b8f20a217a1db6ca972d726bcbaab45b
TMQMG_DIR="$PWD/.local/datasets/tmqmg/$TMQMG_REV"
mkdir -p "$TMQMG_DIR"

curl --fail --location --retry 3 \
  "https://raw.githubusercontent.com/uiocompcat/tmQMg/$TMQMG_REV/data/tmQMg_properties_and_targets.csv" \
  --output "$TMQMG_DIR/tmQMg_properties_and_targets.csv"
curl --fail --location --retry 3 \
  "https://raw.githubusercontent.com/uiocompcat/tmQMg/$TMQMG_REV/data/tmQMg_xyz.zip" \
  --output "$TMQMG_DIR/tmQMg_xyz.zip"

(
cd "$TMQMG_DIR"
sha256sum --check <<'EOF'
3920c1c8f4ec81bc8e44b8d0256a7da1e36c8805c3c0adfd47e50c46e633f473  tmQMg_properties_and_targets.csv
e0d15a70bcba294717cd9f9792e7fac99ef0c5c61c3a6e08dcc8a8643f53660a  tmQMg_xyz.zip
EOF
)

unzip -q "$TMQMG_DIR/tmQMg_xyz.zip" -d "$TMQMG_DIR/tmQMg_xyz"
```

固定版本应包含 74,547 个 XYZ 文件：

```bash
find "$TMQMG_DIR/tmQMg_xyz/xyz" -maxdepth 1 -name '*.xyz' | wc -l
```

## 生成 tmQMg 审核队列

`prepare_tmqmg_queue.py` 遍历 tmQMg、执行 MolGR benchmark，并将需要人工判断的
case 写为审核队列：

```bash
uv run python tools/molgr_review/prepare_tmqmg_queue.py \
  --csv "$TMQMG_DIR/tmQMg_properties_and_targets.csv" \
  --xyz-dir "$TMQMG_DIR/tmQMg_xyz/xyz" \
  --cpp-accelerations all
```

默认队列路径为 `.local/molgr_review/tmqmg/tmqmg_cases.csv`。该命令只生成队列，
不会读取或修改审核数据库。

局部验证可以使用 `--ids`、`--limit`、`--start-row` 和 `--end-row`。完整参数见：

```bash
uv run python tools/molgr_review/prepare_tmqmg_queue.py --help
```

## 导入或更新审核队列

```bash
uv run python tools/molgr_review/import_cases.py \
  --input .local/molgr_review/tmqmg/tmqmg_cases.csv \
  --db .local/molgr_review/review.sqlite
```

重复导入会更新 case 数据。仍存在于新队列中的审核状态会保留；不再出现的 case
及其本地审核状态会被移除。需要跨机器迁移尚未固化的审核状态时，可以使用
`export_reviews.py` 导出 JSONL，再通过 `--reviews-jsonl` 恢复。

## 启动审核服务

```bash
uv run python tools/molgr_review/server.py \
  --db .local/molgr_review/review.sqlite \
  --xyz-dir "$TMQMG_DIR/tmQMg_xyz/xyz" \
  --fixtures-dir tests/data/reviewed/tmqmg \
  --port 8765
```

服务地址为 `http://127.0.0.1:8765`。Trace 在独立页面展示当前 case 的完整重建
过程，默认最多渲染 1000 张 `rdkit-dof` 图像。

## 审核状态与 fixture

| 审核决定 | 持久化结果 |
| --- | --- |
| `accept_candidate` | 将当前 MolGR 图保存为 `approved_graph/<id>.sdf` |
| `accept_reference` | 保存原始 XYZ、电子态和参考 SMILES |
| `manual_reference` | 保存原始 XYZ、电子态和人工指定的 SMILES |
| `needs_followup` | 保留本地审核状态，不生成 fixture |
| `skip` | 保留本地审核状态，不生成 fixture |

前三类决定会立即更新 `tests/data/reviewed/<数据集>/manifest.json`。manifest 只允许
`approved_graph`、`reference_graph` 和 `manual_reference`，并且不保存审核人、
审核时间、临时备注或未确认结论。

将 case 改为 `needs_followup` 或 `skip` 时，工具会移除该 case 已有的 fixture，
避免本地审核状态与测试答案不一致。

## 审核其他数据集

审核服务不要求输入来自 tmQMg。CSV 至少需要以下字段：

| 字段 | 含义 |
| --- | --- |
| `case_id` | case 的唯一标识 |
| `xyz_path` | XYZ 文件路径 |

建议同时提供 `total_charge`、`total_radical_electrons` 和
`spin_multiplicity`。`reference_smiles` 用于提供参考答案，`candidate_smiles` 用于
记录生成队列时的候选快照。其他字段会保存在 `metadata_json` 中，供诊断页面使用。

tmQMg 的字段解析与 benchmark 汇总仅存在于 `prepare_tmqmg_queue.py`。接入其他
数据集时，应单独生成上述 CSV，而不是修改审核数据库 schema。

## 验证与提交

```bash
uv run ruff check tools/molgr_review
uv run pyright tools/molgr_review tests/test_tmqmg_review_fixtures.py
uv run pytest -q tools/molgr_review/tests
uv run pytest -q tests/test_tmqmg_review_fixtures.py
```

提交范围包括审核工具源码、开发文档和已确认的 fixture。`.local/` 中的审核状态、
数据集、队列和运行产物不应进入版本控制。
