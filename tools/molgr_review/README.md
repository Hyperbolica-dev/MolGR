# MolGR molecule review

`tools/molgr_review` provides human validation of MolGR reconstruction results.
It executes the current checkout, renders the input, candidate, reference, and
trace, and persists confirmed answers as regression fixtures.

## tmQMg workflow

Generate the review queue:

```bash
uv run python tools/molgr_review/prepare_tmqmg_queue.py \
  --csv /path/to/tmQMg_properties_and_targets.csv \
  --xyz-dir /path/to/tmQMg_xyz/xyz
```

The preparation command updates both `tmqmg_cases.csv` and `review.sqlite` by
default. `--ids`, row ranges, and limits merge only the refreshed scope into the
existing queue. Pass `--no-sync-review-db` only when producing artifacts for a
database that must remain unchanged.

For an external or manually generated complete queue, import explicitly:

```bash
uv run python tools/molgr_review/import_cases.py \
  --input .local/molgr_review/tmqmg/tmqmg_cases.csv \
  --db .local/molgr_review/review.sqlite
```

Start the service:

```bash
uv run python tools/molgr_review/server.py \
  --db .local/molgr_review/review.sqlite \
  --xyz-dir /path/to/tmQMg_xyz/xyz \
  --fixtures-dir tests/data/reviewed/tmqmg \
  --port 8765
```

The review database and generated queue remain under `.local/`. Confirmed
answers are written to `tests/data/reviewed/` and are the only review artifacts
intended for version control.

The input CSV requires `case_id` and `xyz_path`. Electronic-state fields and
reference/candidate SMILES should be supplied when available; additional fields
are retained as diagnostic metadata.

See the [Chinese development guide](../../docs/development/MOLECULE_REVIEW_TOOL.zh-CN.md)
for data provenance, runtime checks, decision semantics, and validation commands.
