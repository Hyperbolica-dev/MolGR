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

## Isolated Review UI preview

Use the Docker preview launcher to compare the current UI checkout with an
existing Review service without sharing writable review state:

```bash
scripts/review-ui-preview.sh
```

The launcher auto-detects these paths only when they exist:

- `.local/molgr_review/review.sqlite`
- `.local/tmQMg/data/xyz`
- `tests/data/reviewed/tmqmg`

Otherwise, provide them explicitly:

```bash
scripts/review-ui-preview.sh \
  --source-db /path/to/review.sqlite \
  --xyz-dir /path/to/tmQMg/xyz \
  --fixtures-dir /path/to/reviewed/tmqmg \
  --port 8766
```

The source SQLite database is copied with SQLite's backup API, and the fixture
tree is copied into `.local/molgr_review_preview/`. Only those preview copies
are mounted writable in the container; the source XYZ directory is mounted
read-only. Saving reviews in preview therefore cannot update the source DB or
source fixture manifest.

An existing preview is never overwritten implicitly. Stop and resume the same
preview with:

```bash
docker stop molgr-review-ui-preview
docker start molgr-review-ui-preview
```

Use `--fresh` to explicitly discard preview changes, copy the current source
state again, rebuild, and start. Remove the preview container and copied state
with:

```bash
scripts/review-ui-preview.sh --clean
```

The same settings can be supplied with
`MOLGR_REVIEW_PREVIEW_SOURCE_DB`, `MOLGR_REVIEW_PREVIEW_XYZ_DIR`,
`MOLGR_REVIEW_PREVIEW_SOURCE_FIXTURES`, `MOLGR_REVIEW_PREVIEW_PORT`, and
`MOLGR_REVIEW_PREVIEW_DIR`.

This is a local preview workflow, not a production deployment script. A
one-command deployment on another machine still requires confirmation of:

- Docker Engine/Desktop version, daemon access, CPU architecture, and available build resources;
- checkout/update method and the exact revision or release to deploy;
- host locations and ownership/permissions for the production DB, XYZ data, and fixtures;
- production bind address, port, firewall/reverse-proxy, TLS, and authentication requirements;
- backup, restore, retention, and concurrent-writer policy for SQLite and fixtures;
- required package/apt mirrors, outbound network restrictions, and image registry policy;
- service supervision, restart policy, logging, health checks, and rollback procedure.
