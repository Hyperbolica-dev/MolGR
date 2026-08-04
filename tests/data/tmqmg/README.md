# tmQMg fixtures

This directory freezes reviewed tmQMg inputs for offline regression tests.

The source dataset is the official
[`uiocompcat/tmQMg`](https://github.com/uiocompcat/tmQMg) release. Use the
[download and checksum instructions](../../../docs/development/MOLECULE_REVIEW_TOOL.zh-CN.md#获取-tmqmg-数据)
before regenerating fixtures.

- `reconstruction/` contains valid source inputs with published reference graphs.
  Each backend is checked only for equivalence to the final reference result;
  search phases, recovery tiers, candidate counts, and other internal paths are
  deliberately not fixture assertions.
- Boron-cluster systems are intentionally outside this reconstruction fixture
  contract; MolGR makes no support or backend-behavior guarantee for them.
- `source_issues/` contains cases where the tmQMg source or published reference
  is known to be wrong or incomplete. They are tested only for the documented
  issue and must not be used as expected-reference-equivalence successes.
- `fixture_sources.json` is the reviewed selection and provenance source. The
  generated `manifest.csv` files retain tmQMg row indices, charges, reference
  SMILES, and historical classifications, but classifications are not acceptance
  criteria.
- `../reviewed/tmqmg/` is written directly by the manual review client. Its
  `approved_graph/*.sdf` files are authoritative reviewed graphs with original
  coordinates and electronic state; `../reviewed/tmqmg/reference_graph/*.xyz` files use tmQMg
  SMILES and reviewed electronic-state metadata from `manifest.json`.

Reviewed fixture tests normalize only the answer graph before comparison: they
remove RDKit `DATIVE` bonds and clear metal-atom stereochemistry. The stored
evidence remains unchanged, and reconstruction succeeds when
`check_equivalence(..., use_chirality=False)` accepts the normalized answer.

Regenerate the frozen XYZ files and manifests from a local tmQMg checkout with:

```bash
uv run python scripts/freeze_tmqmg_fixtures.py \
  --csv /mnt/e/download/tmQMg_properties_and_targets.csv \
  --xyz-dir /mnt/e/download/tmQMg_xyz/xyz
```

Rebuild the review-driven fixtures from the local review database with:

```bash
uv run python tools/molgr_review/build_fixtures.py \
  --db .local/molgr_review/review.sqlite \
  --fixtures-dir tests/data/reviewed/tmqmg \
  --xyz-dir /mnt/e/download/tmQMg_xyz/xyz
```

Only confirmed review decisions are written to the repository fixture manifest.
Pending and skipped cases remain local review state and are not test fixtures.

Trace the current reviewed fixtures directly from the manifest; no separate
trace case list is maintained:

```bash
uv run python scripts/reconstruction_trace.py \
  --review-fixtures-manifest tests/data/reviewed/tmqmg/manifest.json \
  --fixture-id ABEGOD \
  --format html --out /tmp/abegod-review-trace.html
```

The trace output includes a `review_fixture` synchronization check comparing
the selected trace graph to the approved answer with
`check_equivalence(..., use_chirality=False)`.
