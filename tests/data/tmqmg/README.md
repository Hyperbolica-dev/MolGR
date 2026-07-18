# tmQMg fixtures

This directory freezes reviewed tmQMg inputs for offline regression tests.

- `reconstruction/` contains valid source inputs selected to exercise distinct
  no-metal reconstruction paths. These cases are eligible for backend parity
  and reconstruction assertions.
- `source_issues/` contains cases where the tmQMg source or published reference
  is known to be wrong or incomplete. They are tested only for the documented
  issue and must not be used as expected-reference-equivalence successes.
- `fixture_sources.json` is the reviewed selection and classification source of
  truth. The generated `manifest.csv` files retain tmQMg row indices, charges,
  reference SMILES, classifications, and reasons.
- `reviewed/` is written directly by the manual review client. Its
  `approved_graph/*.sdf` files are authoritative reviewed graphs with original
  coordinates and electronic state; `tmqmg_reference/*.xyz` files use tmQMg
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
uv run python .local/tmqmg_review/build_fixtures.py \
  --xyz-dir /mnt/e/download/tmQMg_xyz/xyz
```

Trace the current reviewed fixtures directly from the manifest; no separate
trace case list is maintained:

```bash
uv run python scripts/reconstruction_trace.py \
  --review-fixtures-manifest tests/data/tmqmg/reviewed/manifest.json \
  --fixture-id ABEGOD \
  --format html --out /tmp/abegod-review-trace.html
```

The trace output includes a `review_fixture` synchronization check comparing
the selected trace graph to the approved answer with
`check_equivalence(..., use_chirality=False)`.
