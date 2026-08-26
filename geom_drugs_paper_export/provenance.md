# GEOM-Drugs formal benchmark provenance

## Dataset

- Release: GEOM Harvard Dataverse v4, 2022-02-11
- DOI: `10.7910/DVN/JNGTDF`
- Dataverse file ID: `4360331`
- Source archive: `drugs_crude.msgpack.tar.gz`
- Source archive MD5: `7778e84c50b7cde755cca670d1f75091` (verified)
- Source records: 292,035
- Eligible unique molecules: 291,709
- Exclusions: 15 non-neutral/open-shell references; 311 duplicate canonical references; zero
  parse, fragment, XYZ, finite-relative-energy, or reference/XYZ atom-count exclusions
- `xyz2mol_smiles` was not used as ground truth and the featurized release was not acquired.

## Frozen protocol

- Primary unit: one molecule, one deterministic conformer
- Conformer: minimum finite `relativeenergy`, original source index tie-break
- Reference: source/canonical SMILES molecular graph
- Electronic-state scope: neutral, zero-radical, single fragment; charge 0, multiplicity 1
- Reconstruction: `molgr_cpp`, git SHA `e2d3b19afd2f459bc84cd6f40061f709584bd6a8`
- Evaluator: evaluator v1; Candidate = predicted, Reference = reference,
  `use_chirality=False`
- Decisions: `equivalent`, `not_equivalent`, `inconclusive`
- Exact SMILES and chirality are diagnostics only.
- No permissive sulfoxide `S=O` / `[S+][O-]` rule was added.

Exact command:

```bash
.venv/bin/python benchmarks/geom_xyz_benchmark/formal_run.py --archive benchmarks/_data/geom/drugs_crude.msgpack.tar.gz --out benchmarks/geom_xyz_benchmark/_runs/drugs_formal_full --expected-eligible 291709 --case-timeout-seconds 10
```

## Environment and execution

- Python: 3.8.20
- RDKit: 2024.03.5
- Platform: Linux-7.0.11-76070011-generic-x86_64-with-glibc2.2.5
- Wall time: 8138.696427 seconds

`results.csv.gz`, `failures.csv`, `review_cases.csv`, and `summary.json` are byte-for-byte copies
of the frozen evaluated run. Creating this export did not invoke reconstruction or the evaluator.
