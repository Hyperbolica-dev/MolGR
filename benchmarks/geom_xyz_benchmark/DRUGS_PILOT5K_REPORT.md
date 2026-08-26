# GEOM-Drugs deterministic 5k pilot

This pilot ran only 5,000 deterministically sampled GEOM-Drugs molecules with `molgr_cpp`. It is
not a full-dataset run. MolGR reconstruction semantics and evaluator v1 policy were unchanged.

## Acquisition and sampling

- Official source: GEOM Harvard Dataverse v4 (2022-02-11), DOI `10.7910/DVN/JNGTDF`, file ID
  `4360331`, `drugs_crude.msgpack.tar.gz`, verified MD5
  `7778e84c50b7cde755cca670d1f75091`.
- Scanned 292,035 source molecule records. There were 291,709 eligible unique references, 15
  electronic-state exclusions, and 311 duplicate canonical references. All other loader exclusion
  counts were zero.
- The sample is the lowest 5,000 SHA-256 ranks of `seed:source_smiles`, seed 20260825, over all
  eligible unique molecules. Size strata do not enter selection or overall weighting.
- Each molecule uses the minimum finite `relativeenergy` conformer and source index as tie-break.
  All selected minima were 0.0 kcal/mol; 889/5,000 selected conformer indices were nonzero.
- Reference is the source/canonical SMILES graph. All cases use charge 0 and multiplicity 1.
  `xyz2mol_smiles` remains unavailable in the crude artifact and was not downloaded or used.

## Overall results

| Measure | Count |
|---|---:|
| Total | 5,000 |
| Reconstruction success / failure | 5,000 / 0 |
| Equivalent / non-equivalent / inconclusive | 4,942 / 46 / 12 |
| Normalized graph identity / resonance equivalence | 4,940 / 2 |
| Exact non-isomeric SMILES match / mismatch | 4,070 / 930 |
| Chirality-aware equivalent / not equivalent | 3,206 / 1,794 |
| Charge consistent / inconsistent | 5,000 / 0 |
| Radical consistent / inconsistent | 5,000 / 0 |
| Timeout / exception | 0 / 0 |

Runtime in milliseconds: p50 19.294, p95 26.779, p99 30.317, max 239.372. Total measured case
runtime was 99.116 seconds.

Primary evaluation is Candidate = predicted, Reference = reference, evaluator v1,
`use_chirality=False`. Exact SMILES and chirality remain diagnostics only.

## Heavy-atom-size analysis

| Heavy atoms | N | Reconstruction failure | Equivalent | Non-equivalent | Inconclusive |
|---|---:|---:|---:|---:|---:|
| 1–15 | 164 | 0 | 163 | 1 | 0 |
| 16–25 | 2,652 | 0 | 2,624 | 23 | 5 |
| 26–35 | 2,070 | 0 | 2,043 | 20 | 7 |
| 36–50 | 110 | 0 | 109 | 1 | 0 |
| 51+ | 4 | 0 | 3 | 1 | 0 |

The natural population sample contains only four 51+ molecules, so it is valid for overall pilot
weighting but insufficient for a precise large-molecule-specific rate estimate.

## Review queue and protocol assessment

The deterministic review queue contains 65 cases: all 46 non-equivalent, all 12 inconclusive,
the one 51+ failure, 27 cases whose RDKit tautomer-canonical diagnostic agrees, and the top 20
runtime cases (categories overlap). The tautomer diagnostic is prioritization only and does not
override evaluator v1. All 58 primary non-pass cases remain `manual_blocked`; no permissive rule
was added.

All 12 inconclusive cases form a visible sulfoxide representation pattern: source neutral `S=O`
versus predicted charge-separated `[S+][O-]`. Evaluator v1 lacks sufficient approved evidence to
accept this class and correctly leaves it inconclusive. This is a bounded evaluator coverage
limitation, not a loader/protocol failure; formal runs must preserve and report it rather than
silently normalize predictions.

The 46 strict mismatches include a substantial suspected tautomer/protomer subset and one
72-heavy-atom zwitterionic/protonation mismatch already seen in smoke100. No systematic loader,
conformer-selection, reconstruction-exception, timeout, charge, or radical problem was found.
The max-runtime cases are dominated by evaluator-inconclusive representation cases, but remain far
below the 10-second timeout.

STATUS: READY_FOR_FORMAL_RUN
