# GEOM-Drugs deterministic 100-molecule smoke report

This is a protocol smoke test, not a benchmark accuracy claim. It ran only `molgr_cpp`; no 5k or
full run was started and no MolGR reconstruction semantics were changed.

## Acquisition and sampling

- Source: official GEOM Harvard Dataverse v4 (2022-02-11), DOI
  `10.7910/DVN/JNGTDF`, file ID `4360331`, `drugs_crude.msgpack.tar.gz`, verified MD5
  `7778e84c50b7cde755cca670d1f75091`.
- Full archive index scanned: 292,035 source molecule records.
- Eligible unique neutral, zero-radical, single-fragment references: 291,709.
- Excluded: 15 electronic-state exclusions and 311 duplicate canonical references. Zero records
  were excluded for reference parse failure, fragmentation, missing XYZ, missing/invalid relative
  energy, or reference/XYZ atom-count mismatch.
- Primary unit: molecule. Sampling uses the lowest SHA-256 score of
  `seed:source_smiles` within fixed heavy-atom strata, with seed 20260825; it is independent of
  archive traversal order and is not the first 100 records.
- Heavy-atom distribution: 1–15: 15; 16–25: 25; 26–35: 25; 36–50: 20; 51+: 15. Observed range:
  11–77 heavy atoms.
- One conformer per molecule was chosen by minimum finite `relativeenergy`, with original source
  index as the deterministic tie-break. All selected minima were 0.0 kcal/mol; nonzero selected
  indices (including conformer 252) confirm that source order was not assumed to be energy order.

The source/canonical SMILES is the reference graph. `xyz2mol_smiles` is unavailable in the crude
artifact and was neither reconstructed nor substituted; its provenance/integrity diagnostic is
therefore recorded as unavailable. Acquiring the much larger featurized artifact would be needed
to add that diagnostic.

## Smoke results

| Measure | Count |
|---|---:|
| Reconstruction success / failure | 100 / 0 |
| Equivalent / non-equivalent / inconclusive | 97 / 3 / 0 |
| Normalized graph identity / resonance equivalence | 97 / 0 |
| Exact non-isomeric SMILES match / mismatch | 85 / 15 |
| Chirality-aware equivalent / not equivalent | 56 / 44 |
| Charge-consistent / radical-consistent | 100 / 100 |
| Timeout / exception | 0 / 0 |

Primary evaluation used Candidate = prediction, Reference = source SMILES graph, evaluator v1,
and `use_chirality=False`. Exact SMILES and chirality are secondary diagnostics. The 44
chirality-aware mismatches support retaining chirality outside the primary metric: source SMILES
may omit stereo while a 3D conformer realizes it. Exact-SMILES mismatches include representation
differences that evaluator v1 correctly classifies as normalized graph identity.

The three strict primary failures are preserved without repair:

1. `geom-drugs-830922aa8d367f23-c0` (22 heavy atoms): tautomer/proton-position and carbonyl bonding
   differ from the source reference.
2. `geom-drugs-c1de119f3ab63c9f-c252` (72 heavy atoms): a large peptide-like molecule reconstructs
   with zwitterionic/protonation differences despite matching total charge and radicals.
3. `geom-drugs-63d1d0fe86197812-c0` (75 heavy atoms): a large peptide-like molecule has localized
   protonation/bonding differences, again with matching total charge and radicals.

## Drugs-specific findings and disposition

The crude format was readable across the complete molecule index and supplied source SMILES,
conformer coordinates, relative energies, and total energies consistently for all otherwise
eligible unique references. No systemic loader, reference, conformer-selection, timeout, or C++
runner issue was observed. The relevant Drugs-specific risk is localized protomer/tautomer
behavior in flexible, functional-group-rich molecules, especially very large peptide-like cases.
These must remain strict failures unless independent provenance supports a reference/conformer
mismatch; evaluator permissiveness must not be introduced.

Before a 5k pilot, the only useful human decision is whether to acquire the featurized release to
audit `xyz2mol_smiles` for integrity/provenance on failures. This is not a blocker for the stated
source-SMILES primary protocol.

STATUS: READY_FOR_5K
