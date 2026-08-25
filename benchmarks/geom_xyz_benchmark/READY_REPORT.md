# GEOM benchmark READY audit

## Source and version

- Official paper: Axelrod and Gómez-Bombarelli, *Scientific Data* 9, 185 (2022),
  <https://doi.org/10.1038/s41597-022-01288-4>.
- Official repository: <https://github.com/learningmatter-mit/geom>.
- Official release: Harvard Dataverse DOI `10.7910/DVN/JNGTDF`, dataset version 4,
  released 2022-02-11.  The audit checksummed the immutable MessagePack pair
  `qm9_crude.msgpack.tar.gz` (MD5 `aad0081ed5d9b8c93c2bd0235987573b`) and
  `qm9_featurized.msgpack.tar.gz` (MD5 `09655f470f438e3a7a0dfd20f40f6f22`).
  The fixture coordinates come from `qm9_crude`; `qm9_featurized` was used to confirm the graph
  fields and their xyz2mol provenance.  The release README says MessagePack is frozen while the Python-specific RDKit folder may be
  updated; these are therefore distinct dataset tracks, not interchangeable versions.

The paper reports about 133,000 QM9 species and 317,000 drug-like species, and about 37 million
conformations over more than 450,000 molecules.  QM9 consists of small C/N/O/F molecules from the
QM9 source set.  Drugs combines experimentally annotated drug-like sources; its preprocessing
canonicalized SMILES and treated clusters/salts (3.9%, 11,886 records), generally retaining the
heaviest component and applying documented acid/base transformations.  Original uncleaned SMILES
are retained for Drugs.  Exact per-release molecule/conformer counts beyond the paper totals are
to be computed from each chosen official artifact before a 5k run; mixing updated RDKit counts
with frozen MessagePack counts is prohibited.

## Graph and geometry provenance

The reference identity is the input/source canonical SMILES used to generate initial conformers.
For Drugs, RDKit embedded 50 starting conformers, GFN2-xTB optimized them, and the lowest-energy
one seeded CREST.  CREST/GFN2-xTB generated the ensemble.  GEOM then used Jensen-group `xyz2mol`
to *re-identify* a graph for every generated geometry, primarily to capture realized stereo and
possible reactions.  The paper reports source-graph recovery for 88.4% of QM9 species and 94.7%
of Drugs; it explicitly notes xyz2mol resonance mistakes, tautomerization, dissociation, and ring
formation.  Hence MolGR truth uses the source SMILES, while `xyz2mol_smiles` is provenance and a
source-integrity diagnostic only.  Ground truth is independent of xyz2mol in this protocol.

Formal charge exists on featurized atoms and is encoded in source SMILES, but no audited,
authoritative species-level spin multiplicity field was found.  Multiplicity is therefore
`UNKNOWN` in the release.  READY v1 is restricted to neutral, zero-RDKit-radical source graphs and
passes charge 0 / multiplicity 1 under an explicit closed-shell inclusion policy.  It does not
infer multiplicity for other records.

Stereo in source SMILES comes from upstream records when specified; conformer-side stereo was
assigned after 3D generation.  Since these are not uniformly equivalent truth sources, stereo is
not primary.

## Corrected / kekulized protocol audit

Direct RDKit graph or canonical-SMILES equality can disagree for aromatic bonds represented as
aromatic versus alternating single/double Kekule bonds, and for legitimate resonance forms.
Sanitization can also fail on strained/reacted GEOM coordinates or nonstandard valence forms.
Reference-side normalization may parse and canonicalize the published source representation and,
when explicitly declared, use a symmetric aromatic/Kekule representation.  It must not regenerate
the reference graph from XYZ.

Prediction-only sanitization, kekulization, valence/charge correction, fragment deletion, bond
rewriting, or choosing whichever resonance form happens to match is unfair: each changes the
method output after seeing the target.  The 2025 “GEOM-Drugs Revisited” correction concerns a 3D
*generative-model stability* pipeline (disconnected structures, sanitization/kekulization and
valency tables); it is useful evidence of dataset hazards but is not itself a graph-reconstruction
equivalence definition.

Evaluator v1 is the primary comparator: it gates formula, atom count, charge, radicals and graph
invariants, then records the decision/relation and bounded-search state.  No GEOM-specific
prediction normalization is needed.  GEOM-specific preprocessing is solely source-integrity
screening and the one-conformer-per-molecule selection, kept outside the evaluator.

## Recommended benchmark design

Primary: one deterministic conformer per molecule, chosen by minimum published relative energy
with source order as tie-breaker.  This tests molecular graph recovery and prevents a flexible
molecule with many near-duplicate conformers from dominating accuracy.  A later multi-conformer
robustness analysis should sample a fixed number per molecule and aggregate within molecule before
dataset-level reporting; it is secondary and is not part of READY.

For a 5k pilot, use GEOM-Drugs as the general-organic target because it has broader size,
heteroatom, aromatic and conformational coverage.  Run QM9 separately as a small-molecule/strained
stress slice, never merge their denominators.  Before Drugs pilot, explicitly decide whether to
use frozen MessagePack v4 or the updated RDKit release and implement the published source-integrity
filter without silently repairing fractured/tautomerized conformers.

## Human confirmations before 5k

1. Approve frozen Dataverse v4 MessagePack as the benchmark version, or choose and checksum an
   updated RDKit-folder snapshot.
2. Approve the neutral closed-shell-only scope; otherwise obtain authoritative multiplicities.
3. Approve source-SMILES truth and exclusion/reporting of conformers whose re-identified graph has
   changed connectivity (rather than treating xyz2mol output as truth).
4. Approve minimum-energy conformer selection and molecule-level denominator.
5. Approve chirality as secondary only.
6. For Drugs, approve the handling of source desalting/protonation provenance and the fracture
   exclusion list before sampling.

## 100-molecule smoke result

The deterministic fixture contains 100 distinct QM9 molecules and one minimum-relative-energy
conformer each: 20 stereo-containing, 21 aromatic, 43 with at least three heteroatoms, and 28 with
at least two RDKit rotatable bonds (categories overlap); heavy-atom counts span 3--9.  This is a
deliberately stratified engineering sample, not “the first 100” and not a paper-accuracy estimate.
The acquisition audit scanned 133,471 source molecule keys: 120,691 met the declared neutral
zero-radical parse policy, 727 failed reference parsing, and 12,053 were excluded by the electronic
state policy.  These counts are recorded rather than repaired or silently dropped.

`molgr_cpp` reconstructed all 100 without timeout or exception.  Evaluator v1 returned
`equivalent / normalized_graph_identity` for all 100, with charge and radical consistency in all
cases.  Non-isomeric exact canonical SMILES agreed for 97; the three diagnostic disagreements
were representation differences accepted by the evaluator, not prediction repairs.  The
secondary chirality-aware diagnostic agreed for all 100.  Total measured per-case wall time was
693.2 ms in the final verification run in this environment.  Machine-readable outputs are in
`_runs/smoke100/`.

These smoke results are an engineering validation, not paper accuracy.

STATUS: READY

Blockers to READY: none.  The confirmations above gate the 5k pilot, not this completed READY
stage.
