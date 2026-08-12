# Benchmark Execution Notes

## BDE-db

- Branch: `benchmark/bde-db-ready`
- The DOI `10.6084/m9.figshare.12609365` returned HTTP 403 during both automated and
  manual access. This record is the approximately 7 kB machine-readable metadata record,
  not the molecular structure dataset, so the access failure is not a readiness blocker.
- Primary provenance path: the Springer Nature Figshare collection associated with
  *Quantum chemical calculations for over 200,000 organic radical species and 40,000
  associated closed-shell molecules* -> `10.6084/m9.figshare.12158646` ->
  `20200415_radical_database.sdf.gz` (approximately 264 MB), Figshare file id
  `22357962`.
- The Springer Nature Figshare article page currently returns HTTP 403 during manual and
  automated access.
- Direct download attempt on 2026-08-10:
  `https://springernature.figshare.com/ndownloader/files/22357962` returned HTTP 403 on
  four attempts (initial request plus three retries). The same official endpoint returned
  HTTP 403 on four attempts through the configured SOCKS5 proxy. Each error response was
  118 bytes; no partial SDF was retained.
- `https://api.figshare.com/v2/articles/12158646` returned HTTP 403, so its `.files[]`
  metadata, API-provided `download_url`, and supplied/computed MD5 values could not be read.
- `https://ndownloader.figshare.com/files/22357962` and
  `https://api.figshare.com/v2/file/download/22357962` each returned HTTP 403 on the initial
  request and two retries.
- `DATA_DOWNLOAD_STATUS: RESOLVED`. The official endpoints above were temporarily blocked,
  but the user supplied the official SDF and metadata locally. The official repository
  exposes no separate download script, alternative official URL, mirror, or published
  checksum.
- The approximately 39.7 GB raw Gaussian log archive is intentionally excluded from the
  initial download. A small number of logs will only be considered if the processed SDF and
  official `pstjohn/bde` code cannot establish charge, radical, or multiplicity semantics.
- Raw files are read-only inputs under `.local/bde_db/` and are not committed.
- Multiplicity derivation was initially `UNKNOWN`; actual-file verification below resolves
  it for the released zero/one-radical-electron scope.

### Official-code evidence

- Official repository: `https://github.com/pstjohn/bde`, inspected at commit
  `5677af8dcbb992c7888746aa018302e6fb04e67d` (2022-06-30). The commit immediately
  following the SDF publication is `2c1af8f` (2020-04-16, `updating for new properties`).
- `bde/fragment.py` builds radicals by breaking an RDKit bond and calling `SanitizeMol` to
  assign the resulting formal radical representation.
- `bde/gaussian.py` creates the reference molecule from the input SMILES, adds explicit
  hydrogens, and replaces its conformer coordinates with the optimized Gaussian geometry.
  The resulting SDF graph therefore comes from the source SMILES/RDKit path, independently
  of xyz2mol.
- For Gaussian input generation, the official code writes the RDKit molecule to a temporary
  SDF and delegates conversion to Open Babel. The fragment calculation reads the generated
  charge/multiplicity line back from that file for its second Gaussian job, but does not
  explicitly calculate or persist multiplicity in Python.
- The official loading notebook expects 289,639 SDF records. It classifies records using
  RDKit `NumRadicalElectrons`, reads `AtomSpins` only when present, and shows both closed-shell
  `molecule` and radical `fragment` records.
- Actual-file verification below confirms that every released record encodes exactly zero or
  one radical electron; no multi-radical record is present.

### Source and version evidence

- Primary paper: P. C. St. John et al., *Quantum chemical calculations for over 200,000
  organic radical species and 40,000 associated closed-shell molecules*, Scientific Data 7,
  244 (2020), DOI `10.1038/s41597-020-00588-x`.
- Structure release: Springer Nature Figshare DOI `10.6084/m9.figshare.12158646`, file
  `20200415_radical_database.sdf.gz`, file id `22357962`.
- BDE relation release: BDE-db v1, DOI `10.6084/m9.figshare.10248932.v1`, described as
  290,664 homolytic BDE rows. This relation table is distinct from the 289,639-record
  optimized-structure SDF.
- The paper reports 43,276 closed-shell molecules and 246,363 radicals, totaling 289,639
  optimized structures. The official loading notebooks use the same expected total.
- The paper reports that parent SMILES were selected from PubChem, restricted to connected
  C/H/N/O structures with no atom formal charges. Radicals were produced by breaking single,
  non-ring bonds, then canonicalized and deduplicated with RDKit.

### Field semantics and benchmark implications

- Published SDF properties are `SMILES`, `Enthalpy`, `FreeEnergy`, `SCFEnergy`,
  `AtomCharges`, `AtomSpins` for radicals, `VibFreqs`, `RotConstants`, and `IRIntensity`.
  `AtomCharges` and `AtomSpins` are Mulliken values aligned to coordinate atom order; they
  are not formal charges or a multiplicity field.
- No explicit multiplicity SDF property is documented in the paper or used by the official
  notebooks. The official generation code labels jobs as `molecule` or `fragment`, creates
  monovalent formal radicals, and delegates the Gaussian charge/multiplicity line to Open
  Babel from the RDKit SDF representation.
- Protocol policy: derive singlet only from exactly zero encoded radical electrons and
  doublet only from exactly one encoded radical electron. Reject records with more than one
  encoded radical electron because radical count alone does not uniquely determine their
  spin multiplicity.
- Formal total charge is read from SDF atom formal-charge records, not from Mulliken
  `AtomCharges`. Actual-file verification identifies one charged released fragment.
- Radical site is read from the atom carrying RDKit formal radical electrons. `AtomSpins`
  is retained as source metadata and may diagnose delocalization, but does not replace the
  formal radical-site reference.
- The official radical-type notebook returns the first atom with positive formal radical
  electrons and does not test for multiple radical centers. Its published analysis therefore
  supports the intended monovalent-radical design but is not proof that every SDF record has
  exactly one radical electron.
- SDF bond connectivity is independent of xyz2mol: `gaussian.py` reconstructs the molecule
  from source SMILES with RDKit, adds explicit hydrogens, and replaces only conformer
  coordinates with the final Gaussian geometry before writing the mol block.
- SDF/SMILES correspondence, RDKit radical/formal-charge reading, charged-species counts,
  and multi-radical presence are resolved by the actual-file verification below.

### Pre-data readiness snapshot

- `PROVENANCE_READY`: YES. The authoritative paper, structure DOI/file id, BDE relation DOI,
  official repository, generation path, and documented field meanings are identified. All
  claims that require the unavailable SDF remain explicitly `UNKNOWN`.
- `ADAPTER_READY`: YES. A minimal read-only `molgr_cpp` adapter supports deterministic
  stratified sampling, `--limit`, inclusive `--start`/`--end`, `--seed`, explicit loader and
  reconstruction failures, existing MolGR equivalence, diagnostic charge/radical metrics,
  and review-package output. Synthetic SDF tests pass; this is not a dataset smoke test.
- `DATA_LOCAL_READY`: NO (`TEMPORARY_NETWORK_BLOCKER`).
- `SMOKE_TEST_READY`: NO (requires the official SDF).

### Actual-file verification

- Local file: `.local/bde_db/raw/20200415_radical_database.sdf.gz`.
- Size: 276,696,272 bytes compressed; gzip integrity check passed.
- SHA-256: `9de7dc389c941c0d804a9fca468f43951f5987bb7299db18841e35a8458eec60`.
- RDKit strictly parsed all 289,639 records with no null records or property-parse failures.
- Every record has `SMILES`, `AtomCharges`, `Enthalpy`, `FreeEnergy`, `SCFEnergy`, and
  `RotConstants`. `VibFreqs` and `IRIntensity` are each absent on two records.
- `AtomCharges` length matches atom count for all records. `AtomSpins` length matches atom
  count wherever present.
- Formal radical encoding is exhaustive and unambiguous for this release: 43,276 records
  have zero radical electrons and zero radical centers; 246,363 records have exactly one
  radical electron on exactly one atom. No record encodes two or more radical electrons.
- Radical centers: C 205,572; N 25,434; O 15,356; H 1.
- Formal charge: 289,638 neutral records and one charge -1 record. The charged record is
  `508590_dc4e99`, SMILES `[O-]`, with one radical electron and multiplicity 2.
- The isolated hydrogen radical is `12_0`, SMILES `[H]`, with one radical electron and
  multiplicity 2.
- No explicit multiplicity property exists. For the actual release, zero/one formal radical
  electron maps uniquely to the intended singlet/doublet calculation classes documented by
  the official generation pipeline. The dataset contains no triplet or higher-spin cases.
- `AtomSpins` occurs on all 246,363 formal radicals and six additional closed-shell records:
  `9584_0`, `57156_07e644`, `78821_f5e25c`, `150308_dbfd3c`, `192348_5f9968`, and
  `218537_c16069`. Therefore `AtomSpins` presence is not a safe radical classifier.
- After normalizing explicit versus implicit hydrogens for graph comparison, all 289,639 SDF
  reference graphs agree with their `SMILES` property in heavy-atom connectivity, bond order,
  formal charge, and radical representation.
- 6,208 records disagree in isomeric/stereo SMILES. Primary graph equivalence must therefore
  use `use_chirality=False`; exact isomeric SMILES is retained only as a diagnostic metric.
- Deterministic 20-case inspection output:
  `.local/bde_db/inspection/inspected_cases.csv` and `inspection_summary.json`.

### 100-case smoke test

- Output: `.local/bde_db/runs/pilot-100-seed0/`.
- Seed: 0. Selected strata: 20 closed-shell, 49 C radicals, 15 N radicals, 15 O radicals,
  and the single charged radical.
- Reconstruction success: 100/100.
- Reference graph equivalence: 100/100. This is a pilot diagnostic, not a final accuracy
  estimate.
- Exact isomeric SMILES: 86/100.
- Charge consistency: 100/100.
- Radical-electron consistency: 100/100.
- Exact radical-site index consistency: 91/100; the remaining nine are graph-equivalent
  radical localization/resonance differences for manual review.
- Runtime per case: mean 1.35 ms, median 1.14 ms, p95 1.63 ms, maximum 19.14 ms in this run.
- Reconstruction failures, exceptions, and loader failures: zero.

### Final readiness

- `PROVENANCE_READY`: YES.
- `ADAPTER_READY`: YES.
- `DATA_LOCAL_READY`: YES.
- `SMOKE_TEST_READY`: YES.
- Overall BDE-db benchmark status: `READY` for manual review before any 2k-5k expansion.

### Approved 5k protocol

- Primary metric remains resonance-aware graph equivalence with
  `use_chirality=False` and `max_resonance=100`.
- Exact formal-radical localization is named `formal_radical_atom_index_match` and is only
  evaluated after an element/order/coordinate guard maps prediction atoms to original
  XYZ/reference atom identities.
- `MolGRCppMethod` is instantiated once per benchmark run. Initialization and a closed-shell
  warm-up are recorded separately from formal per-molecule runtime.
- Ordinary method exceptions and timeouts are isolated per case and recorded without
  aborting the run.
- The deterministic sample follows observed closed-shell/C/N/O population proportions and
  explicitly includes `[H]` and `[O-]`; no artificial charged quota is used.
- Manual review prioritizes reconstruction failures, non-equivalent cases,
  resonance-equivalent cases, formal-radical atom-index mismatches, charge/radical-electron
  mismatches, and exact-SMILES mismatches.
- Before an eventual 289,639-record run, loader and output handling must be redesigned for
  chunked/streaming execution so all cases, results, XYZ strings, bond dumps, and metadata
  are not retained in memory.
