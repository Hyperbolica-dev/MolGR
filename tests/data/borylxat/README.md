# BorylXAT-DB boron-radical test inputs

Thirty SDF fixtures extracted from [Zenodo record 21330198, version
v2](https://doi.org/10.5281/zenodo.21330198), by Jia-Qi Li, Shuo-Qing Zhang,
Qiang Zhao, De-Guang Liu, Xi Lu, Yi-Feng Wang, Yao Fu, and Xin Hong.
The data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
the publisher's license is preserved in `LICENSE`.

| Category | Count | Coverage |
| --- | ---: | --- |
| `B` | 20 | BH2, fluoro, alkyl, aryl, cyano, O/N/S/Si substituents, rings, conjugation |
| `complex_r` | 5 | BH2 with trimethylamine, trimethylphosphine, PF3, and two carbene donors |
| `complex_p` | 5 | Matched chloroborane products for closed-shell controls |

`reference_graph/*.sdf` is the only stored structure format. Each V3000 SDF
contains the original atom order and double-precision coordinates (angstrom,
17-digit precision), corrected source-derived bonds, formal charges, explicit
radicals, total charge, multiplicity, and electronic-state provenance. Source
SMILES and available atom-indexed spin densities are preserved as properties.
No optimization, hydrogen addition, or coordinate reordering was applied.

`selection.json` pins source row IDs, structure keys, categories, and selection
reasons. `manifest.json` records source metadata, SDF checksums, corrected SMILES,
source-to-geometry atom mappings, and corrections. Its `source_geometry_sha256`
hashes the original SQLite element and position blobs (little-endian int32 atomic
numbers followed by float64 XYZ coordinates), allowing precision and atom order
to be verified without retaining separate XYZ files. Atom indices are zero-based.

The correction procedure is independent of MolGR reconstruction:

1. Preserve the source connectivity, bond orders, formal charges, and stereochemistry.
   RDKit may normalize an equivalent aromatic/Kekule representation.
2. Keep all source hydrogen atoms, disable implicit hydrogen addition, and remove
   extra bracket hydrogen counts on B. All 20 standalone source SMILES otherwise
   imply one more hydrogen than the coordinates contain.
3. Explicitly assign one radical electron to B for the 25 radical references,
   and zero radical electrons for the five product controls. Standalone B is
   neutral; complexes retain the source's covalent B-/donor+ Lewis convention.
4. Map the answer to source coordinate atom order. Source atom maps establish the standalone
   mappings; element-labelled graph isomorphism maps complexes, with equivalent
   hydrogen permutations allowed. Check full connectivity independently using
   covalent radii (factor 1.2), and check supplied stereo against the geometry.
5. Verify SDF and SMILES round trips preserve composition, charge, and radical
   count. Preserve the source spin-density array separately from formal radicals.

Examples of corrected SMILES:

| Structure | Corrected reference |
| --- | --- |
| `B_00388` (BH2 radical) | `[BH2]` |
| `B_00409` (BF2 radical) | `F[B]F` |
| `B_00388_LB_00001_r` (amine complex) | `[BH2-][N+](C)(C)C` |
| `B_00388_LB_00141_r` (PF3 complex) | `[BH2-][P+](F)(F)F` |

These are one chemically consistent **Lewis representation**, not unique
atom-localized spin distributions or experimentally reviewed answer graphs.
The B-centred radical form is retained for conjugated species even when another
atom has higher spin density. Evaluate resonance-equivalent outputs accordingly;
a different dative/charge-separated representation needs explicit normalization,
not automatic rejection by literal SMILES comparison.

Two simple species (BH2 and BF2) additionally check supported radical
reconstruction on both backends. The trimethylamine, trimethylphosphine, and PF3
radical complexes check that assigning negative charge to Bv3 also creates its
unpaired electron without breaking C-H or P-F bonds. The two carbene complexes
check the optional B=X to B-X bond-order/charge adjustment before resonance
traversal, retaining the original seed and respecting charge/radical budgets.
X is selected by the local electronic model: the neutral closed-shell starting
state must remain a closed-shell +1 center after reduction, while B- gains
exactly one radical. This is not restricted to carbon and does not permit
arbitrary elements or additional active electrons.
`B_00401` also accepts its reconstructed C-centred radical through an explicit
two-step resonance witness along source atom indices 0–4–3–2–1: the radical
moves from B0 to C3, then to C1, with alternating bond-order changes. Each step
is independently checked by the equivalence evaluator. The endpoint pair can
remain inconclusive in the general evaluator, which only proves a single such
move; the stored B-centred reference and general equivalence rules are unchanged.
Separate idealized BH3/BF3
geometry probes cover radical anions and neutral closed-shell controls, and
BH2-SiH2/BH2-GeH2 probes exercise non-carbon bond reductions on both backends.
These synthetic probes are not database records. General boron-radical reconstruction
support is not implied by these targeted checks.

- Multiplicity is not an explicit field in these source rows. For complexes,
  the author's [production workflow configuration](https://github.com/jackie-illiilli/BorylXAT-DB/blob/main/1_Calc_Reactant.ipynb)
  specifies `B_N_r` as doublet and `B_N_p` as singlet. The author's
  [category documentation](https://github.com/jackie-illiilli/BorylXAT-DB/blob/main/Database_Structure.md)
  also identifies product complexes as closed-shell. Standalone B doublets are
  inferred from the author's open-shell category, odd electron count, and
  spin-density sums within 1e-4 of 1. These distinct bases and their source URLs
  are preserved in each SDF; individual raw Gaussian logs have not been verified.
- Spin density is a continuous, potentially negative population, not an integer
  atom radical count. `B_00440`, for example, has B spin density 0.530277 and a C
  spin density of 0.720354. Preserve this evidence when reviewing localization.
- Source SMILES are retained verbatim and must not be used directly as reference
  answers. `B_00409` has BF2 coordinates, while its source SMILES
  `F[BH:1][F:2]` adds a hydrogen absent from the geometry. Even the parent BH2
  SMILES gains an implicit hydrogen when parsed normally. Do not regenerate
  coordinates or infer the electronic state from these uncorrected strings.

The record page and bundled documentation disagree on dataset totals and some
old release checksums. The manifest pins the **record-specific SQLite file** by
size and its published Zenodo MD5. Initial extraction used read-only HTTPS byte
ranges; the full database checksum was **not** verified. SDF hashes verify the
frozen outputs, not the entire source database. Full downloads stay outside Git.

Regenerate from a complete download (verifies its size and published MD5):

```bash
curl -L -o /tmp/BorylXAT-DB.db \
  'https://zenodo.org/records/21330198/files/BorylXAT-DB.db?download=1'
.venv/bin/python scripts/freeze_borylxat_fixtures.py --db /tmp/BorylXAT-DB.db
```

The extraction script corrects and validates the source graphs in memory before
writing SDF files and the manifest; it does not create intermediate XYZ files.

Load a reference and reconstruct from its coordinates alone:

```python
from pathlib import Path
from rdkit import Chem
from molgr.interface import xyz_to_rdmol

path = Path("tests/data/borylxat/reference_graph/B_00388.sdf")
reference = next(Chem.SDMolSupplier(str(path), removeHs=False))
conf = reference.GetConformer()
lines = [str(reference.GetNumAtoms()), reference.GetProp("_Name")]
for atom in reference.GetAtoms():
    coords = conf.GetAtomPosition(atom.GetIdx())
    lines.append(atom.GetSymbol() + " " + " ".join(format(x, ".17g") for x in coords))
mol = xyz_to_rdmol(
    "\n".join(lines) + "\n",
    total_charge=reference.GetIntProp("TOTAL_CHARGE"),
    spin_multiplicity=reference.GetIntProp("SPIN_MULTIPLICITY"),
    backend="cpp",
)
```

Do not pass reference bonds or atom-localized radicals into reconstruction tests.
Reconstruction may fail while boron-radical support is incomplete; a successful
return must still be compared with the corrected answer, allowing the agreed
resonance/coordination representation conventions.
