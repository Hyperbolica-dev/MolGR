# MolGR operational evaluator contract

The evaluator measures graph/cheminformatics-level equivalence for MolGR reconstruction. It does not prove identical quantum electronic states or exact localization of radical electrons in a transition-metal Lewis representation.

Policy designation: `molgr-evaluator-v2` (an evaluator semantic-version change, not a reconstruction or accuracy-improvement release).

Transition-metal coordination and Lewis-form representations can preserve structural molecular-graph equivalence, formula, connectivity, and total formal charge while differing in RDKit-local radical bookkeeping. Therefore component-electron mismatch remains diagnostic evidence but no longer independently establishes graph non-equivalence.

## Hard invariants

The following independently establish `NOT_EQUIVALENT` when they fail:

- heavy-atom formula;
- explicit-hydrogen formula;
- atom count;
- metal identity/formal-valence signature;
- non-metal connectivity;
- raw whole-molecule total formal charge.

These checks are not weakened by the component-electron policy.

The raw charge gate compares the sum of atom formal charges on each input before
metal-bond standardization, metal removal, or resonance normalization. Its values
and status are recorded in `invariants["raw_total_formal_charge"]`; the existing
`formal_charge` check remains a prepared-organic diagnostic. Charge mismatch
returns `NOT_EQUIVALENT` before any positive comparison path: resonance cannot
equate species that differ in total electron count. The local regression uses
`[S]C(=NN=Cc1ccccc1)N` versus `[S-]C(=NN=Cc1ccccc1)N` (charges 0 and -1),
which previously matched after thiosemicarbazone normalization.

## Positive evidence

After the hard structural checks pass, the evaluator may return `EQUIVALENT` for:

- normalized graph identity;
- established carbene/zwitterion normalization;
- independent full InChIKey agreement after the hard structural checks, including when it supplies the positive fallback for a component-electron diagnostic mismatch;
- successful resonance equivalence.

For open-shell forms, a mapped elementary three-center radical/bond-order
shift is explicit positive evidence. The mapping search is bounded and its
attempted count and limit status remain diagnostic; reaching the bound returns
`INCONCLUSIVE` unless stronger accepted identifier evidence is already present.

The witness requires equal, positive total RDKit radical labels on its two
normalized inputs. If this necessary condition fails, only this witness is
skipped: its search metadata remains unattempted, and all existing comparison
paths continue. Radical-label mismatch is not a top-level `NOT_EQUIVALENT` gate.

Prepared-organic InChIKey agreement without another accepted pathway remains `INCONCLUSIVE`.

With `use_chirality=False`, InChIKey generation removes stereochemistry and
conformers only from its private molecule copy, preventing InChI from inferring
stereochemistry from 3D coordinates. Input molecules remain unchanged.
With `use_chirality=True`, coordinate-derived stereochemistry remains enabled.

## Diagnostic signals

The per-component electron signature remains computed and persisted in `invariants["component_electrons"]`, together with aggregate radical counts, raw diagnostics, normalization effects, and contradiction metadata. A mismatch is diagnostic by default. RDKit radical bookkeeping, metal–ligand Lewis representation, and component-local electron accounting do not independently establish `NOT_EQUIVALENT`.

No representation-independent global electronic-state contradiction is currently implemented. The evaluator therefore does not synthesize an `electronic_state_conflict` from RDKit atom radical flags alone.

An accepted result can intentionally contain both:

```text
decision = equivalent
relation = identifier_equivalence
component_electrons.status = failed
```

## Inconclusive results

The evaluator returns `INCONCLUSIVE` when evidence is insufficient, including:

- prepared-organic identifier-only agreement;
- a bounded resonance search that reaches `max_resonance` without a match;
- a bounded open-shell topology-mapping search that reaches its limit without a witness;
- other comparison failures where no safe negative conclusion is available.

The resonance bound remains configuration-controlled; this contract does not increase it.

Ordinary resonance-search exceptions are not negative evidence. Failed or partial
enumerations are not cached as completed results. The evaluator continues to the
existing independent positive comparison paths; if none succeeds, it returns
`INCONCLUSIVE`. Search failure sets `bounded_search.exhaustive = False`, and
`reason` includes the exception type and message even when an independent path
returns `EQUIVALENT`. A later evaluation after recovery needs no manual cache
clear. `TimeoutError` still propagates to the caller; no retry or larger search
bound is introduced.

## Reproducibility

Benchmark releases should pin the RDKit major/minor version and persist the evaluator source hash, runtime inventory, input manifest hash, `use_chirality`, and `max_resonance`. RDKit-dependent resonance and sanitization behavior must be reported as version-sensitive rather than treated as proof of physical electronic-state identity.

Current authoritative validation uses RDKit `2026.03.6`. Resonance classification is not fully cross-version stable, so results involving resonance enumeration or normalization must be interpreted against the pinned runtime and revalidated when the RDKit version changes.

When two generated resonance sets intersect, the recorded hit is their
lexicographically smallest common SMILES, unless the existing direct-target
match applies. This makes hit selection independent of Python's hash seed
without changing which representations match.

## Benchmark scoring and reference eligibility

The three benchmark callers preserve `evaluator_decision`: EQ maps to
`equivalent=True`, NEQ to `False`, and INC to `None`. CSV output keeps the
decision even when the boolean cell is empty. tmQMg scores the original RDKit
candidate graph; SMILES roundtrip results remain separate diagnostics.

Reference eligibility is independent of evaluator and method outcomes:
`reference_status` is `valid`, `excluded`, or `unresolved`, with the supporting
`reference_reason`. Explicit confirmed-error metadata with a reason, or the
existing independent reference/XYZ formula contradiction, permits exclusion.
Missing/unparseable references and unresolved review records are not automatically
confirmed errors. The historical manual whitelist is read-only reference-quality
evidence, never an equivalence rule. `accepted_reference_issue` requires explicit
`reference_answer_wrong=true` and a reason to establish exclusion;
`accepted_ambiguous` and unconfirmed entries remain unresolved.

`total_count` includes every input. `count` aliases `reference_valid_count`;
confirmed exclusions and unresolved references have separate counts and do not
enter that denominator. Within valid references, `success_count`, `fail_count`,
`inconclusive_count`, and `runtime_failure_count` partition the denominator.
INC and failed/skipped methods retain their valid-reference denominator slots.
Runtime failures do not become NEQ or reference exclusions.

`tmqmg_regression.py` preserves its raw `equivalence_decision` and no longer lets
manual review or formula diagnostics overwrite `equivalent`, `strict_equivalent`,
or `effective_equivalent`. The last field is compatibility data, not the formal
score. Official `equivalent_fraction` is EQ divided by the valid-reference count.
Freeze the same reference-status membership and graph inputs before comparing
the original evaluator, merged v1, and v2; do not derive eligibility from their
decisions. Review records and existing artifacts are not rewritten.

An explicit frozen reference-status membership is reused rather than reopened
by an unconfirmed legacy whitelist entry. Independent reference/XYZ composition
contradictions are still recorded; evaluator outcomes never select membership.

## Frozen three-version replay

`scripts/evaluator_replay.py` consumes a fixed reference-validity manifest,
candidate byte inventory, and three standalone evaluator snapshots. Frozen
membership uses `valid`, `confirmed_excluded`, and `unresolved`; existing curated
adjudications remain evidence without requiring newly introduced review fields.
Original published references are not replaced with accepted candidate graphs.

Every selected case gets raw per-version records before validity aggregation.
INC, evaluator errors, timeouts, input errors, and unavailable candidates remain
distinct. Empty reference SMILES are input errors, not empty-molecule NEQ.
Lossless MolBlocks take precedence over V2000, with no SMILES fallback.
Calls follow benchmark order, reference then candidate, retaining conformers.

SQLite commits each `(case_id, version)` independently. Identical invocations
resume missing rows; changed inputs, source, runtime binaries, or parameters
require a new output directory. Use `--limit` for smoke only; the default exact
RDKit requirement is `2026.03.6`. The full run requires separate authorization.

The replay runner supports Linux and macOS environments with the multiprocessing
`fork` start method. Windows is not supported; replay integration tests are skipped
there. Summary/statistics imports do not require matplotlib; plotting loads the
optional dependency only when requested.
