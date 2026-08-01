# Molecule reconstruction review and trace

`tools/molgr_review/` provides development-time validation of MolGR
reconstruction results. It executes the current checkout, presents XYZ,
candidate, reference, and trace views, and persists confirmed answers as
regression fixtures.

The tmQMg workflow has four explicit stages:

1. `prepare_tmqmg_queue.py` runs the dataset benchmark, updates the review CSV,
   and synchronizes the local review database by default. Partial runs merge
   only their selected scope into the existing queue.
2. `import_cases.py` imports an external or manually generated complete CSV.
3. `server.py` serves the review interface using the current MolGR runtime.
4. Confirmed answers are written to `tests/data/reviewed/tmqmg/`.

Review databases, generated queues, benchmark and trace output, and JSONL
backups remain under `.local/`. Pending and skipped decisions are local state;
only confirmed answers are version-controlled fixtures.

See the [Chinese guide](MOLECULE_REVIEW_TOOL.zh-CN.md) for the pinned tmQMg
source, runtime consistency checks, commands, decision semantics, queue schema,
and validation procedure.
