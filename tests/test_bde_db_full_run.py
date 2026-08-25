from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from benchmarks.bde_db_benchmark.full_run import run_full
from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize


def _write_sdf(path: Path, smiles: list[str]) -> None:
    plain_path = path.with_suffix("")
    writer = Chem.SDWriter(str(plain_path))
    for index, value in enumerate(smiles, start=1):
        mol = Chem.AddHs(Chem.MolFromSmiles(value))
        assert AllChem.EmbedMolecule(mol, randomSeed=index) == 0
        mol.SetProp("_Name", f"case-{index}")
        mol.SetProp("SMILES", value)
        writer.write(mol)
    writer.close()
    with plain_path.open("rb") as source, gzip.open(path, "wb") as destination:
        shutil.copyfileobj(source, destination)


class _IdentityMethod(BenchmarkMethod):
    def run(self, case: dict) -> MethodRunOutput:
        mol = remove_hs_without_sanitize(case["ground_truth_rdmol"])
        return MethodRunOutput(
            status="ok",
            predicted_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            rdkit_mol=mol,
        )


def test_full_runner_streams_compact_results_and_freezes_bundle(
    tmp_path: Path, monkeypatch
) -> None:
    input_path = tmp_path / "test.sdf.gz"
    out_dir = tmp_path / "run"
    _write_sdf(input_path, ["CC", "[CH3]", "[OH]", "[H]"])
    monkeypatch.setattr("benchmarks.bde_db_benchmark.full_run._git_clean", lambda: True)
    monkeypatch.setattr("benchmarks.bde_db_benchmark.full_run._git_sha", lambda: "test-sha")
    monkeypatch.setattr(
        "benchmarks.bde_db_benchmark.full_run.MolGRCppMethod",
        lambda: _IdentityMethod("identity"),
    )

    summary = run_full(input_path, out_dir, expected_records=4, seed=7, timeout_seconds=None)

    with gzip.open(out_dir / "results.csv.gz", "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 4
    assert "xyz" not in rows[0]
    assert "reference_bonds" not in rows[0]
    assert "source_metadata" not in rows[0]
    assert summary["metrics"]["overall"]["equivalent"] == 4
    bundle = out_dir / "bde_db_paper_export"
    assert (bundle / "results.csv.gz").read_bytes() == (out_dir / "results.csv.gz").read_bytes()
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["molgr_git_sha"] == "test-sha"
    assert set(manifest["files"]) == {
        "failures.csv",
        "provenance.md",
        "results.csv.gz",
        "review_cases.csv",
        "summary.json",
    }
