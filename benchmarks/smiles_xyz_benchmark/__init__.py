from __future__ import annotations


def main() -> int:
    from benchmarks.smiles_xyz_benchmark.run import main as _main

    return _main()


__all__ = ["main"]
