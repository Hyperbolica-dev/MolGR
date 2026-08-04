"""Shared metal-ligand coordination distance helpers."""

from openbabel import openbabel as ob


def coordination_distance_cutoff(
    metal_atomic_num: int,
    ligand_atomic_num: int,
    *,
    radius_scale: float,
    extra_tolerance_angstrom: float,
) -> float:
    """Return the covalent-radius based metal-ligand coordination cutoff."""

    metal_radius = float(ob.GetCovalentRad(int(metal_atomic_num)))
    ligand_radius = float(ob.GetCovalentRad(int(ligand_atomic_num)))
    return float(radius_scale) * (metal_radius + ligand_radius) + float(extra_tolerance_angstrom)


__all__ = ["coordination_distance_cutoff"]
