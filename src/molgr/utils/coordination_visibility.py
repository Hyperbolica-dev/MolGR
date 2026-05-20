"""Geometry helpers for metal coordination visibility checks."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Tuple

import numpy as np


Point3D = Tuple[float, float, float]
CoordinationBlocker = Tuple[int, Point3D, float]
CoordinationBlockerArrays = Tuple[np.ndarray, np.ndarray, np.ndarray]


def empty_coordination_blocker_arrays() -> CoordinationBlockerArrays:
    return (
        np.empty((0,), dtype=np.int64),
        np.empty((0, 3), dtype=np.float64),
        np.empty((0,), dtype=np.float64),
    )


def coordination_blockers_to_arrays(
    blockers: Sequence[CoordinationBlocker],
) -> CoordinationBlockerArrays:
    if not blockers:
        return empty_coordination_blocker_arrays()

    blocker_indices, blocker_coordinates, blocker_radii = zip(*blockers)
    return (
        np.asarray(blocker_indices, dtype=np.int64),
        np.asarray(blocker_coordinates, dtype=np.float64),
        np.asarray(blocker_radii, dtype=np.float64),
    )


def distance_point_to_segment(
    point: Point3D,
    segment_start: Point3D,
    segment_end: Point3D,
) -> float:
    point_array = np.asarray(point, dtype=np.float64)
    segment_start_array = np.asarray(segment_start, dtype=np.float64)
    segment_end_array = np.asarray(segment_end, dtype=np.float64)
    segment_vector = segment_end_array - segment_start_array
    seg_len_sq = float(np.dot(segment_vector, segment_vector))
    if seg_len_sq <= 1.0e-12:
        return float(np.linalg.norm(point_array - segment_start_array))

    projection = np.clip(
        np.dot(point_array - segment_start_array, segment_vector) / seg_len_sq,
        0.0,
        1.0,
    )
    closest = segment_start_array + projection * segment_vector
    return float(np.linalg.norm(point_array - closest))


def coordination_visibility_mask(
    atom_indices: Sequence[int] | np.ndarray,
    atom_coordinates: Sequence[Point3D] | np.ndarray,
    metal_coordinates: Point3D | np.ndarray,
    blocker_arrays: CoordinationBlockerArrays,
) -> np.ndarray:
    atom_indices_array = np.asarray(atom_indices, dtype=np.int64)
    atom_coordinates_array = np.asarray(atom_coordinates, dtype=np.float64)
    if atom_indices_array.size == 0:
        return np.empty((0,), dtype=bool)
    if atom_coordinates_array.ndim == 1:
        atom_coordinates_array = atom_coordinates_array.reshape(1, 3)

    metal_coordinates_array = np.asarray(metal_coordinates, dtype=np.float64)
    segment_vectors = atom_coordinates_array - metal_coordinates_array
    segment_lengths_sq = np.einsum("ij,ij->i", segment_vectors, segment_vectors)
    visible = segment_lengths_sq > 1.0e-12

    blocker_indices, blocker_coordinates, blocker_radii = blocker_arrays
    if blocker_indices.size == 0:
        return visible

    safe_segment_lengths_sq = np.where(visible, segment_lengths_sq, 1.0)
    blocker_offsets = blocker_coordinates[:, np.newaxis, :] - metal_coordinates_array
    projections = np.einsum("bij,ij->bi", blocker_offsets, segment_vectors)
    projections = np.clip(projections / safe_segment_lengths_sq[np.newaxis, :], 0.0, 1.0)
    closest_points = (
        metal_coordinates_array[np.newaxis, np.newaxis, :]
        + projections[:, :, np.newaxis] * segment_vectors[np.newaxis, :, :]
    )
    delta = blocker_coordinates[:, np.newaxis, :] - closest_points
    distance_sq = np.einsum("bij,bij->bi", delta, delta)
    blocker_radius_sq = blocker_radii[:, np.newaxis] * blocker_radii[:, np.newaxis]
    self_mask = blocker_indices[:, np.newaxis] == atom_indices_array[np.newaxis, :]
    blocked = (~self_mask) & (distance_sq < blocker_radius_sq)
    return visible & ~np.any(blocked, axis=0)


def has_unobstructed_coordination_path(
    atom_idx: int,
    atom_coordinates: Point3D,
    metal_coordinates: Point3D,
    blockers: Sequence[CoordinationBlocker],
) -> bool:
    blocker_arrays = coordination_blockers_to_arrays(blockers)
    return bool(
        coordination_visibility_mask(
            np.asarray([atom_idx], dtype=np.int64),
            np.asarray([atom_coordinates], dtype=np.float64),
            np.asarray(metal_coordinates, dtype=np.float64),
            blocker_arrays,
        )[0]
    )


__all__ = [
    "CoordinationBlocker",
    "CoordinationBlockerArrays",
    "Point3D",
    "coordination_blockers_to_arrays",
    "coordination_visibility_mask",
    "distance_point_to_segment",
    "empty_coordination_blocker_arrays",
    "has_unobstructed_coordination_path",
]
