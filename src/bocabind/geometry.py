from __future__ import annotations

from collections import defaultdict
from math import pi, sqrt

import numpy as np
from scipy.spatial import cKDTree

from .models import Atom, PathResult, StepMetrics

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "B": 1.92, "SI": 2.10,
    "FE": 1.80, "ZN": 1.39, "MG": 1.73, "CA": 1.94, "NA": 2.27, "K": 2.75,
}


def fibonacci_directions(count: int) -> np.ndarray:
    if count < 6:
        raise ValueError("directions must be at least 6")
    golden = pi * (3.0 - sqrt(5.0))
    indices = np.arange(count, dtype=float)
    y = 1.0 - 2.0 * (indices + 0.5) / count
    radius = np.sqrt(1.0 - y * y)
    theta = golden * indices
    return np.column_stack((np.cos(theta) * radius, y, np.sin(theta) * radius))


def binding_site_residues(protein: list[Atom], ligand: list[Atom], cutoff: float) -> list[str]:
    ligand_coords = np.array([a.coordinate for a in ligand])
    tree = cKDTree(np.array([a.coordinate for a in protein]))
    residues = set()
    for point in ligand_coords:
        residues.update(protein[i].residue_id for i in tree.query_ball_point(point, cutoff))
    return sorted(residues)


def _step_metrics(
    moved: np.ndarray,
    ligand: list[Atom],
    protein: list[Atom],
    tree: cKDTree,
    protein_radii: np.ndarray,
    tolerance: float,
    distance: float,
) -> StepMetrics:
    penetrations = []
    ligand_atoms_blocked = set()
    residues: dict[str, float] = defaultdict(float)
    max_query = tolerance * (max(VDW_RADII.values()) + max(VDW_RADII.values()))
    for i, point in enumerate(moved):
        ligand_radius = VDW_RADII.get(ligand[i].element, 1.70)
        for j in tree.query_ball_point(point, max_query):
            threshold = tolerance * (ligand_radius + protein_radii[j])
            separation = float(np.linalg.norm(point - protein[j].coordinate))
            penetration = threshold - separation
            if penetration > 0:
                penetrations.append(penetration)
                ligand_atoms_blocked.add(i)
                residues[protein[j].residue_id] += penetration
    return StepMetrics(
        distance=float(distance),
        clash_pairs=len(penetrations),
        max_penetration=max(penetrations, default=0.0),
        total_penetration=float(sum(penetrations)),
        blocked_ligand_fraction=len(ligand_atoms_blocked) / max(len(ligand), 1),
        residues=dict(residues),
    )


def scan_paths(
    protein: list[Atom],
    ligand: list[Atom],
    directions: int = 256,
    step_size: float = 0.5,
    path_length: float | None = None,
    tolerance: float = 0.80,
) -> list[PathResult]:
    if step_size <= 0 or tolerance <= 0:
        raise ValueError("step_size and clash_tolerance must be positive")
    protein_coords = np.array([a.coordinate for a in protein])
    ligand_coords = np.array([a.coordinate for a in ligand])
    tree = cKDTree(protein_coords)
    protein_radii = np.array([VDW_RADII.get(a.element, 1.70) for a in protein])
    center = ligand_coords.mean(axis=0)
    if path_length is None:
        nearest_surface = np.linalg.norm(protein_coords - center, axis=1).max()
        ligand_radius = np.linalg.norm(ligand_coords - center, axis=1).max()
        path_length = min(max(nearest_surface + ligand_radius + 4.0, 8.0), 45.0)
    distances = np.arange(0.0, path_length + step_size * 0.5, step_size)
    results = []
    for direction in fibonacci_directions(directions):
        profile = [
            _step_metrics(ligand_coords + distance * direction, ligand, protein, tree,
                          protein_radii, tolerance, distance)
            for distance in distances
        ]
        max_penetration = max(s.max_penetration for s in profile)
        cumulative = float(np.trapezoid([s.total_penetration for s in profile], distances))
        obstructed = [s.distance for s in profile if s.max_penetration > 0.10]
        obstructed_length = (max(obstructed) - min(obstructed) + step_size) if obstructed else 0.0
        results.append(PathResult(
            direction=direction.tolist(), profile=profile, max_penetration=max_penetration,
            cumulative_penetration=cumulative,
            peak_clash_pairs=max(s.clash_pairs for s in profile),
            max_blocked_fraction=max(s.blocked_ligand_fraction for s in profile),
            obstructed_length=float(obstructed_length), clear=bool(max_penetration <= 0.10),
        ))
    return sorted(results, key=lambda p: (not p.clear, p.max_penetration, p.max_blocked_fraction,
                                         p.cumulative_penetration, p.obstructed_length))

