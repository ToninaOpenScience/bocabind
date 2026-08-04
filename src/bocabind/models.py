from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    element: str
    residue_name: str
    chain: str
    residue_number: int
    insertion_code: str
    coordinate: np.ndarray
    record: str
    occupancy: float = 1.0
    altloc: str = ""

    @property
    def residue_id(self) -> str:
        suffix = self.insertion_code.strip()
        return f"{self.residue_name}:{self.chain or '_'}:{self.residue_number}{suffix}"


@dataclass
class Structure:
    path: Path
    atoms: list[Atom]
    warnings: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class LigandCandidate:
    identifier: str
    residue_name: str
    chain: str
    residue_number: int
    heavy_atoms: int
    protein_contacts: int
    buried_fraction: float
    score: float
    confidence: float
    likely_role: str


@dataclass
class InspectionResult:
    input_path: str
    ligand_candidates: list[LigandCandidate]
    warnings: list[str]
    automatic_selection: str | None


@dataclass
class StepMetrics:
    distance: float
    clash_pairs: int
    max_penetration: float
    total_penetration: float
    blocked_ligand_fraction: float
    residues: dict[str, float]


@dataclass
class PathResult:
    direction: list[float]
    profile: list[StepMetrics]
    max_penetration: float
    cumulative_penetration: float
    peak_clash_pairs: int
    max_blocked_fraction: float
    obstructed_length: float
    clear: bool


@dataclass
class AnalysisResult:
    input_path: str
    ligand: LigandCandidate
    binding_site_residues: list[str]
    classification: str
    association_compatibility_score: float
    evidence_strength: str
    best_path: PathResult
    directions_tested: int
    clear_direction_fraction: float
    obstruction_location: str
    primary_obstructing_residues: list[dict[str, Any]]
    interpretation: str
    limitations: list[str]
    warnings: list[str]
    parameters: dict[str, Any]
    _structure: Structure | None = field(default=None, repr=False)
    _ligand_atoms: list[Atom] | None = field(default=None, repr=False)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data.pop("_structure", None)
        data.pop("_ligand_atoms", None)
        return data

    def write(self, output: str | Path) -> Path:
        from .reporting import write_outputs

        return write_outputs(self, Path(output))

