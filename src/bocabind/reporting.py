from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .models import AnalysisResult, Atom


def _pdb_line(atom: Atom, coordinate: np.ndarray) -> str:
    return (f"{atom.record:<6}{atom.serial:>5} {atom.name:<4}{atom.altloc:1}{atom.residue_name:>3} "
            f"{atom.chain:1}{atom.residue_number:>4}{atom.insertion_code:1}   "
            f"{coordinate[0]:>8.3f}{coordinate[1]:>8.3f}{coordinate[2]:>8.3f}"
            f"{atom.occupancy:>6.2f}{0.0:>6.2f}          {atom.element:>2}\n")


def _summary(result: AnalysisResult) -> str:
    best = result.best_path
    residues = ", ".join(x["residue"] for x in result.primary_obstructing_residues[:5]) or "None"
    return f"""Tonina BocaBind v0.1.0

Ligand: {result.ligand.identifier}
Structural classification: {result.classification}
Association Compatibility Score: {result.association_compatibility_score:.1f}/100 (heuristic)
Evidence strength: {result.evidence_strength}

Best rigid extraction route
- Direction: {', '.join(f'{x:.4f}' for x in best.direction)}
- Maximum penetration: {best.max_penetration:.3f} A
- Cumulative penetration: {best.cumulative_penetration:.3f} A^2
- Peak clashing pairs: {best.peak_clash_pairs}
- Maximum blocked ligand fraction: {best.max_blocked_fraction:.1%}
- Obstructed path length: {best.obstructed_length:.2f} A
- Clear tested directions: {result.clear_direction_fraction:.1%}
- Obstruction location: {result.obstruction_location}
- Primary obstructing residues: {residues}

Interpretation
{result.interpretation}

Important: this result is a rigid-body geometric classification, not a probability of induced fit.
"""


def write_outputs(result: AnalysisResult, output: Path) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    source = Path(result.input_path)
    data = result.to_dict()
    data["run_metadata"] = {
        "bocabind_version": "0.1.0",
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "classification_rule_version": "heuristic-v0.1",
    }
    (output / "report.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    (output / "summary.txt").write_text(_summary(result), encoding="utf-8")
    with (output / "path_profile.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["distance_A", "clash_pairs", "max_penetration_A", "total_penetration_A", "blocked_ligand_fraction"])
        for step in result.best_path.profile:
            writer.writerow([step.distance, step.clash_pairs, step.max_penetration,
                             step.total_penetration, step.blocked_ligand_fraction])
    with (output / "obstructing_residues.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["residue", "penetration_sum", "fraction_of_obstruction"])
        writer.writeheader()
        writer.writerows(result.primary_obstructing_residues)
    if result._structure is not None and result._ligand_atoms is not None:
        ligand_ids = {id(a) for a in result._ligand_atoms}
        fixed = [a for a in result._structure.atoms if id(a) not in ligand_ids]
        direction = np.array(result.best_path.direction)
        with (output / "best_path.pdb").open("w", encoding="utf-8") as handle:
            for model, step in enumerate(result.best_path.profile, 1):
                handle.write(f"MODEL     {model:4d}\n")
                for atom in fixed:
                    handle.write(_pdb_line(atom, atom.coordinate))
                for atom in result._ligand_atoms:
                    handle.write(_pdb_line(atom, atom.coordinate + step.distance * direction))
                handle.write("ENDMDL\n")
            handle.write("END\n")
    return output

