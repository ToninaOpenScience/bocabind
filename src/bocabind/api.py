from __future__ import annotations

from pathlib import Path

from .classification import classify, summarize_obstruction
from .geometry import binding_site_residues, scan_paths
from .ligands import atoms_for_candidate, choose_candidate, detect_ligands, protein_atoms
from .models import AnalysisResult, InspectionResult
from .pdb import read_pdb


def inspect_structure(path: str | Path) -> InspectionResult:
    structure = read_pdb(path)
    candidates = detect_ligands(structure)
    automatic = None
    if candidates:
        if len(candidates) == 1 or (
            candidates[0].confidence >= 0.80
            and candidates[0].score - candidates[1].score >= 2.0
        ):
            automatic = candidates[0].identifier
    return InspectionResult(str(structure.path), candidates, structure.warnings, automatic)


def analyze(
    path: str | Path,
    ligand: str | None = None,
    *,
    directions: int = 256,
    step_size: float = 0.5,
    site_cutoff: float = 5.0,
    clash_tolerance: float = 0.80,
    path_length: float | None = None,
) -> AnalysisResult:
    structure = read_pdb(path)
    candidates = detect_ligands(structure)
    selected = choose_candidate(candidates, ligand)
    ligand_atoms = atoms_for_candidate(structure, selected)
    proteins = protein_atoms(structure)
    site = binding_site_residues(proteins, ligand_atoms, site_cutoff)
    if not site:
        raise ValueError("The selected ligand has no protein atoms within the binding-site cutoff.")
    paths = scan_paths(proteins, ligand_atoms, directions, step_size, path_length, clash_tolerance)
    best = paths[0]
    clear_fraction = sum(p.clear for p in paths) / len(paths)
    label, score, evidence = classify(best, clear_fraction)
    location, obstructors = summarize_obstruction(best)
    if best.clear:
        interpretation = (
            "At least one tested rigid translation reaches the exterior without meaningful van der Waals overlap. "
            "The supplied conformation is geometrically compatible with rigid entry or withdrawal along a tested route."
        )
    else:
        names = ", ".join(item["residue"] for item in obstructors[:3]) or "local protein atoms"
        interpretation = (
            f"The least-obstructed rigid route is blocked primarily by {names}. Protein motion, ligand rotation, "
            "ligand conformational change, or another route would be required under this model. This is consistent "
            "with a ligand-accommodated endpoint but does not distinguish induced fit from conformational selection."
        )
    return AnalysisResult(
        input_path=str(structure.path), ligand=selected, binding_site_residues=site,
        classification=label, association_compatibility_score=score,
        evidence_strength=evidence, best_path=best, directions_tested=directions,
        clear_direction_fraction=clear_fraction, obstruction_location=location,
        primary_obstructing_residues=obstructors, interpretation=interpretation,
        limitations=[
            "Protein and ligand are rigid; ligand rotation and conformational changes are not sampled.",
            "Only straight-line paths are tested.",
            "The association compatibility score and classification thresholds are heuristic and not calibrated probabilities.",
            "A single holo structure cannot distinguish induced fit from conformational selection or determine a kinetic barrier.",
            "Hydrogens, electrostatics, solvent, and energetic effects are not modeled.",
        ], warnings=structure.warnings,
        parameters={"directions": directions, "step_size_angstrom": step_size,
                    "site_cutoff_angstrom": site_cutoff, "clash_tolerance": clash_tolerance,
                    "path_length_angstrom": best.profile[-1].distance},
        _structure=structure, _ligand_atoms=ligand_atoms,
    )
