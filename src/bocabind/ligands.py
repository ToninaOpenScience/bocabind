from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from .models import Atom, LigandCandidate, Structure
from .pdb import InputError

STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "SEC", "PYL", "DA", "DC", "DG", "DT", "A", "C", "G", "U",
}
WATER = {"HOH", "WAT", "DOD"}
ADDITIVES = {
    "GOL", "EDO", "PEG", "PG4", "MPD", "DMS", "ACT", "FMT", "ACE", "TRS",
    "MES", "HEP", "BME", "SO4", "PO4", "NO3", "EOH", "IPA",
}
IONS = {
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "CD",
    "HG", "BR", "I",
}


def protein_atoms(structure: Structure) -> list[Atom]:
    return [a for a in structure.atoms if a.record == "ATOM" and a.element != "H"]


def ligand_groups(structure: Structure) -> dict[str, list[Atom]]:
    groups: dict[str, list[Atom]] = defaultdict(list)
    for atom in structure.atoms:
        if atom.record != "HETATM" or atom.element == "H" or atom.residue_name in WATER:
            continue
        groups[atom.residue_id].append(atom)
    return dict(groups)


def detect_ligands(structure: Structure) -> list[LigandCandidate]:
    proteins = protein_atoms(structure)
    if not proteins:
        raise InputError("No protein heavy atoms were found.")
    tree = cKDTree(np.array([a.coordinate for a in proteins]))
    raw = []
    for identifier, atoms in ligand_groups(structure).items():
        resname = atoms[0].residue_name
        if resname in STANDARD or resname in IONS:
            continue
        coords = np.array([a.coordinate for a in atoms])
        contact_indices = set()
        close_atoms = 0
        for point in coords:
            near = tree.query_ball_point(point, 4.5)
            contact_indices.update(near)
            if near:
                close_atoms += 1
        heavy = len(atoms)
        contacts = len(contact_indices)
        buried = close_atoms / max(heavy, 1)
        additive_penalty = 3.0 if resname in ADDITIVES else 0.0
        tiny_penalty = max(0, 6 - heavy) * 0.5
        score = 0.12 * min(heavy, 40) + 0.10 * min(contacts, 40) + 2.0 * buried - additive_penalty - tiny_penalty
        role = "likely crystallization additive" if resname in ADDITIVES else (
            "bound small-molecule candidate" if heavy >= 8 and contacts >= 4 else "low-confidence heterogen"
        )
        raw.append((identifier, atoms, heavy, contacts, buried, score, role))
    if not raw:
        return []
    values = np.array([r[5] for r in raw], dtype=float)
    # Independent heuristic strength, not a calibrated probability and not
    # forced to 100% merely because a file contains only one heterogen.
    probabilities = 1.0 / (1.0 + np.exp(-(values - 3.0)))
    candidates = []
    for item, confidence in zip(raw, probabilities):
        identifier, atoms, heavy, contacts, buried, score, role = item
        a = atoms[0]
        candidates.append(LigandCandidate(identifier, a.residue_name, a.chain, a.residue_number,
                                          heavy, contacts, buried, float(score), float(confidence), role))
    return sorted(candidates, key=lambda candidate: candidate.score, reverse=True)


def choose_candidate(candidates: list[LigandCandidate], selector: str | None = None) -> LigandCandidate:
    if not candidates:
        raise InputError("No supported bound small-molecule ligand was detected. Apo structures are not supported in BocaBind v0.1.")
    if selector:
        normalized = selector.upper().replace(" ", "")
        matches = [c for c in candidates if c.identifier.upper() == normalized or c.residue_name == normalized]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise InputError(f"Ligand selector '{selector}' matches multiple residues; use RESNAME:CHAIN:NUMBER.")
        raise InputError(f"Ligand selector '{selector}' did not match any candidate.")
    if len(candidates) == 1 or (
        candidates[0].confidence >= 0.80 and candidates[0].score - candidates[1].score >= 2.0
    ):
        return candidates[0]
    raise InputError("Ligand selection is ambiguous; specify --ligand RESNAME:CHAIN:NUMBER.")


def atoms_for_candidate(structure: Structure, candidate: LigandCandidate) -> list[Atom]:
    return ligand_groups(structure)[candidate.identifier]
