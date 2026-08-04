from __future__ import annotations

from collections import defaultdict

from .models import PathResult


def classify(best: PathResult, clear_fraction: float) -> tuple[str, float, str]:
    severity = (
        0.45 * min(best.max_penetration / 1.5, 1.0)
        + 0.25 * min(best.max_blocked_fraction / 0.5, 1.0)
        + 0.20 * min(best.cumulative_penetration / 25.0, 1.0)
        + 0.10 * min(best.obstructed_length / 8.0, 1.0)
    )
    score = round(100.0 * max(0.0, min(1.0, 1.0 - severity)), 1)
    if best.clear:
        label = "Open binding site(extraction-compatible)"
    elif best.max_penetration <= 0.35 and best.obstructed_length <= 2.0:
        label = "Narrowly accessible binding site"
    elif best.max_penetration <= 1.0 and best.obstructed_length <= 6.0:
        label = "Partially closed binding site"
    else:
        label = "Closed binding site"
    evidence = "high" if clear_fraction > 0.02 or best.max_penetration > 0.75 else "moderate"
    return label, score, evidence


def summarize_obstruction(best: PathResult) -> tuple[str, list[dict]]:
    totals: dict[str, float] = defaultdict(float)
    obstructed_steps = [s for s in best.profile if s.max_penetration > 0.10]
    for step in obstructed_steps:
        for residue, value in step.residues.items():
            totals[residue] += value
    total = sum(totals.values())
    residues = [
        {"residue": name, "penetration_sum": round(value, 3),
         "fraction_of_obstruction": round(value / total, 4) if total else 0.0}
        for name, value in sorted(totals.items(), key=lambda item: item[1], reverse=True)[:10]
    ]
    if not obstructed_steps:
        return "no meaningful obstruction", residues
    peak = max(best.profile, key=lambda s: s.total_penetration)
    path_end = best.profile[-1].distance
    fraction = peak.distance / max(path_end, 0.001)
    if peak.distance <= 2.0:
        location = "deep-pocket obstruction"
    elif fraction < 0.45:
        location = "internal-gate obstruction"
    elif fraction < 0.75:
        location = "pocket-mouth obstruction"
    else:
        location = "extended-channel obstruction"
    return location, residues

