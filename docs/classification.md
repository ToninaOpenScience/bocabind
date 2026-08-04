# Structural classification

Classification is based on the least-obstructed sampled rigid path.

| Classification | Current heuristic rule | Structural meaning |
|---|---|---|
| Open / extraction-compatible | Maximum penetration ≤ 0.10 Å | At least one sampled rigid path has no meaningful overlap |
| Narrowly accessible | Maximum penetration ≤ 0.35 Å and obstruction ≤ 2 Å | Only brief, minor overlap occurs |
| Partially closed | Maximum penetration ≤ 1.0 Å and obstruction ≤ 6 Å | Local rearrangement would be required |
| Ligand-accommodated / closed | More severe or extended overlap | No tested rigid path is geometrically compatible |

These thresholds are versioned experimental rules. They require validation against curated structural systems before use as calibrated predictors.

“Ligand-accommodated / closed” is a structural description. “Induced fit” and “conformational selection” are mechanistic hypotheses that cannot be distinguished from one holo structure.

