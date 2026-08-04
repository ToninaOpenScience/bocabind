# Methodology

## Input and ligand identification

BocaBind reads PDB `ATOM` and `HETATM` coordinate records and retains the first/main alternate location. Candidate ligands are non-water, non-ion heterogen residues. Candidate ranking combines heavy-atom count, number of nearby protein atoms, fraction of ligand atoms contacting protein, and penalties for common crystallization additives or very small heterogens.

The reported candidate confidence is the softmax-normalized ranking score among heterogens in that file. It is a heuristic selection confidence, not a probability that the candidate is biologically relevant.

## Extraction geometry

For a ligand atom at bound coordinate \(\mathbf{x}_i(0)\) and sampled unit direction \(\mathbf{u}\), a path is:

\[
\mathbf{x}_i(s) = \mathbf{x}_i(0) + s\mathbf{u}.
\]

Directions are distributed on a Fibonacci sphere. At each translation step, BocaBind queries nearby protein atoms and computes penetration:

\[
\delta_{ij}(s) = \max\{0,\alpha(r_i + r_j) - d_{ij}(s)\},
\]

where \(r_i\) and \(r_j\) are element-specific van der Waals radii, \(d_{ij}\) is atomic separation, and \(\alpha\) is the clash-tolerance scale (default 0.80).

Each path records maximum penetration, integrated total penetration, peak clashing pairs, maximum fraction of ligand atoms obstructed, and the length of the obstructed interval. Paths are ranked lexicographically: clear paths first, then penetration severity, blocked fraction, cumulative obstruction, and obstruction length.

## Reproducibility

The direction set is deterministic. `report.json` records input checksum, parameters, software version, and rule version. The same input and parameters should produce the same result on supported NumPy/SciPy versions, apart from insignificant floating-point differences.

