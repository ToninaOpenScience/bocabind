# bocabind

Determines if a binding site is open or closed for small molecule binding.

Classifies binding-site accessibility by quantifying steric obstruction along rigid ligand-extraction paths from protein–ligand structures.

BocaBind tests a deliberately narrow first-principles question: can the ligand in a supplied holo structure be rigidly withdrawn along a straight path without intersecting the protein? It automatically ranks likely bound ligands, defines the local binding site, samples extraction directions, quantifies van der Waals overlap along each path, identifies blocking residues, and explains the structural classification.

> **Status: research MVP (`v0.1.1`).** The output is a rigid-body geometric classification—not a binding pathway, energy barrier, or probability of induced fit.

## What it produces

- automatic ligand candidate ranking and a heuristic selection confidence;
- binding-site residues within a configurable cutoff;
- 256 evenly distributed rigid extraction directions by default;
- maximum and cumulative penetration, clashing pairs, blocked ligand fraction, and obstruction length;
- classification as **open / extraction-compatible**, **narrowly accessible**, **partially closed**, or **ligand-accommodated / closed**;
- the location and residue composition of the best path's bottleneck;
- a JSON report, readable summary, CSV profiles, and multi-model PDB trajectory.

## Requirements

- Python 3.10 or newer
- Internet access for PDB-ID downloads, or a local PDB file containing a protein and at least one bound, noncovalent small molecule

Apo structures, covalent ligands, mmCIF inputs, curved paths, protein flexibility, and ligand rotation are outside `v0.1.1`.

## Install locally

### macOS or Linux

```bash
unzip bocabind-v0.1.1.zip
cd bocabind

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

### Windows PowerShell

```powershell
Expand-Archive bocabind-v0.1.1.zip
cd bocabind-v0.1.1\bocabind

py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

Confirm installation:

```bash
bocabind --version
python -c "import bocabind; print(bocabind.__version__)"
```

## Test the repository

Run the full test suite:

```bash
python -m pytest
```

The suite uses synthetic structures and mocks the RCSB download, so it does not need internet access.

Run a complete smoke test with the included synthetic holo structure:

```bash
bocabind inspect tests/fixtures/synthetic_open_complex.pdb
bocabind analyze tests/fixtures/synthetic_open_complex.pdb \
  --directions 32 \
  --path-length 7 \
  --output demo_results

cat demo_results/summary.txt
```

The fixture exists only to verify installation and file generation; do not use it for scientific validation.

Build the installable source and wheel packages:

```bash
python -m pip install build
python -m build
```

## Analyze a real holo PDB

Give BocaBind either a four-character RCSB PDB ID:

```bash
bocabind inspect 1HVR
bocabind analyze 1HVR --output results
```

or a local experimental or prepared holo PDB:

```bash
bocabind inspect path/to/complex.pdb
```

Then run the analysis:

```bash
bocabind analyze path/to/complex.pdb --output results
```

If selection is ambiguous, BocaBind asks interactively. For a reproducible run, provide the candidate explicitly:

```bash
bocabind analyze path/to/complex.pdb \
  --ligand X77:A:401 \
  --directions 256 \
  --step-size 0.5 \
  --output results
```

For automated workflows, add `--no-interactive`; an ambiguous ligand then produces an error instead of pausing.

The output directory contains:

```text
results/
├── report.json
├── summary.txt
├── path_profile.csv
├── obstructing_residues.csv
└── best_path.pdb
```

Open `best_path.pdb` as a trajectory in PyMOL or VMD to inspect the least-obstructed tested route.

## Python API

```python
from bocabind import analyze, inspect_structure

inspection = inspect_structure("complex.pdb")
for candidate in inspection.ligand_candidates:
    print(candidate.identifier, candidate.confidence, candidate.likely_role)

result = analyze(
    "complex.pdb",
    ligand="X77:A:401",
    directions=256,
    step_size=0.5,
)

print(result.classification)
print(result.association_compatibility_score)
print(result.primary_obstructing_residues)
result.write("results")
```

The Python API never requests terminal input. Callers must give `ligand=` when automatic selection is ambiguous.

Downloaded structures are cached in `~/.cache/bocabind/structures` and reused. Set
`BOCABIND_CACHE_DIR` to choose a different cache directory.

## Main parameters

| Option              |   Default | Meaning                                                         |
| ------------------- | --------: | --------------------------------------------------------------- |
| `--directions`      |       256 | Straight-line directions sampled on a Fibonacci sphere          |
| `--step-size`       |     0.5 Å | Translation increment                                           |
| `--site-cutoff`     |     5.0 Å | Ligand–protein distance used to define site residues            |
| `--clash-tolerance` |      0.80 | Scale factor applied to summed van der Waals radii              |
| `--path-length`     | automatic | Distance followed toward the structure exterior, capped at 45 Å |

## Scientific interpretation

The **Association Compatibility Score** is a transparent heuristic derived from the least-obstructed path's penetration depth, blocked ligand fraction, cumulative overlap, and obstructed length. A score of 90 means low geometric obstruction under the tested rigid model; it does **not** mean a 90% probability that the deposited structure is biologically open.

Likewise, a ligand-accommodated/closed structure is compatible with several explanations, including induced fit, conformational selection, ligand reorientation, side-chain gating, or an untested curved route. A single holo structure cannot distinguish among them.

See [docs/methodology.md](docs/methodology.md), [docs/classification.md](docs/classification.md), and [docs/limitations.md](docs/limitations.md).

## Development

```bash
python -m pip install -e ".[dev]"
ruff check src tests
python -m pytest --cov=bocabind
```

Contributions should include a test and should not strengthen the scientific claim beyond what the implemented geometry supports. See [CONTRIBUTING.md](CONTRIBUTING.md).

## License and citation

BocaBind is licensed under the [Mozilla Public License 2.0](LICENSE). See [CITATION.cff](CITATION.cff) for citation metadata.
