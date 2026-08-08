# BocaBind

[![PyPI version](https://img.shields.io/pypi/v/bocabind.svg)](https://pypi.org/project/bocabind/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Status: Research MVP](https://img.shields.io/badge/status-research%20MVP%20v0.1.1-orange.svg)](#status)
[![Tests](https://img.shields.io/badge/tests-pytest-informational.svg)](#test-the-repository)

From the Tonina Open Science initiative.

**Is a protein's binding site open or closed?** BocaBind answers a narrow, well defined version of that question for protein-ligand structures: given a holo structure, can the bound ligand be rigidly withdrawn along a straight path without colliding with the protein? This is useful groundwork for structure based drug design, pocket comparison across homologs, and triaging crystal structures before running more expensive simulations.

BocaBind does not model energy, induced fit, or ligand flexibility. It answers one geometric question, quickly and transparently, and tells you exactly what it did.

---

## Table of contents

- [Quick start](#quick-start)
- [What it produces](#what-it-produces)
- [Example output](#example-output)
- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Analyzing a real PDB](#analyzing-a-real-pdb)
- [Python API](#python-api)
- [Main parameters](#main-parameters)
- [Classification categories](#classification-categories)
- [Scientific interpretation](#scientific-interpretation)
- [How this differs from other pocket and tunnel tools](#how-this-differs-from-other-pocket-and-tunnel-tools)
- [Performance](#performance)
- [Limitations](#limitations)
- [Roadmap](#roadmap)
- [Status](#status)
- [Development and contributing](#development-and-contributing)
- [Getting help](#getting-help)
- [License and citation](#license-and-citation)

---

## Quick start

```bash
pip install bocabind

# Analyze a structure directly from the PDB
bocabind analyze 1HVR --output results

# Look at the result
cat results/summary.txt
```

That's it. For a deeper walkthrough, including local files and reproducible runs, see [Analyzing a real PDB](#analyzing-a-real-pdb).

---

## What it produces

- Automatic ligand candidate ranking with a heuristic selection confidence
- Binding site residues within a configurable cutoff
- 256 evenly distributed rigid extraction directions by default
- Maximum and cumulative penetration, clashing pairs, blocked ligand fraction, and obstruction length
- Classification as open/extraction-compatible, narrowly accessible, partially closed, or ligand-accommodated/closed
- The location and residue composition of the best path's bottleneck
- A JSON report, readable summary, CSV profiles, and multi-model PDB trajectory

---

## Example output

A trimmed `summary.txt` for a well opened pocket might look like this:

```
BocaBind v0.1.1 — Analysis Summary
Structure: 1HVR
Ligand: XK2:A:401 (selection confidence: high)

Classification: open / extraction-compatible
Association compatibility score: 91.4

Best path direction: (0.12, -0.84, 0.53)
Max penetration: 0.31 Å
Blocked ligand fraction: 0.04
Bottleneck: none identified above clash tolerance

Site residues (5.0 Å cutoff): 12 residues
Obstructing residues (best path): none
```

And a closed or accommodated site:

```
Classification: ligand-accommodated / closed
Association compatibility score: 22.7

Best path direction: (0.67, 0.11, -0.73)
Max penetration: 2.84 Å
Blocked ligand fraction: 0.61
Bottleneck: Phe82, Ile84 (side chain occlusion)
```

Full machine readable results are written to `report.json`; the trajectory of the least obstructed path is written to `best_path.pdb` for visual inspection in PyMOL or VMD.

---

## How it works

1. Identify the bound ligand (automatically ranked, or specified explicitly).
2. Define the local binding site using a configurable distance cutoff.
3. Sample a set of straight, evenly distributed extraction directions on a Fibonacci sphere.
4. For each direction, rigidly translate the ligand outward and measure van der Waals overlap with the protein at each step.
5. Record penetration depth, clashing residue pairs, and blocked fraction along each path.
6. Classify the site based on the least obstructed path found, and report the residues responsible for any bottleneck.

The method deliberately does not move the protein, rotate the ligand, or follow curved routes. It is a fast, first-principles geometric test, not a simulation of unbinding.

---

## Requirements

- Python 3.10 or newer
- Internet access for PDB-ID downloads, or a local PDB file containing a protein and at least one bound, noncovalent small molecule

Apo structures, covalent ligands, mmCIF inputs, curved paths, protein flexibility, and ligand rotation are outside the scope of v0.1.1. See [Limitations](#limitations).

---

## Installation

### Pip install

```bash
pip install bocabind
```

### Install locally

**macOS or Linux**

```bash
unzip bocabind-v0.1.1.zip
cd bocabind

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

**Windows PowerShell**

```powershell
Expand-Archive bocabind-v0.1.1.zip
cd bocabind-v0.1.1\bocabind

py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

### Confirm installation

```bash
bocabind --version
python -c "import bocabind; print(bocabind.__version__)"
```

### Test the repository

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

The fixture exists only to verify installation and file generation. Do not use it for scientific validation.

Build the installable source and wheel packages:

```bash
python -m pip install build
python -m build
```

---

## Analyzing a real PDB

Give BocaBind either a four character RCSB PDB ID:

```bash
bocabind inspect 1HVR
bocabind analyze 1HVR --output results
```

or a local experimental or prepared holo PDB:

```bash
bocabind inspect path/to/complex.pdb
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

For automated workflows, add `--no-interactive`. An ambiguous ligand then produces an error instead of pausing.

The output directory contains:

```
results/
├── report.json
├── summary.txt
├── path_profile.csv
├── obstructing_residues.csv
└── best_path.pdb
```

Open `best_path.pdb` as a trajectory in PyMOL or VMD to inspect the least obstructed tested route.

---

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

The Python API never requests terminal input. Callers must supply `ligand=` when automatic selection is ambiguous.

Downloaded structures are cached in `~/.cache/bocabind/structures` and reused. Set `BOCABIND_CACHE_DIR` to choose a different cache directory.

---

## Main parameters

| Option | Default | Meaning |
|---|---|---|
| `--directions` | 256 | Straight line directions sampled on a Fibonacci sphere |
| `--step-size` | 0.5 Å | Translation increment |
| `--site-cutoff` | 5.0 Å | Ligand-protein distance used to define site residues |
| `--clash-tolerance` | 0.80 | Scale factor applied to summed van der Waals radii |
| `--path-length` | automatic | Distance followed toward the structure exterior, capped at 45 Å |

---

## Classification categories

| Classification | What it means geometrically |
|---|---|
| Open / extraction-compatible | At least one sampled path shows minimal to no obstruction; the ligand can be rigidly withdrawn essentially unimpeded |
| Narrowly accessible | A path exists but passes close to protein atoms, with limited clearance |
| Partially closed | All paths show meaningful obstruction, but the best path clears with moderate penetration |
| Ligand-accommodated / closed | Every sampled direction shows substantial steric overlap; no rigid straight line path escapes without significant clash |

These categories describe the tested rigid model only. They are not claims about the true, dynamic behavior of the protein.

---

## Scientific interpretation

The Association Compatibility Score is a transparent heuristic derived from the least obstructed path's penetration depth, blocked ligand fraction, cumulative overlap, and obstructed length. A score of 90 means low geometric obstruction under the tested rigid model. It does not mean a 90 percent probability that the deposited structure is biologically open.

Likewise, a ligand-accommodated/closed structure is compatible with several explanations, including induced fit, conformational selection, ligand reorientation, side chain gating, or an untested curved route. A single holo structure cannot distinguish among them.

See `docs/methodology.md`, `docs/classification.md`, and `docs/limitations.md` for full detail.

---

## How this differs from other pocket and tunnel tools

Tools like CAVER, MOLE, and POVME map tunnels and cavities, often accounting for flexible or curved routes and sometimes incorporating molecular dynamics. BocaBind takes a deliberately narrower approach:

- It tests one specific, already bound ligand pose rather than searching for generic cavities.
- It uses straight, rigid extraction paths only, with no protein flexibility or ligand rotation.
- It is fast and dependency light, intended as a quick first-principles triage step, not a replacement for tunnel detection or MD based unbinding studies.

If you need curved path tunnel networks, flexible protein sampling, or free energy estimates, a dedicated tunnel detection or MD tool is likely a better fit. BocaBind is a good fit when you have a specific holo structure and want a fast, interpretable, geometry only answer.

---

## Performance

Runtime scales roughly linearly with `--directions` and step resolution. On a typical single domain protein-ligand complex with default settings (256 directions, 0.5 Å step size), analysis generally completes in a few seconds to low tens of seconds on a standard laptop CPU. Larger structures, smaller step sizes, or higher direction counts will increase runtime accordingly. No GPU is required or used.

---

## Limitations

Beyond what falls outside v0.1.1 scope (apo structures, covalent ligands, mmCIF inputs, curved paths, protein flexibility, ligand rotation), a few practical points to keep in mind:

- Results depend on the quality and completeness of the input structure; missing side chains or poor resolution near the binding site can distort the classification.
- The straight line assumption means genuinely tortuous but biologically real escape routes will be missed and may be misclassified as closed.
- The clash tolerance and step size parameters can shift a structure between adjacent classification categories, so results near a boundary should be treated as such and inspected manually.

See `docs/limitations.md` for the complete list.

---

## Roadmap

v0.1.1 is a research MVP focused on establishing the core rigid-path methodology cleanly and correctly. Areas under consideration for future versions include:

- Curved path sampling
- Basic side chain flexibility
- Apo structure support
- mmCIF input support
- Covalent ligand handling

None of these are committed timelines. Interest and contributions on any of them are welcome; see [Development and contributing](#development-and-contributing).

---

## Status

Research MVP (v0.1.1). The output is a rigid body geometric classification, not a binding pathway, energy barrier, or probability of induced fit. This project is under active development and interfaces may change between minor versions.

---

## Development and contributing

```bash
python -m pip install -e ".[dev]"
ruff check src tests
python -m pytest --cov=bocabind
```

Contributions should include a test and should not strengthen the scientific claim beyond what the implemented geometry supports. See `CONTRIBUTING.md` for details.

---

## Getting help

- Bug reports and feature requests: open an issue on the project's GitHub repository
- Questions about methodology or interpretation: see `docs/methodology.md` and `docs/limitations.md` first, then open a discussion thread
- General questions about the Tonina Open Science initiative: see the initiative's main site

---

## License and citation

BocaBind is licensed under the Apache-2.0 license. See `LICENSE` for details.

If you use BocaBind in published work, please cite the project. A formal citation entry (BibTeX and DOI) will be added here once available; in the meantime, cite the repository URL and version number (v0.1.1).
