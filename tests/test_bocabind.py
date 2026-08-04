from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from bocabind import analyze, inspect_structure
from bocabind.pdb import InputError, resolve_structure


def pdb_line(record, serial, atom, residue, chain, number, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {atom:<4} {residue:>3} {chain:1}{number:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}          {element:>2}\n"
    )


def write_open_complex(path: Path, include_additive: bool = False) -> Path:
    lines = []
    serial = 1
    protein_points = [(-3.4, 0, 0), (-3.2, 2.8, 0), (-3.2, -2.8, 0), (-3.0, 0, 2.8)]
    for residue, point in enumerate(protein_points, 1):
        lines.append(pdb_line("ATOM", serial, "CA", "ALA", "A", residue, *point, "C"))
        serial += 1
    for atom, point, element in [("C1", (0, 0, 0), "C"), ("O1", (0, 1.2, 0), "O"), ("N1", (0, 0, 1.2), "N"),
                                 ("C2", (1.2, 0, 0), "C"), ("C3", (1.2, 1.2, 0), "C"), ("C4", (1.2, 0, 1.2), "C"),
                                 ("C5", (0, 1.2, 1.2), "C"), ("C6", (1.2, 1.2, 1.2), "C")]:
        lines.append(pdb_line("HETATM", serial, atom, "LIG", "A", 401, *point, element))
        serial += 1
    if include_additive:
        for i in range(4):
            lines.append(pdb_line("HETATM", serial, f"C{i+1}", "GOL", "A", 500, 15 + i, 0, 0, "C"))
            serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def write_caged_complex(path: Path) -> Path:
    lines = []
    serial = 1
    residue = 1
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for z in (-1, 0, 1):
                if (x, y, z) == (0, 0, 0):
                    continue
                vector = np.array([x, y, z], dtype=float)
                vector = vector / np.linalg.norm(vector) * 4.0
                lines.append(pdb_line("ATOM", serial, "CA", "ALA", "A", residue, *vector, "C"))
                serial += 1
                residue += 1
    ligand_points = [(0, 0, 0), (0.8, 0, 0), (-0.8, 0, 0), (0, 0.8, 0),
                     (0, -0.8, 0), (0, 0, 0.8), (0, 0, -0.8), (0.6, 0.6, 0.6)]
    for i, point in enumerate(ligand_points, 1):
        lines.append(pdb_line("HETATM", serial, f"C{i}", "LIG", "A", 401, *point, "C"))
        serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_inspect_ranks_bound_ligand_over_additive(tmp_path):
    pdb = write_open_complex(tmp_path / "open.pdb", include_additive=True)
    result = inspect_structure(pdb)
    assert result.ligand_candidates[0].identifier == "LIG:A:401"
    assert result.ligand_candidates[0].confidence > result.ligand_candidates[1].confidence


def test_open_structure_has_clear_path(tmp_path):
    pdb = write_open_complex(tmp_path / "open.pdb")
    result = analyze(pdb, directions=32, step_size=0.5, path_length=7.0)
    assert result.classification == "open / extraction-compatible"
    assert result.best_path.clear
    assert result.association_compatibility_score >= 95


def test_cage_is_classified_as_obstructed(tmp_path):
    pdb = write_caged_complex(tmp_path / "caged.pdb")
    result = analyze(pdb, directions=48, step_size=0.5, path_length=7.0)
    assert not result.best_path.clear
    assert result.classification in {"partially closed", "ligand-accommodated / closed"}
    assert result.primary_obstructing_residues


def test_apo_input_is_rejected(tmp_path):
    pdb = tmp_path / "apo.pdb"
    pdb.write_text(pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0, 0, 0, "C") + "END\n")
    with pytest.raises(InputError, match="Apo structures"):
        analyze(pdb, directions=12, path_length=4.0)


def test_outputs_and_cli(tmp_path):
    pdb = write_open_complex(tmp_path / "open.pdb")
    result = analyze(pdb, directions=16, step_size=1.0, path_length=5.0)
    output = result.write(tmp_path / "results")
    for name in ["report.json", "summary.txt", "path_profile.csv", "obstructing_residues.csv", "best_path.pdb"]:
        assert (output / name).is_file()
    report = json.loads((output / "report.json").read_text())
    assert report["ligand"]["identifier"] == "LIG:A:401"
    completed = subprocess.run(
        [sys.executable, "-m", "bocabind.cli", "inspect", str(pdb)],
        text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0
    assert "LIG:A:401" in completed.stdout


def test_pdb_id_is_downloaded_and_cached(tmp_path):
    remote_pdb = write_open_complex(tmp_path / "remote.pdb").read_bytes()
    response = MagicMock()
    response.__enter__.return_value.read.return_value = remote_pdb
    response.__exit__.return_value = False
    cache = tmp_path / "cache"

    with patch("bocabind.pdb.urlopen", return_value=response) as mocked_urlopen:
        downloaded = resolve_structure("1hvr", cache_dir=cache)

    assert downloaded == (cache / "1HVR.pdb").resolve()
    assert downloaded.read_bytes() == remote_pdb
    request = mocked_urlopen.call_args.args[0]
    assert request.full_url == "https://files.rcsb.org/download/1HVR.pdb"

    with patch("bocabind.pdb.urlopen") as mocked_urlopen:
        assert resolve_structure("1HVR", cache_dir=cache) == downloaded
        mocked_urlopen.assert_not_called()


def test_missing_input_explains_supported_sources(tmp_path):
    with pytest.raises(InputError, match="existing PDB file or a four-character PDB ID"):
        resolve_structure(tmp_path / "missing.pdb")
