from __future__ import annotations

import os
import re
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np

from .models import Atom, Structure


class InputError(ValueError):
    """Raised when an input structure cannot be analyzed defensibly."""


PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def _default_cache_dir() -> Path:
    configured = os.environ.get("BOCABIND_CACHE_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".cache" / "bocabind" / "structures"


def resolve_structure(source: str | Path, cache_dir: str | Path | None = None) -> Path:
    """Resolve a local PDB path or download a four-character PDB ID from RCSB."""
    path = Path(source).expanduser()
    if path.is_file():
        return path.resolve()

    source_text = str(source).strip()
    if not PDB_ID_PATTERN.fullmatch(source_text):
        raise InputError(
            f"Structure not found: {source}. Provide an existing PDB file or a "
            "four-character PDB ID such as 1HVR."
        )

    pdb_id = source_text.upper()
    destination_dir = Path(cache_dir).expanduser() if cache_dir else _default_cache_dir()
    destination = destination_dir / f"{pdb_id}.pdb"
    if destination.is_file() and destination.stat().st_size > 0:
        return destination.resolve()

    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    request = Request(url, headers={"User-Agent": "BocaBind/0.1.1"})
    try:
        with urlopen(request, timeout=30) as response:
            payload = response.read()
    except HTTPError as exc:
        if exc.code == 404:
            raise InputError(f"PDB ID not found in RCSB: {pdb_id}") from exc
        raise InputError(f"RCSB returned HTTP {exc.code} while downloading {pdb_id}.") from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise InputError(
            f"Could not download {pdb_id} from RCSB. Check your internet connection "
            "or provide a local PDB file."
        ) from exc

    if not payload.strip():
        raise InputError(f"RCSB returned an empty structure for {pdb_id}.")

    try:
        destination_dir.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(".pdb.part")
        temporary.write_bytes(payload)
        temporary.replace(destination)
    except OSError as exc:
        raise InputError(f"Could not cache downloaded structure at {destination}.") from exc
    return destination.resolve()


def _element(line: str, atom_name: str) -> str:
    value = line[76:78].strip().upper() if len(line) >= 78 else ""
    if value:
        return value
    letters = "".join(c for c in atom_name if c.isalpha()).upper()
    if not letters:
        return "C"
    two_letter = {"CL", "BR", "ZN", "FE", "MG", "MN", "CA", "NA", "CU", "NI", "CO"}
    return letters[:2] if letters[:2] in two_letter else letters[0]


def read_pdb(path: str | Path) -> Structure:
    path = resolve_structure(path)
    atoms: list[Atom] = []
    warnings: list[str] = []
    seen_altloc: set[tuple] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                atom_name = line[12:16].strip()
                altloc = line[16:17].strip()
                residue_name = line[17:20].strip().upper()
                chain = line[21:22].strip()
                residue_number = int(line[22:26])
                insertion_code = line[26:27].strip()
                coordinate = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float
                )
                occupancy = float(line[54:60].strip() or 1.0)
                serial = int(line[6:11])
            except (ValueError, IndexError) as exc:
                raise InputError(f"Malformed atom record on line {line_number}") from exc
            key = (record, chain, residue_number, insertion_code, residue_name, atom_name)
            if altloc not in {"", "A", "1"}:
                if key not in seen_altloc:
                    warnings.append(f"Ignored alternate location {altloc} for {residue_name}:{chain}:{residue_number}.")
                continue
            seen_altloc.add(key)
            atoms.append(
                Atom(serial, atom_name, _element(line, atom_name), residue_name, chain,
                     residue_number, insertion_code, coordinate, record, occupancy, altloc)
            )
    if not atoms:
        raise InputError("No ATOM or HETATM coordinate records were found.")
    if not any(a.record == "ATOM" for a in atoms):
        raise InputError("No protein ATOM records were found.")
    return Structure(path=path, atoms=atoms, warnings=warnings)
