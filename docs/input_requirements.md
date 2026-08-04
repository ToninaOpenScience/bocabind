# Input requirements

- a four-character RCSB PDB ID, or a local PDB file with standard fixed-width coordinate records
- at least one protein chain represented by `ATOM` records
- at least one noncovalent bound ligand represented by `HETATM` records
- element symbols in columns 77–78, or conventional atom names from which elements can be inferred

PDB-ID inputs are downloaded from RCSB and cached in `~/.cache/bocabind/structures` by default.
The cache location can be changed with the `BOCABIND_CACHE_DIR` environment variable.

Not supported in `v0.1.1`: apo structures, mmCIF, covalent ligands, multi-residue ligands, and automatic biological-assembly reconstruction.
