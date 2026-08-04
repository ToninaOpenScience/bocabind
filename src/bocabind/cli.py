from __future__ import annotations

import argparse
import sys

from .api import analyze, inspect_structure
from .pdb import InputError, PDB_ID_PATTERN, resolve_structure


def _candidate_table(inspection) -> None:
    if not inspection.ligand_candidates:
        print("No supported ligand candidates detected.")
        return
    print("Ligand candidates\n")
    for index, candidate in enumerate(inspection.ligand_candidates, 1):
        print(f"[{index}] {candidate.identifier}")
        print(f"    Heavy atoms: {candidate.heavy_atoms}")
        print(f"    Protein contacts: {candidate.protein_contacts}")
        print(f"    Candidate confidence: {candidate.confidence:.1%} (heuristic)")
        print(f"    Estimated role: {candidate.likely_role}")
    if inspection.automatic_selection:
        print(f"\nAutomatic selection: {inspection.automatic_selection}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bocabind", description="Rigid ligand-extraction accessibility analysis")
    parser.add_argument("--version", action="version", version="BocaBind 0.1.1")
    commands = parser.add_subparsers(dest="command", required=True)
    inspect_cmd = commands.add_parser("inspect", help="detect ligand candidates")
    inspect_cmd.add_argument("pdb", help="local PDB file or four-character RCSB PDB ID")
    analyze_cmd = commands.add_parser("analyze", help="run a rigid extraction scan")
    analyze_cmd.add_argument("pdb", help="local PDB file or four-character RCSB PDB ID")
    analyze_cmd.add_argument("--ligand", help="RESNAME or RESNAME:CHAIN:NUMBER")
    analyze_cmd.add_argument("--directions", type=int, default=256)
    analyze_cmd.add_argument("--step-size", type=float, default=0.5)
    analyze_cmd.add_argument("--site-cutoff", type=float, default=5.0)
    analyze_cmd.add_argument("--clash-tolerance", type=float, default=0.80)
    analyze_cmd.add_argument("--path-length", type=float)
    analyze_cmd.add_argument("--output", default="bocabind_results")
    analyze_cmd.add_argument("--no-interactive", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        source = args.pdb
        if PDB_ID_PATTERN.fullmatch(str(source).strip()):
            print(f"Resolving PDB ID {str(source).upper()} from RCSB...")
        resolved_source = resolve_structure(source)
        inspection = inspect_structure(resolved_source)
        if args.command == "inspect":
            _candidate_table(inspection)
            return 0 if inspection.ligand_candidates else 2
        selector = args.ligand
        if selector is None and inspection.automatic_selection is None and not args.no_interactive:
            _candidate_table(inspection)
            selection = input("\nSelect ligand number or identifier: ").strip()
            if selection.isdigit() and 1 <= int(selection) <= len(inspection.ligand_candidates):
                selector = inspection.ligand_candidates[int(selection) - 1].identifier
            else:
                selector = selection
        result = analyze(resolved_source, selector, directions=args.directions, step_size=args.step_size,
                         site_cutoff=args.site_cutoff, clash_tolerance=args.clash_tolerance,
                         path_length=args.path_length)
        result.write(args.output)
        print(f"Selected ligand: {result.ligand.identifier} ({result.ligand.confidence:.1%} heuristic confidence)")
        print(f"Structural classification: {result.classification}")
        print(f"Association Compatibility Score: {result.association_compatibility_score:.1f}/100 (heuristic)")
        print(f"Results written to: {args.output}")
        return 0
    except (InputError, ValueError) as exc:
        print(f"BocaBind stopped: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
