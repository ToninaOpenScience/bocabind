from bocabind import analyze, inspect_structure


inspection = inspect_structure("complex.pdb")
for candidate in inspection.ligand_candidates:
    print(candidate)

result = analyze("complex.pdb", ligand=inspection.ligand_candidates[0].identifier)
result.write("bocabind_results")
print(result.classification)

