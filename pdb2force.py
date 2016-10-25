from PdbEditor import *
import json

mol = GoProtein("1SRL.pdb")

output = {"_comment" : "This file is generated from mashmol 0.9"}
output["resID"] = list(mol.residues.keys())
output["position"] = mol.positions
output["bond"] = mol.nativeBond
output["angle"] = mol.nativeAngle
output["dihedral"] = mol.nativeDihedral
output["contact"] = mol.nativeContact

print(json.dumps(output, indent=4, sort_keys=True))

