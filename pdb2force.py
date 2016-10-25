from PdbEditor import *
from collections import OrderedDict
import json
from datetime import datetime

VERSION = 0.9
inputfile = "1SRL.pdb"
mol = GoProtein("1SRL.pdb")
date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

output = OrderedDict()
output["_comment"] = {
	"generated_by": "mashmol {}".format(VERSION),
	"generated_from": inputfile,
	"generated_at": date,
	"unit_length" : "angstrom",
	"unit_angle" : "degree"
	}
output["resID"] = list(mol.residues.keys())
output["position"] = mol.positions
output["bond"] = mol.nativeBond
output["angle"] = mol.nativeAngle
output["dihedral"] = mol.nativeDihedral
output["contact"] = mol.nativeContact

print(json.dumps(output, indent=4))

