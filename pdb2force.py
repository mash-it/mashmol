from PdbEditor import *
from collections import OrderedDict
import json
from datetime import datetime

VERSION = 0.9
inputfile = "1SRL.pdb"
mol = GoProtein("1SRL.pdb")
date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

output = OrderedDict()

output["_comment"] = OrderedDict()
output["_comment"]["generated_by"] = "mashmol {}".format(VERSION)
output["_comment"]["generated_from"] = inputfile
output["_comment"]["generated_at"] = date
output["_comment"]["unit_length"] = "angstrom"
output["_comment"]["unit_angle"] = "degree"

output["summary"] = OrderedDict()
output["summary"]["N_atoms"] = len(mol.positions)
output["summary"]["N_bonds"] = len(mol.nativeBond)
output["summary"]["N_angle"] = len(mol.nativeAngle)
output["summary"]["N_dihedral"] = len(mol.nativeDihedral)
output["summary"]["N_contact"] = len(mol.nativeContact)

output["parameters"] = OrderedDict()
output["parameters"]["Temperature"] = 300.0
output["parameters"]["TimePerStepInPs"] = 0.02
output["parameters"]["SimulationSteps"] = 10000
output["parameters"]["NStepSave"] = 100

output["resSeq"] = list(mol.residues.keys())
output["position"] = mol.positions
output["bond"] = mol.nativeBond
output["angle"] = mol.nativeAngle
output["dihedral"] = mol.nativeDihedral
output["contact"] = mol.nativeContact

print(json.dumps(output, indent=4))

