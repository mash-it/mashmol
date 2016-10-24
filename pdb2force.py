from PdbEditor import *
import json

mol = GoProtein("1SRL.pdb")

output = {"_comment" : "This file is generated from mashmol 0.9"}

print(json.dumps(output))

