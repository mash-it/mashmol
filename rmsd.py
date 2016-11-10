import sys
import warnings
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser()
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	structure = parser.get_structure("traj", sys.argv[1])

reference = list(structure[0].get_atoms())

for model in structure:
	sup = Superimposer()
	sup.set_atoms(reference, list(model.get_atoms()))
	print(sup.rms)

