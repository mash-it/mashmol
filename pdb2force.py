from PdbEditor import *

mol = GoProtein("1SRL.pdb")


for k in list(mol.residues.keys())[:-1]:
	print(mol.nativeBondLength(k, k+1))
