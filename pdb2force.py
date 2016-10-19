from PdbEditor import Molecule

mol = Molecule("1SRL.pdb")

for i in mol.residues.keys():
	for j in mol.residues.keys():
		if i < j-3:
			print(i, j, mol.is_native_contact(i, j, 6.0))

