from PdbEditor import Molecule

mol = Molecule("1SRL.pdb")

for resSeq, res in mol.residues.items():
	print(res.get_ca())

