from PdbEditor import Protein

mol = Protein("1SRL.pdb")

for resSeq, res in mol.residues.items():
	print(res.getCa())

