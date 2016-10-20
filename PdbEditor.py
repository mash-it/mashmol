import itertools
import numpy as np

# global constant 
degPerRadian = 180. / np.pi

# angle between two vector
def vecAngle(v1, v2):
	return np.arccos(np.sum(v1*v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

# remove spaces at both sides from string
def strstrip(s):
	return str(s).strip()

class Molecule:
	def __init__(self, filename):
		self.filename = filename
		self.file = open(filename)
		self.readAtoms()
	
	def readAtoms(self):
		self.file.seek(0)
		self.atoms = []
		for line in self.file:
			if len(line) > 4 and line[0:4] == "ATOM":
				self.atoms.append(self.readAtomline(line))
	
	def readAtomline(self, line):
		# to avoid IndexError
		line += " " * 80

		columns = [
			 ("serial", 7,11, int)
			,("atomname", 13,16, strstrip)
			,("altLoc", 17,17, strstrip)
			,("resName", 18,20, strstrip)
			,("chainID", 22,22, strstrip)
			,("resSeq", 23,26, int)
			,("iCode", 27,27, strstrip)
			,("x", 31,38, float)
			,("y", 39,46, float)
			,("z", 47,53, float)
			,("occupancy", 55,60, float)
			,("tempFactor", 61,66, float)
			,("element", 77,78, strstrip)
			,("charge", 79,80, strstrip)
		]

		atom = {}
		for c in columns:
			try:
				atom[c[0]] = c[3](line[c[1]-1: c[2]])
			except ValueError:
				pass
				
		atom['pos'] = np.array((atom['x'], atom['y'], atom['z']))
		return atom
	
	def getBoxSize(self):
		xs = [ atom['x'] for atom in self.atoms ]
		ys = [ atom['y'] for atom in self.atoms ]
		zs = [ atom['z'] for atom in self.atoms ]

		return ((min(xs), max(xs)),(min(ys), max(ys)),(min(zs), max(zs)))
		
	def shift(self, distance, direction):
		for atom in self.atoms:
			atom[distance] += direction

	def getPdbText(self):
		lines = []
		for atom in self.atoms:
			lines.append("ATOM  {serial:5d} {atomname:4s}{altLoc:1s}{resName:3s} {chainID:1s}{resSeq:4d}{iCode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}\n".format(
				 serial = atom['serial']
				,atomname = atom['atomname']
				,altLoc = atom['altLoc']
				,resName = atom['resName']
				,chainID = atom['chainID']
				,resSeq = atom['resSeq']
				,iCode = atom['iCode']
				,x = atom['x']
				,y = atom['y']
				,z = atom['z']
			))

		return "".join(lines)
	
	def output(self, filename):
		f = open(filename, 'w')
		f.write(self.getPdbText())

class Residue(list):
	def getCa(self):
		for atom in self:
			if atom['atomname'] == 'CA':
				return atom
		raise RuntimeError("No CA in this residue")

class Protein(Molecule):
	def __init__(self, filename):
		super().__init__(filename)
		self.classifyResidues()

	def classifyResidues(self):
		self.residues = {} # map from int:resSeq to atoms

		for atom in self.atoms:
			resSeq = atom['resSeq']
			if resSeq in self.residues.keys():
				self.residues[resSeq].append(atom)
			else:
				self.residues[resSeq] = Residue()
				self.residues[resSeq].append(atom)


class GoProtein(Protein):
	""" Protein implemented as an Go model data"""
	def __init__(self, filename):
		super().__init__(filename)
		self.getNativeInfo()

	def getNativeInfo(self):
		# todo: missing of CA atom position
		res = list(self.residues.keys())

		for i in range(len(res)-1):
			print(self.getNativeBondLength(res[i], res[i+1]))

		for i in range(len(res)-2):
			print(self.getNativeAngle(res[i], res[i+1], res[i+2]))

		for i in range(len(res)-3):
			print(self.getNativeDihedral(res[i], res[i+1], res[i+2], res[i+3]))
		
	def getNativeBondLength(self, i, j):
		resI = self.residues[i].getCa()
		resJ = self.residues[j].getCa()
		distance = np.linalg.norm(resI['pos'] - resJ['pos'])
		return distance
	
	def getNativeAngle(self, i, j, k):
		resI = self.residues[i].getCa()
		resJ = self.residues[j].getCa()
		resK = self.residues[k].getCa()
		ji = resI['pos'] - resJ['pos']
		jk = resK['pos'] - resJ['pos']
		angle = vecAngle(ji, jk)
		return angle

	def getNativeDihedral(self, i, j, k, l):
		# 二面角は法線ベクトルの積である
		resI = self.residues[i].getCa()
		resJ = self.residues[j].getCa()
		resK = self.residues[k].getCa()
		resL = self.residues[l].getCa()
		ij = resJ['pos'] - resI['pos']
		jk = resK['pos'] - resJ['pos']
		kl = resL['pos'] - resK['pos']
		ijk = np.cross(ij, jk)
		jkl = np.cross(jk, kl)
		dihedral = vecAngle(ijk, jkl)
		return dihedral 

	def isNativeContact(self, resSeqA, resSeqB, threshold):
		try:
			atomsA = self.residues[resSeqA]
			atomsB = self.residues[resSeqB]
		except KeyError:
			raise RuntimeError(e + "No such residue(s) in this protein")

		for a, b in itertools.product(atomsA, atomsB):
			# remove hydrogen
			if a['element'] == 'H' or b['element'] == 'H':
				continue

			# calc distance
			distance = np.linalg.norm(a['pos'] - b['pos'])

			if distance < threshold:
				return True

		return False

