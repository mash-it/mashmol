import itertools
import numpy as np
import warnings

# global constant 
degPerRadian = 180. / np.pi
NATCONT_THRESHOLD = 6.5

# angle between two vector
def vecAngle(v1, v2):
	return np.arccos(np.dot(v1,v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

# remove spaces at both sides from string
def strstrip(s):
	return str(s).strip()

class PdbMolecule:
	def __init__(self, filename):
		self.filename = filename
		self.file = open(filename)
		self.readAtoms()
		self.pos = np.array([atom['pos'] for atom in self.atoms])
	
	def readAtoms(self):
		self.file.seek(0)
		self.atoms = []
		for line in self.file:
			if len(line) > 4 and line[0:4] == "ATOM":
				self.atoms.append(self.readAtomLine(line))
	
	def readAtomLine(self, line):
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

class Protein(PdbMolecule):
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
		self.getCaPositions()
		self.getNativeInfo()
		self.checkMissingResidues()

	def checkMissingResidues(self):
		prev = False
		for resSeq in self.residues.keys():
			if prev and resSeq - prev > 1:
				warnings.warn("There are some missing residues between {} and {}".format(prev, resSeq))
			prev = resSeq

	def getCaPositions(self):
		self.positions = []
		for seq, res in self.residues.items():
			self.positions.append({"resSeq": seq, "position": res.getCa()['pos'].tolist()})

	def getNativeInfo(self):
		res = list(self.residues.keys())

		self.nativeBond = []
		for i in range(len(res)-1):
			distance = self.getCaDistance(res[i], res[i+1])
			self.nativeBond.append({"resSeq": [res[i], res[i+1]], "length": distance})

		self.nativeAngle = []
		for i in range(len(res)-2):
			angle = self.getCaAngle(res[i], res[i+1], res[i+2]) * degPerRadian
			self.nativeAngle.append({"resSeq": [res[i], res[i+1], res[i+2]], "angle": angle})

		self.nativeDihedral = []
		for i in range(len(res)-3):
			dihedral = self.getCaDihedral(res[i], res[i+1], res[i+2], res[i+3]) * degPerRadian
			self.nativeDihedral.append({"resSeq": [res[i], res[i+1], res[i+2], res[i+3]], "dihedral":dihedral})

		self.nativeContact = []
		for i in range(len(res)):
			for j in range(i+4, len(res)):
				if self.isNativeContact(res[i], res[j], NATCONT_THRESHOLD):
					distance = self.getCaDistance(res[i], res[j])
					self.nativeContact.append({"resSeq": [res[i], res[j]], "length": distance})
		
	def getCaDistance(self, a, b):
		# a, b: integer, resSeq
		posA = self.residues[a].getCa()['pos']
		posB = self.residues[b].getCa()['pos']
		distance = np.linalg.norm(posB - posA)
		return distance
	
	def getCaAngle(self, a, b, c):
		# a, b, c: integer, resSeq
		posA = self.residues[a].getCa()['pos']
		posB = self.residues[b].getCa()['pos']
		posC = self.residues[c].getCa()['pos']
		angle = vecAngle(posA-posB, posC-posB)
		return angle

	def getCaDihedral(self, a, b, c, d):
		# a, b, c, d: integer, resSeq
		# 二面角は二平面の法線ベクトルの積として計算する
		posA = self.residues[a].getCa()['pos']
		posB = self.residues[b].getCa()['pos']
		posC = self.residues[c].getCa()['pos']
		posD = self.residues[d].getCa()['pos']
		ab = posB - posA
		bc = posC - posB
		cd = posD - posC
		abc = np.cross(ab, bc)
		bcd = np.cross(bc, cd)
		dihedral = vecAngle(abc, bcd)
		# 平面の表裏を考えると二面角にはプラスとマイナスがある
		# A-BxB-C の外積と C-D が同じ向きの場合をプラスと定義する
		if np.dot(abc, posD-posC) < 0:
			dihedral *= -1
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

