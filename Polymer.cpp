#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <json.hpp>
#include <map>

void simulate();
void writePDBFrame(int, const OpenMM::State&, std::ofstream&);

using json = nlohmann::json;

// global parameters
// K: spring constant
// E: Energy
// D: Distance

const double K_BondStretch = 40000.0; // amu ps^-2
const double K_BondAngle= 80.0; // amu nm^2 ps^-2
const double K_BondTorsion1 = 4.0; // amu nm^2 ps^-2
const double K_BondTorsion3 = 2.0; // amu nm^2 ps^-2
const double E_ExclusionPair = 0.8; // amu nm^2 ps^-2
const double E_GoContactPair = 1.2; // amu nm^2 ps^-2
const double D_ExclusiveInNm = 0.4; // nm
const double D_ExclusionCutoffInNm = 2.0; // nm
const double LangevinFrictionPerPs = 5.0;
const bool UseConstraints = false;

int main() {
	simulate();
	return 0;
}

void simulate() {

	// Load any shared libraries containing GPU implementations.
	OpenMM::Platform::loadPluginsFromDirectory(
      	OpenMM::Platform::getDefaultPluginsDirectory());

	// force to use CPU platform that is better than CUDA in my simulation scale
	OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("CPU");
	// force to use single thread 
	platform.setPropertyDefaultValue("Threads", "1");
	
	// Create a system with forces.
	OpenMM::System system;

	// load global paremters and atom information from input file.
	json molinfo;
	std::ifstream ifs("input.json");
	if(ifs.fail()) {
		std::cerr << "File do not exist.\n";
		exit(1);
	}
	ifs >> molinfo;

	const int N_particles = molinfo["resSeq"].size();
	std::cout << "# Number of Particles : " << N_particles << '\n';

	// set MD paramters
	const double Temperature = molinfo["parameters"]["Temperature"];
	int SimulationSteps = molinfo["parameters"]["SimulationSteps"];
	const double TimePerStepInPs = molinfo["parameters"]["TimePerStepInPs"];
	const int NStepSave = molinfo["parameters"]["NStepSave"];
	const int RandomSeed = molinfo["parameters"]["RandomSeed"];

	// set other parameters
	const std::string filename = molinfo["output"]["filename"];

	// Create Atoms
	std::vector<OpenMM::Vec3> initPosInNm(N_particles);

	// map from resSeq to particle index(0..N-1)
	std::map<int,int> rs2pi;
	for (int i=0; i<N_particles; ++i) {
		int resSeq = molinfo["resSeq"][i];
		rs2pi.insert(std::pair<int,int>(resSeq, i));
	}

	double x, y, z;
	for (int i=0; i<N_particles; ++i) {
		x = molinfo["position"][i]["position"][0];
		y = molinfo["position"][i]["position"][1];
		z = molinfo["position"][i]["position"][2];
		initPosInNm[i] = OpenMM::Vec3(x, y, z) * OpenMM::NmPerAngstrom;
		system.addParticle(137.0);	// average mass of amino acid residues
	}

	// add non-local repulsive force
	// This pair potential is excluded from bonded pairs 
	OpenMM::CustomNonbondedForce& nonlocalRepulsion = *new OpenMM::CustomNonbondedForce("epsilon*(d_exclusive/r)^12");
	system.addForce(&nonlocalRepulsion);
	nonlocalRepulsion.setNonbondedMethod(OpenMM::CustomNonbondedForce::NonbondedMethod::CutoffNonPeriodic);
	nonlocalRepulsion.setCutoffDistance(D_ExclusionCutoffInNm);
	nonlocalRepulsion.addGlobalParameter("d_exclusive", D_ExclusiveInNm);
	nonlocalRepulsion.addGlobalParameter("epsilon", E_ExclusionPair);
	for (int i=0; i<N_particles; ++i) {
		nonlocalRepulsion.addParticle();
	}

	// add non-local Go contact for native contact pairs
	OpenMM::CustomBondForce goContactForce = *new OpenMM::CustomBondForce("epsilon*(5*(r_native/r)^12 - 6*(r_native/r)^10)");
	system.addForce(&goContactForce);
	goContactForce.addGlobalParameter("epsilon", E_GoContactPair);
	goContactForce.addPerBondParameter("r_native");
	for (int i=0; i<molinfo["contact"].size(); ++i) {
		int p1 = rs2pi[molinfo["contact"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["contact"][i]["resSeq"][1]];
		double r_native = molinfo["contact"][i]["length"];
		std::vector<double> perBondParams = {r_native * OpenMM::NmPerAngstrom};
		goContactForce.addBond(p1, p2, perBondParams);
		nonlocalRepulsion.addExclusion(p1, p2);
	}

	// add bonding force or constraints between two particles
	OpenMM::HarmonicBondForce& bondStretch = *new OpenMM::HarmonicBondForce();
	system.addForce(&bondStretch);
	for (int i=0; i<molinfo["bond"].size(); ++i) {
		int p1 = rs2pi[molinfo["bond"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["bond"][i]["resSeq"][1]];
		double length = molinfo["bond"][i]["length"];
		if (UseConstraints) {
			system.addConstraint(p1, p2, length * OpenMM::NmPerAngstrom);
		} else {
			bondStretch.addBond(p1, p2, length * OpenMM::NmPerAngstrom, K_BondStretch);
		}
		nonlocalRepulsion.addExclusion(p1, p2);
	}

	// add angle force
	OpenMM::HarmonicAngleForce& bondBend = *new OpenMM::HarmonicAngleForce();
	system.addForce(&bondBend);
	for (int i=0; i<molinfo["angle"].size(); ++i) {
		int p1 = rs2pi[molinfo["angle"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["angle"][i]["resSeq"][1]];
		int p3 = rs2pi[molinfo["angle"][i]["resSeq"][2]];
		double angle = molinfo["angle"][i]["angle"];
		bondBend.addAngle(p1, p2, p3, angle * OpenMM::RadiansPerDegree, K_BondAngle);
		nonlocalRepulsion.addExclusion(p1, p3);
	}

	// add dihedral force
	OpenMM::PeriodicTorsionForce& bondTorsion = *new OpenMM::PeriodicTorsionForce();
	system.addForce(&bondTorsion);
	for (int i=0; i<molinfo["dihedral"].size(); ++i) {
		int p1 = rs2pi[molinfo["dihedral"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["dihedral"][i]["resSeq"][1]];
		int p3 = rs2pi[molinfo["dihedral"][i]["resSeq"][2]];
		int p4 = rs2pi[molinfo["dihedral"][i]["resSeq"][3]];
		double dihedral = molinfo["dihedral"][i]["dihedral"];
		bondTorsion.addTorsion(p1, p2, p3, p4, 
			1, dihedral * OpenMM::RadiansPerDegree, K_BondTorsion1);
		bondTorsion.addTorsion(p1, p2, p3, p4, 
			3, dihedral * OpenMM::RadiansPerDegree, K_BondTorsion3);
		nonlocalRepulsion.addExclusion(p1, p4);
	}

	// set langevin integrator
	OpenMM::LangevinIntegrator integrator(Temperature, LangevinFrictionPerPs, TimePerStepInPs);
	integrator.setRandomNumberSeed(RandomSeed);

	// Let OpenMM Context choose best platform
	OpenMM::Context context(system, integrator, platform);

	std::cout << "# Using OpenMM Platform: ";
	std::cout << context.getPlatform().getName().c_str() << std::endl;

	// Set starting positions of the atoms. Leave time and velocity zero.
	context.setPositions(initPosInNm);

	// set output files
	std::ofstream opdb;
	opdb.open(filename + ".pdb", std::ios::out);
	std::cout << "# output: " << filename + ".pdb\n";

	// Simulate.
	for (int frameNum = 1;; ++frameNum) {
		// Output current state information
		OpenMM::State state = context.getState(OpenMM::State::Positions);

		// write PDB frame
		writePDBFrame(frameNum, state, opdb);

		// show progress in stdout
		int steps = frameNum * NStepSave;
		std::cout << '\r';
		std::cout << std::setw(8) << steps << " steps / ";
		std::cout << std::setw(8) << SimulationSteps << std::flush;
		if (frameNum * NStepSave >= SimulationSteps) break;
		integrator.step(NStepSave);
	}
	std::cout << std::endl;
	opdb.close();

}

void writePDBFrame(int frameNum, const OpenMM::State& state, std::ofstream &opdb) {
	// Reference atomic positions in the OpenMM State.
	const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

	opdb << "MODEL " << frameNum << "\n";
	for (int a = 0; a < (int)posInNm.size(); ++a)
	{
		opdb << "ATOM  " << std::setw(5) << a+1 << "  C    C      1    "; // atom number
		opdb << std::setw(8) << std::setprecision(3) << posInNm[a][0]*10;
		opdb << std::setw(8) << std::setprecision(3) << posInNm[a][1]*10;
		opdb << std::setw(8) << std::setprecision(3) << posInNm[a][2]*10;
		opdb << "  1.00  0.00\n";
	}
	opdb << "ENDMDL\n"; // end of frame
}

