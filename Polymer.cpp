#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <json.hpp>
#include <map>

void simulate();
void writePdbFrame(int, const OpenMM::State&);

using json = nlohmann::json;

// global parameters
const double K_BondStretch = 40000.0; // amu ps^-2
const double K_BondAngle= 80.0; // amu nm^2 ps^-2
const double K_BondTorsion1 = 4.0; // amu nm^2 ps^-2
const double K_BondTorsion3 = 2.0; // amu nm^2 ps^-2
const bool UseConstraints = false;

const double Temperature = 300.0; // Kelvin
const double LangevinFrictionPerPs = 5.0;
const double SimulationTimeStepInPs = 0.02;
const double SimulationTimeInPs = 100.0;
const int N_StepSave = 10;

int main() {
	simulate();
	return 0;
}

void simulate() {

	// Load any shared libraries containing GPU implementations.
	OpenMM::Platform::loadPluginsFromDirectory(
      	OpenMM::Platform::getDefaultPluginsDirectory());
	
	// Create a system with forces.
	OpenMM::System system;

	OpenMM::HarmonicBondForce& bondStretch = *new OpenMM::HarmonicBondForce();
	OpenMM::HarmonicAngleForce& bondBend = *new OpenMM::HarmonicAngleForce();
	OpenMM::PeriodicTorsionForce& bondTorsion = *new OpenMM::PeriodicTorsionForce();
	system.addForce(&bondStretch);
	system.addForce(&bondBend);
	system.addForce(&bondTorsion);
	// OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
	// system.addForce(nonbond);


	// load atom information from input file.
	json molinfo;
	std::ifstream ifs("input.json");
	if(ifs.fail()) {
		std::cerr << "File do not exist.\n";
		exit(1);
	}
	ifs >> molinfo;

	const int N_particles = molinfo["resSeq"].size();
	std::cout << "# Number of Particles : " << N_particles << '\n';

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

	// add constraints between two particles
	for (int i=0; i<molinfo["bond"].size(); ++i) {
		int p1 = rs2pi[molinfo["bond"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["bond"][i]["resSeq"][1]];
		double length = molinfo["bond"][i]["length"];
		if (UseConstraints) {
			system.addConstraint(p1, p2, length * OpenMM::NmPerAngstrom);
		} else {
			bondStretch.addBond(p1, p2, length * OpenMM::NmPerAngstrom, K_BondStretch);
		}
	}

	// add angle force
	for (int i=0; i<molinfo["angle"].size(); ++i) {
		int p1 = rs2pi[molinfo["angle"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["angle"][i]["resSeq"][1]];
		int p3 = rs2pi[molinfo["angle"][i]["resSeq"][2]];
		double angle = molinfo["angle"][i]["angle"];
		bondBend.addAngle(p1, p2, p3, angle * OpenMM::RadiansPerDegree, K_BondAngle);
	}

	// add dihedral force
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
	}

	//OpenMM::VerletIntegrator integrator(0.004);	// step size in ps
	OpenMM::LangevinIntegrator integrator(Temperature, LangevinFrictionPerPs, SimulationTimeStepInPs);

	// force to use CPU platform
	OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("CPU");

	//Let OpenMM Context choose best platform
	OpenMM::Context context(system, integrator, platform);
	std::cout << "REMARK Using OpenMM platform ";
	std::cout << context.getPlatform().getName().c_str() << std::endl;
	std::cout << std::endl;

	// Set starting positions of the atoms. Leave time and velocity zero.
	context.setPositions(initPosInNm);

	// Simulate.
	for (int frameNum = 1;; ++frameNum) {
		// Output current state information
		OpenMM::State state = context.getState(OpenMM::State::Positions);
		const double timeInPs = state.getTime();
		writePdbFrame(frameNum, state);

		if (timeInPs >= SimulationTimeInPs) break;
		integrator.step(N_StepSave);
	}

}

// Handy homebrew PDB writer for quick-and-dirty trajectory output.
void writePdbFrame(int frameNum, const OpenMM::State& state)
{
	// Reference atomic positions in the OpenMM State.
	const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

	// Use PDB MODEL cards to number trajectory frames
	printf("MODEL	 %d\n", frameNum); // start of frame
	for (int a = 0; a < (int)posInNm.size(); ++a)
	{
		printf("ATOM  %5d  C    C      1    ", a+1); // atom number
		printf("%8.3f%8.3f%8.3f  1.00  0.00\n",	  // coordinates
			// "*10" converts nanometers to Angstroms
			posInNm[a][0]*10, posInNm[a][1]*10, posInNm[a][2]*10);
	}
	printf("ENDMDL\n"); // end of frame
}

