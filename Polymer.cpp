#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <json.hpp>

void simulate();
void writePdbFrame(int, const OpenMM::State&);

using json = nlohmann::json;

// parameters
const bool UseConstraints = false;

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

	OpenMM::HarmonicBondForce&  bondStretch = *new OpenMM::HarmonicBondForce();
	OpenMM::HarmonicAngleForce& bondBend    = *new OpenMM::HarmonicAngleForce();
	system.addForce(&bondStretch);
	system.addForce(&bondBend);
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
	const int N_atoms = molinfo["resSeq"].size();

	// Create Atoms
	std::vector<OpenMM::Vec3> initPosInNm(N_atoms);

	// map from 

	double x, y, z;
	for (int i=0; i<N_atoms; ++i) {
		x = molinfo["position"][i][0];
		y = molinfo["position"][i][1];
		z = molinfo["position"][i][2];
		initPosInNm[i] = OpenMM::Vec3(x, y, z) * OpenMM::NmPerAngstrom;
		system.addParticle(137.0);	// average mass of amino acid residues
	}

	// add constraints between two particles
	for (int i=0; i<N_atoms-1; ++i) {
		if (UseConstraints) {
			system.addConstraint(i, i+1, bondLength * OpenMM::NmPerAngstrom);
		} else {
			bondStretch.addBond(i, i+1, bondLength * OpenMM::NmPerAngstrom, 100);
		}
	}

	// add angle
	for (int i=0; i<N_atoms-2; ++i) {
		bondBend.addAngle(i, i+1, i+2, 120.0 * OpenMM::RadiansPerDegree, 50);
	}

	//OpenMM::VerletIntegrator integrator(0.004);	// step size in ps
	OpenMM::LangevinIntegrator integrator(3, 0.04, 0.004);	// step size in ps

	// force to use CPU platform
	OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("CPU");

	//Let OpenMM Context choose best platform
	OpenMM::Context context(system, integrator, platform);
	std::cout << "REMARK Using OpenMM platform ";
	std::cout << context.getPlatform().getName().c_str() << std::endl;

	// Set starting positions of the atoms. Leave time and velocity zero.
	context.setPositions(initPosInNm);

	// Simulate.
	for (int frameNum = 1;; ++frameNum) {
		// Output current state information
		OpenMM::State state = context.getState(OpenMM::State::Positions);
		const double timeInPs = state.getTime();
		writePdbFrame(frameNum, state);

		if (timeInPs >= 100.) break;
		integrator.step(100);
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

