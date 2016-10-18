#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

void simulate();

int main() {
	simulate();
	return 0;
}

void simulate() {
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());
	
	// Create a system with nonbonded forces.
	OpenMM::System system;
	OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
	system.addForce(nonbond);
	
	// Create Atoms
	const int Natoms = 3;
	double x, y, z;
	std::vector<OpenMM::Vec3> initPosInNm(Natoms);

	// load atom information from input file.
	std::ifstream ifs("input.dat");
	std::string str;
	if(ifs.fail()) {
		std::cerr << "File do not exist.\n";
		exit(1);
	}
	for (int i=0; i<Natoms; ++i) {
		getline(ifs, str); // TODO: Exception
		sscanf(str.data(), "%lf\t%lf\t%lf", &x, &y, &z);
		initPosInNm[i] = OpenMM::Vec3(x, y, z);

		system.addParticle(39.95);	// mass of Ar, gram per mole
		nonbond->addParticle(0.0, 0.3350, 0.996);	// charge, LJ sigma (nm), well depth (kJ) 
	}

	OpenMM::VerletIntegrator integrator(0.004);	// step size in ps

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
		// std::cout << timeInPs << std::endl;

		if (timeInPs >= 100.) break;
		integrator.step(100);
	}
}
