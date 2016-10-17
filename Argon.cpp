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
		std::cout << x << ',' << y << ',' << z  << std::endl;
	}
}
