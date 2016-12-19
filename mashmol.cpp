#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <json.hpp>
#include <map>
#include <cmath>
#include <chrono>
#include <unistd.h>

using json = nlohmann::json;

void simulate(json&);
void writePDBFrame(int, const OpenMM::State&, json&, std::ofstream&);
void writeTimeSeries(int, const OpenMM::State&, std::ofstream&, double);
double getQscore(const OpenMM::State&, json&);


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
const double QscoreThreshold = 1.44;

// global variables
std::map<int,int> rs2pi; // residue sequence to particle index
bool SilentMode = false; // do not show progress in stdout

int main(int argc, char* argv[]) {

	if (argc == 1) {
		std::cerr << "Usage:\t" << argv[0] << " input.json\n";
		exit(1);
	}

	// load global paremters and atom information from input file.
	json molinfo;
	std::ifstream ifs(argv[1]);
	if(ifs.fail()) {
		std::cerr << argv[1] << " :cannot open input file.\n";
		exit(1);
	}
	ifs >> molinfo;

	// read other options
	int result;
	while ((result = getopt(argc, argv, "s"))!=-1) {
		switch(result) {
			case 's':
				SilentMode = true;
				break;
			case ':':
				break;
			case '?':
				break;
		}
	}

	const auto startTime = std::chrono::system_clock::now();

	// run simulation
	simulate(molinfo);

	const auto endTime = std::chrono::system_clock::now();
	const auto timeSpan = endTime - startTime;
	const int timeSpanInmsec = std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count();
	const double timeSpanInSec = timeSpanInmsec / 1000.0;

	std::cout << "Time: " << timeSpanInSec << '\n';
	std::cout << "Speed: " << (double)molinfo["parameters"]["SimulationSteps"] / timeSpanInSec << "steps/sec" << '\n';

	return 0;
}

void simulate(json &molinfo) {

	// Load any shared libraries containing GPU implementations.
	OpenMM::Platform::loadPluginsFromDirectory(
      	OpenMM::Platform::getDefaultPluginsDirectory());

	// force to use CPU platform that is better than CUDA in my simulation scale
	OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("CPU");
	// force to use single thread 
	platform.setPropertyDefaultValue("Threads", "1");
	
	// Create a system with forces.
	OpenMM::System system;

	const int N_particles = molinfo["resSeq"].size();
	std::cout << "# Number of Particles : " << N_particles << '\n';

	// set MD paramters
	const double Temperature = molinfo["parameters"]["Temperature"];
	const int SimulationSteps = molinfo["parameters"]["SimulationSteps"];
	const double TimePerStepInPs = molinfo["parameters"]["TimePerStepInPs"];
	const int NStepSave = molinfo["parameters"]["NStepSave"];
	const int RandomSeed = molinfo["parameters"]["RandomSeed"];

	// set other parameters
	const std::string filename = molinfo["output"]["filename"];

	// Create Atoms
	std::vector<OpenMM::Vec3> initPosInNm(N_particles);

	// map from resSeq to particle index(0..N-1)
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
	OpenMM::CustomNonbondedForce& nonlocalRepulsion = *new OpenMM::CustomNonbondedForce("epsilon_repul*(d_exclusive/r)^12");
	system.addForce(&nonlocalRepulsion);
	nonlocalRepulsion.setNonbondedMethod(OpenMM::CustomNonbondedForce::NonbondedMethod::CutoffNonPeriodic);
	nonlocalRepulsion.setCutoffDistance(D_ExclusionCutoffInNm);
	nonlocalRepulsion.addGlobalParameter("d_exclusive", D_ExclusiveInNm);
	nonlocalRepulsion.addGlobalParameter("epsilon_repul", E_ExclusionPair);
	for (int i=0; i<N_particles; ++i) {
		nonlocalRepulsion.addParticle();
	}

	// add non-local Go contact for native contact pairs
	OpenMM::CustomBondForce& goContactForce = *new OpenMM::CustomBondForce("epsilon_ngo*(5*(r_native/r)^12 - 6*(r_native/r)^10)");
	system.addForce(&goContactForce);
	goContactForce.addGlobalParameter("epsilon_ngo", E_GoContactPair);
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
			1, (dihedral+180.0) * OpenMM::RadiansPerDegree, K_BondTorsion1);
		bondTorsion.addTorsion(p1, p2, p3, p4, 
			3, (dihedral+180.0) * OpenMM::RadiansPerDegree, K_BondTorsion3);
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
	// set initial velocity
	context.setVelocitiesToTemperature(Temperature, RandomSeed);

	// set output files
	std::ofstream opdb, ots;
	opdb.open(filename + ".pdb", std::ios::out);
	ots.open(filename + ".ts", std::ios::out);
	std::cout << "# PDB: " << filename + ".pdb\n";
	std::cout << "# TimeSereis: " << filename + ".ts\n";

	int infoMask = 0;
	infoMask = OpenMM::State::Positions;
	infoMask += OpenMM::State::Energy;

	for (int frameNum = 0;; ++frameNum) {

		// Output current state information
		OpenMM::State state = context.getState(infoMask);
		int steps = frameNum * NStepSave;
		double qscore = getQscore(state, molinfo);

		// write PDB frame
		writePDBFrame(frameNum, state, molinfo, opdb);
		writeTimeSeries(steps, state, ots, qscore);

		// show progress in stdout
		if (SilentMode == false) {
			std::cout << '\r';
			std::cout << "Q=" << std::setw(5) << std::setprecision(3) << qscore;
			std::cout << "; T=" << std::setw(5) << integrator.getTemperature();
			std::cout << "; E=" << std::setw(5) << state.getPotentialEnergy();
			std::cout << std::setw(8) << steps << " steps / ";
			std::cout << std::setw(8) << SimulationSteps << std::flush;
		}

		if (frameNum * NStepSave >= SimulationSteps) break;

		integrator.step(NStepSave);
	}
	std::cout << std::endl;
	opdb.close();
	ots.close();

}

void writeTimeSeries(int steps, const OpenMM::State& state, std::ofstream &ots, double qscore) {
	ots << std::setw(8) << steps;
	ots << ' ' << std::fixed << std::setw(8) << std::setprecision(3) << qscore;
	ots << ' ' << std::fixed << std::setw(8) << state.getPotentialEnergy();
	ots << ' ' << std::fixed << std::setw(8) << state.getKineticEnergy();
	ots << '\n';
}

void writePDBFrame(int frameNum, const OpenMM::State& state, json &molinfo, std::ofstream &opdb) {
	// Reference atomic positions in the OpenMM State.
	const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

	opdb << "MODEL " << frameNum << "\n";
	for (int a = 0; a < (int)posInNm.size(); ++a)
	{
		opdb << "ATOM  " << std::setw(5) << a+1 << "  C    C  A";
		opdb << std::setw(4) << (int)molinfo["resSeq"][a];
		opdb << "    "; // atom number
		opdb << std::fixed << std::setw(8) << std::setprecision(3) << posInNm[a][0]*10;
		opdb << std::fixed << std::setw(8) << std::setprecision(3) << posInNm[a][1]*10;
		opdb << std::fixed << std::setw(8) << std::setprecision(3) << posInNm[a][2]*10;
		opdb << "  1.00  0.00\n";
	}
	opdb << "ENDMDL\n"; // end of frame
}

double getQscore(const OpenMM::State& state, json &molinfo) {
	const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

	int contact = 0;
	for (int i=0; i<molinfo["contact"].size(); ++i) {
		int p1 = rs2pi[molinfo["contact"][i]["resSeq"][0]];
		int p2 = rs2pi[molinfo["contact"][i]["resSeq"][1]];
		double r_native = molinfo["contact"][i]["length"];
		double r, r2;
		r2  = pow(posInNm[p1][0]-posInNm[p2][0], 2);
		r2 += pow(posInNm[p1][1]-posInNm[p2][1], 2);
		r2 += pow(posInNm[p1][2]-posInNm[p2][2], 2);
		r = sqrt(r2) / OpenMM::NmPerAngstrom;
		if (r < r_native * QscoreThreshold ) contact++;
	}
	return ((double)contact) / molinfo["contact"].size();
}

