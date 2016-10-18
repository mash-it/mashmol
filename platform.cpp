#include <OpenMM.h>
#include <iostream>
using namespace std;

int main() {
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());
	cout << OpenMM::Platform::getNumPlatforms() << endl;	// 3
	cout << OpenMM::Platform::getPlatform(0).getName() << endl;	// Reference
	cout << OpenMM::Platform::getPlatform(1).getName() << endl;	// Reference
	cout << OpenMM::Platform::getPlatform(2).getName() << endl;	// Reference
	return 0;
}
