#include <OpenMM.h>
#include <iostream>
using namespace std;

int main() {
	cout << OpenMM::Platform::getNumPlatforms() << endl;	// 1
	cout << OpenMM::Platform::getPlatform(0).getName() << endl;	// Reference
	return 0;
}
