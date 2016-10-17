#include <OpenMM.h>

void simulate();

int main() {
	simulate();
	return 0;
}

void simulate() {
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

}
