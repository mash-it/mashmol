CC = g++
INLUCDES = -I/usr/loca/openmm/include
LIBS = -L/usr/local/openmm/lib
CO = -std=c++11

Polymer: Polymer.cpp
	$(CC) Polymer.cpp -o Polymer $(INCLUDES) $(LIBS) -lOpenMM $(CO)

