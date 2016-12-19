CC = g++
INLUCDES = -I/usr/loca/openmm/include
LIBS = -L/usr/local/openmm/lib
CO = -std=c++11

mashmol: mashmol.cpp
	$(CC) mashmol.cpp -o mashmol $(INCLUDES) $(LIBS) -lOpenMM $(CO)

