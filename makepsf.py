import sys

N_atoms = int(sys.argv[1])
N_bonds = N_atoms-1
template = "    {:4d} U    {:4d}  MET  CA   C     0.000000       12.0110           0"

print("  {:4d} !NATOM".format(N_atoms))
for i in range(1,N_atoms+1):
	print(template.format(i,i))

print("")
print(" {:4d} !NBOND".format(N_bonds))

count = 0
for i in range(1, N_bonds+1):
	print("{:8d}{:8d}".format(i, i+1), end='')
	count += 1
	if count % 4 == 0:
		print("")

print("")
