#!/usr/bin/python

import sys

if ( len(sys.argv) < 2 ) :
	print("Usage:")
	print(sys.argv[0]+ " .../mesa/chem/data/isotopes.data > .../flash/Simulation/SimulationComposition/Burn/SpeciesList.txt")
	sys.exit(1)



mesa_iso_data_filename = sys.argv[1]

mfile = open(mesa_iso_data_filename)

for line in mfile :
	# skip empty lines
	if (len(line.strip())<1) : continue
	# lines starting with a space are partition function values
	if (line[0]==' ') : continue
	sl = line.split()
	# now get the data
	nuclidename = sl[0].strip()
	# skip isomeric states, names would be confusing for flash
	if ( '-' in nuclidename ) : continue
	# flash only does 4-character names
	if ( len(nuclidename) > 4) : continue
	Z = int(sl[2])
	N = int(sl[3])
	spin = float(sl[4])
	mass_excess = float(sl[5])

	binding_energy = Z*(938.2723+0.5110)+N*939.5656 -((Z+N)*931.494028 + mass_excess)

	spin_factor = 2*spin+1

	# this factor is in the flash file, but I don't know what it is
	if (N+Z==1) :
		unknown = 0.5
	else :
		unknown = 0

	print(nuclidename +'\t'+ ( "%i\t%i\t%i\t%7.3f\t%i\t%3.1f" % (Z, (Z+N), N, binding_energy, spin_factor, unknown ) ))
