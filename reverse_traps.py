'''
This script can be used to create force files with two harmonic traps after force-extension simulations to reverse them.

The easiest way to run it is to copy this source file to the same directory as the raw data and execute
	python reverse_traps.py file1 file2 [tequi] [null]
with
	python	an interpreter for Python 3.x (usually python3 or python if using conda),
	file1	configuration file,
	file2	the old force file,
	[tequi]	the number of time steps for equilibration (defaults to 0);
		pass a negative value if you want to restart the step counter for each run.

The force file must have the following format:

{
type = trap
particle = <int>
pos0 = <float>, <float>, <float>
stiff = <float>
rate = <float>
dir = <float>, <float>, <float>
}

{
type = trap
particle = <int>
pos0 = <float>, <float>, <float>
stiff = <float>
rate = <float>
dir = <float>, <float>, <float>
}

Generated output:
* force file for equilibration
* force file for pulling

This code is cc-by-sa.
-- Sebastian V. Bauer, Walther Lab at the University of Mainz, 19.04.2024
-- sebastian.bauer@uni-mainz.de
'''



import sys
import numpy as np

###############################################
# customized functions to evaluate oxDNA data #
###############################################

# determine center of mass for a single nucleotide from a respective line of an oxDNA trajectory file
def get_oxDNA_cm_site(line):
	line_arr = (np.asarray(line.split())).astype(float)
	return line_arr[0:3]

# determine 1D parameter of type int from a line "key = <int>"
def get_value_int1d(line):
	line_arr = np.asarray(line.split())
	return int(line_arr[2])

# determine 1D parameter of type float from a line "key = <float>"
def get_value_float1d(line):
	line_arr = np.asarray(line.split())
	return float(line_arr[2])

# determine 3D parameter of type float from a line "key = <float>, <float>, <float>"
def get_value_float3d(line):
	temp_arr = np.zeros(3, float)
	line_arr = np.asarray(line.split())
	for i in range(3):
		temp_arr[i] = line_arr[i+2].replace(',', '')
	return temp_arr



###################
# data extraction #
###################

# parse positional arguments
if len(sys.argv) < 3:
	print("execute\n\tpython reverse_traps.py file1 file2 [tequi]\nwith")
	print("""
	file1	configuration file,
	file2	the old force file,
	[tequi]	the number of time steps for equilibration (defaults to 0);
		pass a negative value if you want to restart the step counter for each run.
	""")
	sys.exit(2)
configfile_path = sys.argv[1]
forcefile_path  = sys.argv[2]
try:
	tequi   = int(sys.argv[3])
except:
	tequi   = 0
	print("no number of time steps for equilibration provided, defaulting to 0")
else:
	if tequi < 0:
		print("received a negative number of time steps for equalibration, will assume you want to restart the step counter")
	else:
		print("number of time steps for equilibration is", tequi)
print()


# extract time step from configuration file
configfile = open(configfile_path, 'r')
timestep = get_value_int1d(configfile.readline())
configfile.close()
print(configfile_path, "-> found configuration for time step no.", timestep)
print()



# extract parameters from force file
forcefile = open(forcefile_path, 'r')
for l in range(2):
	line = forcefile.readline()
alpha_index = get_value_int1d(forcefile.readline())
alpha_pos0  = get_value_float3d(forcefile.readline())
alpha_stiff = get_value_float1d(forcefile.readline())
alpha_rate  = get_value_float1d(forcefile.readline())
alpha_dir   = get_value_float3d(forcefile.readline())
for l in range(4):
	line = forcefile.readline()
omega_index = get_value_int1d(forcefile.readline())
omega_pos0  = get_value_float3d(forcefile.readline())
omega_stiff = get_value_float1d(forcefile.readline())
omega_rate  = get_value_float1d(forcefile.readline())
omega_dir   = get_value_float3d(forcefile.readline())
forcefile.close()

print(forcefile_path, "-> found harmonic trap alpha:")
print("\tparticle\t", alpha_index)
print("\tpos0\t",     alpha_pos0)
print("\tstiff\t",    alpha_stiff)
print("\trate\t",     alpha_rate)
print("\tdir\t",      alpha_dir)
print(forcefile_path, "-> found harmonic trap omega:")
print("\tparticle\t", omega_index)
print("\tpos0\t",     omega_pos0)
print("\tstiff\t",    omega_stiff)
print("\trate\t",     omega_rate)
print("\tdir\t",      omega_dir)
print()



###################
# data processing #
###################

# calculate position of traps for equilibration
alpha_pos1 = alpha_pos0 + timestep * alpha_rate * alpha_dir
omega_pos1 = omega_pos0 + timestep * omega_rate * omega_dir
print("current trap positions for equilibration:")
print("\talpha\t", alpha_pos1)
print("\tomega\t", omega_pos1)
print()



# calculate virtual origin of traps for pullback
if tequi < 0:
	alpha_pos2 = alpha_pos1
	omega_pos2 = omega_pos1
else:
	alpha_pos2 = alpha_pos1 + (timestep + tequi) * alpha_rate * alpha_dir
	omega_pos2 = omega_pos1 + (timestep + tequi) * omega_rate * omega_dir
print("virtual trap origins for pullback:")
print("\talpha\t", alpha_pos2)
print("\tomega\t", omega_pos2)
print()



#######################
# force file creation #
#######################

# write force file for equilibration
equifile_path = "forces_equi.conf"
equifile = open(equifile_path, 'w')
equifile.write("{")
equifile.write("\ntype = trap")
equifile.write("\nparticle = " + repr(alpha_index))
equifile.write("\npos0 = "     + repr(alpha_pos1[0]) + ", " + repr(alpha_pos1[1]) + ", " + repr(alpha_pos1[2]))
equifile.write("\nstiff = "    + repr(alpha_stiff))
equifile.write("\nrate = "     + repr(0.0))
equifile.write("\ndir = "      + repr(-alpha_dir[0]) + ", " + repr(-alpha_dir[1]) + ", " + repr(-alpha_dir[2]))
equifile.write("\n}")
equifile.write("\n")
equifile.write("\n{")
equifile.write("\ntype = trap")
equifile.write("\nparticle = " + repr(omega_index))
equifile.write("\npos0 = "     + repr(omega_pos1[0]) + ", " + repr(omega_pos1[1]) + ", " + repr(omega_pos1[2]))
equifile.write("\nstiff = "    + repr(omega_stiff))
equifile.write("\nrate = "     + repr(0.0))
equifile.write("\ndir = "      + repr(-omega_dir[0]) + ", " + repr(-omega_dir[1]) + ", " + repr(-omega_dir[2]))
equifile.write("\n}")
equifile.close()
print(equifile_path, "<- harmonic traps for equilibration written to file")
print()



# write force file for pulling
pullfile_path = "forces_pull.conf"
pullfile = open(pullfile_path, 'w')
pullfile.write("{")
pullfile.write("\ntype = trap")
pullfile.write("\nparticle = " + repr(alpha_index))
pullfile.write("\npos0 = "     + repr(alpha_pos2[0]) + ", " + repr(alpha_pos2[1]) + ", " + repr(alpha_pos2[2]))
pullfile.write("\nstiff = "    + repr(alpha_stiff))
pullfile.write("\nrate = "     + repr(alpha_rate))
pullfile.write("\ndir = "      + repr(-alpha_dir[0]) + ", " + repr(-alpha_dir[1]) + ", " + repr(-alpha_dir[2]))
pullfile.write("\n}")
pullfile.write("\n")
pullfile.write("\n{")
pullfile.write("\ntype = trap")
pullfile.write("\nparticle = " + repr(omega_index))
pullfile.write("\npos0 = "     + repr(omega_pos2[0]) + ", " + repr(omega_pos2[1]) + ", " + repr(omega_pos2[2]))
pullfile.write("\nstiff = "    + repr(omega_stiff))
pullfile.write("\nrate = "     + repr(omega_rate))
pullfile.write("\ndir = "      + repr(-omega_dir[0]) + ", " + repr(-omega_dir[1]) + ", " + repr(-omega_dir[2]))
pullfile.write("\n}")
pullfile.close()
print(pullfile_path, "<- harmonic traps for reverse pulling written to file")
print()



