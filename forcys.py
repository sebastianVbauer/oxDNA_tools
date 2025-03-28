'''
This script can be used to analyze force-extension hysteresis data from oxDNA simulations with two harmonic traps as external forces, each acting on a specific nucleotide.
Force is calculated from nucleotide-trap displacements as projection of the total force for two springs connected in series onto the axis between the springs.

The easiest way to run the analysis is to copy the script to the same directory as the raw data and execute
	python forcys.py file1 file2 file3 file4 file5 [file6] [file7]
with
	python	an interpreter for Python 3.x (usually python3 or python if using conda),
	file1	topology file (only classic oxDNA format, i. e. 3'->5'),
	file2	forward trajectory file,
	file3	force file for the forward trajectory,
	file4	backward trajectory file,
	file5	force file for the backward trajectory,
	[file6]	equilibration trajectory file (optional),
	[file7]	force file for the equilibration trajectory (optional).

The force files must have the following format:

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
* local maxima/minima in force for the forward/backward trajectories (written to standard output)
* force-extension integrals for each trajectory (written to standard output)
* area enclosed by the force-extension hysteresis curve (written to standard output)
* trajectories of nucleotide-nucleotide and both nucleotide-trap distances (written to file forcys_distances.png)
* force-extension hysteresis curve (written to file forcys_hysteresis.png)

This code is cc-by-sa.
-- Sebastian V. Bauer, Walther Lab at the University of Mainz, 04.07.2024
-- sebastian.bauer@uni-mainz.de
'''



import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt



###############################################
# customized functions to evaluate oxDNA data #
###############################################

# determine center of mass for a single nucleotide from a respective line of an oxDNA trajectory file
def get_oxDNA_cm_site(line):
	line_arr = (np.asarray(line.split())).astype(float)
	return line_arr[0:3]

# determine time of a configuration from a respective line of an oxDNA trajectory file
def get_oxDNA_step(line):
	line_arr = np.asarray(line.split())
	return int(line_arr[2])

# extract 3d value (pos0 or dir) from a respective line of an oxDNA force file
def get_oxDNA_ffval(line):
	val3d = line.split()[2:5]
	for i in range(3):
		val3d[i] = val3d[i].replace(',', ' ')
	return (np.asarray(val3d)).astype(float)

# convert length from oxDNA simulation units to nanometer
def su2nm(x):
	return 0.8518 * x

# convert time from oxDNA simulation units to picoseconds
def su2ps(t):
	return 3.031 * t

# convert force from oxDNA simulation units to picoNewton
def su2pN(F):
	return 48.63 * F

# convert force constant from oxDNA simulation units to picoNewton per nanometer
def su2pN_nm(k):
	return 57.09 * k

# convert energy from oxDNA simulation units to Joules
def su2J(E):
	return 4.142e-20 * E

# exponential moving average with pre-smoothing
def ema(data, alpha, pre):
	m = min(pre, len(data))
	s = np.zeros_like(data)
	s[0] = np.mean(data[:m])
	for k in range(1, len(data)):
		s[k] = s[k-1] + alpha*(data[k] - s[k-1])
	return s



###################
# data extraction #
###################

# parse positional arguments
if len(sys.argv) < 6:
	print("execute\n\tpython forcys.py file1 file2 file3 file4 file5 [file6] [file7]\nwith")
	print("""
	file1	topology file (only classic oxDNA format, i. e. 3'->5'),
	file2	forward trajectory file,
	file3	force file for the forward trajectory,
	file4	backward trajectory file,
	file5	force file for the backward trajectory,
	[file6]	equilibration trajectory file (optional),
	[file7]	force file for the equilibration trajectory (optional).
	""")
	sys.exit(2)
topfile_path      = sys.argv[1]
trajfile_path_fw  = sys.argv[2]
forcefile_path_fw = sys.argv[3]
trajfile_path_bw  = sys.argv[4]
forcefile_path_bw = sys.argv[5]
try:
	trajfile_path_eq = sys.argv[6]
except:
	eq = False
else:
	eq = True
try:
	forcefile_path_eq = sys.argv[7]
except:
	eq = False
else:
	eq = True



# evaluate topology file
topfile = open(topfile_path, 'r')
line = topfile.readline()
line_lst = line.split()
number_of_nucleotides = int(line_lst[0])
topfile.close()
print(topfile_path, "-> number of nucleotides:", number_of_nucleotides)
print()



# browse trajectory files to get numbers of configurations
## forward trajectory
num_configs_fw = 0
trajfile = open(trajfile_path_fw, 'r')
line = trajfile.readline()
while line:
	for n in range (number_of_nucleotides + 2):
		line = trajfile.readline()
	num_configs_fw += 1
	line = trajfile.readline()
trajfile.close()
print(trajfile_path_fw, "-> found", num_configs_fw, "forward configurations")

## equilibration trajectory
num_configs_eq = 0
if eq:
	trajfile = open(trajfile_path_eq, 'r')
	line = trajfile.readline()
	while line:
		for n in range (number_of_nucleotides + 2):
			line = trajfile.readline()
		num_configs_eq += 1
		line = trajfile.readline()
	trajfile.close()
	print(trajfile_path_eq, "-> found", num_configs_eq, "equilibration configurations")
else:
	print("no equilibration period defined")

## backward trajectory
num_configs_bw = 0
trajfile = open(trajfile_path_bw, 'r')
line = trajfile.readline()
while line:
	for n in range (number_of_nucleotides + 2):
		line = trajfile.readline()
	num_configs_bw += 1
	line = trajfile.readline()
trajfile.close()
print(trajfile_path_bw, "-> found", num_configs_bw, "backward configurations")
print()



# evaluate external force files
## forward force file
forcefile = open(forcefile_path_fw, 'r')
alpha_check = False
omega_check = False
line = forcefile.readline()
while line:
	if line[0] == '{':
		if alpha_check == False:
			for l in range(2):
				line = forcefile.readline()
			alpha_fw = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_pos0_fw = get_oxDNA_ffval(line)
			line = forcefile.readline()
			alpha_stiff_fw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_rate_fw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_dir_fw = get_oxDNA_ffval(line)
			if (np.linalg.norm(alpha_dir_fw) > 0.0):
				alpha_dir_fw /= np.linalg.norm(alpha_dir_fw)
			alpha_check = True
		else:
			for l in range(2):
				line = forcefile.readline()
			omega_fw = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_pos0_fw = get_oxDNA_ffval(line)
			line = forcefile.readline()
			omega_stiff_fw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_rate_fw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_dir_fw = get_oxDNA_ffval(line)
			if (np.linalg.norm(omega_dir_fw) > 0.0):
				omega_dir_fw /= np.linalg.norm(omega_dir_fw)
			omega_check = True
	line = forcefile.readline()
forcefile.close()
if alpha_check == False:
	print(forcefile_path_fw, "-> unable to determine harmonic trap alpha in forward force file, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path_fw, "-> found harmonic trap alpha in forward force file for nucleotide", alpha_fw)
	print("\tinitial position [s.u.]:", alpha_pos0_fw)
	print("\tstiffness [s.u.]:", alpha_stiff_fw)
	print("\tmovement rate [s.u.]:", alpha_rate_fw)
	print("\tmovement direction [s.u.]:", alpha_dir_fw)
if omega_check == False:
	print(forcefile_path_fw, "-> unable to determine harmonic trap omega in forward force file, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path_fw, "-> found harmonic trap omega in forward force file for nucleotide", omega_fw)
	print("\tinitial position [s.u.]:", omega_pos0_fw)
	print("\tstiffness [s.u.]:", omega_stiff_fw)
	print("\tmovement rate [s.u.]:", omega_rate_fw)
	print("\tmovement direction [s.u.]:", omega_dir_fw)

## equilibration force file
if eq:
	forcefile = open(forcefile_path_eq, 'r')
	alpha_check = False
	omega_check = False
	line = forcefile.readline()
	while line:
		if line[0] == '{':
			if alpha_check == False:
				for l in range(2):
					line = forcefile.readline()
				alpha_eq = int(np.asarray(line.split())[2])
				line = forcefile.readline()
				alpha_pos0_eq = get_oxDNA_ffval(line)
				line = forcefile.readline()
				alpha_stiff_eq = float(np.asarray(line.split())[2])
				line = forcefile.readline()
				alpha_rate_eq = float(np.asarray(line.split())[2])
				line = forcefile.readline()
				alpha_dir_eq = get_oxDNA_ffval(line)
				if (np.linalg.norm(alpha_dir_eq) > 0.0):
					alpha_dir_eq /= np.linalg.norm(alpha_dir_eq)
				alpha_check = True
			else:
				for l in range(2):
					line = forcefile.readline()
				omega_eq = int(np.asarray(line.split())[2])
				line = forcefile.readline()
				omega_pos0_eq = get_oxDNA_ffval(line)
				line = forcefile.readline()
				omega_stiff_eq = float(np.asarray(line.split())[2])
				line = forcefile.readline()
				omega_rate_eq = float(np.asarray(line.split())[2])
				line = forcefile.readline()
				omega_dir_eq = get_oxDNA_ffval(line)
				if (np.linalg.norm(omega_dir_eq) > 0.0):
					omega_dir_eq /= np.linalg.norm(omega_dir_eq)
				omega_check = True
		line = forcefile.readline()
	forcefile.close()
	if alpha_check == False:
		print(forcefile_path_eq, "-> unable to determine harmonic trap alpha in equilibration force file, aborting.\n")
		sys.exit(1)
	else:
		print(forcefile_path_eq, "-> found harmonic trap alpha in equilibration force file for nucleotide", alpha_eq)
		print("\tinitial position [s.u.]:", alpha_pos0_eq)
		print("\tstiffness [s.u.]:", alpha_stiff_eq)
		print("\tmovement rate [s.u.]:", alpha_rate_eq)
		print("\tmovement direction [s.u.]:", alpha_dir_eq)
	if omega_check == False:
		print(forcefile_path_eq, "-> unable to determine harmonic trap omega in equilibration force file, aborting.\n")
		sys.exit(1)
	else:
		print(forcefile_path_eq, "-> found harmonic trap omega in equilibration force file for nucleotide", omega_eq)
		print("\tinitial position [s.u.]:", omega_pos0_eq)
		print("\tstiffness [s.u.]:", omega_stiff_eq)
		print("\tmovement rate [s.u.]:", omega_rate_eq)
		print("\tmovement direction [s.u.]:", omega_dir_eq)
	print()

## backward force file
forcefile = open(forcefile_path_bw, 'r')
alpha_check = False
omega_check = False
line = forcefile.readline()
while line:
	if line[0] == '{':
		if alpha_check == False:
			for l in range(2):
				line = forcefile.readline()
			alpha_bw = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_pos0_bw = get_oxDNA_ffval(line)
			line = forcefile.readline()
			alpha_stiff_bw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_rate_bw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_dir_bw = get_oxDNA_ffval(line)
			if (np.linalg.norm(alpha_dir_bw) > 0.0):
				alpha_dir_bw /= np.linalg.norm(alpha_dir_bw)
			alpha_check = True
		else:
			for l in range(2):
				line = forcefile.readline()
			omega_bw = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_pos0_bw = get_oxDNA_ffval(line)
			line = forcefile.readline()
			omega_stiff_bw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_rate_bw = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_dir_bw = get_oxDNA_ffval(line)
			if (np.linalg.norm(omega_dir_bw) > 0.0):
				omega_dir_bw /= np.linalg.norm(omega_dir_bw)
			omega_check = True
	line = forcefile.readline()
forcefile.close()
if alpha_check == False:
	print(forcefile_path_bw, "-> unable to determine harmonic trap alpha in backward force file, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path_bw, "-> found harmonic trap alpha in backward force file for nucleotide", alpha_bw)
	print("\tinitial position [s.u.]:", alpha_pos0_bw)
	print("\tstiffness [s.u.]:", alpha_stiff_bw)
	print("\tmovement rate [s.u.]:", alpha_rate_bw)
	print("\tmovement direction [s.u.]:", alpha_dir_bw)
if omega_check == False:
	print(forcefile_path_bw, "-> unable to determine harmonic trap omega in backward force file, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path_bw, "-> found harmonic trap omega in backward force file for nucleotide", omega_bw)
	print("\tinitial position [s.u.]:", omega_pos0_bw)
	print("\tstiffness [s.u.]:", omega_stiff_bw)
	print("\tmovement rate [s.u.]:", omega_rate_bw)
	print("\tmovement direction [s.u.]:", omega_dir_bw)
print()

## some checks
if ( eq and not ( (alpha_fw, omega_fw) == (alpha_eq, omega_eq) == (alpha_bw, omega_bw) ) ) or (alpha_fw, omega_fw) != (alpha_bw, omega_bw):
	print("WARNING: nucleotides assigned to traps differ between the forward and backward runs.")
if ( eq and (alpha_rate_eq > 0.0 or omega_rate_eq > 0.0) and np.linalg.norm(alpha_dir_eq) > 0.0 ):
	print("WARNING: harmonic traps are moving during equilibration")



# evaluate configurations in trajectory files
## forward trajectory
r_alpha_fw = np.zeros((num_configs_fw, 3), float)
r_omega_fw = np.zeros((num_configs_fw, 3), float)
steparr_fw = np.zeros(num_configs_fw, int)
trajfile = open(trajfile_path_fw, 'r')
for k in range(num_configs_fw):
	## get time step at which the current configuration was recorded
	line = trajfile.readline()
	steparr_fw[k] = get_oxDNA_step(line)
	
	## skip box size and energy
	for l in range(2):
		line = trajfile.readline()
	
	## scroll down nucleotides
	for n in range(alpha_fw):
		line = trajfile.readline()
	
	## get position of nucleotide alpha
	line = trajfile.readline()
	r_alpha_fw[k] = get_oxDNA_cm_site(line)
	
	## scroll down nucleotides
	for n in range(omega_fw - alpha_fw - 1):
		line = trajfile.readline()
	
	## get position of nucleotide omega
	line = trajfile.readline()
	r_omega_fw[k] = get_oxDNA_cm_site(line)
	
	## scroll down to next configuration
	for n in range(number_of_nucleotides - omega_fw - 1):
		line = trajfile.readline()
trajfile.close()

## equilibration trajectory
if eq:
	r_alpha_eq = np.zeros((num_configs_eq, 3), float)
	r_omega_eq = np.zeros((num_configs_eq, 3), float)
	steparr_eq = np.zeros(num_configs_eq, int)
	trajfile = open(trajfile_path_eq, 'r')
	for k in range(num_configs_eq):
		## get time step at which the current configuration was recorded
		line = trajfile.readline()
		steparr_eq[k] = get_oxDNA_step(line)
		
		## skip box size and energy
		for l in range(2):
			line = trajfile.readline()
		
		## scroll down nucleotides
		for n in range(alpha_eq):
			line = trajfile.readline()
		
		## get position of nucleotide alpha
		line = trajfile.readline()
		r_alpha_eq[k] = get_oxDNA_cm_site(line)
		
		## scroll down nucleotides
		for n in range(omega_eq - alpha_eq - 1):
			line = trajfile.readline()
		
		## get position of nucleotide omega
		line = trajfile.readline()
		r_omega_eq[k] = get_oxDNA_cm_site(line)
		
		## scroll down to next configuration
		for n in range(number_of_nucleotides - omega_eq - 1):
			line = trajfile.readline()
	trajfile.close()

## backward trajectory
r_alpha_bw = np.zeros((num_configs_bw, 3), float)
r_omega_bw = np.zeros((num_configs_bw, 3), float)
steparr_bw = np.zeros(num_configs_bw, int)
trajfile = open(trajfile_path_bw, 'r')
for k in range(num_configs_bw):
	## get time step at which the current configuration was recorded
	line = trajfile.readline()
	steparr_bw[k] = get_oxDNA_step(line)
	
	## skip box size and energy
	for l in range(2):
		line = trajfile.readline()
	
	## scroll down nucleotides
	for n in range(alpha_bw):
		line = trajfile.readline()
	
	## get position of nucleotide alpha
	line = trajfile.readline()
	r_alpha_bw[k] = get_oxDNA_cm_site(line)
	
	## scroll down nucleotides
	for n in range(omega_bw - alpha_bw - 1):
		line = trajfile.readline()
	
	## get position of nucleotide omega
	line = trajfile.readline()
	r_omega_bw[k] = get_oxDNA_cm_site(line)
	
	## scroll down to next configuration
	for n in range(number_of_nucleotides - omega_bw - 1):
		line = trajfile.readline()
trajfile.close()



#################
# data analysis #
#################

# forward trajectory
## calculate distances between nucleotides alpha and omega
dr_fw = r_omega_fw - r_alpha_fw
distances_fw = np.linalg.norm(dr_fw, axis=1)

## calculate trajectories of each trap and of the trap-trap displacement
trap_alpha_fw = np.zeros((num_configs_fw, 3), float)
trap_omega_fw = np.zeros((num_configs_fw, 3), float)
for k in range(num_configs_fw):
	trap_alpha_fw[k] = alpha_pos0_fw + steparr_fw[k] * alpha_rate_fw * alpha_dir_fw
	trap_omega_fw[k] = omega_pos0_fw + steparr_fw[k] * omega_rate_fw * omega_dir_fw
traps_disp_fw = trap_omega_fw - trap_alpha_fw
traps_dist_fw = np.linalg.norm(traps_disp_fw, axis=1)
traps_axis_fw = traps_disp_fw / traps_dist_fw[:, None]

## calculate nucleotide-trap displacements
disp_alpha_fw = r_alpha_fw - trap_alpha_fw
disp_omega_fw = r_omega_fw - trap_omega_fw
dist_alpha_fw = np.linalg.norm(disp_alpha_fw, axis=1)
dist_omega_fw = np.linalg.norm(disp_omega_fw, axis=1)

## calculate the force for two springs connected in series along the trap-trap axis
keff_fw       = alpha_stiff_fw * omega_stiff_fw / (alpha_stiff_fw + omega_stiff_fw)
projection_fw = np.einsum('ij,ij->i', disp_alpha_fw - disp_omega_fw, traps_axis_fw)
F_fw          = keff_fw * projection_fw
print("forward keff [pN/nm]:\n\t", su2pN_nm(keff_fw))



# equilibration trajectory
if eq:
	## calculate distances between nucleotides alpha and omega
	dr_eq = r_omega_eq - r_alpha_eq
	distances_eq = np.linalg.norm(dr_eq, axis=1)

	## calculate trajectories of each trap and of the trap-trap displacement
	trap_alpha_eq = np.zeros((num_configs_eq, 3), float)
	trap_omega_eq = np.zeros((num_configs_eq, 3), float)
	for k in range(num_configs_eq):
		trap_alpha_eq[k] = alpha_pos0_eq + steparr_eq[k] * alpha_rate_eq * alpha_dir_eq
		trap_omega_eq[k] = omega_pos0_eq + steparr_eq[k] * omega_rate_eq * omega_dir_eq
	traps_disp_eq = trap_omega_eq - trap_alpha_eq
	traps_dist_eq = np.linalg.norm(traps_disp_eq, axis=1)
	traps_axis_eq = traps_disp_eq / traps_dist_eq[:, None]

	## calculate nucleotide-trap displacements
	disp_alpha_eq = r_alpha_eq - trap_alpha_eq
	disp_omega_eq = r_omega_eq - trap_omega_eq
	dist_alpha_eq = np.linalg.norm(disp_alpha_eq, axis=1)
	dist_omega_eq = np.linalg.norm(disp_omega_eq, axis=1)

	## calculate the force for two springs connected in series along the trap-trap axis
	keff_eq       = alpha_stiff_eq * omega_stiff_eq / (alpha_stiff_eq + omega_stiff_eq)
	projection_eq = np.einsum('ij,ij->i', disp_alpha_eq - disp_omega_eq, traps_axis_eq)
	F_eq          = keff_eq * projection_eq
	print("equilibration keff [pN/nm]:\n\t", su2pN_nm(keff_eq))



# backward trajectory
## calculate distances between nucleotides alpha and omega
dr_bw = r_omega_bw - r_alpha_bw
distances_bw = np.linalg.norm(dr_bw, axis=1)

## calculate trajectories of each trap and of the trap-trap displacement
trap_alpha_bw = np.zeros((num_configs_bw, 3), float)
trap_omega_bw = np.zeros((num_configs_bw, 3), float)
for k in range(num_configs_bw):
	trap_alpha_bw[k] = alpha_pos0_bw + steparr_bw[k] * alpha_rate_bw * alpha_dir_bw
	trap_omega_bw[k] = omega_pos0_bw + steparr_bw[k] * omega_rate_bw * omega_dir_bw
traps_disp_bw = trap_omega_bw - trap_alpha_bw
traps_dist_bw = np.linalg.norm(traps_disp_bw, axis=1)
traps_axis_bw = traps_disp_bw / traps_dist_bw[:, None]

## calculate nucleotide-trap displacements
disp_alpha_bw = r_alpha_bw - trap_alpha_bw
disp_omega_bw = r_omega_bw - trap_omega_bw
dist_alpha_bw = np.linalg.norm(disp_alpha_bw, axis=1)
dist_omega_bw = np.linalg.norm(disp_omega_bw, axis=1)

## calculate the force for two springs connected in series along the trap-trap axis
keff_bw       = alpha_stiff_bw * omega_stiff_bw / (alpha_stiff_bw + omega_stiff_bw)
projection_bw = np.einsum('ij,ij->i', disp_alpha_bw - disp_omega_bw, traps_axis_bw)
F_bw          = keff_bw * projection_bw
print("backward keff [pN/nm]:\n\t", su2pN_nm(keff_bw))
print()



# smoothing
## concatenation
if eq:
	steparr_cat    = np.concatenate((steparr_fw, steparr_eq, steparr_bw))
	distances_cat  = np.concatenate((distances_fw, distances_eq, distances_bw))
	F_cat          = np.concatenate((F_fw, F_eq, F_bw))
	dist_alpha_cat = np.concatenate((dist_alpha_fw, dist_alpha_eq, dist_alpha_bw))
	dist_omega_cat = np.concatenate((dist_omega_fw, dist_omega_eq, dist_omega_bw))
else:
	steparr_cat    = np.concatenate((steparr_fw, steparr_bw))
	distances_cat  = np.concatenate((distances_fw, distances_bw))
	F_cat          = np.concatenate((F_fw, F_bw))
	dist_alpha_cat = np.concatenate((dist_alpha_fw, dist_alpha_bw))
	dist_omega_cat = np.concatenate((dist_omega_fw, dist_omega_bw))

## smoothing with exponential moving averages
alpha_ema         = 5e-3
m_ema             = 100
dema_cat          = ema(distances_cat, alpha_ema, m_ema)
Fema_cat          = ema(F_cat, alpha_ema, m_ema)

## definition of aliases
dema_fw = dema_cat[:num_configs_fw]
dema_eq = dema_cat[num_configs_fw:num_configs_fw+num_configs_eq]
dema_bw = dema_cat[num_configs_fw+num_configs_eq:]
Fema_fw = Fema_cat[:num_configs_fw]
Fema_eq = Fema_cat[num_configs_fw:num_configs_fw+num_configs_eq]
Fema_bw = Fema_cat[num_configs_fw+num_configs_eq:]

## determining zipping and unzipping forces, respectively
Fmin_fw           = 15.0/48.63
Fmax_bw           = 15.0/48.63
peak_distance     = m_ema
Fema_peaks_fw, _  = sp.signal.find_peaks(Fema_fw, distance=m_ema)
dist_peaks_fw     = dema_fw[Fema_peaks_fw]
selection_prob_fw = []
selection_poss_fw = []
for i in range(1, len(dist_peaks_fw) - 1):
	## lower experimental unzipping forces as lower limit
	if Fema_fw[Fema_peaks_fw[i]] > Fmin_fw:
		## during rupture, distance must increase
		if dist_peaks_fw[i] > dist_peaks_fw[i-1]:
			selection_poss_fw.append(Fema_peaks_fw[i])
			## the quality of this criterion depends on remaining oscillations after smoothing, i. a.
			if Fema_fw[Fema_peaks_fw[i]] > Fema_fw[Fema_peaks_fw[i+1]]:
				selection_prob_fw.append(Fema_peaks_fw[i])
print("forward trajectory: local maxima of force, followed by increasing distance:")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
for i in range(len(selection_poss_fw)):
	print("\t", steparr_fw[selection_poss_fw[i]], "\t", su2nm(dema_fw[selection_poss_fw[i]]), "\t", su2pN(Fema_fw[selection_poss_fw[i]]))
print("... with lower force at next local maximum:")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
for i in range(len(selection_prob_fw)):
	print("\t", steparr_fw[selection_prob_fw[i]], "\t", su2nm(dema_fw[selection_prob_fw[i]]), "\t", su2pN(Fema_fw[selection_prob_fw[i]]))
Fema_peaks_bw, _  = sp.signal.find_peaks(-Fema_bw, distance=m_ema)
dist_peaks_bw     = dema_bw[Fema_peaks_bw]
selection_prob_bw = []
selection_poss_bw = []
for i in range(1, len(dist_peaks_bw) - 1):
	## lower experimental unzipping forces as upper limit
	if Fema_bw[Fema_peaks_bw[i]] < Fmax_bw:
		## during rezipping, distance must decrease
		if dist_peaks_bw[i] < dist_peaks_bw[i-1]:
			selection_poss_bw.append(Fema_peaks_bw[i])
			## the quality of this criterion depends on remaining oscillations after smoothing, i. a.
			if Fema_bw[Fema_peaks_bw[i]] < Fema_bw[Fema_peaks_bw[i+1]]:
				selection_prob_bw.append(Fema_peaks_bw[i])
print("backward trajectory: local minima of force, followed by decreasing distance:")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
for i in range(len(selection_poss_bw)):
	print("\t", steparr_bw[selection_poss_bw[i]], "\t", su2nm(dema_bw[selection_poss_bw[i]]), "\t", su2pN(Fema_bw[selection_poss_bw[i]]))
print("... with higher force at next local minimum:")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
print()



# integration
## integration over noisy raw data
int_raw_fw = np.sum(0.5*(F_fw[1:]+F_fw[:-1]) * np.diff(distances_fw))
int_raw_eq = np.sum(0.5*(F_eq[1:]+F_eq[:-1]) * np.diff(distances_eq)) if eq else 0.0
int_raw_bw = np.sum(0.5*(F_bw[1:]+F_bw[:-1]) * np.diff(distances_bw))
int_raw    = int_raw_fw + int_raw_eq + int_raw_bw

## integration over smoothed data
int_ema_fw = np.sum(0.5*(Fema_fw[1:]+Fema_fw[:-1]) * np.diff(dema_fw))
int_ema_eq = np.sum(0.5*(Fema_eq[1:]+Fema_eq[:-1]) * np.diff(dema_eq)) if eq else 0.0
int_ema_bw = np.sum(0.5*(Fema_bw[1:]+Fema_bw[:-1]) * np.diff(dema_bw))
int_ema    = int_ema_fw + int_ema_eq + int_ema_bw

## output results
print("results of force-extension integration:")
print("[J]", "\t", "raw\t\t\t",      "  \t", "ema")
print("fw",  "\t", su2J(int_raw_fw), "  \t", su2J(int_ema_fw))
print("eq",  "\t", su2J(int_raw_eq), "  \t", su2J(int_ema_eq))
print("bw",  "\t", su2J(int_raw_bw), "  \t", su2J(int_ema_bw))
print("---------------------------------------------------------------")
print("sum", "\t", su2J(int_raw),    "  \t", su2J(int_ema))
print()
print("calculated molar hysteresis energy [kJ/mol]:")
print("\traw:", 6.02214076e20*su2J(int_raw))
print("\tema:", 6.02214076e20*su2J(int_ema))
print()



# plotting

## distances between nucleotides
fig1, ax1 = plt.subplots()
ax1.set_xlabel("simulation steps")
ax1.set_ylabel("distance [nm]")
ax1.plot(steparr_cat, su2nm(distances_cat),  color='black', alpha=0.33, label="distance between nucleotides")
ax1.plot(steparr_cat, su2nm(dema_cat),       color='black', alpha=1.00, label="exp. moving average")
ax1.plot(steparr_cat, su2nm(dist_alpha_cat), color='blue',  alpha=0.33, label="distance to trap (alpha)")
ax1.plot(steparr_cat, su2nm(dist_omega_cat), color='red',   alpha=0.33, label="distance to trap (omega)")
ax1.legend()
plt.grid()
fig1.savefig(f'forcys_distances.png')
plt.cla()
plt.clf()
print("forcys_distances.png <- trajectories plotted")
print()


## force-extension plot
fig2, ax2 = plt.subplots()
ax2.set_xlabel("distance [nm]")
ax2.set_ylabel("force [pN]")
ax2.plot(su2nm(distances_fw), su2pN(F_fw), color='red', alpha=0.33, label="projected force (forward)")
if eq:
	ax2.plot(su2nm(distances_eq), su2pN(F_eq), color='grey', alpha=0.33, label="projected force (equilibration)")
ax2.plot(su2nm(distances_bw), su2pN(F_bw), color='blue', alpha=0.33, label="projected force (backward)")
ax2.plot(su2nm(dema_cat), su2pN(Fema_cat), color='black', alpha=1.00, label="exp. moving average")
if eq:
	dema_temp = np.append(dema_fw, dema_eq)
	Fema_temp = np.append(Fema_fw, Fema_eq)
	plt.fill(su2nm(np.append(dema_temp, dema_bw)), su2pN(np.append(Fema_temp, Fema_bw)), color='grey', alpha=0.33)
else:
	plt.fill(su2nm(np.append(dema_fw, dema_bw)), su2pN(np.append(Fema_fw, Fema_bw)), color='grey', alpha=0.33)
ax2.legend()
plt.grid()
fig2.savefig(f'forcys_hysteresis.png')
plt.cla()
plt.clf()
print("forcys_hysteresis.png <- force-extension hysteresis curve plotted")
print()




