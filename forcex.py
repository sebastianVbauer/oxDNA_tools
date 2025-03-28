'''
This script can be used to analyze the data from oxDNA simulations with two harmonic traps as external forces, each acting on a specific nucleotide.
Force is calculated from nucleotide-trap displacements as projection of the total force for two springs connected in series onto the axis between the springs.

The easiest way to run the analysis is to copy the script to the same directory as the raw data and execute
	python forcex.py file1 file2 file3 [mark] [skip1] [skip2]
with
	python	an interpreter for Python 3.x (usually python3 or python if using conda),
	file1	topology file (only classic oxDNA format, i. e. 3'->5'),
	file2	trajectory file,
	file3	force file (forces need to be derived from exactly two harmonic traps),
	[mark]	keyword to set options for the force-extension plot, could be
		poss	to mark possible rupture events only,
		prob	to mark probable rupture events only,
		both	to mark both possible and probable rupture events,
		none	to not mark anything (default),
	[skip1]	number of configurations to skip at the beginning (defaults to 0),
	[skip2]	number of configurations to skip at the end (defaults to 0).

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
* critical values of force-distance pairs (written to standard output)
* trajectories of nucleotide-nucleotide and both nucleotide-trap distances (written to file distances.png)
* force-extension curve (written to file force_extension.png)

This code is cc-by-sa.
-- Sebastian V. Bauer, Walther Lab at the University of Mainz, 05.04.2024
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
if len(sys.argv) < 4:
	print("execute\n\tpython forcex.py file1 file2 file3 [mark] [skip1] [skip2]\nwith")
	print("""
	file1	topology file (only classic oxDNA format, i. e. 3'->5'),
	file2	trajectory file,
	file3	force file (forces need to be derived from exactly two harmonic traps),
	[mark]	keyword to set options for the force-extension plot, could be
		poss	to mark possible rupture events only,
		prob	to mark probable rupture events only,
		both	to mark both possible and probable rupture events,
		none	to not mark anything (default),
	[skip1]	number of configurations to skip at the beginning (defaults to 0),
	[skip2]	number of configurations to skip at the end (defaults to 0).
	""")
	sys.exit(2)
topfile_path   = sys.argv[1]
trajfile_path  = sys.argv[2]
forcefile_path = sys.argv[3]
print("evaluating nucleotides at their center of mass")
try:
	mark = sys.argv[4]
except:
	mark = "none"
	print("no events will be marked within the force-extension plot")
else:
	if mark == "poss":
		print("only possible rupture events will be marked within the force-extension plot")
	elif mark == "prob":
		print("only probable rupture events will be marked within the force-extension plot")
	elif mark == "both":
		print("both possible and probable rupture events will be marked within the force-extension plot")
	else:
		print("no events will be marked within the force-extension plot")
try:
	skip_first = int(sys.argv[5])
except:
	skip_first = 0
	print("trajectory file will be evaluated from the beginning")
else:
	print("the first", skip_first, "configurations in the trajectory file will be skipped")
try:
	skip_last = int(sys.argv[6])
except:
	skip_last = 0
	print("trajectory file will be evaluated until the end")
else:
	print("the last", skip_last, "configurations in the trajectory file will be skipped")	
print()



# evaluate topology file
topfile = open(topfile_path, 'r')
line = topfile.readline()
line_lst = line.split()
number_of_nucleotides = int(line_lst[0])
topfile.close()
print(topfile_path, "-> number of nucleotides:", number_of_nucleotides)
print()



# browse whole trajectory file to get number of configurations
all_configs = 0
trajfile = open(trajfile_path, 'r')
line = trajfile.readline()
while line:
	for n in range (number_of_nucleotides + 2):
		line = trajfile.readline()
	all_configs += 1
	line = trajfile.readline()
trajfile.close()
num_configs = all_configs - (skip_first + skip_last)
print(trajfile_path, "-> found", all_configs, "configurations")
if num_configs < 1:
	print(trajfile_path, "-> no configurations left after skipping, aborting.")
	sys.exit(1)
else:
	print(trajfile_path, "-> evaluating", num_configs, "configurations")
print()



# evaluate external forces file
forcefile = open(forcefile_path, 'r')
alpha_check = False
omega_check = False
line = forcefile.readline()
while line:
	if line[0] == '{':
		## get values for trap at nucleotide alpha
		if alpha_check == False:
			for l in range(2):
				line = forcefile.readline()
			alpha = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_pos0 = get_oxDNA_ffval(line)
			line = forcefile.readline()
			alpha_stiff = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_rate = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			alpha_dir = get_oxDNA_ffval(line)
			if (np.linalg.norm(alpha_dir) > 0.0):
				alpha_dir /= np.linalg.norm(alpha_dir)
			alpha_check = True
		## get values for trap at nucleotide omega
		else:
			for l in range(2):
				line = forcefile.readline()
			omega = int(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_pos0 = get_oxDNA_ffval(line)
			line = forcefile.readline()
			omega_stiff = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_rate = float(np.asarray(line.split())[2])
			line = forcefile.readline()
			omega_dir = get_oxDNA_ffval(line)
			if (np.linalg.norm(omega_dir) > 0.0):
				omega_dir /= np.linalg.norm(omega_dir)
			omega_check = True
	line = forcefile.readline()
forcefile.close()
if alpha_check == False:
	print(forcefile_path, "-> unable to determine harmonic trap alpha, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path, "-> found harmonic trap alpha for nucleotide", alpha)
	print("\tinitial position [s.u.]:", alpha_pos0)
	print("\tstiffness [s.u.]:", alpha_stiff)
	print("\tmovement rate [s.u.]:", alpha_rate)
	print("\tmovement direction [s.u.]:", alpha_dir)
if omega_check == False:
	print(forcefile_path, "-> unable to determine harmonic trap omega, aborting.\n")
	sys.exit(1)
else:
	print(forcefile_path, "-> found harmonic trap omega for nucleotide", omega)
	print("\tinitial position [s.u.]:", omega_pos0)
	print("\tstiffness [s.u.]:", omega_stiff)
	print("\tmovement rate [s.u.]:", omega_rate)
	print("\tmovement direction [s.u.]:", omega_dir)
print()



# evaluate configurations in trajectory file
r_alpha = np.zeros((num_configs, 3), float)
r_omega = np.zeros((num_configs, 3), float)
steparr = np.zeros(num_configs, int)
trajfile = open(trajfile_path, 'r')
for k in range(skip_first):
	## skip unwanted configurations at the beginning
	for l in range(number_of_nucleotides + 3):
		line = trajfile.readline()
for k in range(num_configs):
	## get time step at which the current configuration was recorded
	line = trajfile.readline()
	steparr[k] = get_oxDNA_step(line)
	
	## skip box size and energy
	for l in range(2):
		line = trajfile.readline()
	
	## scroll down nucleotides
	for n in range(alpha):
		line = trajfile.readline()
	
	## get position of nucleotide alpha
	line = trajfile.readline()
	r_alpha[k] = get_oxDNA_cm_site(line)
	
	## scroll down nucleotides
	for n in range(omega - alpha - 1):
		line = trajfile.readline()
	
	## get position of nucleotide omega
	line = trajfile.readline()
	r_omega[k] = get_oxDNA_cm_site(line)
	
	## scroll down to next configuration
	for n in range(number_of_nucleotides - omega - 1):
		line = trajfile.readline()
	## don't evaluate remaining configurations
trajfile.close()



#################
# data analysis #
#################

# calculate distances between nucleotides alpha and omega
dr = r_omega - r_alpha
distances = np.linalg.norm(dr, axis=1)



# calculate trajectories of each trap and of the trap-trap displacement
trap_alpha = np.zeros((num_configs, 3), float)
trap_omega = np.zeros((num_configs, 3), float)
for k in range(num_configs):
	trap_alpha[k] = alpha_pos0 + steparr[k] * alpha_rate * alpha_dir
	trap_omega[k] = omega_pos0 + steparr[k] * omega_rate * omega_dir
traps_disp = trap_omega - trap_alpha
traps_dist = np.linalg.norm(traps_disp, axis=1)
traps_axis = traps_disp / traps_dist[:, None]



# calculate nucleotide-trap displacements
disp_alpha = r_alpha - trap_alpha
disp_omega = r_omega - trap_omega
dist_alpha = np.linalg.norm(disp_alpha, axis=1)
dist_omega = np.linalg.norm(disp_omega, axis=1)



# calculate the force for two springs connected in series along the trap-trap axis
eff_stiff  = alpha_stiff * omega_stiff / (alpha_stiff + omega_stiff)
projection = np.einsum('ij,ij->i', disp_alpha - disp_omega, traps_axis)
F          = eff_stiff * projection
Fmax_index = np.argmax(F)
Fmax_value = F[Fmax_index]
Fmax_step  = steparr[Fmax_index]
print("keff [pN/nm]:\n\t", su2pN_nm(eff_stiff))
print()



# smoothing with exponential moving averages
alpha_ema     = 5e-3
m_ema         = 100
dema          = ema(distances, alpha_ema, m_ema)
Fema          = ema(F, alpha_ema, m_ema)
Femamax_index = np.argmax(Fema)
Femamax_value = Fema[Femamax_index]
Femamax_step  = steparr[Femamax_index]



# global maxima in projected force (single and e.m.a.)
print("Global force maxima:")
print("\t", "time step no.", "\t", "distance [nm]       ", "\t", "projected force [pN]")
print("\t", Fmax_step, "\t", su2nm(distances[Fmax_index]), "\t", su2pN(Fmax_value))
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
print("\t", Femamax_step, "\t", su2nm(distances[Femamax_index]), "\t", su2pN(Femamax_value))



# search for rupture events
Fmin           = 15.0/48.63
peak_distance  = m_ema
Fema_peaks, _  = sp.signal.find_peaks(Fema, distance=m_ema)
dist_peaks     = dema[Fema_peaks]
selection_prob = []
selection_poss = []
for i in range(1, len(dist_peaks) - 1):
	## lower experimental unzipping forces as lower limit
	if Fema[Fema_peaks[i]] > Fmin:
		## during rupture, distance must increase
		if dist_peaks[i] > dist_peaks[i-1]:
			selection_poss.append(Fema_peaks[i])
			## the quality of this criterion depends on remaining oscillations after smoothing, i. a.
			if Fema[Fema_peaks[i]] > Fema[Fema_peaks[i+1]]:
				selection_prob.append(Fema_peaks[i])
print("Possible rupture events detected (local maxima of force, followed by increasing distance):")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
for i in range(len(selection_poss)):
	print("\t", steparr[selection_poss[i]], "\t", su2nm(dema[selection_poss[i]]), "\t", su2pN(Fema[selection_poss[i]]))
print("Probable rupture events detected (as before, and with lower force at next local maximum):")
print("\t", "time step no.", "\t", "e.m.a. distance [nm]", "\t", "e.m.a. projected force [pN]")
for i in range(len(selection_prob)):
	print("\t", steparr[selection_prob[i]], "\t", su2nm(dema[selection_prob[i]]), "\t", su2pN(Fema[selection_prob[i]]))
print()



# plotting

## distances between nucleotides
fig1, ax1 = plt.subplots()
ax1.set_xlabel("simulation steps")
ax1.set_ylabel("distance [nm]")
ax1.plot(steparr, su2nm(distances),  color='black', alpha=0.33, label="distance between nucleotides")
ax1.plot(steparr, su2nm(dema),       color='black', alpha=1.00, label="exp. moving average")
ax1.plot(steparr, su2nm(dist_alpha), color='blue',  alpha=0.33, label="distance to trap (alpha)")
ax1.plot(steparr, su2nm(dist_omega), color='red',   alpha=0.33, label="distance to trap (omega)")
ax1.legend()
plt.grid()
fig1.savefig(f'distances.png')
plt.cla()
plt.clf()
print("distances.png <- trajectories plotted")
print()

## force-extension plot
fig2, ax2 = plt.subplots()
ax2.set_xlabel("distance [nm]")
ax2.set_ylabel("force [pN]")
ax2.plot(su2nm(distances),            su2pN(F),                     color='green', alpha=0.33, label="projected force")
ax2.plot(su2nm(dema),                 su2pN(Fema),                  color='black', alpha=1.00, label="exp. moving average")
if mark in ["poss", "both"]:
	ax2.plot(su2nm(dema[selection_poss]), su2pN(Fema[selection_poss]), 'bx', label="possible rupture event")
if mark in ["prob", "both"]:
	ax2.plot(su2nm(dema[selection_prob]), su2pN(Fema[selection_prob]), 'rx', label="probable rupture event")
ax2.legend()
plt.grid()
fig2.savefig(f'force_extension.png')
plt.cla()
plt.clf()
print("force_extension.png <- force-extension curve plotted")
print()





