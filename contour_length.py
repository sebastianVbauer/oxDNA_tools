'''
You can run the analysis by placing it into the same directory as the data and then executing
	python contour_length_x005.py xxx.top yyy.dat route site X
with
	xxx.top	the name of the topology file in the 'classic' oxDNA format, i. e. 3'->5' (can be a path as well),
	yyy.dat	the name of the configuration or trajectory file (can be a path as well),
	route	the route along which nucleotide-to-nucleotide distances are summed up to give the contour length:
		This argument takes a string of nucleotide indices, e. g. "51 52 3 4 5 6 7 81 80 79 78 77".
		Successive nucleotides along the same strand can be abbreviated as ..., e. g. "51 52 3 ... 7 81 ... 77".
	site	the site used to define nucleotide positions, can be:
		bb	backbone repulsion site (default),
		cm	center-of-mass position,
		st	stacking site,
		hb	hydrogen-bonding/repulsion site,
	X	an optional argument; if set, the first X configurations in the trajectory file will not be evaluated (defaults to 0).

This code is cc-by-sa.
-- Sebastian V. Bauer, University of Mainz, 21.03.2025
-- sebastian.bauer@uni-mainz.de
'''

import sys
import numpy as np
import time

###############################################
# customized functions to evaluate oxDNA data #
###############################################

# get vectors r, b, n for a single nucleotide from an oxDNA configuration file
def eval_oxDNA_config_line(line):
	line_arr = (np.asarray(line.split())).astype(float)
	r = line_arr[0:3]
	b = line_arr[3:6]
	n = line_arr[6:9]
	return (r, b, n)

# determine position for a single nucleotide from an oxDNA2 configuration file
def get_oxDNA2_site(line, site="bb"):
	r, b, n = eval_oxDNA_config_line(line)
	if site == "bb":
		y = np.cross(n, b)
		return r - 0.34*b + 0.3408*y
	elif site == "cm":
		return r
	elif site == "st":
		return r + 0.34*b
	elif site == "hb":
		return r + 0.4*b

# convert length from oxDNA simulation units to nanometer
def su2nm(x):
	return 0.8518 * x



###################
# data extraction #
###################

# parse positional arguments
if len(sys.argv) < 5:
	print("python contour_length_x005.py topology.top trajectory.dat route site X")
	print("You must define the route along which nucleotide-to-nucleotide distances are summed up to give the contour length.")
	print("\tThis argument takes a string of nucleotide indices, e. g. \"51 52 3 4 5 6 7 81 80 79 78 77\".")
	print("\tSuccessive nucleotides along the same strand can be abbreviated as ..., e. g. \"51 52 3 ... 7 81 ... 77\".")
	print("You must select which positions are to be used as reference positions for nucleotides by setting the argument site:")
	print("\tbb for backbone repulsion site (default)")
	print("\tcm for center-of-mass position")
	print("\tst for stacking site")
	print("\thb for hydrogen-bonding/repulsion site")
	print("You can choose whether to skip the first X configurations by setting the last positional argument (defaults to 0).")
	sys.exit(2)
topfile_path  = sys.argv[1]
trajfile_path = sys.argv[2]
route         = sys.argv[3]
site          = sys.argv[4]
if site == "bb":
	print("evaluating nucleotides at backbone repulsion site")
elif site == "cm":
	print("evaluating nucleotides at center-of-mass position")
elif site == "st":
	print("evaluating nucleotides at stacking site")
elif site == "hb":
	print("evaluating nucleotides at hydrogen-bonding/repulsion site")
else:
	print("no valid argument provided for site, evaluating nucleotides at backbone repulsion sites (default)")
	site = "bb"
try:
	skip = int(sys.argv[5])
except:
	skip = 0
	print("all configurations found in the trajectory file will be evaluated\n")
else:
	print("the first", skip, "configurations found in the trajectory file will be skipped\n")



# determine chain length of first strand in nt from topology file
topfile = open(topfile_path, "r")
line = topfile.readline()
line_lst = line.split()
number_of_nucleotides = int(line_lst[0])
number_of_strands = int(line_lst[1])
chain_length = 0
if number_of_strands == 1:
	chain_length = number_of_nucleotides
else:
	while line:
		line = topfile.readline()
		if int(line.split()[0]) > 1:
			break
		else:
			chain_length += 1
topfile.close()
print(topfile_path, "-> chain length of evaluated strand [nt]:", chain_length)



# browse whole trajectory file to get number of configurations
all_configs = 0
trajfile = open(trajfile_path, "r")
line = trajfile.readline()
while line:
	for n in range (number_of_nucleotides + 2):
		line = trajfile.readline()
	all_configs += 1
	line = trajfile.readline()
trajfile.close()
num_configs = all_configs - skip
print(trajfile_path, "-> found", all_configs, "configurations")
if num_configs < 1:
	print(trajfile_path, "-> no configurations left after skipping, aborting.")
	sys.exit(1)
else:
	print(trajfile_path, "-> evaluating", num_configs, "configurations\n")



# definition of sequence of interest
seq_int_raw = route.split()
seq_int_tmp = []
for s in seq_int_raw:
	try:
		seq_int_tmp.append(int(s))
	except:
		seq_int_tmp.append("...")
seq_int_lst = []
for i in range(len(seq_int_tmp)):
	if seq_int_tmp[i] != "...":
		seq_int_lst.append(seq_int_tmp[i])
	else:
		pre = seq_int_tmp[i-1]
		suc = seq_int_tmp[i+1]
		p2s = 1 if pre < suc else -1
		for nt in range(pre+p2s, suc, p2s):
			seq_int_lst.append(nt)
seq_int_srt = seq_int_lst.copy()
seq_int_srt.sort()
seq_int_len = len(seq_int_lst)
seq_int_ind = []
for n in range(seq_int_len):
	seq_int_ind.append(seq_int_srt.index(seq_int_lst[n]))
print("seq_int_raw =", seq_int_raw)
print("seq_int_tmp =", seq_int_tmp)
print("sequence of interest:\n\t", seq_int_lst)
print("-> sorted:\n\t", seq_int_srt)
print("-> cross-indexing:\n\t", seq_int_ind)
print()



# extract positions from trajectory file
seq_int_arr = np.zeros((num_configs, seq_int_len, 3), float)
trajfile = open(trajfile_path, "r")
for k in range(skip):
	## skip unwanted configurations
	for l in range(number_of_nucleotides + 3):
		line = trajfile.readline()
for k in range(num_configs):
	## scroll down to 3'-end
	for l in range(4):
		line = trajfile.readline()
	
	## scroll down to first nucleotide within the sorted sequence of interest and collect position
	for l in range(seq_int_srt[0]):
		line = trajfile.readline()
		seq_int_arr[k, 0] = get_oxDNA2_site(line, site)
	
	## proceed to the next nucleotides within the sorted sequence of interest and collect positions
	for n in range(1, seq_int_len):
		for l in range(seq_int_srt[n] - seq_int_srt[n-1]):
			line = trajfile.readline()
			seq_int_arr[k, n] = get_oxDNA2_site(line, site)
	
	## scroll down to next configuration
	for l in range (number_of_nucleotides - seq_int_srt[-1] - 1):
		line = trajfile.readline()
trajfile.close()



###############
# do the math #
###############

seq_dst_arr = np.zeros((num_configs, seq_int_len-1), float)
for n in range(seq_int_len-1):
	this = seq_int_ind[n]
	next = seq_int_ind[n+1]
	seq_dst_arr[:, n] = np.linalg.norm(seq_int_arr[:, next] - seq_int_arr[:, this], axis=1)

seq_dst_avg_per_n2n = np.mean(seq_dst_arr, axis=0)
seq_dst_std_per_n2n = np.std(seq_dst_arr, axis=0)
print("seq_dst_avg_per_n2n [nm] =", su2nm(seq_dst_avg_per_n2n))
print("seq_dst_std_per_n2n [nm] =", su2nm(seq_dst_std_per_n2n))
print()

contour_length_over_time = np.sum(seq_dst_arr, axis=1)
contour_length_avg = np.mean(contour_length_over_time)
contour_length_std = np.std(contour_length_over_time)
print("contour_length_over_time [nm] =", su2nm(contour_length_over_time))
print("contour_length_avg [nm] =", su2nm(contour_length_avg))
print("contour_length_std [nm] =", su2nm(contour_length_std))






