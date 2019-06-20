# coding: utf-8
import numpy as np
import mdtraj as md
import sys
from sys import argv
import math
import os
#get_ipython().magic(u'pylab inline')

script, input_system = argv
input_system = int(input_system)

# system_list with folder names of different systems. 
system_list = ["0_system_A","1_system_B"]
# Input one or more trajectory names into filname list. 
filename_list = ["short_traj_aligned.dcd"]

topology_filename = "first_frame.pdb"

md_data = ["0_TRAJECTORIES"]

filename_list_1_traj = []
filename_list_1_pdb = []

# Make a list with all file locations of trajectory data
all_files_list = []

for i in range(0,len(system_list)):
    for j in range(0,len(filename_list)):
        filename_traj = "%s/%s/%s" % (md_data[0],system_list[i],filename_list[j])
        filename_list_1_traj.append(filename_traj)
        filename_pdb = "%s/%s/%s" % (md_data[0],system_list[i],topology_filename)
        filename_list_1_pdb.append(filename_pdb)


# Make a list of lists to separate file locations for each simulation.
input_files = []
for i in range(0,len(system_list)):
    inside_list = []
    for j in range(0,len(filename_list)):
        filenames = "%s/%s/%s" % (md_data[0],system_list[i],filename_list[j])
        inside_list.append(filenames)
    input_files.append(inside_list)

for i in system_list:
    if not os.path.exists("4_CA_DISTANCES/%s/OUTPUT/CA_dist" % i):
        filename = "4_CA_DISTANCES/%s/OUTPUT/CA_dist" % i
        cmd = "mkdir -p %s" % filename
        os.system(cmd)
        
for i in system_list:
    if not os.path.exists("4_CA_DISTANCES/%s/OUTPUT/CA_raw_data" % i):
        filename = "4_CA_DISTANCES/%s/OUTPUT/CA_raw_data" % i
        cmd = "mkdir -p %s" % filename
        os.system(cmd)

outfile = "4_CA_DISTANCES/%s/OUTPUT" % system_list[input_system]
traj_input = filename_list_1_traj[input_system]
pdb_input = filename_list_1_pdb[input_system]

print ("Trajectory input: ") , traj_input
print ("Topology input: ") , pdb_input

test = md.load_pdb(pdb_input)
top = test.topology

traj = md.load_dcd(traj_input,top=pdb_input)
print ("")
print ("Computing CA contacts...")
CA_contacts = md.compute_contacts(traj, contacts='all', scheme="ca")

distance_per_snapshot = CA_contacts[0]
indices_per_snapshot = CA_contacts[1]

print ("")
print ("Loading min max values from 4_CA_DISTANCES/global_min_max_array.dat")

min_max_arr_margin_int  = np.loadtxt("4_CA_DISTANCES/global_min_max_array.dat")

# bin values already in angstrom
MINBIN = min_max_arr_margin_int[:,0]
MAXBIN = min_max_arr_margin_int[:,1]

print ("Min bin :"), (MINBIN)
print ("Max bin :"), (MAXBIN)

for i in range(0,len(distance_per_snapshot[0,:])):
    dist_angstrom = distance_per_snapshot[:,i] * 10
    # load bin ranges from min_max file
    min_bin = MINBIN[i]
    max_bin = MAXBIN[i]
    (n, bins) = np.histogram(dist_angstrom, bins = 60, range = (min_bin, max_bin), normed=True)
    n = n / (sum(n))
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    total_bin_addition = 0.000001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    if zero_bins != 0:
        bin_addition = total_bin_addition/float(zero_bins)
        # Adds the bin_addition amount into all zero-count bins
        for j in xrange(len(n)):
            if n[j] == 0.0:
                n[j] = bin_addition
    data = np.vstack((index,n)).T
    np.savetxt("%s/CA_dist/distance_%s_distribution.dat" % (outfile,i) , data, fmt=['%d','%.20f'])   

print ("Distribution data saved to %s/CA_dist" % outfile)

