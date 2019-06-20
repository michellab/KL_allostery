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

print ("All input files: ") , (input_files)

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
print (top)

traj = md.load_dcd(traj_input,top=pdb_input)
print (traj)

print ("")
print ("Computing CA contacts...")
CA_contacts = md.compute_contacts(traj, contacts='all', scheme="ca")

distance_per_snapshot = CA_contacts[0]
indices_per_snapshot = CA_contacts[1]

# Output files with all atom pairs, and with all distances vs all snapshots
print ("")
np.savetxt("%s/ALL_atom_pairs.dat" % outfile , CA_contacts[1], fmt='%s')
print ("Saved atom pairs to %s/ALL_atom_pairs.dat" % outfile)
min_list = []
max_list = []
for i in range(0,len(distance_per_snapshot[:][0])):
    min_list.append(min(distance_per_snapshot[:,i]))
    max_list.append(max(distance_per_snapshot[:,i]))

min_array = np.array(min_list)
max_array = np.array(max_list)
min_max_array = np.vstack((min_list,max_list)).T

print ("")
min_max_array_angstrom = min_max_array * 10
np.savetxt("%s/min_max_rawdata.dat" % outfile,min_max_array_angstrom)
print ("Saved min max values to %s/min_max_rawdata.dat" % outfile)

# Also need to output a file with the atom numbers of CA atoms instead of the residue numbers 
# col1 and col2 are the two columns of RESIDUE numbers - which make each pair
# Want to make the same array but with atom numbers of CA atoms

column1 = CA_contacts[1][:,0]
column2 = CA_contacts[1][:,1]

column_1_atom_num = []
column_2_atom_num = []

for i in column1:
    column_1_atom_num.append(int(top.select("name CA and resid %s"%i)))

for i in column2:
    column_2_atom_num.append(int(top.select("name CA and resid %s"%i)))

column_1_atom_array = np.array(column_1_atom_num)
column_2_atom_array = np.array(column_2_atom_num)

atom_number_pair_array = np.vstack((column_1_atom_array,column_2_atom_array)).T
print ("")
np.savetxt("%s/atom_pairs_ATOM_NUMBERS.dat" % outfile, atom_number_pair_array, fmt='%s')
print ("Saved atom indices to %s/atom_pairs_ATOM_NUMBERS.dat" % outfile)

