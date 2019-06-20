# coding: utf-8
import scipy as sp
import numpy as np
import mdtraj as md
import sys
import os
from sys import argv
import math
#get_ipython().magic(u'pylab inline')

script, system , number_of_bins = argv
# *system* is a number (zero indexed) corresponding to the list of systems you input into system_list. 
# *number_of_bins* is for histogram data. Aim for around 15 bins per 1000 snapshots. (Ideally from a trajectory with around 100k snapshots for 1 microsecond simulation).

input_system = int(system)

# edit system_list with folder names of different systems.
system_list = ["0_system_A","1_system_B"]

# Input one or more trajectory names into filname list.
filename_list = ["traj0001.dcd","traj0002.dcd"]

topology_filename = "topology.parm7"

md_data = ["0_TRAJECTORIES"]

print ("2.1: input system is: %s" % system_list[input_system])
print ("MD data located in directory: %s" % md_data[0])
print ("Topology filename is: %s" % topology_filename)
print ("Trajectories are %s" % filename_list)
	
# Make a list with all file locations of trajectory data
all_files_list = []

for i in range(0,len(system_list)):
    for j in range(0,len(filename_list)):
        filenames = "%s/%s/%s" % (md_data[0],system_list[i],filename_list[j])
        all_files_list.append(filenames)
# Make a list of lists to separate file locations for each simulation.
input_files = []
for i in range(0,len(system_list)):
    inside_list = []
    for j in range(0,len(filename_list)):
        filenames = "%s/%s/%s" % (md_data[0],system_list[i],filename_list[j])
        inside_list.append(filenames)
    input_files.append(inside_list)

print ("")
print ("2.2: list of all input file locations.") 
print (input_files)

print ("")
print ("2.3: Creating output folders for all systems")

for i in system_list:
    if not os.path.exists("1_DIHEDRALS/%s/OUTPUT/PSI" % i):
        filename = "1_DIHEDRALS/%s/OUTPUT/PSI" % i
        cmd = "mkdir -p %s" % filename
	print ("Folder created for %s" % filename)
        os.system(cmd)
        
for i in system_list:
    if not os.path.exists("1_DIHEDRALS/%s/OUTPUT/PHI" % i):
        filename = "1_DIHEDRALS/%s/OUTPUT/PHI" % i
        cmd = "mkdir -p %s" % filename
	print ("Folder created for %s" % filename)
        os.system(cmd)
        
for i in system_list:
    if not os.path.exists("1_DIHEDRALS/%s/OUTPUT/CHI1" % i):
        filename = "1_DIHEDRALS/%s/OUTPUT/CHI1" % i
        cmd = "mkdir -p %s" % filename
	print ("Folder created for %s" % filename)
        os.system(cmd)
        
for i in system_list:
    if not os.path.exists("1_DIHEDRALS/%s/OUTPUT/CHI2" % i):
        filename = "1_DIHEDRALS/%s/OUTPUT/CHI2" % i
        cmd = "mkdir -p %s" % filename
	print ("Folder created for %s" % filename)
        os.system(cmd)  

print ("")
print ("2.4: Loading trajectory and calculating for system: ") 
print ("")
print (system_list[input_system])
print ("")
# Start by loading the first trajectory

traj = md.load_dcd(input_files[input_system][0] , top="%s/%s/%s" % (md_data[0],system_list[input_system],topology_filename))
top = traj.topology
print ("One trajectory file traj0001.dcd is: ") 
print (traj)

# Then loop through the remaining trajectories to add them to the first.

for i in range(1,len(filename_list)):
    traj_ADD = md.load_dcd(input_files[input_system][i] , top="%s/%s/%s" % (md_data[0],system_list[input_system],topology_filename))
    traj = traj.join(traj_ADD,check_topology=True, discard_overlapping_frames=False)
    
print ("Combined traj0001.dcd and traj0002.dcd is: ") 
print (traj)

## Uncomment to save pdb.
## traj[0].save_pdb("%s/%s/first_frame.pdb" % (md_data[0],system_list[input_system]))

print ("")
print ("2.5: Computing Dihedrals")

print ("computing psi")
psi_list = md.compute_psi(traj, periodic=True, opt=True)

print ("computing phi")
phi_list = md.compute_phi(traj, periodic=True, opt=True)

print ("computing chi1")
chi1_list = md.compute_chi1(traj, periodic=True, opt=True)

print ("computing chi2")
chi2_list = md.compute_chi2(traj, periodic=True, opt=True)

number_angles_list = []
number_angles_list.append(len(psi_list[0][:,0]))
number_angles_list.append(len(phi_list[0][:,0]))
number_angles_list.append(len(chi1_list[0][:,0]))
number_angles_list.append(len(chi2_list[0][:,0]))
number_angles_array = np.vstack((np.array(['Number of psi','Number of phi','Number of chi1','Number of chi2']),np.array(number_angles_list))).T

print ("The system has the following numbers of psi, phi, chi1 and chi2: ") 
print (number_angles_array)

print ("Saving number of angles to file: 1_DIHEDRALS/number_angles.dat")
np.savetxt("1_DIHEDRALS/number_angles.dat",number_angles_array,fmt='%s')

print ("")
print ("2.6: Output file with CA index and torsions: ")
c_alphas = top.select("name CA")
c_alphas_list = []
for i in c_alphas:
    c_alphas_list.append(i)

psi_ca_list = []
for i in psi_list[0][:,1]:
    psi_ca_list.append(i)
phi_ca_list = []
for i in phi_list[0][:,2]:
    phi_ca_list.append(i)
chi1_ca_list = []
for i in chi1_list[0][:,1]:
    chi1_ca_list.append(i)
chi2_ca_list = []
for i in chi2_list[0][:,0]:
    chi2_ca_list.append(i)

psi_all_list = []
phi_all_list = []
chi1_all_list = []
chi2_all_list = []

for i in xrange(0, len(c_alphas)):
    CA_index = c_alphas[i]
    if CA_index in psi_list[0][:,1]:
        psi_all_list.append(psi_ca_list.index(CA_index))
    else:
        psi_all_list.append("none")
    if CA_index in phi_list[0][:,2]:
        phi_all_list.append(phi_ca_list.index(CA_index))
    else:
        phi_all_list.append("none")
    if CA_index in chi1_list[0][:,1]:
        chi1_all_list.append(chi1_ca_list.index(CA_index))
    else:
        chi1_all_list.append("none")
    if CA_index in chi2_list[0][:,0]:
        chi2_all_list.append(chi2_ca_list.index(CA_index))
    else:
        chi2_all_list.append("none")

psi_all_arr = np.array(psi_all_list)
phi_all_arr = np.array(phi_all_list)
chi1_all_arr = np.array(chi1_all_list)
chi2_all_arr = np.array(chi2_all_list)

CA_psi_phi_chi1_chi2_indexes = np.vstack((c_alphas,psi_all_arr,phi_all_arr,chi1_all_arr,chi2_all_arr))
CA_psi_phi_chi1_chi2_indexes = CA_psi_phi_chi1_chi2_indexes.T
np.savetxt("1_DIHEDRALS/%s/DIHEDRALS_CA_angles_indexes.dat" % system_list[input_system],CA_psi_phi_chi1_chi2_indexes,fmt='%s')
print ("File saved with CA indices and torision numbers to: 1_DIHEDRALS/%s/DIHEDRALS_CA_angles_indexes.dat" % system_list[input_system])

print ("")
print ("2.7: Calculating dihedral distributions: ")
print ("Computing psi distributions")
for i in xrange(psi_list[1].shape[1]):
    dihedral_traj = psi_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    # Save torsion to output
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/PSI/raw_data_psi_%d.dat' % (system_list[input_system],(i+1)),dihedral_traj_deg)
    # Histogram
    (n, bins) = np.histogram(dihedral_traj_deg, bins=int(number_of_bins), range=(-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Total amount to be split over empty bins only
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    # Adds the bin_addition amount into all zero-count bins
    if zero_bins != 0:
        bin_addition = total_bin_addition/float(zero_bins)
        # Adds the bin_addition amount into all zero-count bins
        for j in xrange(len(n)):
            if n[j] == 0.0:
                n[j] = bin_addition
    #normalise
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/PSI/psi_hist_%d.dat' % (system_list[input_system],(i+1)), data, fmt=['%d', '%.30f'])

index = np.linspace(1,(psi_list[1].shape[1]),num=(psi_list[1].shape[1]))
index = index.reshape(len(index),1)
psi_atoms_array = np.array(psi_list[0])
psi_atoms_array =  np.hstack((index , psi_atoms_array))

np.savetxt('1_DIHEDRALS/%s/OUTPUT/PSI/1_psi_indices_list.dat' % system_list[input_system] , psi_atoms_array , fmt='%d')
print ("Psi distributions saved to: '1_DIHEDRALS/%s/OUTPUT/PSI'" % system_list[input_system])


print ("Computing phi distributions")
for i in xrange(phi_list[1].shape[1]):
    dihedral_traj = phi_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    # Save torsion to output
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/PHI/raw_data_phi_%d.dat' % (system_list[input_system],(i+1)),dihedral_traj_deg)
    # Histogram
    (n, bins) = np.histogram(dihedral_traj_deg, bins=int(number_of_bins), range = (-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    if zero_bins != 0:
        bin_addition = total_bin_addition/float(zero_bins)
        # Adds the bin_addition amount into all zero-count bins
        for j in xrange(len(n)):
            if n[j] == 0.0:
                n[j] = bin_addition
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/PHI/phi_hist_%d.dat' % (system_list[input_system],(i+1)), data, fmt=['%d', '%.30f'])

index = np.linspace(1,(phi_list[1].shape[1]),num=(phi_list[1].shape[1]))
index = index.reshape(len(index),1)
phi_atoms_array = np.array(phi_list[0])
phi_atoms_array =  np.hstack((index , phi_atoms_array))
np.savetxt('1_DIHEDRALS/%s/OUTPUT/PHI/1_phi_indices_list.dat' % system_list[input_system] , phi_atoms_array , fmt='%d')
print ("Phi distributions saved to: '1_DIHEDRALS/%s/OUTPUT/PHI'" % system_list[input_system])

print ("Computing chi1 distributions")
for i in xrange(chi1_list[1].shape[1]):
    dihedral_traj = chi1_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    # Save torsion to output
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI1/raw_data_chi1_%d.dat' % (system_list[input_system],(i+1)),dihedral_traj_deg)
    # Histogram
    (n, bins) = np.histogram(dihedral_traj_deg, bins=int(number_of_bins), range = (-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Total amount to be split over empty bins only
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    # Adds the bin_addition amount into all zero-count bins
    if zero_bins != 0:
        bin_addition = total_bin_addition/float(zero_bins)
        # Adds the bin_addition amount into all zero-count bins
        for j in xrange(len(n)):
            if n[j] == 0.0:
                n[j] = bin_addition
    #normalise
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI1/chi1_hist_%d.dat' % (system_list[input_system],(i+1)), data, fmt=['%d', '%.30f'])
    
index = np.linspace(1,(chi1_list[1].shape[1]),num=(chi1_list[1].shape[1]))
index = index.reshape(len(index),1)
chi1_atoms_array = np.array(chi1_list[0])
chi1_atoms_array =  np.hstack((index , chi1_atoms_array))
np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI1/1_chi1_indices_list.dat' % system_list[input_system] , chi1_atoms_array , fmt='%d')
print ("Chi1 distributions saved to: '1_DIHEDRALS/%s/OUTPUT/CHI1'" % system_list[input_system])


print ("Computing chi2 distributions")
for i in xrange(chi2_list[1].shape[1]):
    dihedral_traj = chi2_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    # Save torsion to output
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI2/raw_data_chi2_%d.dat' % (system_list[input_system],(i+1)),dihedral_traj_deg)
    # Histogram
    (n, bins) = np.histogram(dihedral_traj_deg, bins=int(number_of_bins), range = (-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Total amount to be split over empty bins only
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    # Adds the bin_addition amount into all zero-count bins
    if zero_bins != 0:
        bin_addition = total_bin_addition/float(zero_bins)
        # Adds the bin_addition amount into all zero-count bins
        for j in xrange(len(n)):
            if n[j] == 0.0:
                n[j] = bin_addition
    #normalise
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI2/chi2_hist_%d.dat' % (system_list[input_system],(i+1)), data, fmt=['%d', '%.30f'])
    
index = np.linspace(1,(chi2_list[1].shape[1]),num=(chi2_list[1].shape[1]))
index = index.reshape(len(index),1)
chi2_atoms_array = np.array(chi2_list[0])
chi2_atoms_array =  np.hstack((index , chi2_atoms_array))
np.savetxt('1_DIHEDRALS/%s/OUTPUT/CHI2/1_chi2_indices_list.dat'% system_list[input_system] , chi2_atoms_array , fmt='%d')
print ("Chi2 distributions saved to: '1_DIHEDRALS/%s/OUTPUT/CHI2'" % system_list[input_system])

