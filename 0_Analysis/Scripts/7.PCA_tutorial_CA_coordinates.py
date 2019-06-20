# coding: utf-8
# get_ipython().magic(u'pylab inline')
import pyemma
import pyemma.coordinates as coor
from pyemma.coordinates import pca
import mdtraj as md
from pyemma.coordinates import load
from pyemma.coordinates import source
import numpy as np
import os
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#

# edit system_list with folder names of different systems.
system_list = ["0_system_A","1_system_B","2_system_C","3_system_D"]
md_data = "0_TRAJECTORIES"

# Select only CA and also not all residues - C-terminal region is large. 
# This is approx focus on residues around active site, around allosteric site, and between
atom_sel = "resid 20 to 180 and name CA"

# Select trajectory:
#traj_filename = "longtraj_aligned_PCA.dcd"
traj_filename = "traj0001.dcd"

number_of_bins = 60
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#

top = "%s/%s/topology.parm7" % (md_data,system_list[0])

topology = md.load_prmtop(top)
print (topology)

# all CA indices:
CA_topo = topology.select("name CA")

# selected CA indices
atom_selection = topology.select(atom_sel)


# Make a list with all file locations of trajectory data
all_files_list = []

for i in range(0,len(system_list)):
    filenames = "%s/%s/%s" % (md_data,system_list[i],traj_filename)
    all_files_list.append(filenames)

print ("")
print ("Input trajectory locations: ")
print (all_files_list)

filename = "2_PCA/CA_COOR_OUTPUT" 
if not os.path.exists(filename):
    cmd = "mkdir -p %s" % filename
    os.system(cmd)

print ("")
feat = coor.featurizer(top)

feat.add_selection(atom_selection)
inp = coor.source([all_files_list[0],all_files_list[1],all_files_list[2],all_files_list[3]], features=feat)

print ("")
print ("Data input to PCA as follows: ")
print ('trajectory length = '),(inp.trajectory_length(0))
print ('number of dimensions = '),(inp.dimension())
print ('number of trajectories ='),(inp.number_of_trajectories())
print ('total number of frames = '),(inp.n_frames_total())

print ("")
print ("Running PCA")
pca_obj = coor.pca(inp, dim=10)

print ("")
print ("Generating per system output")
system_A_pca_out = pca_obj.get_output()[0]
system_B_pca_out = pca_obj.get_output()[1]
system_C_pca_out = pca_obj.get_output()[2]
system_D_pca_out = pca_obj.get_output()[3]

# for pc1, min and max values are these frames:
A_pc1_max = (np.argmax(system_A_pca_out[:,0]))
A_pc1_min = (np.argmin(system_A_pca_out[:,0]))

B_pc1_max = (np.argmax(system_B_pca_out[:,0]))
B_pc1_min = (np.argmin(system_B_pca_out[:,0]))

C_pc1_max = (np.argmax(system_C_pca_out[:,0]))
C_pc1_min = (np.argmin(system_C_pca_out[:,0]))

D_pc1_max = (np.argmax(system_D_pca_out[:,0]))
D_pc1_min = (np.argmin(system_D_pca_out[:,0]))

system_A_traj = md.load(all_files_list[0],top=topology)
system_B_traj = md.load(all_files_list[1],top=topology)
system_C_traj = md.load(all_files_list[2],top=topology)
system_D_traj = md.load(all_files_list[3],top=topology)

print ("")
print ("")
print ("Saving frames corresponding to the minimum and maximum values of PCA for each system to folder %s" % filename)

system_A_traj[A_pc1_max].save_pdb("%s/A_PC1_max_frame%s.pdb" % (filename,A_pc1_max))
system_A_traj[A_pc1_min].save_pdb("%s/A_PC1_min_frame%s.pdb" % (filename,A_pc1_min))

system_B_traj[B_pc1_max].save_pdb("%s/B_PC1_max_frame%s.pdb" % (filename,B_pc1_max))
system_B_traj[B_pc1_min].save_pdb("%s/B_PC1_min_frame%s.pdb" % (filename,B_pc1_min))

system_C_traj[C_pc1_max].save_pdb("%s/C_PC1_max_frame%s.pdb" % (filename,C_pc1_max))
system_C_traj[C_pc1_min].save_pdb("%s/C_PC1_min_frame%s.pdb" % (filename,C_pc1_min))

system_D_traj[D_pc1_max].save_pdb("%s/D_PC1_max_frame%s.pdb" % (filename,D_pc1_max))
system_D_traj[D_pc1_min].save_pdb("%s/D_PC1_min_frame%s.pdb" % (filename,D_pc1_min))

# pc1 - residues 1-264 xyz
pc1_ca_contributions = []
for i in pca_obj.feature_PC_correlation[:,0]:
    pc1_ca_contributions.append(np.absolute(i))

pc1_ca_X = np.array(pc1_ca_contributions[0:len(pc1_ca_contributions):3])
pc1_ca_Y = np.array(pc1_ca_contributions[1:len(pc1_ca_contributions):3])
pc1_ca_Z = np.array(pc1_ca_contributions[2:len(pc1_ca_contributions):3])

# Since we don't include all residues in the PCA, we need to assign the per atom contribution to the 
# correct residue in the topology. To do this, rearrange the data by adding zeros for residues not included in the PCA

def rearrange_data(data):
    residue_and_per_atom_contribution = np.vstack((atom_selection,data)).T
    new_atoms = []
    for i in CA_topo:
        if i in atom_selection: 
            new_atoms.append(i)
        else:
            new_atoms.append("--")
    new_data = []
    j = 0
    for i in range (0,len(CA_topo)):
        if CA_topo[i] in residue_and_per_atom_contribution[:,0]:
            new_data.append(residue_and_per_atom_contribution[j][1])
            j = j + 1 
        else:
            new_data.append(0.0)
    return new_data

# Sum x y and z components 
PC1_all_ca_contrib = np.vstack((pc1_ca_X,pc1_ca_Y,pc1_ca_Z)).T
summed_cont_PC1 = PC1_all_ca_contrib.sum(axis=1)
summed_cont_PC1 = rearrange_data(summed_cont_PC1)
np.savetxt("%s/PC1_atom_contribution.dat" % filename , summed_cont_PC1)

plt.plot(summed_cont_PC1)
plt.xlabel("Residue number")
plt.ylabel("Contribution to PC1")
plt.savefig("%s/Per_res_PC1.png" % filename)
plt.clf()

# PC2
pc2_ca_contributions = []
for i in pca_obj.feature_PC_correlation[:,1]:
    pc2_ca_contributions.append(np.absolute(i))


pc2_ca_X = np.array(pc2_ca_contributions[0:len(pc2_ca_contributions):3])
pc2_ca_Y = np.array(pc2_ca_contributions[1:len(pc2_ca_contributions):3])
pc2_ca_Z = np.array(pc2_ca_contributions[2:len(pc2_ca_contributions):3])

PC2_all_ca_contrib = np.vstack((pc2_ca_X,pc2_ca_Y,pc2_ca_Z)).T

summed_cont_PC2 = PC2_all_ca_contrib.sum(axis=1)
summed_cont_PC2 = rearrange_data(summed_cont_PC2)
np.savetxt("%s/PC2_atom_contribution.dat" % filename , summed_cont_PC2)

plt.plot(summed_cont_PC2)
plt.xlabel("Residue number")
plt.ylabel("Contribution to PC2")
plt.savefig("%s/Per_res_PC2.png" % filename)
plt.clf()

#print (pca_obj.eigenvalues.shape)
#print (pca_obj.eigenvectors)

# 264 Ca atoms x (atom selection may vary)- each has 3 dimensions (x,y,z)
plt.subplot2grid((2,1),(0,0))
plt.plot(pca_obj.eigenvectors[0])
plt.subplot2grid((2,1),(1,0))
plt.plot(pca_obj.eigenvectors[1])
plt.xlabel("")
plt.clf()

#plt.plot(pca_obj.eigenvectors[0])
first_eigenvector = []
for i in pca_obj.eigenvectors[0]:
    first_eigenvector.append(i)

#print ("")    
#print (len(first_eigenvector))    
#print ("Min value")
#print ("index", first_eigenvector.index(min(first_eigenvector)))
#print (min(first_eigenvector))
#print ("Max value")
#print ("index" , first_eigenvector.index(max(first_eigenvector)))
#print (max(first_eigenvector))

#print ("First 10 PCs with 40k ss each for first sim" ),( system_A_pca_out.shape)
#print ("First PC for first input sim"),( system_A_pca_out[:,0])
#print ("Second PC for first input sim"),( system_A_pca_out[:,1])


# plot of pc1 and pc2
plt.plot(pca_obj.eigenvectors[0],pca_obj.eigenvectors[1], marker='.', lw=0)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.xlim(-0.3,0.3)
plt.ylim(-0.3,0.3)
plt.savefig('%s/1st_2nd_eigenvec_2d_plot.png' % filename)
plt.clf()

print ("")
print ("Calculating variance for calculated PCs")
percentage_variance = []
total = sum(pca_obj.eigenvalues)

for i in xrange(0,10):
    first = pca_obj.eigenvalues[i]
    x = first/total * 100
    percentage_variance.append(x)
    print ("Percentage variance PC%s: " % (i+1)) , (x)

print ("Sum of variance for first 10 PCs: "),(sum(percentage_variance))


index = np.linspace(1 , len(percentage_variance), num=len(percentage_variance))
plt.plot(index,percentage_variance)
plt.ylabel("Percentage variance")
plt.xlabel("Principal component")
plt.savefig("%s/variance.png" % filename)
plt.clf()

# Plots of PC1 and PC2 vs. snapshot
plt.subplot2grid((2,1),(0,0))
plt.plot(system_A_pca_out[:,0])
plt.title("PC1 for each system")
plt.ylabel('PC1 system A')
plt.ylim((-5 , 5))
plt.subplot2grid((2,1),(1,0))
plt.plot(system_B_pca_out[:,0])
plt.ylabel('PC1 system B')
plt.ylim((-5 , 5))
plt.savefig('%s/PC1_sysA_and_B_vs_ss.png' % filename)
plt.clf()

# Plots of PC1 and PC2 vs. snapshot
plt.subplot2grid((2,1),(0,0))
plt.plot(system_C_pca_out[:,0])
plt.title("PC1 for each system")
plt.ylabel('PC1 system C')
plt.ylim((-5 , 5))
plt.subplot2grid((2,1),(1,0))
plt.plot(system_D_pca_out[:,0])
plt.ylabel('PC1 system D')
plt.ylim((-5 , 5))
plt.savefig('%s/PC1_sysC_and_D_vs_ss.png' % filename)
plt.clf()

# Plots of PC1 and PC2 vs. snapshot
plt.subplot2grid((2,1),(0,0))
plt.plot(system_A_pca_out[:,1])
plt.title("PC2 for each system")
plt.ylabel('PC2 system A')
plt.ylim((-5 , 5))
plt.subplot2grid((2,1),(1,0))
plt.plot(system_B_pca_out[:,1])
plt.ylabel('PC2 system B')
plt.ylim((-5 , 5))
plt.savefig('%s/PC2_sysA_and_B_vs_ss.png' % filename)
plt.clf()

# Plots of PC1 and PC2 vs. snapshot
plt.subplot2grid((2,1),(0,0))
plt.plot(system_C_pca_out[:,1])
plt.title("PC2 for each system")
plt.ylabel('PC2 system C')
plt.ylim((-5 , 5))
plt.subplot2grid((2,1),(1,0))
plt.plot(system_D_pca_out[:,1])
plt.ylabel('PC2 system D')
plt.ylim((-5 , 5))
plt.savefig('%s/PC2_sysC_and_D_vs_ss.png' % filename)
plt.clf()


#print (np.vstack((system_A_pca_out[:,0],system_B_pca_out[:,0],system_C_pca_out[:,0],system_D_pca_out[:,0])).shape)

pc1_all = np.vstack((system_A_pca_out[:,0],system_B_pca_out[:,0],system_C_pca_out[:,0],system_D_pca_out[:,0])).reshape(inp.n_frames_total(),)
pc2_all = np.vstack((system_A_pca_out[:,1],system_B_pca_out[:,1],system_C_pca_out[:,1],system_D_pca_out[:,1])).reshape(inp.n_frames_total(),)

z_,x_,y_ = np.histogram2d(pc1_all,pc2_all, bins=number_of_bins, range=[[-5, 5], [-5, 5]])
plot_surface = [x_[0], x_[-1], y_[0], y_[-1]]
plt.contourf(z_.T, 100, extent=plot_surface)
plt.xlim(-4,3)
plt.ylim(-3,4)
plt.savefig('%s/2dhistogram.png' % filename)
plt.clf()

z_,x_,y_ = np.histogram2d(pc1_all,pc2_all, bins=number_of_bins, range=[[-5, 5], [-5, 5]])
plot_surface = [x_[0], x_[-1], y_[0], y_[-1]]
plt.contourf(z_.T, 100, extent=plot_surface)

plt.plot(system_A_pca_out[:,0][0::500],system_A_pca_out[:,1][0::500], marker='_', color='r')
plt.plot(system_B_pca_out[:,0][0::500],system_B_pca_out[:,1][0::500], marker='_', color='w')
plt.plot(system_C_pca_out[:,0][0::500],system_C_pca_out[:,1][0::500], marker='_', color='y')
plt.plot(system_D_pca_out[:,0][0::500],system_D_pca_out[:,1][0::500], marker='_', color='g')
plt.xlim(-4,3)
plt.ylim(-3,4)
plt.savefig('%s/2dhistogram_TRAJ.png' % filename)
plt.clf()

print ("")
print ("Calculating bin ranges to histogram PC1 and PC2.")
min_max_list_pc1 = []
min_max_list_pc2 = []

min_max_list_pc1.append(min(system_A_pca_out[:,0]))
min_max_list_pc1.append(max(system_A_pca_out[:,0]))
min_max_list_pc1.append(min(system_B_pca_out[:,0]))
min_max_list_pc1.append(max(system_B_pca_out[:,0]))

min_max_list_pc1.append(min(system_C_pca_out[:,0]))
min_max_list_pc1.append(max(system_C_pca_out[:,0]))
min_max_list_pc1.append(min(system_D_pca_out[:,0]))
min_max_list_pc1.append(max(system_D_pca_out[:,0]))

min_max_list_pc2.append(min(system_A_pca_out[:,1]))
min_max_list_pc2.append(max(system_A_pca_out[:,1]))
min_max_list_pc2.append(min(system_B_pca_out[:,1]))
min_max_list_pc2.append(max(system_B_pca_out[:,1]))

min_max_list_pc2.append(min(system_C_pca_out[:,1]))
min_max_list_pc2.append(max(system_D_pca_out[:,1]))
min_max_list_pc2.append(min(system_D_pca_out[:,1]))
min_max_list_pc2.append(max(system_D_pca_out[:,1]))

bin_max_pc1 = int((max(min_max_list_pc1))+1)
bin_min_pc1 = int((min(min_max_list_pc1))-1)
bin_max_pc2 = int((max(min_max_list_pc2))+1)
bin_min_pc2 = int((min(min_max_list_pc2))-1)

print ("Bin range will be from %s to %s for PC1" % (bin_min_pc1,bin_max_pc1))
print ("Bin range will be from %s to %s for PC2" % (bin_min_pc2,bin_max_pc2))

#Can then output distributions of PC1 and PC2 for kL/MI etc. 

sysA = plt.hist(system_A_pca_out[:,0], bins=number_of_bins, range=(bin_min_pc1,bin_max_pc1), normed=True, histtype='step', color='r', label='A')
sysB = plt.hist(system_B_pca_out[:,0], bins=number_of_bins, range=(bin_min_pc1,bin_max_pc1), normed=True, histtype='step', color='b', label='B')
sysB = plt.hist(system_C_pca_out[:,0], bins=number_of_bins, range=(bin_min_pc1,bin_max_pc1), normed=True, histtype='step', color='y', label='C')
sysB = plt.hist(system_D_pca_out[:,0], bins=number_of_bins, range=(bin_min_pc1,bin_max_pc1), normed=True, histtype='step', color='g', label='D')

plt.ylim(0,3)
plt.xlim(bin_min_pc1,bin_max_pc1)
#pylab.legend(loc='upper left')
plt.savefig('%s/PC1_histogram.png' % filename)
plt.clf()

sysA = plt.hist(system_A_pca_out[:,1], bins=number_of_bins, range=(bin_min_pc2,bin_max_pc2), normed=True, histtype='step', color='r', label='A')
sysB = plt.hist(system_B_pca_out[:,1], bins=number_of_bins, range=(bin_min_pc2,bin_max_pc2), normed=True, histtype='step', color='b', label='B')
sysB = plt.hist(system_C_pca_out[:,1], bins=number_of_bins, range=(bin_min_pc2,bin_max_pc2), normed=True, histtype='step', color='y', label='C')
sysB = plt.hist(system_D_pca_out[:,1], bins=number_of_bins, range=(bin_min_pc2,bin_max_pc2), normed=True, histtype='step', color='g', label='D')

plt.ylim(0,3)
plt.xlim(bin_min_pc2,bin_max_pc2)
#pylab.legend(loc='upper left')
plt.savefig('%s/PC2_histogram.png' % filename)
plt.clf()

list_of_sys = [system_A_pca_out , system_B_pca_out , system_C_pca_out , system_D_pca_out]    
list_of_names = ["system_A" , "system_B" , "system_C" , "system_D"]

for p in range(0, len(list_of_sys)):
    np.savetxt("%s/PC1_raw_data_%s.dat" % (filename,list_of_names[p]) , list_of_sys[p][:,0])
    (n, bins) = np.histogram(list_of_sys[p][:,0], bins = number_of_bins, range=(bin_min_pc1,bin_max_pc1), normed=True)
    n = n / (sum(n))
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
    total_bin_addition = 0.000001
    all_bins = len(bincentre)
    # To count the number of populated and non populated bins, to allow dividion of the total bin addition
    non_zero = np.count_nonzero(n)
    print ("Number of populated bins:"), (non_zero)
    zero_bins = all_bins - non_zero
    print ("Number of zero bins:"), (zero_bins)
    bin_addition = total_bin_addition/float(zero_bins)
    print ("Amount added to empty bins:"), (bin_addition)
    for i in xrange(len(n)):
        if n[i]==0.0:
            n[i] = bin_addition
    data = np.vstack((index, n)).T
    np.savetxt("%s/PC1_hist_%s.dat" % (filename,list_of_names[p]) , data, fmt=['%d', '%.20f'])

for p in range(0, len(list_of_sys)):
    np.savetxt("%s/PC2_raw_data_%s.dat" % (filename,list_of_names[p]) , list_of_sys[p][:,1])
    (n, bins) = np.histogram(list_of_sys[p][:,1], bins = number_of_bins, range=(bin_min_pc2,bin_max_pc2), normed=True)
    n = n / (sum(n))
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
    total_bin_addition = 0.000001
    all_bins = len(bincentre)
    # To count the number of populated and non populated bins, to allow dividion of the total bin addition
    non_zero = np.count_nonzero(n)
    print ("Number of populated bins:"), (non_zero)
    zero_bins = all_bins - non_zero
    print ("Number of zero bins:"), (zero_bins)
    bin_addition = total_bin_addition/float(zero_bins)
    print ("Amount added to empty bins: "), (bin_addition)
    for i in xrange(len(n)):
        if n[i]==0.0:
            n[i] = bin_addition
    data = np.vstack((index, n)).T
    np.savetxt("%s/PC2_hist_%s.dat" % (filename,list_of_names[p]) , data, fmt=['%d', '%.20f'])

print ("")
print ("_________________________________________")
print ("All output saved to %s." % filename)
print ("_________________________________________")
