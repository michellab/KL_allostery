# coding: utf-8
import numpy as np
from sys import argv
#get_ipython().magic(u'pylab inline')

script, output_folder = argv

psi_KL = np.loadtxt("%s/PSI/ALL_KL_psi.dat" % output_folder)
phi_KL = np.loadtxt("%s/PHI/ALL_KL_phi.dat" % output_folder)
chi1_KL = np.loadtxt("%s/CHI1/ALL_KL_chi1.dat" % output_folder)
chi2_KL = np.loadtxt("%s/CHI2/ALL_KL_chi2.dat" % output_folder)

ca_indexes = np.loadtxt("DIHEDRALS_CA_angles_indexes.dat",dtype=str)

#print ca_indexes
#print len(ca_indexes[:,0])

index_list_ca = np.linspace(1,len(ca_indexes),num=len(ca_indexes))
index_psi = np.linspace(0,len(psi_KL)-1,num=len(psi_KL))
index_phi = np.linspace(0,len(phi_KL)-1,num=len(phi_KL))
index_chi1 = np.linspace(0,len(chi1_KL)-1,num=len(chi1_KL))
index_chi2 = np.linspace(0,len(chi2_KL)-1,num=len(chi2_KL))

#print psi_KL.shape
#print index_psi.shape
psi_KL = np.vstack((index_psi,psi_KL)).T
phi_KL = np.vstack((index_phi,phi_KL)).T
chi1_KL = np.vstack((index_chi1,chi1_KL)).T
chi2_KL = np.vstack((index_chi2,chi2_KL)).T

gg_psi = []
for i in (ca_indexes[:,1]):
    try:
        gg_psi.append(float(i))
    except ValueError:
        gg_psi.append("-")
#print len(gg_psi)

gg_phi = []
for i in (ca_indexes[:,2]):
    try:
        gg_phi.append(float(i))
    except ValueError:
        gg_phi.append("-")
#print len(gg_phi)

gg_chi1 = []
for i in (ca_indexes[:,3]):
    try:
        gg_chi1.append(float(i))
    except ValueError:
        gg_chi1.append("-")
#print len(gg_chi1)

gg_chi2 = []
for i in (ca_indexes[:,4]):
    try:
        gg_chi2.append(float(i))
    except ValueError:
        gg_chi2.append("-")
#print len(gg_chi2)

#  gg_psi  gg_phi  gg_chi1  gg_chi2 are a list - one entry for each CA with either angle number 
#  or "-" if no angle for that CA

# TORSIONAL ANGLE NUMBER IS psi_KL[:,0]
# KL VALUE IS psi_KL[:,1]

#print (gg_psi)
#print (psi_KL)
for i in xrange(len(gg_psi)):
    for j in psi_KL[:,0]:
        if gg_psi[int(i)] == psi_KL[int(j),0]:
            gg_psi[int(i)] = psi_KL[int(j),1]
for i in xrange(len(gg_phi)):
    for j in phi_KL[:,0]:
        if gg_phi[int(i)] == phi_KL[int(j),0]:
            gg_phi[int(i)] = phi_KL[int(j),1]
for i in xrange(len(gg_chi1)):
    for j in chi1_KL[:,0]:
        if gg_chi1[int(i)] == chi1_KL[int(j),0]:
            gg_chi1[int(i)] = chi1_KL[int(j),1]
for i in xrange(len(gg_chi2)):
    for j in chi2_KL[:,0]:
        if gg_chi2[int(i)] == chi2_KL[int(j),0]:
            gg_chi2[int(i)] = chi2_KL[int(j),1]

gg_psi_arr = np.array(gg_psi)
gg_phi_arr = np.array(gg_phi)
gg_chi1_arr = np.array(gg_chi1)
gg_chi2_arr = np.array(gg_chi2)

ee = np.vstack((index_list_ca , ca_indexes[:,0].astype(int) , gg_psi_arr, gg_phi_arr, gg_chi1_arr, gg_chi2_arr))

# angle numbers are chi1_KL[:,0] and corresponding KL is chi1_KL[:,1]
all_data_ = ee.T
#print all_data_

new_data = [] 
for d in all_data_:    
    small_list = []    
    for i in d:        
        try:            
            small_list.append(float(i))        
        except:            
            small_list.append(0.0)    
    new_data.append(small_list)

new_data_array = np.array(new_data)

# First column is residue. Second is CA atom number. Rest are psi phi chi chi KL
np.savetxt("KL_residue_caindex_psi_phi_chi_chi.dat", new_data_array , fmt=['%.d','%.d','%.6f','%.6f','%.6f','%.6f'])

psi_ALL_KL = new_data_array[:,2]
phi_ALL_KL = new_data_array[:,3]
chi1_ALL_KL = new_data_array[:,4]
chi2_ALL_KL = new_data_array[:,5]

np.savetxt("KL_psi.dat", psi_ALL_KL, fmt='%.4f')
np.savetxt("KL_phi.dat", phi_ALL_KL, fmt='%.4f')
np.savetxt("KL_chi1.dat", chi1_ALL_KL, fmt='%.4f')
np.savetxt("KL_chi2.dat", chi2_ALL_KL, fmt='%.4f')

backbone_KL = np.vstack((psi_ALL_KL,phi_ALL_KL)).T
backbone_KL2 = backbone_KL.sum(axis=1)
np.savetxt("KL_backbone.dat",backbone_KL2,fmt='%.6f')

sidechain_KL = np.vstack((chi1_ALL_KL,chi2_ALL_KL)).T
sidechain_KL2 = sidechain_KL.sum(axis=1)
np.savetxt("KL_sidechain.dat", sidechain_KL2, fmt='%.6f')
