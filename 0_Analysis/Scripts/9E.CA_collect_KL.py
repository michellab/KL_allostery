# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

atom_number_pair_array = np.loadtxt("atom_pairs_ATOM_NUMBERS.dat")

KL_values = np.loadtxt("ALL_KL.dat")

KL_values = KL_values.reshape(len(KL_values),1)

atom_indices_KLs = np.concatenate((atom_number_pair_array,KL_values),axis=1)

print atom_indices_KLs
print max(KL_values)
print min(KL_values)

np.savetxt("KL_OUTPUT/0A_1B_atoms_indices_KLs.dat",atom_indices_KLs, fmt='%d  %d  %0.10f')
