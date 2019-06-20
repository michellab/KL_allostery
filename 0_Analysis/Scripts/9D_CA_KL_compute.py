# get_ipython().magic(u'pylab inline')
# -*- coding: utf-8 -*-

# Make soft links to folders for output
import scipy as sp
from scipy import stats
import numpy as np
from sys import argv

script, reference, target, distances = argv
distances = int(distances)
KLs = []
for i in range(0,distances):
    #Load two text files as np arrays
    P = np.loadtxt("%s/distance_%s_distribution.dat" % (reference,i))
    Q = np.loadtxt("%s/distance_%s_distribution.dat" % (target,i))
    KL = sp.stats.entropy(pk=P, qk=Q, base=None)
    KLs.append(KL[1])

KLs_arr = np.array(KLs)

np.savetxt("ALL_KL.dat", KLs_arr, fmt='%.10f')


