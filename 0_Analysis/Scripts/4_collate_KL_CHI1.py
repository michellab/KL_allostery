# -*- coding: utf-8 -*-
import numpy as np
from sys import argv

script, var1 = argv

KL_list = []
for i in xrange(1, int(var1)+1):
    KL = np.loadtxt("chi1_out_%s.dat" % i)
    KL_list.append(KL[1])

KL_list_6dp = [ '%.6f' % elem for elem in KL_list ]
#print KL_list_6dp

np.savetxt('ALL_KL_chi1.dat', KL_list_6dp, fmt="%s")
