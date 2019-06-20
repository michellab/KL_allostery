# coding: utf-8
import numpy as np
import sys
from sys import argv

script, filepathA, filepathB = argv

print ("Loading min max values from %s and %s" % (filepathA,filepathB))
_0_system_A_min_max = np.loadtxt("%s/min_max_rawdata.dat" % filepathA)
_1_system_B_min_max = np.loadtxt("%s/min_max_rawdata.dat" % filepathB)

min_col_sys_A = _0_system_A_min_max[:,0]
max_col_sys_A = _0_system_A_min_max[:,1]

min_col_sys_B = _1_system_B_min_max[:,0]
max_col_sys_B = _1_system_B_min_max[:,1]

MIN_cols = np.vstack((min_col_sys_A,min_col_sys_B)).T
MAX_cols = np.vstack((max_col_sys_A,max_col_sys_B)).T

print ("")
print ("Calculating overall min and max for each CA-CA distance... ")
global_mins = []
for i in MIN_cols:
    if i[0] < i[1]:
        global_mins.append(i[0])
    else: 
        global_mins.append(i[1])
        
global_maxs = []
for i in MAX_cols:
    if i[0] > i[1]:
        global_maxs.append(i[0])
    else: 
        global_maxs.append(i[1])

global_mins_int = []
for i in range(0,len(global_mins)):
    global_mins_int.append(int(global_mins[i]))
global_mins_arr = np.array(global_mins_int).clip(min=0)

global_maxs_int = []
for i in range(0,len(global_maxs)):
    global_maxs_int.append(int(global_maxs[i]))
global_maxs_arr = np.array(global_maxs_int).clip(min=0)

col1 = (global_mins_arr - 2).clip(min=0)
col2 = (global_maxs_arr + 2).clip(min=0)
min_max_arr_margin_int = np.vstack((col1,col2)).T
print ("")
print ("Saving overall min max to 4_CA_DISTANCES/global_min_max_array.dat")
np.savetxt("4_CA_DISTANCES/global_min_max_array.dat" ,min_max_arr_margin_int,fmt='%d')

