#get_ipython().magic(u'pylab inline')
import mdtraj as md
import numpy as np
from sys import argv

def entropy(data):
    data_norm = data / float(np.sum(data))
    data_norm = data_norm[np.nonzero(data_norm)]
    H = -sum(data_norm * np.log(data_norm))  
    return H

def mutual_information(x , y , x_bins , y_bins):
    
    bins = (x_bins , y_bins)
    
    joint_dist = np.histogram2d(x,y,bins=bins)[0]
    x_dist = np.histogram(x,x_bins)[0]
    y_dist = np.histogram(y,y_bins)[0]
    
    H_joint = entropy(joint_dist)
    H_x = entropy(x_dist)
    H_y = entropy(y_dist)
    
    MutInfo = H_x + H_y - H_joint
    return MutInfo

# load some data: 
# Example calculating one comparison, but scripts will be written to loop through all possible combinations. 
# I.e. calculate MI of all torsions wrt all torsions (find correlated motions)
# Or calculate MI of all torsions wrt pairwise interaction 
# energies of ligand (find motions that correlate to interactions of ligand)

data_1 = md.loadtxt(variableA)
data_2 = md.loadtxt(variableB)

# compute MI
print mutual_information(data_1, data_2, binsA, binsB)
