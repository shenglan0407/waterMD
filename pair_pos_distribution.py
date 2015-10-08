##############################################################################
# Copyright 2015 Stanford University and the Authors
#
# Author: Shenglan Qiao
# 
#############################################################################


##############################################################################
# Imports
##############################################################################

import mdtraj as md

import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

##############################################################################
# Code
##############################################################################

data_path = '/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj

#atom.index for all Oxygen of the water molecules, get pariwise distances
water_inds = traj.topology.select_atom_indices(selection='water')
water_pairs = np.array(list(combinations(sorted(water_inds),2)))
water_dist = md.compute_distances(traj,water_pairs)

print water_dist[0].shape

plt.hist(water_dist[0])
plt.show()