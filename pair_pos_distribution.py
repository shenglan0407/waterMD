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

data_path = '/Users/shenglanqiao/zauber/MD_simulations/water_box/cubic_1.2nm'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj

time_step=traj.timestep # in ps

#atom.index for all Oxygen of the water molecules, get pariwise distances
water_inds = traj.topology.select_atom_indices(selection='water')
water_pairs = np.array(list(combinations(sorted(water_inds),2)))
water_dist = md.compute_distances(traj,water_pairs)

#examine statisitics for every frame
mean_dist = np.mean(water_dist,axis=1)
sd_dist = np.std(water_dist,axis=1)
traj_time = np.array(range(len(mean_dist)))*time_step

#visualize
fig1, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(traj_time, mean_dist)
axarr[0].set_ylim([1.2325,1.2335])
axarr[0].set_ylabel('mean pairwise distances (nm)')
axarr[1].scatter(traj_time, sd_dist)
axarr[1].set_ylim([0.353,0.354])
axarr[1].set_ylabel('std pairwise distances (nm)')
axarr[1].set_xlabel('time(ps)')

fig1.savefig('/Users/shenglanqiao/Documents/Github/waterMD/output/pairwise_timeseries.png')
plt.close(fig1)


fig2 = plt.figure()
for this_frame in water_dist:
    plt.hist(this_frame,alpha=0.05)
plt.title('Distribution of pairwise water distances superimposed by frame')
plt.xlabel('Pairwise distances (nm)')
fig2.savefig('/Users/shenglanqiao/Documents/Github/waterMD/output/pairwise_histogram.png')
plt.close(fig2)