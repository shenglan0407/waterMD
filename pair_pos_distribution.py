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

# data_path = '/Users/shenglanqiao/zauber/MD_simulations/water_box/cubic_5nm'
data_path='/Users/shenglanqiao/Documents/Github/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr-5nm.trr', top = data_path+'/water-sol-5nm.gro', frame = 1)
print ('here is some info about the trajectory we are looking at:')
print traj

time_step = 1 # in ps
# time_step=traj.timestep # in ps

#atom.index for all Oxygen of the water molecules, get pariwise distances
water_inds = traj.topology.select_atom_indices(selection='water')
water_pairs = np.array(list(combinations(sorted(water_inds),2)))
water_dist = md.compute_distances(traj,water_pairs) # unit in nm

#examine statisitics for every frame
mean_dist = np.mean(water_dist,axis=1)
sd_dist = np.std(water_dist,axis=1)
traj_time = np.array(range(len(mean_dist)))*time_step

#visualize
# fig1, axarr = plt.subplots(2, sharex=True)
# axarr[0].scatter(traj_time, mean_dist)
# axarr[0].set_xlim([0,100])
# axarr[0].set_ylabel('mean pairwise distances (nm)')
# axarr[1].scatter(traj_time, sd_dist)
# axarr[1].set_ylim([0.353,0.354])
# axarr[1].set_ylabel('std pairwise distances (nm)')
# axarr[1].set_xlabel('time(ps)')
# 
# fig1.savefig('/Users/shenglanqiao/Documents/Github/waterMD/output/pairwise_timeseries.png')
# plt.close(fig1)
# 
# 
# fig2 = plt.figure()
# for this_frame in water_dist:
#     plt.hist(this_frame,alpha=0.05)
# plt.title('Distribution of pairwise water distances superimposed by frame')
# plt.xlabel('Pairwise distances (nm)')
# fig2.savefig('/Users/shenglanqiao/Documents/Github/waterMD/output/pairwise_histogram.png')
# plt.close(fig2)
# 
# fig3 = plt.figure()
# plt.plot(traj_time, mean_dist,'bo')
# plt.errorbar(traj_time, mean_dist,yerr=sd_dist)
# plt.xlabel('time(ps)')
# plt.ylabel('mean pairwise distances (nm)')
# fig3.savefig('/Users/shenglanqiao/Documents/Github/waterMD/output/pairwise_mean_dist.png')
# plt.close(fig3)

test_frame = water_dist[0]
n_bin = 20
bin_edges=np.linspace(np.min(test_frame),np.max(test_frame),n_bin+1)
pressure = 10**5.0 # Pa
temp = 300 # K
k_b= 1.3806e-23 # SI unit
rho_id = pressure/(k_b*temp)*1e-27 # number of atomers per nm^3
n_ideal = 4.* np.pi*rho_id/3*(bin_edges[1:]**3.0-bin_edges[:-1]**3.0)
print len(n_ideal)
n_hist, _ = np.histogram(test_frame, bins = bin_edges,density=True)
n_pairs = len(test_frame)
plt.plot(bin_edges[:-1],n_hist/n_ideal)
plt.show()

print min(test_frame)