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
from scipy.fftpack import fft
from itertools import combinations
import matplotlib.pyplot as plt

##############################################################################
# Code
##############################################################################

data_path = '/Users/shenglanqiao/Documents/GitHub/waterMD/data'
# data_path='/home/shenglan/MD_simulations/water_box/cubic_2nm'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj

# time_step = 1 # in ps
time_step=traj.timestep # in ps

#atom.index for all Oxygen of the water molecules, get pariwise distances
water_inds = traj.topology.select_atom_indices(selection='water')
water_pairs = np.array(list(combinations(sorted(water_inds),2)))
water_dist = md.compute_distances(traj,water_pairs) # unit in nm

#examine statisitics for every frame
mean_dist = np.mean(water_dist,axis=1)
#sd_dist = np.std(water_dist,axis=1)
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

n_bin = 100

# pressure = 10**5.0 # Pa
# temp = 300 # K
# k_b= 1.3806e-23 # SI unit
# rho_id = pressure/(k_b*temp)*1e-27 # number of atomers per nm^3

# print len(water_dist[0])
# print traj.n_residues
# print traj.unitcell_volumes
# rho_id = traj.n_residues/np.mean(traj.unitcell_volumes)
rho_id = 1.0


# look at 0 to 1.0 nm

g2s = []
bin_edges = np.linspace(0.0,1.0,n_bin)
for this_frame in water_dist:
    
    n_ideal = 4.* np.pi*rho_id/3.*(bin_edges[1:]**3.0-bin_edges[:-1]**3.0)
    
    n_hist, _ = np.histogram(this_frame, bins = bin_edges,density=True)
    this_g2 = n_hist/n_ideal
    g2s.append(this_g2)

g2s = np.array(g2s)
bin_edges = 0.5*(bin_edges[1:]+bin_edges[:-1])
ave_g2 = np.mean(g2s,axis=0)
g2_err = np.std(g2s,axis=0)/(ave_g2[-1])
ave_g2 = ave_g2/(ave_g2[-1])

r,grs = md.compute_rdf(traj, pairs = water_pairs,bin_width = 0.01)

fig4 = plt.figure()
plt.plot(r,grs,color='Green')
plt.errorbar(bin_edges,ave_g2,yerr=g2_err)
plt.xlabel('pairwise distances (nm)')
plt.ylabel('normalized g(r)')
plt.xlim(0,1.0)
plt.ylim(0,3.0)
plt.title('pair distribution function for simulated water')
#plt.axvline(x = 1.0,ymin= 0,ymax = 1000,linestyle='--')
fig4.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/norm_g_cubic_2nm.png')
plt.close(fig4)

# cut_off=np.where(mean_bin_edges >=1.0)[0][0]
# intensity = np.absolute(fft(ave_g2[:cut_off]))**2.0
# fig5 = plt.figure()
# plt.plot(intensity)
# fig5.savefig('/home/shenglan/GitHub/waterMD/output/intensity_cubic_2nm.png')
# plt.close(fig5)
