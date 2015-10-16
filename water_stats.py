##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# The computation done here is based on the following papers:
# 
# [1] Salacuse et al. Finite-szize effects of molecular dynamics simulations: static structure 
# factor and compressibility. I. theoretical method
#
# [2] Salacuse et al. Finite-szize effects of molecular dynamics simulations: static structure 
# factor and compressibility. II. application to a model krypton fluid
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import mdtraj as md

import numpy as np
from itertools import combinations, product
import matplotlib.pyplot as plt
import random

##############################################################################
# Code
##############################################################################

class WaterStats:
    def __init__(self,traj):
        self.traj = traj
        self.n_waters = traj.topology.n_residues
        self.water_inds = traj.topology.select_atom_indices(selection='water')
        self.time_step=traj.timestep # in ps
        self.n_frames = self.traj.n_frames
        self.total_time = self.time_step*self.n_frames # in ps
        
        self.rho = np.mean(self.n_waters/traj.unitcell_volumes) # in nm^-3

        

    def single_frame_N_QR(self,Q,R,pair_dist,time_dependent = False):
        """See equations (1a) and (8) in ref [2]
        """
        if time_dependent:
            return 1.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
            for this_pair in pair_dist if this_pair <= R])
        else:
            return 2.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
            for this_pair in pair_dist if this_pair <= R]) \
            +self.n_waters
    

    def struct_factor(self,Q,R,dt):
        """
        See equations (6) and (4) in ref [2]
        """
        water_pairs = np.array(list(combinations(sorted(self.water_inds),2)))
        step = int(dt/self.time_step)
        frame_steps= [step+random.randint(-int(step/2),int(step/2)) \
        for ii in range(int(self.total_time/dt))]
        
        while sum(frame_steps)>self.n_frames:
            frame_steps = frame_steps[:-1]
            
        print ('the average time between each configuration is %f ps' % (np.mean(frame_steps)*self.time_step))
        
        N_QR = []
        current_frame = 0
        for this_step in frame_steps:
            current_frame += this_step
            water_dist = md.compute_distances(self.traj[current_frame],water_pairs)[0] # unit in nm
            print water_dist
            N_QR.append(self.single_frame_N_QR(Q,R,water_dist))
        
        # equation (4) in ref [2]
        err_Sn_QR = np.std(N_QR)/len(N_QR)
        
        # equations (4a-d) in ref [1]
        Sn_QR = np.mean(N_QR)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return Sn_QR, err_Sn_QR
            
    

    def scat_func(self,Q,R,t,dt):
        water_pairs = np.array(list(product(sorted(self.water_inds),repeat = 2)))
        
##############################################################################
# test
##############################################################################

data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj)
Sn_QR = test.struct_factor(2*np.pi/0.5,0.5,2)