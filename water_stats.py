##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# The computation done here is based on the following papers:
# 
# [1] Salacuse et al. Finite-size effects of molecular dynamics simulations: static structure 
# factor and compressibility. I. theoretical method
#
# [2] Salacuse et al. Finite-size effects of molecular dynamics simulations: static structure 
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
        """
        See equations (1a) and (8) in ref [2]
        """
        if Q == 0:
            # np.sin(Q*this_pair)/(Q*this_pair) = 1.0 in the limit of Q -> 0
            if time_dependent:
                return 1.0/self.n_waters*sum([ 1.0 \
                for this_pair in pair_dist if this_pair <= R])
            else:
                # I disgree with equation (1a) and think there must be a typo
                return 2.0/self.n_waters*sum([1.0 \
                for this_pair in pair_dist if this_pair <= R]) \
                + 1
                # + self.n_waters -> this is what equation (1a) in ref [2] thinks
        else:
            # See equations (1a) and (8) in ref [2]
            if time_dependent:
                return 1.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
                for this_pair in pair_dist if this_pair <= R])
            else:
                # I disgree with equation (1a) and think there must be a typo
                return 2.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
                for this_pair in pair_dist if this_pair <= R]) \
                + 1
                # + self.n_waters -> this is what equation (1a) in ref [2] thinks
    

    def struct_factor(self,Q,R,dt):
        """
        See equations (6) and (4) in ref [2]
        """
        water_pairs = np.array(list(combinations(sorted(self.water_inds),2)))
        step = int(dt/self.time_step)
        frame_steps= [step+random.randint(-int(step/2),int(step/2)) \
        for ii in range(int(self.total_time/dt))]
        
        while sum(frame_steps)>self.n_frames-1:
            frame_steps = frame_steps[:-1]
            
        # print ('the average time between each configuration is %f ps' % (np.mean(frame_steps)*self.time_step))
        
        N_QR = []
        current_frame = 0
        for this_step in frame_steps:
            current_frame += this_step
            water_dist = md.compute_distances(self.traj[current_frame],water_pairs)[0] # unit in nm
            N_QR.append(self.single_frame_N_QR(Q,R,water_dist))
        
        # equation (4) in ref [2]
        err_Sn_QR = np.std(N_QR)/len(N_QR)
        
        if Q == 0:
            # equation (7) in ref [2]
            Sn_QR = np.mean(N_QR)-4./3.*np.pi*self.rho*R**3.0
            
        else:
            # equations (4a-d) in ref [1]
            Sn_QR = np.mean(N_QR)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return Sn_QR, err_Sn_QR
            
    

    def scat_func(self,Q,R,t):
        """
        Very approximate
        """
        water_pairs = np.array(list(product(sorted(self.water_inds),repeat = 2)))
        
        frame_step = int(t/self.time_step)
        xyz_pos=self.traj.xyz
        
        N_QRt = []
        
        # I using all frames here in the trajectory. 
        # Caution: they may not be statistically independent
        for ii in range(self.n_frames):
            if (ii+frame_step)<self.n_frames:
                water_dist = []
                for this_pair in water_pairs:
                    water_dist.append(np.sqrt(np.sum((xyz_pos[ii,this_pair[0],:]- \
                                        xyz_pos[ii+frame_step,this_pair[1],:])**2.0)))
            
                N_QRt.append(self.single_frame_N_QR(Q,R,water_dist,time_dependent = True))
            else:
                print "Frame number %d and beyond not include in I_n(Q,R,t) calculations." % (ii+1)
                break
                
        
        if Q == 0:
            # equation (7) in ref [2]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*R**3.0
            
        else:
            # equations (28a-b) in ref [1]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return In_QRt
        
    def radial_dist(self, r_range, n_bins = 200):
        """
        Computer Simulation of Liquids book
        """
        water_pairs = np.array(list(combinations(sorted(self.water_inds),2)))
        water_dist = md.compute_distances(self.traj,water_pairs) # unit in nm
        
        g2s = []
        bin_edges = np.linspace(r_range[0],r_range[1],n_bins)
        for this_frame in water_dist:
            n_ideal = 4.* np.pi/3.*(bin_edges[1:]**3.0-bin_edges[:-1]**3.0)
            n_hist, _ = np.histogram(this_frame, bins = bin_edges,density=True)
            this_g2 = n_hist/n_ideal
            g2s.append(this_g2)

        g2s = np.array(g2s)
        bin_edges = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ave_g2 = np.mean(g2s,axis=0)
        g2_err = np.std(g2s,axis=0)/(ave_g2[-1]) # normalize assuming g2(r) converges to 1.0.
        ave_g2 = ave_g2/(ave_g2[-1])
        
        self.rdf = [bin_edges,ave_g2,g2_err]
    
    def estimate_struct_factor(self,Qs,R_inf,dt_ind):
        """
        Equations (22) and (23) in ref [1]
        """
        S_Q = []
        S_Qerr = []
        
        for this_Q in Qs:
            this_Sn_Q = self.struct_factor(this_Q,R_inf,dt_ind)
            this_Sn_0 = self.struct_factor(0,R_inf,dt_ind)
            if this_Q == 0:
                S_Q.append(this_Sn_Q[0])
                S_Qerr.append(this_Sn_Q[1])
            else:    
                this_S = this_Sn_Q[0] \
                + this_Sn_0[0]/(1-1/self.n_waters*4./3.*np.pi*self.rho*R_inf**3.0) \
                * 1./self.n_waters*4./3.*np.pi*self.rho \
                *3./this_Q**3.*(np.sin(this_Q*R_inf)-this_Q*R_inf*np.cos(this_Q*R_inf))
                S_Q.append(this_S)
                S_Qerr.append(np.sqrt(this_Sn_Q[1]**2.0 + this_Sn_0[1]**2.0)) # estimate of error
        
        self.ssf = [np.array(Qs), np.array(S_Q), np.array(S_Qerr)]
        
    
        
##############################################################################
# test
##############################################################################
# 
# data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
# traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
# print ('here is some info about the trajectory we are looking at:')
# print traj
# test = WaterStats(traj)
