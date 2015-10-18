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
            
        print ('the average time between each configuration is %f ps' % (np.mean(frame_steps)*self.time_step))
        
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
        
        frame_step = int(self.time_step/t)
        xyz_pos=self.traj.xyz
        
        N_QRt = []
        
        # I using all frames here in the trajectory. 
        # Caution: they may not be statistically independent
        for ii in range(self.n_frames):
            water_dist = []
            for this_pair in water_pairs:
                if (ii+frame_step)<self.n_frames:
                    water_dist.append(np.sqrt(np.sum((xyz_pos[ii,this_pair[0],:]- \
                                        xyz_pos[ii+frame_step,this_pair[1],:])**2.0)))
                else:
                    print "Frame number %d and beyond not include in I_n(Q,R,t) calculations." % (ii+1)
                    break
             
            N_QRt.append(self.single_frame_N_QR(Q,R,water_dist,time_dependent = True))
                
        
        if Q == 0:
            # equation (7) in ref [2]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*R**3.0
            
        else:
            # equations (28a-b) in ref [1]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return In_QRt
        
    def radial_dist(self):
        return md.compute_rdf(self.traj, pairs = self.water_ind,bin_width = 0.05)
              
        
        
##############################################################################
# test
##############################################################################

data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj)

rs, g_R = test.radial_dist()
fig3 = plt.figure()
plt.plot(rs,g_R)
plt.title('gn(r)')
plt.xlabel('r (nm)')
plt.ylabel('gn(r)')
fig3.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/gn_r.png')
plt.close(fig3)


R = 0.5
Rs = np.linspace(0.1,0.95,10)
Qs = 2.*np.pi/np.linspace(0.1,R,10)


Sn_0R = []
for R in Rs:
    Sn_0R.append(test.struct_factor(0,R,1)[0])
fig1 = plt.figure()
plt.plot(Rs, Sn_0R)
plt.title("Sn(0,R) with dt = 1.0 ps")
plt.xlabel("R (nm)")
plt.ylabel("Sn(0,R)") 
fig1.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/Sn_0R.png')
plt.close(fig1)

ts = np.linspace(1,10,10)
In_0tR1 = []
R1 = 0.1 # nm

for tt in ts:
    In_0tR1.append(test.scat_func(0,R1,tt))
fig2 = plt.figure()
plt.plot(ts, In_0TR1)
plt.title("In(0,t,R)")
plt.xlabel("t (pd)")
plt.ylabel("In(0,t,R)") 
fig2.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/In_0tR.png')
plt.close(fig2)


# Sn_QR=[]
# for Q in Qs:
#     Sn_QR.append(test.struct_factor(Q,R,1)[0])
# 
# print Qs
# print Sn_QR
# plt.plot(Qs,Sn_QR,'o')
# plt.show()