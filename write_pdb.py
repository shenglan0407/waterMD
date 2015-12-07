##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# write pdb
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import mdtraj as md

import h5py
from water_stats import WaterStats
import numpy as np
import os

import time

##############################################################################
# Code
##############################################################################

run_name = 'run5'
data_path = os.getcwd()+'/data'
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
ws = WaterStats(traj,run_name)
cut_off = 0.5
frame_ind = 1

half_box = ws.traj.unitcell_lengths[0][0]/2.*10

inds = range(1001)[101:]
count = 33897
for frame_ind in inds:
    nbs = ws.find_nearest_nbs(cut_off,frame_ind,3)
    xyz_pos = ws.traj[frame_ind].xyz
    
    with open(os.getcwd()+'/output_data/tthd_pdb_1000.pdb','a') as f:
        for this_nb in nbs:
            
            r1 = xyz_pos[0,this_nb[0],:]*10
            r2 = xyz_pos[0,this_nb[1],:]*10
            r3 = xyz_pos[0,this_nb[2],:]*10
            r4 = xyz_pos[0,this_nb[3],:]*10
        
            dr = r1-r2
            for ii in range(3):
                if np.abs(dr[ii]) > half_box:
                    if r2[ii]<r1[ii]:
                        r2[ii] = r2[ii]+half_box*2.
                    else:
                        r2[ii] = r2[ii]-half_box*2.
        
            dr = r1-r3
            for ii in range(3):
                if np.abs(dr[ii]) > half_box:
                    if r3[ii]<r1[ii]:
                        r3[ii] = r3[ii]+half_box*2.
                    else:
                        r3[ii] = r3[ii]-half_box*2.
        
            dr = r1-r4
            for ii in range(3):
                if np.abs(dr[ii]) > half_box:
                    if r4[ii]<r1[ii]:
                        r4[ii] = r4[ii]+half_box*2.
                    else:
                        r4[ii] = r4[ii]-half_box*2.
            centeroid = (r1+r2+r3+r4)/4.
            r1=r1-centeroid
            r2=r2-centeroid
            r3=r3-centeroid
            r4=r4-centeroid
            count+=1
            f.write('%s        %d \n' % ('MODEL' , count))
            f.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 1,'O','HOH','I',1,' '
            ,r1[0],r1[1],r1[2], 1.0, 0.0,'O'))
            f.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 2,'O','HOH','I',2,' '
            ,r2[0],r2[1],r2[2], 1.0, 0.0,'O'))
            f.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 3,'O','HOH','I',3,' '
            ,r3[0],r3[1],r3[2], 1.0, 0.0,'O'))
            f.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 4,'O','HOH','I',4,' '
            ,r4[0],r4[1],r4[2], 1.0, 0.0,'O'))
            f.write('%s \n' % 'ENDMDL')
            
            f.flush()
        