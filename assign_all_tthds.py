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

# this is a change
##############################################################################
# Imports
##############################################################################

import mdtraj as md

import h5py
from water_stats import WaterStats
import numpy as np
import os

##############################################################################
# Code
##############################################################################

run_name = 'run4'
data_path = os.getcwd()+'/data'
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
ws = WaterStats(traj,run_name)
cut_off = 0.5

for this_frame in range(ws.n_frames):
    if str(this_frame) in ws.all_tthds:
        pass
    else:
        tthds = []
        for this_water in ws.water_inds:
            tthds.extend(ws.make_tthd(this_water,cut_off,this_frame)) 
        ws.all_tthds.create_dataset(str(this_frame),data = tthds)

# print len(tthds)
# print tthds[0][0]


#ws.all_tthds['1'][:]= tthds

ws.all_tthds.close()
