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

import h5py
from water_stats import WaterStats
import numpy as np
import os

import time

##############################################################################
# Code
##############################################################################

run_name = 'run13'
data_path = '/home/shenglan/MD_simulations/water_box/cubic_1nm_%s'%run_name
traj = md.load_trr(data_path+'/water-md.trr', top = data_path+'/water-md.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
ws = WaterStats(traj,run_name)
cut_off = 0.5

tic = time.clock()
for this_frame in range(ws.n_frames):
    tic_loop = time.clock()
    print 'finding tthds for frame %d' % this_frame
    ws.make_nearest_nb_tthds(cut_off,this_frame)
    toc_loop = time.clock()
    print('Total time for this frame: %.2f' %(toc_loop-tic_loop)) 

toc = time.clock()
print('Total time: %.2f' %(toc-tic))
