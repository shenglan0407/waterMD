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

##############################################################################
# Code
##############################################################################

data_path = '/home/shenglan/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
ws = WaterStats(traj)
cut_off = 0.5

for this_frame in range(ws.n_frames):
    tthds = []
    for this_water in ws.water_inds:
        tthds.extend(ws.make_tthd(this_water,cut_off,this_frame)) 
    ws.all_tthds.create_dataset(str(this_frame),data = tthds)
# 
# print len(tthds)
# print tthds[0][0]


#ws.all_tthds['1'][:]= tthds

ws.all_tthds.close()