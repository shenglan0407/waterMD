##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
#
# Testing class WaterStats
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

from water_stats import WaterStats

import mdtraj as md
import numpy as np

import os
import time

##############################################################################
# Code
##############################################################################

frames = np.arange(101)[1:50]
run_name = 'run1'

data_path = os.getcwd()+'/data'
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
run = WaterStats(traj,run_name,read_mod='r')

q = 1/0.3*np.pi*2.0
theta_1 = np.pi/12.
phi = np.linspace(-np.pi,np.pi,10)
dt = 1.0 # ps


tic = time.clock()
run.correlator(q,theta_1,frames,phi,cut_off = 0.5)
toc = time.clock()

print("Correlator process time: %.2f" %(toc-tic))

run.all_tthds.close()
run.nearest_tthds.close()
