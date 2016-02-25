#! /usr/bin/env python

##############################################################################
# Copyright 2016 Stanford University and the Author
#
# Author: Shenglan Qiao
#
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import numpy as np
import mdtraj
import h5py
import os
import time

from thor import xray
import matplotlib.pyplot as plt
##############################################################################
# Code
##############################################################################

run_name = 'run7'
t = mdtraj.load('data/nvt-pr_'+run_name+'.trr',top='data/water-sol_'+run_name+'.gro')

n_frames=t.n_frames
# n_frames =10

output_path = 'output/thor_run7_1.hdf5'
while os.path.isfile(output_path):
    print "Will not overwrite old file %s. Please enter new name:" % output_path
    output_path = raw_input()

# simulation parameters
n_shots = 1                      # total number of shots to do
n_molecules = 1                     # the number of molecules to include per shot
q_values = np.arange(1.5,2.0,0.01)
# q_values = np.arange(1.0,1.01,0.01)  # the |q| values of the rings to sim
n_phi = 360                         # number of pts around the rings


PI = []
tic=time.clock()
for ii in range(n_frames):
    rings = xray.Rings.simulate(t[ii],n_molecules,q_values,n_phi,1)
    if ii == 0:
        phis=rings.phi_values
        qs=rings.q_values
    PI.append(rings.polar_intensities[0])
toc=time.clock()
print toc-tic
     
output=h5py.File(output_path,'a')
output.create_dataset('polar_intensities',data=PI)
output.create_dataset('phi_values',data=phis)
output.create_dataset('q_values',data=qs)
output['polar_intensities'].attrs.create('num_molecules',n_molecules)

output.close()