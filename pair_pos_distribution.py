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

##############################################################################
# Code
##############################################################################

data_path = '/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj

#atom.index for all Oxygen of the water molecules
water_inds = traj.topology.select_atom_indices(selection='water')


