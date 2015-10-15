##############################################################################
# Copyright 2015 Stanford University and the Authors
#
# Author: Shenglan Qiao
# 
# The computation done here is based on the two following papers:
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
from itertools import combinations
import matplotlib.pyplot as plt

##############################################################################
# Code
##############################################################################

class WaterStat:
    def __init__(traj):

    def single_frame_N(Q,R,pair_dist,time_dependent = False):
    """
    See equations (1a) and (8) in ref [2]
    """
        if time_dependent:
        else:
            2*sum([this_pair for this_pair in pair_dist if this_pair <= R])+self.Nwaters
    

    def struct_factor(traj,Q,R,dt_min):
    

    def scat_func(traj,Q,R,t,dt_min):