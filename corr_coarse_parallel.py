##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# compute the correlator using using coarse-grained multiprocessing
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import mdtraj as md

import os
import h5py

import numpy as np
from itertools import combinations, product
import matplotlib.pyplot as plt
import random
from joblib import Parallel, delayed

##############################################################################
# Code
##############################################################################