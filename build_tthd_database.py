#! /usr/bin/env python

##############################################################################
# Copyright 2016 Stanford University and the Author
#
# Author: Shenglan Qiao
#
# This script take positions of tthd vertices, stored in pdb format and builds 
# a consolidated database (hdf5)
#
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
import sys
import getopt
##############################################################################
# Code
##############################################################################



run_name = 'run8'

tthd_vertices = []
with open(os.getcwd()+'/output_data/tthds_'+run_name+'.pdb','r') as pdb_file:
    for line in pdb_file:
        line_element = line.split()
        if line_element[0] == 'MODEL':
            this_tthd = []
        elif line_element[0] == 'HETATM':
            this_tthd.append([float(line_element[6]),float(line_element[7]),float(line_element[8])])
        elif line_element[0] == 'ENDMDL':
            tthd_vertices.append(this_tthd)