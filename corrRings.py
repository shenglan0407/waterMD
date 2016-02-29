#! /usr/bin/env python

##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
#
# Functions to compute auto and cross correlator directly from simulated scattering
# intensities
# 
#############################################################################


##############################################################################
# Imports
##############################################################################

import numpy as np

##############################################################################
# Code
##############################################################################
def calc_autocorr(ring, delta):
    
    return np.mean([ring[ii]*ring[ii+delta] for ii in range(len(ring)-delta)])

def calc_crosscorr(ring1,ring2, delta):
    try:
        assert len(ring1)==len(ring2)
        return np.mean([ring1[ii]*ring2[ii+delta] for ii in range(len(ring1)-delta)])
    except AssertionError:
        print "Error! The two Bragg rings you are trying to cross-correlate have different number of data points."
        return