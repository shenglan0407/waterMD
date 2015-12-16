##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# Functions for plotting computed results for 4-point correlator
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import os
import h5py
import mdtraj as md
from water_stats import WaterStats

from pylab import *
import numpy as np
from scipy.interpolate import griddata

##############################################################################
# Code
##############################################################################

def load_data(path,stride,partial=None,with_frame=False,debug = False):
    """Load computed correlator results from .csv file. return the average correlator, 
    its uncertainty, and list of phi used to compute the correlators
    
    Parameters
    ----------
    path : str
        path to the csv file
    stride : int
        if 1 load data computed from every frame. if >1 skip every stride-th frame
    partial : tuple-like, (i, j)
        range of frames for which computed data is to be loaded, from the ith frame to the jth frame
    with_frame : bool
        if True, the csv file has a line that records the indices of frames used to compute. default False
    debug : bool
        if True, return the full data set loaded to computed the average correlator and the magnitude of q1 as well
    """
    all_data = []
    with open(path,'r') as f:
            csvreader = csv.reader(f,delimiter=' ',quotechar='|')
            q1=(next(csvreader))
            q2=(next(csvreader))
            phi=(next(csvreader))
            if with_frame:
                frames = (next(csvreader))
            count = 0
            for row in csvreader:
                count += 1
                if count%stride==0:
                    all_data.append(row)
    if partial == None:
        all_data = np.array(all_data,dtype=float)
    else:
        all_data = np.array(all_data[partial[0]:partial[1]]
                            ,dtype=float)
    phi = np.array(phi,dtype = float)/np.pi*180
    q_mag = np.linalg.norm(np.array(q1,dtype=float))

    Corr = [np.mean(all_data[:,ii]) for ii in range(len(phi))]
    C_err = [np.std(np.array(all_data[:,ii])) \
             /np.sqrt(len(all_data[:,ii])) for ii in range(len(phi))]
    
    print len(all_data[:,ii])
    print 'The average error is %.2g. ' % np.mean(C_err)
    
    if debug:
        return Corr,C_err,phi,all_data,q_mag
    else:
        return Corr,C_err,phi

def load_q1_q2(path):
    """returns the q1 and the list of q2s used for computed correlator
    
    Parameters
    ----------
    path : str
        path to the csv file
    """
    with open(path,'r') as f:
            csvreader = csv.reader(f,delimiter=' ',quotechar='|')
            q1=(next(csvreader))
            q2=(next(csvreader))
    new_q2=[]
    for q2_item in q2:
        q2_float=[]
        for item in q2_item.split('[')[-1].split(']')[0].split(' '):
            try:
                q2_float.append(float(item))
            except ValueError:
                pass
        new_q2.append(q2_float)
    return np.array(q1,dtype=float),new_q2