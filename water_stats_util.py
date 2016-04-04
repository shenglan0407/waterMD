##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# Functions for plotting, handling, and analyzing computed results for 4-point correlator
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import h5py
import sys

from pylab import *
import numpy as np
from scipy.interpolate import griddata

##############################################################################
# Code
##############################################################################

def load_data(path,debug = False):
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
    debug : bool
        if True, return the full data set loaded to computed the average correlator and the magnitude of q1 as well
    """
    corr_data = []
    corr_err = []
    with open(path,'r') as f:
            csvreader = csv.reader(f,delimiter=' ',quotechar='|')
            q1=(next(csvreader))
            q2=(next(csvreader))
            phi=(next(csvreader))
            sets = (next(csvreader))
            even = True
            for row in csvreader:
                if even:
                    corr_data.append(row)
                    even=False
                else:
                    corr_err.append(row)
                    even = True
    corr_data = np.array(corr_data,dtype=float)
    corr_err = np.array(corr_err,dtype=float)
    
    phi = np.array(phi,dtype = float)/np.pi*180
    q_mag = np.linalg.norm(np.array(q1,dtype=float))
    
    q1,new_q2 = load_q1_q2(path)
    
    cos_psi = [np.dot(q1,this_q2)/(q_mag*np.linalg.norm(np.array(this_q2,dtype=float))) for this_q2 in new_q2]

    
    if debug:
        return corr_data,corr_err,phi,cos_psi,sets,q_mag
    else:
        return corr_data,corr_err,phi,cos_psi,sets

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
    
def check_convergence(path,increment):
    """produce correlator average over increasing number of tthds. 
    return a list of differences between consecutive correlator curves as more tthds from 
    more simulation frames are added to the overall average. Also return the magnitude of qs
    for which the convergence test is done (assuming autocorrelator at this point)
    
    Parameters
    ----------
    path : str
        path to the csv file
    increment : int
        number of frames to increase by when computing the next correlator curve
        
    Returns
    -------
    diffs : list
        difference in consecutive correlator curves averaged over increasing number of simulation frames
    q_mag : float
        magnitude of q's for which the convergence test is done (assuming auto-correlation for now) 
    """
    
    _,_,phi,_,all_data,q_mag =load_data(path
                        ,1
                        ,debug = True
                        )
    increment = 1
    partials = [[0,(i+1)*increment] for i in range(len(all_data)/increment)]
    accu_corr = []
    for this_partial in partials:
        data = all_data[this_partial[0]:this_partial[1]]
        Corr = [np.mean(data[:,ii]) for ii in range(len(phi))]
        accu_corr.append(Corr)
    
    diffs = []
    for ii in range(len(accu_corr)-1):
        diff = np.sqrt(np.linalg.norm(np.array(accu_corr[ii])\
                       -np.array(accu_corr[ii+1])))/len(all_data[0])
        diffs.append(diff)
    return diffs, q_mag

def make_file_paths(q_invs,run_name,n_points,tag):
    """returns a list of file names that are consistent with correlator naming conventions
    
    Parameters
    ----------
    q_invs : list of str
        list of q_invs in nm for which correlator has been computed
    run_name : str
        name of the simulation run by naming convention

    Returns
    -------
    file_path: list
        list of str, paths to computed data files
    """
    file_paths = [os.getcwd()+\
                        '/computed_results/corr_'+run_name+'_'\
                                    +this_q+'q_'+str(n_points)+'p_'+tag+'.csv'\
             for this_q in q_invs]
    return file_paths

def make_polar_plot_coordinates(file_paths,q_s,this_set,align = True,phi_or_psi = True):
    """Return the polar coordinates needed to make a polar plot of correlator
    
    Parameters
    ----------
    file_paths : list of str
        list of paths to file whose data we want to combine to make a polar plot
    q_s : list of float
        values of used for computing all the correlators. must be in the same order
        as paths in file_path

    Returns
    -------
    angle_rad : list of float
        angles in the polar plot where there is data, radian
    data : array-like
        data points to be plotted on the polar plot. shape is len(file_paths) by len(angle_rad) 

    """
    all_corr = []
    for this_file in file_paths:
        corr_data,_,phi,cos_psi,sets =load_data(this_file
                        ,1
                        )
        if this_set in sets:
            Corr = corr_data[np.where(np.array(sets)==this_set)[0][0]]
        else:
            print "the set you want does not exist in this file:"
            print this_file
            sys.exit(2)

        if align:
            all_corr.append(Corr-min(Corr))
        else:
            all_corr.append(Corr)
    all_corr = np.array(all_corr)
    
    if phi_or_psi:
        angle_rad = [this_phi/360*2*np.pi for this_phi in phi]
    else: 
        angle_rad = [np.arccos(this_cos_psi) for this_cos_psi in cos_psi]
    points = [[q_s[i],angle_rad[j]] for i in range(len(q_s)) for j in range(len(phi))] 
    
    # create grid values and interpolate
    grid_r, grid_phi = np.meshgrid(q_s, angle_rad)
    data = griddata(points, all_corr.ravel(), (grid_r, grid_phi), method='cubic',fill_value=0)
    data = data.T
    
    return angle_rad, data