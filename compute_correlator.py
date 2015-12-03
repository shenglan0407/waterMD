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
import sys
import getopt
##############################################################################
# Code
##############################################################################
def usage():
    print 'compute_correlator.py -i <runname> -o <outputfile> -q <q_inv> -p <nphi> -s <fstart> -e <fend>'

def main(argv):
    # default values for options
    run_name = None
    outputfile = None
    number_qs = 10
    frame_start = 1
    frame_end = None
    q_inverse = 0.3
    
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:q:p:s:e:",["ifile=","ofile=","q_inv=","n_phi=","fstart=","fend="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            run_name = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-q", "--q_inv"):
            q_inverse = float(arg)
        elif opt in ("-p","--n_phi"):
            number_qs = int(arg)
        elif opt in ("-s","--fstart"):
            frame_start = int(arg)
        elif opt in ("-e","--fend"):
            frame_end = int(arg)
    print 'Input run is %s' % run_name
    print 'Output file is %s'% outputfile
    print 'Number of phi used is %d'%number_qs
    print 'Inverse of q is %.2f nm'%q_inverse
    
    if run_name == None:
        print '<runname> must be provided.'
        usage()
        sys.exit(2)
        
    data_path = os.getcwd()+'/data'
    traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
    print ('here is some info about the trajectory we are looking at:')
    print traj
    run = WaterStats(traj,run_name,read_mod='r')
    if frame_start>=run.n_frames:
        print 'Starting frame cannot be greater than the number of frames in simulation.'
        usage()
        sys.exit(2)
    elif frame_end == None:
        frames = np.arange(run.n_frames)[frame_start:]
    else:
        frames = np.arange(run.n_frames)[frame_start:frame_end]

    q = 1/q_inverse*np.pi*2.0
    # wavelength of laser
    wavelength = 0.1
    phi = np.linspace(-np.pi,np.pi,number_qs)
    dt = 1.0 # ps


    tic = time.clock()
    run.correlator(q,wavelength,frames,phi,cut_off = 0.5,output=outputfile)
    toc = time.clock()

    print("Correlator process time: %.2f" %(toc-tic))

    run.all_tthds.close()
    run.nearest_tthds.close()

if __name__ == "__main__":
   main(sys.argv[1:])


