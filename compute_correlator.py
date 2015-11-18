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

def main(argv):
    run_name = None
    outputfile = None
    number_qs = 10
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:",["ifile=","ofile=","n_phi="])
    except getopt.GetoptError:
        print 'compute_correlator.py -i <runname> -o <outputfile> -p <nphi>'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            
            print 'compute_correlator.py -i <runname> -o <outputfile> -p <nphi>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            run_name = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-p","--n_phi"):
            number_qs = arg
    print 'Input run is %s' % run_name
    print 'Output file is %s'% outputfile
    print 'Number of phi used is %s'%number_qs
    
    if run_name == None:
        print '<runname> must be provided.'
        print 'compute_correlator.py -i <runname> -o <outputfile> -p <nphi>'
        sys.exit(2)
        
    data_path = os.getcwd()+'/data'
    traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
    print ('here is some info about the trajectory we are looking at:')
    print traj
    run = WaterStats(traj,run_name,read_mod='r')
    frames = np.arange(1001)[500:]

    q = 1/0.3*np.pi*2.0
    theta_1 = np.pi/12.
    phi = np.linspace(-np.pi/2.,np.pi/2.,number_qs)
    dt = 1.0 # ps


    tic = time.clock()
    run.correlator(q,theta_1,frames,phi,cut_off = 0.5,output=outputfile)
    toc = time.clock()

    print("Correlator process time: %.2f" %(toc-tic))

    run.all_tthds.close()
    run.nearest_tthds.close()

if __name__ == "__main__":
   main(sys.argv[1:])


