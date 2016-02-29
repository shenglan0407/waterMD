#! /usr/bin/env python

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
import platform
##############################################################################
# Code
##############################################################################
def usage():
    print './compute_correlator.py -i <runname> -o <outputfile> -q <q_inv> -p <nphi> -r <phi_range> -s <sets> -c'

def main(argv):
    # default values for options
    run_name = None
    outputfile = None
    number_qs = 10
    sets = ["all"]
#     q_inverse = 0.3
    q_range = [0.3]
    phi_range = [0,1.0]
    autocorr = True
    
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:q:r:p:s:c",["ifile=","ofile=","q_inv=","n_phi=","phi_range=","sets="])
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
            qs = arg.split('/')
            print qs
            if len(qs) > 1:
                try:
                    q_range = [float(this_q) for this_q in qs]
                except ValueError:
                    print 'Enter a single value of q of a range of values separated by /'
                    sys.exit(2)
            else:
                try:
                    q_range = [float(arg)]
                except ValueError:
                    print 'Enter a single value of q of a range of values separated by space'
                    sys.exit(2)
                
        elif opt in ("-p","--n_phi"):
            number_qs = int(arg)
        elif opt in ("-r","--phi_range"):
            phis = arg.split('/')
            if len(phis) == 2:
                try:
                    phi_range = [float(this_phi) for this_phi in phis]
                except ValueError:
                    print 'Enter two values (in units of pi) for starting and ending phi separated by /. '
                    sys.exit(2)
            else:
                print 'Enter two values (in units of pi) for starting and ending phi separated by /. '
                sys.exit(2)
        elif opt in ("-s","--sets"):
            sets = arg.split('/')
            if set(sets).issubset(set(["tthdSet1","tthdSet2","tthdSet3","all"])):
                pass
            else:
                print 'Enter sets of tthds for computation separated by /. '
                print 'You entered: '
                print sets
                print 'Allowed values are [\"tthdSet1\",\"tthdSet2\",\"tthdSet3\",\"all\"]'
                sys.exit(2)
        elif opt in ("-c","--cross"):
            autocorr=False
            print "Computing cross correlator between two q's..."
    print 'Input run is %s' % run_name
    print 'Output file is %s'% outputfile
    print 'Number of phi used is %d from phi = %.3g pi to phi = %.3g pi'%(number_qs,phi_range[0],phi_range[1])
#     print 'Inverse of q is %.2f nm'%q_inverse
    print 'Computing correlators for the following q_inverse values in nm:'
    print q_range
    
    print "Sets of tthds used are: "
    print sets   
    
    if run_name == None:
        print '<runname> must be provided.'
        usage()
        sys.exit(2)
        

    if platform.node()=='DN0a22c83d.SUNet':
        # running this on shenglan's mac
        data_path = '/Users/shenglanqiao/zauber/MD_simulations/water_box/cubic_1nm_'+run_name
    elif platform.node()=='Zauber':
        # running in on zauber
        data_path = '/home/shenglan/MD_simulations/water_box/cubic_1nm_'+run_name
    else:
        data_path =os.getcwd()+'/data'
    traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
    
    print ('here is some info about the trajectory we are looking at:')
    print traj
    run = WaterStats(traj,run_name,read_mod='r')
    
    # wavelength of laser
    wavelength = 0.1
    phi = np.linspace(phi_range[0]*np.pi,phi_range[1]*np.pi,number_qs)
    dt = 1.0 # ps
    
#########################################################################################
# Auto-correlator
    if autocorr:
        for q_inverse in q_range:
            print('computing for q_invers = %.3g nm' % q_inverse)
            q = 1/q_inverse*np.pi*2.0
        
            if outputfile == None:
                if len(sets)==1:
                    outputfile = 'corr_'+run_name+\
                    '_'+str(q_inverse)+'q_'+str(number_qs)+'p_'+\
                    str(sets[0])+\
                    '.csv'
                elif "all" in sets:
                    outputfile = 'corr_'+run_name+\
                    '_'+str(q_inverse)+'q_'+str(number_qs)+'p_all.csv'
                else:
                    outputfile = 'corr_'+run_name+\
                    '_'+str(q_inverse)+'q_'+str(number_qs)+'p_compare.csv'

            tic = time.clock()
            run.correlator(q,wavelength,sets,phi,cut_off = 0.5,output=outputfile)
            toc = time.clock()

            print("Correlator process time for %.3g nm: %.2f" % (q_inverse,(toc-tic)))
            outputfile = None
#########################################################################################
# cross-correlator
    else:
        try:
            assert len(q_range)==2
        except AssertionError:
            print "need exactly two q values for cross-correlator"
            sys.exit(2)
        q=[1/q_inverse*np.pi*2.0 for q_inverse in q_range]
        print "computing cross correlator between %.3f and %.3f"%(q[0],q[1])
        
        if outputfile == None:
            if len(sets)==1:
                outputfile = 'crosscorr_'+run_name+\
                '_'+str(q[0])+'_'+str(q[1])+'q_'+str(number_qs)+'p_'+\
                str(sets[0])+\
                '.csv'
            elif "all" in sets:
                outputfile = 'crosscorr_'+run_name+\
                '_'+str(q[0])+'_'+str(q[1])+'q_'+str(number_qs)+'p_all.csv'
            else:
                outputfile = 'crosscorr_'+run_name+\
                '_'+str(q[0])+'_'+str(q[1])+'q_'+str(number_qs)+'p_compare.csv'
        tic = time.clock()
        run.correlator(q,wavelength,sets,phi,cut_off = 0.5,output=outputfile)
        toc = time.clock()

        print("Correlator process time: %.2f" % (toc-tic))

    

if __name__ == "__main__":
   main(sys.argv[1:])


