#! /usr/bin/env python

##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
#
# Script to compute cross and auto correlator from scattering data
# 
#############################################################################


##############################################################################
# Imports
##############################################################################

import numpy as np
from corrRings import calc_autocorr
import getopt

import os
import sys
import h5py
##############################################################################
# Code
##############################################################################
def usage():
    print './thor_compute_autocorr.py -i <input> -o <outputfile> -s <start_ind> -e <end_ind> -p <phi_stride>'

def main(argv):
    inputfile = None
    outputfile = None
    s_ind = 0
    e_ind = 5
    stride = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:e:p:",["ifile=","ofile=","start_ind=","end_ind=","phi_stride="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s", "--start_ind"):
            s_ind = int(arg)
        elif opt in ("-e","--end_ind"):
            e_ind = int(arg)
        elif opt in ("-p" or "--phi_stride"):
            stride = int(arg)
    
    if inputfile == None:
        print '<input> must be provided.'
        usage()
        sys.exit(2)
    if outputfile == None:
        print '<outputfile> must be provided.'
        usage()
        sys.exit(2)
    
    outputpath='computed_results/'+outputfile+'.hdf5'
    inputpath='computed_results/'+inputfile+'.hdf5'
    print 'Input file is %s' % inputpath
    print 'Output file is %s'% outputpath
    
    try:
        fin = h5py.File(inputpath,'r')
    except IOError:
        print "input does not exist, double check path %s"%inputpath
        sys.exit(2)
        
    fout = h5py.File(outputpath,'a')
    
    try:
        PI=fin['polar_intensities'][:,s_ind:e_ind]
        phi=fin['phi_values'][:]
        q_values = fin['q_values'][s_ind:e_ind]
    except KeyError:
        print "Key error in input file!"
        sys.exit(2)
    
    n_shots=len(PI)
    print 'Data set has %d snapshots'%n_shots
    
    
    deltas = np.arange(0,len(phi),stride)
    for ii in range(len(q_values)):

        corr = np.zeros(len(deltas))
        for this_ring in PI[:,ii]:
            corr =corr+[calc_autocorr(this_ring,delta) for delta in deltas]
        corr=corr/n_shots
        
        try:
            fout.create_dataset(str(q_values[ii]),data=corr)
        except RuntimeError:
            print "Somehow autocorrelator for q=%.3f already exists in this output file %s"% (q_values[ii], outputpath)
            print "Terminating this run..."
            print opts
            print args
            sys.exit(2)
    
    if 'delta_phi' in fout.keys():
        pass
    else:
        delta_phi=[phi[jj] for jj in deltas]
        try:
            fout.create_dataset('delta_phi',data=delta_phi)
        except RuntimeError:
            pass
        
            
    
    fin.close()
    fout.close()

if __name__ == "__main__":
   main(sys.argv[1:])