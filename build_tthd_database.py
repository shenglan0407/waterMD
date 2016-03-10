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

import h5py
import numpy as np
import os
import sys
import time

import mdtraj as md
##############################################################################
# Code
##############################################################################



run_name = 'run13'

database = h5py.File(os.getcwd()+'/output_data/test_tthd_data.hdf5','a')

if run_name in database:
    if len(database[run_name].items())==0:
        del database[run_name]
        print 'deleted empty group...'
        
        this_run = database.create_group(run_name)
        data_path = '/home/shenglan/MD_simulations/water_box/cubic_1nm_'+run_name
        traj = md.load_trr(data_path+'/water-md.trr', top = data_path+'/water-md.gro')
        print ('here is some info about the trajectory we are looking at:')
        print traj
        n_waters = traj.topology.n_residues
        n_frames = traj.n_frames
        box_length = traj.unitcell_lengths[0][0]
    
        this_run.attrs.create('runName',run_name)
        this_run.attrs.create('boxSize',box_length)
        this_run.attrs.create('nFrames',n_frames)
        this_run.attrs.create('nWaters',n_waters)
    else:   
        print "Data for "+run_name+" already exist in database!"
        this_run = database[run_name]

else:
    this_run = database.create_group(run_name)
    data_path = '/home/shenglan/MD_simulations/water_box/cubic_1nm_'+run_name
    traj = md.load_trr(data_path+'/water-md.trr', top = data_path+'/water-md.gro')
    print ('here is some info about the trajectory we are looking at:')
    print traj
    n_waters = traj.topology.n_residues
    n_frames = traj.n_frames
    box_length = traj.unitcell_lengths[0][0]
    
    this_run.attrs.create('runName',run_name)
    this_run.attrs.create('boxSize',box_length)
    this_run.attrs.create('nFrames',n_frames)
    this_run.attrs.create('nWaters',n_waters)
    
if 'tthdVertices' in this_run:
    print 'tthd vertices exist'
else:
    tic=time.clock()
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
    toc=time.clock()
    print 'Time taken to load pdb is %g sec' % (toc-tic)
    
    tic = time.clock()
    this_run.create_dataset('tthdVertices',data=tthd_vertices)
    toc=time.clock()
    print 'Time taken to store pdb is %g sec' % (toc-tic)
    

if set(['tthdSet1','tthdSet2','tthdSet3']).issubset(set(this_run.keys())):
    print 'all tthd sets already exist. exit without doing anything.'
    database.close()
    sys.exit()
else:
    tic = time.clock()
    pairings = [[0,1,2,3],[0,2,3,1],[0,3,1,2]]

    all_tthds=[[],[],[]]
    for this_tthd in tthd_vertices:
        for ii in range(len(pairings)):
            v_ij = [x[0]-x[1] for x in zip(this_tthd[pairings[ii][0]],this_tthd[pairings[ii][1]])]
            v_lm = [x[0]-x[1] for x in zip(this_tthd[pairings[ii][2]],this_tthd[pairings[ii][3]])]
            all_tthds[ii].append([v_ij,v_lm])
    toc = time.clock()
    toc=time.clock()
    print 'Time taken to permutate vertices is %g sec' % (toc-tic)
    
    tic = time.clock()
    this_run.create_dataset('tthdSet1',data=all_tthds[0])
    this_run.create_dataset('tthdSet2',data=all_tthds[1])
    this_run.create_dataset('tthdSet3',data=all_tthds[2])
    
    this_run['tthdSet1'].attrs.create('verticesPairing',[0,1,2,3])
    this_run['tthdSet2'].attrs.create('verticesPairing',[0,2,3,1])
    this_run['tthdSet3'].attrs.create('verticesPairing',[0,3,1,2])
    toc=time.clock()
    print 'Time taken to store tthd sets is %g sec' % (toc-tic)


database.close()