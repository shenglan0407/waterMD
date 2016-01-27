##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
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

import mdtraj as md

import sys
import os
import pickle
import h5py

import numpy as np
from itertools import combinations, product
import matplotlib.pyplot as plt
import random
import csv
import time

##############################################################################
# Code
##############################################################################

class WaterStats:
    def __init__(self,traj,run_name,read_mod='a'):
        """WaterStats object gives access to functions that compute statistical properties
        of water. 
        
        Parameters
        ----------
        traj : mdtraj.Trajetory object
            MD simulation trajectory of water
        run_name : str
            unique identifier of simulation run
        read_mod : {'a','r','w'}
            mode for opening relevant files. Default is 'a' for append
        
        Attributes
        ----------
        traj : mdtraj.Trajetory object
            MD simulation trajectory of water
        n_waters : int
            number of water molecules in simulation
        water_inds : list
            atom indices of oxygens of water molecules
        time_step : float
            time between simulation frames in ps
        n_frames : int
            total number of simulation frames
        total_time : flost
            total run time of simulation 
        run_name : str
            unique identifier of simulation run, used in naming convention of output files
        rho : float
            number density of water molecules
        nearest_tthds : h5py dataset (dict)
            database to store tetrahedrons made with only the nearest three neighbors of 
            water molecules. They are represented by two vectors r_ij (r_i-r_j) and r_kl
            (r_k-r_l) where i, j, k, l are the indices of the four waters.
        all_tthds : h5py dataset (dict)
            database to store all tetrahedrons formed by water molecules in simulation.
        pdb_tthds : csv file
            records coordinates of the four vertices of all nearest-neighbor tetrahedrons 
            in pdb format
        tthd_counter : int
            keeps count of the total number of nearest tetrahedrons
        """
        self.traj = traj
        self.n_waters = traj.topology.n_residues
        self.water_inds = traj.topology.select_atom_indices(selection='water')
        self.time_step=traj.timestep # in ps
        self.n_frames = traj.n_frames
        self.total_time = self.time_step*self.n_frames # in ps
        self.run_name = run_name
        
        self.rho = np.mean(self.n_waters/traj.unitcell_volumes) # in nm^-3
        
        if os.path.isdir(os.getcwd()+"/output_data"):
            pass
        else:
            os.mkdir(os.getcwd()+"/output_data")
            
        pdb_tthds_path = os.getcwd()+'/output_data/tthds_'+run_name+'.pdb'
        if os.path.isfile(pdb_tthds_path) and read_mod !='r':
            input = raw_input('pdb_tthds pdb file already exists. are you sure you want to append to/write it? [y or n]')
            if input =='y':
                self.pdb_tthds = open(pdb_tthds_path,'a')
                self.tthd_counter = 0
            else:
                self.pdb_tthds = None
        else:
            self.pdb_tthds = open(pdb_tthds_path,'a')
            self.tthd_counter = 0
        
        # test database, future way of organizing
        if os.path.isfile(os.getcwd()+'/output_data/test_tthd_data.hdf5'):
            self.test_datapath = os.getcwd()+'/output_data/test_tthd_data.hdf5'

    def make_tthd(self,vertex_ind,cut_off,frame_ind):
        """
        Given index of one vertex, find three other vertices within radius cut_off to form
        a unique tetrahedron in one single frame of traj
        returns list a array of shape (3,), indices of the three other vertices
        """
        nbs = md.compute_neighbors(self.traj[frame_ind],cut_off,[vertex_ind],haystack_indices = self.water_inds)[0]
        # print nbs.shape
        tthd_inds = np.array(list(combinations(nbs,3))) # I should not need to sort nbs, consider getting rid of the sorted part for all other instances of combinations
        
        tthds = []
        xyz_pos = self.traj[frame_ind].xyz
        for this_tthd in tthd_inds:
            # representing a tetrahedron with just two vectors, is this even right?
            r_ij = xyz_pos[0,vertex_ind,:]- xyz_pos[0,this_tthd[0],:]
            r_kl = xyz_pos[0,this_tthd[1],:]- xyz_pos[0,this_tthd[2],:] # nm
            
            tthds.append([r_ij,r_kl])
        return tthds 
            
    def four_point_struct_factor(self,qs,cut_off,set):
        """Computes the average correlator of tthds from a single simulation frame
        
        Parameters
        ----------
        qs : np.array, (n,2,3)
            pairs of (q1,q2) for computing correlator. n is the number of different pairs
            to compute for.
        cut_off : float
            cuttoff distance for looking for nearest neighbors. Should be set such that at
            least 3 waters can be found within that distance
        frame_ind: int
            index of single frame to compute for
        nearest_nb : bool
            if True, only use the nearest neighbors to form tthds
        """
        # this is a testing phase
        # if test_dataset is true, frame_ind (string) is interpreted as which set of tthd vectors to use.
        # there are 3 sets in the test dataset
        
        
        print set
        database = h5py.File(self.test_datapath,'r')
        if self.run_name in database:
            if set == 'all':
                all_sets = list(database[self.run_name].keys())
                all_sets.pop(all_sets.index('tthdVertices'))
                this_tthds = []
                for this_set in all_sets:
                    this_tthds.extend(database[self.run_name][this_set][:][:]) 
            else:
                this_tthds = database[self.run_name][set][:][:]
            database.close()
        else:        
            print "data for this run do not exist"
            database.close()
            sys.exit()
        
        n_tthds = len(this_tthds)
        
            #To-do: ask user if she wants to create the tthd data for run
        
        
        corr_single_set = []
        err_single_set=[]
        aa = [3.0485,2.2868,1.5463,0.867]
        bb = [13.2771,5.7011,0.3239,32.9089]
        cc = 0.2580
        form_factor = cc
        for this_a,this_b in zip(aa,bb):
            form_factor += this_a * np.exp(-this_b*(np.linalg.norm(qs[0])/(4*np.pi))**2.0)
        
        
        this_I1 = []
        for tt in this_tthds:
            this_I1.append((1+np.cos(np.dot(qs[0][0],tt[0])))*form_factor**2.0*2.0)
        
        tic =time.clock()
        I1_err = np.std(this_I1)/np.sqrt(len(this_I1))
        I1_mean = np.mean(this_I1)
        toc = time.clock()
        print "time taken is %g sec."%(toc-tic)
        assert len(this_I1) == n_tthds
        
        for this_q in qs:
            this_I2 = []
            
            for tt in this_tthds:             
                this_I2.append((1+np.cos(np.dot(this_q[1],tt[1])))*form_factor**2.0*2.0)
            
            assert n_tthds == len(this_I2)
            
            this_I1I2 = [this_I1[ii]*this_I2[ii] for ii in range(len(this_I1))]
            err = np.std(this_I1I2)/np.sqrt(n_tthds)
            I1I2_mean = np.mean(this_I1I2)
            
            I2_err = np.std(this_I2)/np.sqrt(len(this_I2))
            I2_mean = np.mean(this_I2)
            
            this_result = I1I2_mean-I2_mean*I1_mean
            this_err = np.sqrt(err**2.0+((I1_err/I1_mean)**2.0+(I2_err/I2_mean)**2.0)*(I1_mean*I2_mean)**2.0)
            corr_single_set.append(this_result)
            err_single_set.append(this_err)
            
            
        
        return corr_single_set,err_single_set

    
    def atomic_form_factor(self,element):
        """computes the atomic form factor givent the element
        """
        pass

    def correlator(self,q,wavelength,sets,phi,cut_off = 0.5,output = None):
        """Computes 4-point correlator and saves results from each simulation frames in .csv file
        Assume incident beam is along the z axis and water box sample is at origin
        
        Parameters
        ----------
        q : float
            magnitude of the q vectors, assumed to be same for now (auto-correlator) 
        wavelength : float
            wavelength in nm of the incident beam
        frames : list
            indices of frames in which to sample tthds and compute correlator
        phi : list
            angles in radian between q1 and q2 in the xy-plane (detector plane). It's used
            to generate q2 from q1, which is fixed.
        cut_off: float
            distance cutoff for making tthds. default 0.5 nm
        nearest_nb : bool
            If True only use the nearest three nbs to form tthds, default True
        output : str
            File name to save the results in, default None.
        """
 
        
        q_beam=2.0*np.pi/wavelength
        
        # generate q1 and many q2s
        q1 = np.array([q*np.sqrt(1-(q/(2.*q_beam))**2.0),0,-q**2.0/(2.0*q_beam)])
        q2 = np.array([np.array([q*np.sqrt(1-(q/(2.*q_beam))**2.0)*np.cos(this_phi),q*np.sqrt(1-(q/(2.*q_beam))**2.0)*np.sin(this_phi),q1[2]]) for this_phi in phi])
        # make pairs of q1 and q2
        qs = np.array([[q1,this_q2] for this_q2 in q2])
        
        if output == None:
            output_path = os.getcwd()+'/computed_results/corr_'+self.run_name+'.csv'
        else:
            output_path =  os.getcwd()+'/computed_results/'+output
        while os.path.isfile(output_path):
            print 'Computed results exist, will not overwrite file'
            input = raw_input('enter new file name here:')
            output_path = os.getcwd()+'/computed_results/'+input
            
        with open(output_path,'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            csvwriter.writerow(q1)
            csvwriter.writerow(q2)
            csvwriter.writerow(phi)
            csvwriter.writerow(sets)
        csvfile.close()

        for this_set in sets:
            this_row,this_err = self.four_point_struct_factor(qs,cut_off,this_set)
            with open(output_path,'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=' ',
                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                csvwriter.writerow(this_row)
                csvwriter.writerow(this_err)
                csvfile.flush()
    
    
    def find_nearest_nbs(self,cut_off,frame_ind,N_nbs):
        """Finds the N nearest neighbors of all waters in a single simulation frame and returns
        all unique nearest-neighbor tthds as an (n,4) array. n is the number of tthds found in 
        this frame. Every tthd is represented by 4 indices of the water molecules (oxygens) that
        are the for vertices. More generally, this function be used to form ensembles with 
        (N_nbs+1) vertices. The numpy array returned will have shape (n, N_nbs+1).
        
        Parameters
        ----------
        cut_off : float
            cuttoff distance for looking for nearest neighbors. Should be set such that at
            least N_nbs waters can be found within that distance
        frame_ind : int
            the index of the frame in which to find nearest neighbors
        N_nbs : int
            the number of nearest neighbors to return. set to 3 in this case since we are 
            interest in 4-point correlators.
        """
        nearest_nbs = []
        frame = self.traj[frame_ind]
        
        for this_ind in self.water_inds:
            # find all the neighbors with the cut_off distance for each water molecule
            nbs = md.compute_neighbors(frame,
        cut_off,[this_ind],haystack_indices = self.water_inds)[0]
            
            while len(nbs)!= N_nbs:
                if len(nbs) > N_nbs:
                    # compute all the distance between neighbors and only keep the closest
                    # N_nbs ones
                    pairs = [[this_ind,this_nb] for this_nb in nbs]
                    distances=md.compute_distances(frame,pairs)[0]
                    min_three = np.argsort(distances)[:][:N_nbs]
                    nbs = [nbs[ii] for ii in min_three]
                
                else:
                    print 'increase cut_off!'
            # keep nbs a list so append can happen
            nbs = [nbs[ii] for ii in range(len(nbs))]
            # add the ind of the water for whom neighbors have been found to the list
            nbs.append(this_ind)
            # sort the list so it can be comapred to set of neighbors we already found
            nbs.sort()
            
            if nbs in nearest_nbs:
                # if this set of nearest neighbors already exist, don't added it to the list
                print "not unique"
            else:
                # only add to list if unique
                nearest_nbs.append(nbs)
        # return list as an numpy array
        return np.array(nearest_nbs)
    
    def make_nearest_nb_tthds(self,cut_off,frame_ind):
        """Returns a list of nearest_nb tthds from a single simulation frame. Each tthd  is 
        represented by two vectors, r_ij and r_kl. r_ij = (r_i - r_j) and r_kl = (r_k-r_l) 
        while r_i is the position vector of indexed-i water molecule etc. Also write the four
        position vecotrs of the four vertices to a pdb format file.
        
        Parameters
        ----------
        cut_off : float
            cuttoff distance for looking for nearest neighbors. Should be set such that at
            least 3 waters can be found within that distance
        frame_ind : int
            the index of the frame in which to find nearest neighbor tthds
        """
        nbs = self.find_nearest_nbs(cut_off,frame_ind,3)
        
        xyz_pos = self.traj[frame_ind].xyz
        half_box = self.traj.unitcell_lengths[0][0]/2.
        for this_nb in nbs:
            # correcting for periodic boundary condition
            
            r1 = xyz_pos[0,this_nb[0],:]
            r2 = xyz_pos[0,this_nb[1],:]
            r3 = xyz_pos[0,this_nb[2],:]
            r4 = xyz_pos[0,this_nb[3],:]
            
            # check if any of the nearest neighbors are on the other side of the simulation box
            dr = r1-r2
            for ii in range(3):
                # if x, y or z coordinate is greater than half the simulation box size, correct
                # for periodic boundary condition by 'folding' the position vector over
                if np.abs(dr[ii]) > half_box:
                    if r2[ii]<r1[ii]:
                        r2[ii] = r2[ii]+half_box*2.
                    else:
                        r2[ii] = r2[ii]-half_box*2.
        
            dr = r1-r3
            for ii in range(3):
                if np.abs(dr[ii]) > half_box:
                    if r3[ii]<r1[ii]:
                        r3[ii] = r3[ii]+half_box*2.
                    else:
                        r3[ii] = r3[ii]-half_box*2.
        
            dr = r1-r4
            for ii in range(3):
                if np.abs(dr[ii]) > half_box:
                    if r4[ii]<r1[ii]:
                        r4[ii] = r4[ii]+half_box*2.
                    else:
                        r4[ii] = r4[ii]-half_box*2.
            
            if self.pdb_tthds == None:
                pass
            else:
                # store the coordinates the four vertices in pdb format. The center of tthd
                # always set to origin
                self.tthd_counter += 1
                centeroid = (r1+r2+r3+r4)/4.
                r1=r1-centeroid
                r2=r2-centeroid
                r3=r3-centeroid
                r4=r4-centeroid
                
                self.pdb_tthds.write('%s        %d \n' % ('MODEL' , self.tthd_counter))
                self.pdb_tthds.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 1,'O','HOH','I',1,' '
                ,r1[0],r1[1],r1[2], 1.0, 0.0,'O'))
                self.pdb_tthds.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 2,'O','HOH','I',2,' '
                ,r2[0],r2[1],r2[2], 1.0, 0.0,'O'))
                self.pdb_tthds.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 3,'O','HOH','I',3,' '
                ,r3[0],r3[1],r3[2], 1.0, 0.0,'O'))
                self.pdb_tthds.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % ('HETATM', 4,'O','HOH','I',4,' '
                ,r4[0],r4[1],r4[2], 1.0, 0.0,'O'))
                self.pdb_tthds.write('%s \n' % 'ENDMDL')
            
                self.pdb_tthds.flush()