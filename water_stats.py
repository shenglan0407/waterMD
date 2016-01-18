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
        
        # dictionary to store all tthd vectors, keys are frame numbers (str), 0-indexed
        tthds_path = os.getcwd()+'/output_data/all_tthds_'+run_name+'.hdf5'
        if os.path.isfile(tthds_path) and read_mod !='r':
            input = raw_input('all_tthds database already exists. are you sure you want to append to/write it? [y or n]')
            if input =='y':
                self.all_tthds = h5py.File(tthds_path,read_mod)
            else:
                self.all_tthds = h5py.File(tthds_path,'r')
        else:
            self.all_tthds = h5py.File(tthds_path,read_mod)
        
        # dictionary to store all nearest tthd vectors
        nearest_tthds_path =  os.getcwd()+'/output_data/nearest_tthds_'+run_name+'.hdf5'
        if os.path.isfile(nearest_tthds_path) and read_mod !='r':
            input = raw_input('nearest_tthds database already exists. are you sure you want to append to/write it? [y or n]')
            if input =='y':
                self.nearest_tthds = h5py.File(nearest_tthds_path,read_mod)
            else:
                self.nearest_tthds = h5py.File(nearest_tthds_path,'r')
        else:
            self.nearest_tthds = h5py.File(nearest_tthds_path,read_mod)

        # pdb file to store all nearest tthd coordinates
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
        
        
        

    def single_frame_N_QR(self,Q,R,pair_dist,time_dependent = False):
        """
        See equations (1a) and (8) in ref [2]
        """
        if Q == 0:
            # np.sin(Q*this_pair)/(Q*this_pair) = 1.0 in the limit of Q -> 0
            if time_dependent:
                return 1.0/self.n_waters*sum([ 1.0 \
                for this_pair in pair_dist if this_pair <= R])
            else:
                # I disgree with equation (1a) and think there must be a typo
                return 2.0/self.n_waters*sum([1.0 \
                for this_pair in pair_dist if this_pair <= R]) \
                + 1
                # + self.n_waters -> this is what equation (1a) in ref [2] thinks
        else:
            # See equations (1a) and (8) in ref [2]
            if time_dependent:
                return 1.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
                for this_pair in pair_dist if this_pair <= R])
            else:
                # I disgree with equation (1a) and think there must be a typo
                return 2.0/self.n_waters*sum([np.sin(Q*this_pair)/(Q*this_pair) \
                for this_pair in pair_dist if this_pair <= R]) \
                + 1
                # + self.n_waters -> this is what equation (1a) in ref [2] thinks
    

    def struct_factor(self,Q,R,dt):
        """
        See equations (6) and (4) in ref [2]
        """
        water_pairs = np.array(list(combinations(sorted(self.water_inds),2)))
        step = int(dt/self.time_step)
        frame_steps= [step+random.randint(-int(step/2),int(step/2)) \
        for ii in range(int(self.total_time/dt))]
        
        while sum(frame_steps)>self.n_frames-1:
            frame_steps = frame_steps[:-1]
            
        # print ('the average time between each configuration is %f ps' % (np.mean(frame_steps)*self.time_step))
        
        N_QR = []
        current_frame = 0
        for this_step in frame_steps:
            current_frame += this_step
            water_dist = md.compute_distances(self.traj[current_frame],water_pairs)[0] # unit in nm
            N_QR.append(self.single_frame_N_QR(Q,R,water_dist))
        
        # equation (4) in ref [2]
        err_Sn_QR = np.std(N_QR)/len(N_QR)
        
        if Q == 0:
            # equation (7) in ref [2]
            Sn_QR = np.mean(N_QR)-4./3.*np.pi*self.rho*R**3.0
            
        else:
            # equations (4a-d) in ref [1]
            Sn_QR = np.mean(N_QR)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return Sn_QR, err_Sn_QR
            
    

    def scat_func(self,Q,R,t):
        """
        Very approximate
        """
        water_pairs = np.array(list(product(sorted(self.water_inds),repeat = 2)))
        
        frame_step = int(t/self.time_step)
        xyz_pos=self.traj.xyz
        
        N_QRt = []
        
        # I using all frames here in the trajectory. 
        # Caution: they may not be statistically independent
        for ii in range(self.n_frames):
            if (ii+frame_step)<self.n_frames:
                water_dist = []
                for this_pair in water_pairs:
                    water_dist.append(np.sqrt(np.sum((xyz_pos[ii,this_pair[0],:]- \
                                        xyz_pos[ii+frame_step,this_pair[1],:])**2.0)))
            
                N_QRt.append(self.single_frame_N_QR(Q,R,water_dist,time_dependent = True))
            else:
                print "Frame number %d and beyond not include in I_n(Q,R,t) calculations." % (ii+1)
                break
                
        
        if Q == 0:
            # equation (7) in ref [2]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*R**3.0
            
        else:
            # equations (28a-b) in ref [1]
            In_QRt = np.mean(N_QRt)-4./3.*np.pi*self.rho*3./Q**3.*(np.sin(Q*R)-Q*R*np.cos(Q*R))
        
        return In_QRt
        
    def radial_dist(self, r_range, n_bins = 200):
        """
        Computer Simulation of Liquids book
        """
        water_pairs = np.array(list(combinations(sorted(self.water_inds),2)))
        water_dist = md.compute_distances(self.traj,water_pairs) # unit in nm
        
        g2s = []
        bin_edges = np.linspace(r_range[0],r_range[1],n_bins)
        for this_frame in water_dist:
            n_ideal = 4.* np.pi/3.*(bin_edges[1:]**3.0-bin_edges[:-1]**3.0)
            n_hist, _ = np.histogram(this_frame, bins = bin_edges,density=True)
            this_g2 = n_hist/n_ideal
            g2s.append(this_g2)

        g2s = np.array(g2s)
        bin_edges = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ave_g2 = np.mean(g2s,axis=0)
        g2_err = np.std(g2s,axis=0)/(ave_g2[-1]) # normalize assuming g2(r) converges to 1.0.
        ave_g2 = ave_g2/(ave_g2[-1])
        
        self.rdf = [bin_edges,ave_g2,g2_err]
    
    def estimate_struct_factor(self,Qs,R_inf,dt_ind):
        """
        Equations (22) and (23) in ref [1]
        """
        S_Q = []
        S_Qerr = []
        
        for this_Q in Qs:
            this_Sn_Q = self.struct_factor(this_Q,R_inf,dt_ind)
            this_Sn_0 = self.struct_factor(0,R_inf,dt_ind)
            if this_Q == 0:
                S_Q.append(this_Sn_Q[0])
                S_Qerr.append(this_Sn_Q[1])
            else:    
                this_S = this_Sn_Q[0] \
                + this_Sn_0[0]/(1-1/self.n_waters*4./3.*np.pi*self.rho*R_inf**3.0) \
                * 1./self.n_waters*4./3.*np.pi*self.rho \
                *3./this_Q**3.*(np.sin(this_Q*R_inf)-this_Q*R_inf*np.cos(this_Q*R_inf))
                S_Q.append(this_S)
                S_Qerr.append(np.sqrt(this_Sn_Q[1]**2.0 + this_Sn_0[1]**2.0)) # estimate of error
        
        self.ssf = [np.array(Qs), np.array(S_Q), np.array(S_Qerr)]

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
    
    def make_pairs(self,vertex_ind,cut_off,frame_ind):
        """Returns a list of vectors that connect pairs of water molecules that are 
        within a cutoff distances of one water molecules
        
        Parameters
        ----------
        vertex_ind : int
            index of the water molecule whose neighbors one want to find and form pairs with
        cut_off : float
            distance in nm within which to search for neighbors
        frame_ind : int
            index of simulation in which to make pairs
        
        """
        nbs = md.compute_neighbors(self.traj[frame_ind],cut_off,[vertex_ind],haystack_indices = self.water_inds)[0]

        pairs = []
        xyz_pos = self.traj[frame_ind].xyz
        for this_nb in nbs:
            
            r_ij = xyz_pos[0,vertex_ind,:]- xyz_pos[0,this_nb,:]
            pairs.append(r_ij)
        return pairs 
        
    def two_point_struct_factor(self,q1,*args):
        sum = 0
        sin_sum = 0
        
        q_norm = np.linalg.norm(q1)
        for this_water in self.water_inds:
            pairs = self.make_pairs(this_water,*args)
            for pp in pairs:
                r_norm = np.linalg.norm(pp)
                if q_norm == 0:
                    sin_sum += 1
                else:
                    sin_sum += np.sin(q_norm*r_norm)/(q_norm*r_norm)
                sum += np.exp(1j*np.sum(pp*q1))
         
        return sum/self.n_waters+1, sin_sum/self.n_waters + 1
        
    def make_frame_inds(self,dt):
        """
        returns indices of frames for average a quantity over
        """
        step = int(dt/self.time_step)
        frame_steps= [step+random.randint(-int(step/2),int(step/2)) \
        for ii in range(int(self.total_time/dt))]
        
        while sum(frame_steps)>self.n_frames-1:
            frame_steps = frame_steps[:-1]
        
        return np.cumsum(frame_steps)

        
    def compute_term_four_point(self,q_pair,tthd_pair):
        return (1+np.cos(np.dot(q_pair[0],tthd_pair[0])))*(1+np.cos(np.dot(q_pair[1],tthd_pair[1])))
    
    def four_point_struct_factor(self,qs,cut_off,frame_ind,nearest_nb,test_dataset=True):
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
        if test_dataset:
            print frame_ind
            database = h5py.File(self.test_datapath,'r')
            this_tthds = database[self.run_name][frame_ind][:][:]
            database.close()
        else:
            if nearest_nb:
                if str(frame_ind) in self.nearest_tthds:
    #                 print "recycling for frame %d!" % frame_ind
                    this_tthds = self.nearest_tthds[str(frame_ind)][:][:]
                else:
                    this_tthds = ws.make_nearest_nb_tthds(cut_off,frame_ind)
                    self.nearest_tthds.create_dataset(str(frame_ind),data = this_tthds)
            else:
                if str(frame_ind) in self.all_tthds:
    #                 print "recycling for frame %d!" % frame_ind
                    this_tthds = self.all_tthds[str(frame_ind)][:][:]
                else:
                    this_tthds = []
                    for this_water in self.water_inds:
                        this_tthds.extend(self.make_tthd(this_water,cut_off,frame_ind))
                    self.all_tthds.create_dataset(str(frame_ind),data = this_tthds)
        
        corr_single_frame = []
        aa = [3.0485,2.2868,1.5463,0.867]
        bb = [13.2771,5.7011,0.3239,32.9089]
        cc = 0.2580
        form_factor = cc
        for this_a,this_b in zip(aa,bb):
            form_factor += this_a * np.exp(-this_b*(np.linalg.norm(qs[0])/(4*np.pi))**2.0)
        
        for this_q in qs:
            this_sum = 0
            
            for tt in this_tthds:             
                # derived new formula, ingnoring form factor for now for constant q
                this_sum += self.compute_term_four_point(this_q,tt)
            
            n_tthds = len(this_tthds)
            
            
            
            corr_single_frame.append(this_sum/n_tthds*form_factor**4.0*4.0)
        
        return corr_single_frame
        
        
     #    sum = 0
#         if str(frame_ind) in self.all_tthds:
#             print "recycling for frame %d!" % frame_ind
#             this_tthds = self.all_tthds[str(frame_ind)]
#             for tt in this_tthds:                
#                 # derived new formula, ingnoring form factor for now for constant q
#                 
#                 sum += (1+np.cos(np.dot(q1,tt[0])))*(1+np.cos(np.dot(q2,tt[1])))
#         else:
#             this_tthds = []
#             for this_water in self.water_inds:
#                 this_tthds.extend(self.make_tthd(this_water,cut_off,frame_ind))
#             self.all_tthds.create_dataset(str(frame_ind),data = this_tthds)
# 
#             for tt in this_tthds:
#                 # derived new formula, ingnoring form factor for now for constant q
#                 sum += (1+np.cos(np.dot(q1,tt[0])))*(1+np.cos(np.dot(q2,tt[1])))
        # aa = [3.0485,2.2868,1.5463,0.867]
#         bb = [13.2771,5.7011,0.3239,32.9089]
#         cc = 0.2580
#         form_factor = cc
#         for this_a,this_b in zip(aa,bb):
#             form_factor += this_a * np.exp(-this_b*(np.linalg.norm(q1)/(4*np.pi))**2.0)
#          
#         return 1/self.n_waters**2.0*sum #*form_factor**4.0*4.
    
    def atomic_form_factor(self,element):
        """computes the atomic form factor givent the element
        """
        pass

    def correlator(self,q,wavelength,frames,phi,cut_off = 0.5,nearest_nb=True,output = None,test_dataset = True):
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
 
        # print "frames used for averaging..."
#         print frames
        
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
            csvwriter.writerow(frames)
        csvfile.close()

        
        for this_frame in frames:
#             print('computing for frame number %d...' % this_frame)
            print this_frame
            this_row = self.four_point_struct_factor(qs,cut_off,this_frame,nearest_nb,test_dataset=test_dataset)
            with open(output_path,'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=' ',
                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                csvwriter.writerow(this_row)
                csvfile.flush()
        
 ################################################################################       
        
#          outfile = open('C(psi).txt','a')
# 
#         
#         q1 = np.array([np.sin(2*theta_1),0,np.cos(2*theta_1)])*q
#         
#         S_q = []
#         S_qerr = []
#         psi = []
#         
#         phi = np.linspace(-np.pi,np.pi,10)
#         
#         for this_phi in phi:
#             print "calculating for phi = %.2f" % this_phi
#             q2 = np.array([q1[0]*np.cos(this_phi),q1[0]*np.sin(this_phi),q1[2]])
#             sf = [self.four_point_struct_factor(q1,q2,cut_off,this_fr) for this_fr in frames]
#             
#             this_Sqerr = np.std(sf)/np.sqrt(len(sf))
#             S_qerr.append(this_Sqerr)
#             
#             this_Sq = np.mean(sf)
#             S_q.append(this_Sq)
#             
#             this_psi = np.arccos(np.dot(q1/q,q2/q))
#             psi.append(this_psi)
#             
#             outfile.write("%g,%g,%g,%g" % (this_Sq,this_Sqerr,this_psi,this_phi)+"\n")
#             outfile.flush()
#         
#         outfile.close()
#         
#         return np.array(S_q),np.array(S_qerr),np.array(psi),phi
    
    
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
        
        tthds = []
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
            r_ij = r1-r2
            r_kl = r3-r4 # nm
            
            tthds.append([r_ij,r_kl])
            
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
        return tthds 
        
                
##############################################################################
# test
##############################################################################

# data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
# data_path = '/home/shenglan/GitHub/waterMD/data'
# traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
# print ('here is some info about the trajectory we are looking at:')
# print traj
# test = WaterStats(traj)
# 
# 
# this_phi = np.pi/4.
# q = 1/0.3*np.pi*2.0
# theta_1 = np.pi/12.
# q1 = np.array([np.sin(2*theta_1),0,np.cos(2*theta_1)])*q
# q2 = np.array([q1[0]*np.cos(this_phi),q1[0]*np.sin(this_phi),q1[2]])
# print np.arccos(np.dot(q1/q,q2/q))
# 
# print test.four_point_struct_factor(q1,-q2,0.5,10)
# 
# this_phi = np.pi/3
# q2 = np.array([q1[0]*np.cos(this_phi),q1[0]*np.sin(this_phi),q1[2]])
# print test.four_point_struct_factor(q1,q2,0.5,10)

# 
# # q1 = []
# print test.four_point_struct_factor(q1,2*q1,0.5,10)
# 
# vec2 = test.four_point_struct_factor(q1,-q1,0.5,10)
# print vec2
# print "it's magnitude is %g" % np.abs(vec2)
# vec3 = test.four_point_struct_factor(q1,-q1,0.5,10)
# print vec3
# print "it's magnitude is %g" % np.abs(vec3)

# frames = test.make_frame_inds(8.0)
# 
# for this_frame in frames
# print test.two_point_struct_factor(q1,0.5,10)
# print test.two_point_struct_factor(q1,0.5,20)