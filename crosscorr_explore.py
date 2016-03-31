##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# this script explores the cross-correlator
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import matplotlib as mpl
from matplotlib import mlab
import h5py

import time
##############################################################################
# Code
##############################################################################

def get_angle(r1,r2):
    norm1=np.linalg.norm(r1)
    norm2=np.linalg.norm(r2)
    return np.arccos(np.dot(r1,r2)/norm1/norm2)/np.pi*180

# load data
f= h5py.File('output_data/test_tthd_data.hdf5','r')
run13_set1 = f["run13"]["tthdSet1"][:]
run13_set2 = f["run13"]["tthdSet2"][:]
run13_set3 = f["run13"]["tthdSet3"][:]
f.close()

# generate q-vectors

wavelength = 0.1544
q1_mag = 17.0
q_beam = 2.0*np.pi/wavelength
phi = np.linspace(0*np.pi,1.0*np.pi,30)
#         
# # generate q1 and many q2s
# q1 = np.array([q1_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0),0,-q1_mag**2.0/(2.0*q_beam)])
# q2 = np.array([np.array([q2_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0)*np.cos(this_phi),
#                          q2_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0)*np.sin(this_phi),
#                          -q1_mag*q2_mag/(2.0*q_beam)]) for this_phi in phi])
# # make pairs of q1 and q2
# qs = np.array([[q1,this_q2] for this_q2 in q2])
# psi = [get_angle(q_pair[0],q_pair[1]) for q_pair in qs]
# 


###########################################
q2_mag_range =np.arange(16,16.5,0.5) 

c_3p=np.zeros((len(q2_mag_range),len(phi)))
c_4p=np.zeros((len(q2_mag_range),len(phi)))
qind=0

tic=time.clock()
for q2_mag in q2_mag_range:
    
    # generate q1 and many q2s
    q1 = np.array([q1_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0),0,-q1_mag**2.0/(2.0*q_beam)])
    q2 = np.array([np.array([q2_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0)*np.cos(this_phi),
                             q2_mag*np.sqrt(1-(q1_mag/(2.*q_beam))**2.0)*np.sin(this_phi),
                             -q1_mag*q2_mag/(2.0*q_beam)]) for this_phi in phi])
    # make pairs of q1 and q2
    qs = np.array([[q1,this_q2] for this_q2 in q2])
    psi = [get_angle(q_pair[0],q_pair[1]) for q_pair in qs]
    c_tri=[]
    c_tet=[]
    for q_pair in qs:
        cc3=[]
        cc4=[]
        for ind in range(len(run13_set1)):
            r1 = run13_set1[ind][0]
            r2 = run13_set2[ind][1]
            cos1=np.cos(np.dot(q_pair[0],r1))
            cos2=np.cos(np.dot(q_pair[1],r2))
            cc3.append(cos1*cos2)
            
            
            r3 = run13_set1[ind][1]
            cos3=np.cos(np.dot(q_pair[1],r3))
            cc4.append(cos1*cos3)
        c_tri.append(np.mean(cc3))
        c_tet.append(np.mean(cc4))
    c_3p[qind,:]=c_tri
    c_4p[qind,:]=c_tet
    
    qind+=1

toc=time.clock()
print "Time spend in loop is %g sec"%(toc-tic)

# plot stuff
psi=np.nan_to_num(psi)
for ii in range(len(q2_mag_range)):
    this_q=q2_mag_range[ii]
    this_c3=c_3p[ii,:]
    line_color = cm.jet((max(q2_mag_range)-this_q)/(max(q2_mag_range)-min(q2_mag_range)))
    
    plt.plot(psi,this_c3,color=line_color,label=str(this_q/10))

for ii in range(len(q2_mag_range)):
    this_q=q2_mag_range[ii]
    this_c4=c_4p[ii,:]
    line_color = cm.jet((max(q2_mag_range)-this_q)/(max(q2_mag_range)-min(q2_mag_range)))
    
    plt.plot(psi,this_c4,color=line_color,label=str(this_q/10),linestyle='--')
    
plt.legend()