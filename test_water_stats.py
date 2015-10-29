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
import matplotlib.pyplot as plt

##############################################################################
# Code
##############################################################################

# data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
data_path = '/home/shenglan/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj)

R_water = 0.3

# output_path = '/Users/shenglanqiao/Documents/GitHub/waterMD/output'
output_path = '/home/shenglan/GitHub/waterMD/output'

def test_rdf(r_range):
    test.radial_dist(r_range)
    rs, g_R, g_err = test.rdf[0],test.rdf[1],test.rdf[2]
    fig = plt.figure()
    plt.errorbar(rs,g_R, yerr=g_err)
    plt.title('gn(r)')
    plt.xlabel('r (nm)')
    plt.ylabel('gn(r)')
    plt.ylim(0,max(g_R)*1.2)
    fig.savefig(output_path + '/gn_r.png')
    plt.close(fig)

def test_Sn_0R(Rs,dt):
    Sn_0R = []
    for R in Rs:
        Sn_0R.append(test.struct_factor(0,R,1)[0])
    fig = plt.figure()
    plt.plot(Rs, Sn_0R)
    plt.title("Sn(0,R) with dt = 1.0 ps")
    plt.xlabel("R (nm)")
    plt.ylabel("Sn(0,R)") 
    fig.savefig(output_path+'/Sn_0R.png')
    plt.close(fig)
    
def test_In_0tR(R1,R2,tt):
    In_0tR1 = []
    In_0tR2 = []
    R1 = 0.3 # nm
    R2 = 0.5 # nm

    for tt in ts:
        In_0tR1.append(test.scat_func(0,R1,tt))
        In_0tR2.append(test.scat_func(0,R2,tt))
    
    fig = plt.figure()
    plt.plot(ts, In_0tR1-min(In_0tR1),label = "R = %.2f" % R1)
    plt.plot(ts, In_0tR2-min(In_0tR2),label = "R = %.2f" % R2)
    plt.legend()
    plt.title("In(0,t,R)")
    plt.xlabel("t (ps)")
    plt.ylabel("In(0,t,R)") 
    fig.savefig(output_path+'/In_0tR.png')
    plt.close(fig)

def test_estimate_sf(Qs,R_max,dt,plot_Sn_QR=False):
    test.estimate_struct_factor(Qs,R_max,dt)
    
    Sn_QR = []
    Qs,S_Q,S_Qerr=np.array(test.ssf[0]),test.ssf[1],test.ssf[2]
    fig = plt.figure()
    plt.errorbar(Qs*R_water/(2*np.pi),S_Q,yerr=S_Qerr,label='S(Q)')
    
    
    
    plt.title("Estimate of S(Q) with dt = %.1f ps and R_inf = %.1f" % (dt,R_max))
    plt.xlabel("Q*R_water")
    plt.ylabel("S(Q)")
    plt.xlim(0,max(Qs*R_water/(2*np.pi)))
    if plot_Sn_QR:
        for Q in Qs:
            Sn_QR.append(test.struct_factor(Q,R_max,dt)[0])
        plt.plot(Qs*R_water/(2*np.pi),Sn_QR,'go',label='Sn(Q) (no finite size effect)')
        plt.legend(loc=4)
    fig.savefig(output_path+'/S_Q.png')
    plt.close(fig)

def test_Nbar(Qs,R_max):
    """
    The second term contribution to S(Q) that contains the bessel function
    See equations (22) and (23) of ref [1]
    """
    fig = plt.figure()
    N_bar=[]
    for Q in Qs:
        if Q == 0:
            N_bar.append(0)
        else:
            N_bar.append(-4./3.*np.pi*test.rho*3./Q**3.*(np.sin(Q*R_max)-Q*R_max*np.cos(Q*R_max)))
    plt.plot(np.array(Qs)*R_water/(2.*np.pi),N_bar)
    plt.xlabel('Q*R_water')
    fig.savefig(output_path+'S_Q_BessTerm.png')
    plt.close(fig)

def test_Sn_QR(Qs,R_max,dt):
    Sn_QR=[]
    for Q in Qs:
        Sn_QR.append(test.struct_factor(Q,R_max,dt)[0])
    
    fig = plt.figure()
    plt.plot(Qs*R_water/(2*np.pi),Sn_QR,'o')
    plt.title("Sn(Q) with dt = %.1f ps and R_inf = %.1f" % (dt,R_max))
    plt.xlabel("Q*R_water")
    plt.ylabel("Sn(Q)")
    plt.xlim(0,max(Qs*R_water/(2*np.pi)))
    fig.savefig(output_path+'/Sn_QR.png')
    plt.close(fig)

def test_two_point_ft(Qs,R_max,dt=20.0):
    # I want to check I get the same result for the 2-point strcture factor if I try to do the 
    # fourier transform explicitly. So far I have observed that at least the imaginary part
    # is effectively zero. Next I want to check that for a particular length q only the 
    # magnitude of the vector q matters but not the direction, i.e. two q vectors with
    # the same length should return the same answer for the summation term in Sn(Q) 
    frames = test.make_frame_inds(dt)
    print frames
    
    unit_vec = np.array([1,1,1])/np.sqrt(3.)
    Sn_Q1 = []
    Sn_Q2 = []
    for Q in Qs:
        print "calculating for %.2f" % Q
        this_q = Q * unit_vec
        sf2 = [test.two_point_struct_factor(this_q,R_max,this_frame)[0] for this_frame in frames]
        Sn_Q2.append(np.mean(sf2))
        
        this_q = Q * np.array([1,-1,1])/np.sqrt(3.)
        sf1 = [test.two_point_struct_factor(this_q,R_max,this_frame)[0] for this_frame in frames]
        Sn_Q1.append(np.mean(sf2))
        
    Sn_Q2 = np.array(Sn_Q2)
    Sn_Q1 = np.array(Sn_Q1)
    fig = plt.figure()
    plt.plot(Qs,np.abs(Sn_Q2),label='q')
    plt.plot(Qs,np.abs(Sn_Q1),label='-q')
    plt.legend(loc=4)
    plt.xlabel('Q')
    plt.ylabel('sum of fourier terms')
    fig.savefig(output_path+'/test_ft.png')
    plt.close(fig)
    
def test2_two_point_ft(Qs, R_max, dt = 20.0):
    frames = test.make_frame_inds(dt)
    sin_term = []
    
    unit_vec = np.array([1,1,1])/np.sqrt(3.)
    Sn_Q = []
    for Q in Qs:
        print "calculating for %.2f" % Q
        this_q = Q * unit_vec
        sf = [test.two_point_struct_factor(this_q,R_max,this_frame) for this_frame in frames]
        this_mean = np.mean(sf,axis = 0)
        Sn_Q.append(this_mean[0])
        sin_term.append(this_mean[1])
        
    fig, (ax1,ax2) = plt.subplots(1,2,sharex = True,sharey = True)
    ax1.plot(Qs,np.abs(Sn_Q))
    
    ax2.plot(Qs,sin_term)
    ax1.set_title('explicit ft')
    ax2.set_title('Debye formula')
    # plt.legend(loc=4)
#     plt.xlabel('Q')
#     plt.ylabel('sum of fourier terms')
    fig.savefig(output_path+'/test2_ft.png')
    plt.close(fig)
        
def test_corr(q,theta_1,dt,cut_off = 0.5,return_three=False):
    S_q,S_qerr,psi,phi = test.correlator(q,theta_1,dt,cut_off = 0.5,return_three=False)
    
    fig = plt.figure()
    plt.errorbar(phi,S_q,yerr=S_qerr)
    plt.plot(phi,S_q,'--')
    plt.xlabel("phi (rad)")
    plt.ylabel("C(psi)")
    plt.title("C(q1,q2,psi) with q1 = q2 = q")
    # fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)
#     ax1.plot(psi,S_q,'o')
#     ax1.set_title('Intensity vs psi')
#     ax1.set_xlabel('psi')
#     ax1.set_ylabel('Intensity (a.u.)')
#     ax2.plot(phi,S_q,'o')
#     ax2.set_title('Intensity vs phi')
#     ax2.set_xlabel('phi')
    fig.savefig(output_path+'/corr.png')
    plt.close(fig)

##############################################################################
# test
##############################################################################

Rs = np.linspace(0.2,0.4,10)

R_max = 0.5 # nm
dt = 8.0 # ps
Qs = 2.*np.pi*np.linspace(0.0,1.5/R_water,10)

ts = np.linspace(1,10,3)

q = 1/0.3*np.pi*2.0
theta_1 = np.pi/12.


test_corr(q,theta_1,dt)
