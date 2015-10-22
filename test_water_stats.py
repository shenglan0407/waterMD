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

data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj)

R_water = 0.3

output_path = '/Users/shenglanqiao/Documents/GitHub/waterMD/output'

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

def test_estimate_sf(Qs,R_max,dt,plot_Sn_QR=True):
    test.estimate_struct_factor(Qs,R_max,dt)
    
    Sn_QR = []
    Qs,S_Q,S_Qerr=np.array(test.ssf[0]),test.ssf[1],test.ssf[2]
    fig = plt.figure()
    plt.errorbar(Qs*R_water/(2*np.pi),S_Q,yerr=S_Qerr,label='S(Q)')
    
    
    
    plt.title("Estimate of S(Q) with dt = %.1f ps and R_inf = %.1f" % (dt,R_max))
    plt.xlabel("Q*R_water (rad)")
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


##############################################################################
# test
##############################################################################

Rs = np.linspace(0.2,0.4,10)

R_max = 0.5 # nm
dt = 8.0 # ps
Qs = 2.*np.pi*np.linspace(0.0,3/R_water,50)

ts = np.linspace(1,10,10)

test_estimate_sf(Qs,R_max,dt)
 


