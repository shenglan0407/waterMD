import mdtraj as md
import numpy as np
from water_stats import WaterStats

import time

run_name = 'run7'
data_path ='/home/shenglan/MD_simulations/water_box/cubic_1nm_'+run_name
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj,'run7',read_mod = 'r')

q=1/0.465*np.pi*2.0
phis = np.linspace(0,np.pi,30)

tic = time.clock()
test.correlator(q,0.1,['tthdSet1','tthdSet2','tthdSet3'],phis,output='run7_corr_test_q0.465.csv',test_dataset=True)
toc=time.clock()
print 'time taken %g second' % (toc-tic)