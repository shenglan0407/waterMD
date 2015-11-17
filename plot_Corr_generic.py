##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import os
import csv

from water_stats import WaterStats
import mdtraj as md

##############################################################################
# Code
##############################################################################

run_name = 'run1'
data_path = os.getcwd()+'/data'
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj,run_name,read_mod='r')
test.all_tthds.close()
test.nearest_tthds.close()

run_name = 'all'
all_data = []
with open(os.getcwd()+'/computed_results/combined_corr_all.csv','r') as f:
        csvreader = csv.reader(f,delimiter=' ',quotechar='|')
        q1=(next(csvreader))
        q2=(next(csvreader))
        psi=(next(csvreader))
        count = 0
        for row in csvreader:
            count += 1
            if count%1==0:
                all_data.append(row)

all_data = np.array(all_data,dtype=float)



Corr = [np.mean(all_data[:,ii]) for ii in range(len(psi))]
C_err = [np.std(all_data[:,ii])/np.sqrt(len(all_data[:,ii])) for ii in range(len(psi))]


fig = plt.figure()
plt.errorbar(psi,Corr,yerr = C_err,fmt='-.')
plt.title('4-point correlator')
plt.xlabel('psi')
plt.ylabel('C')
fig.savefig('Corr_fine_'+run_name+'_all.png')
plt.close(fig)