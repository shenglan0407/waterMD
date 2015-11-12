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

run_name = 'run3'
data_path = os.getcwd()+'/data'
traj = md.load_trr(data_path+'/nvt-pr_'+run_name+'.trr', top = data_path+'/water-sol_'+run_name+'.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj,run_name,read_mod='r')
test.all_tthds.close()

# run 3 stuff
all_data = []
with open(os.getcwd()+'/computed_results/combined_corr_'+run_name+'.csv','r') as f:
        csvreader = csv.reader(f,delimiter=' ',quotechar='|')
        q1=(next(csvreader))
        q2=(next(csvreader))
        psi=(next(csvreader))
        
        for row in csvreader:
            all_data.append(row)

#run 2 data
run_name = 'run2'

with open(os.getcwd()+'/computed_results/combined_corr_'+run_name+'.csv','r') as f:
        csvreader = csv.reader(f,delimiter=' ',quotechar='|')
        q1=(next(csvreader))
        q2=(next(csvreader))
        psi=(next(csvreader))
        
        for row in csvreader:
            all_data.append(row)

all_data = np.array(all_data,dtype=float)


Corr = [np.mean(all_data[:,ii])/test.n_waters**2.0 for ii in range(len(psi))]
C_err = [np.std(all_data[:,ii])/test.n_waters**2.0/np.sqrt(len(all_data[:,ii])) for ii in range(len(psi))]

#run 1 data

CorrMat = np.loadtxt(os.getcwd()+'/computed_results/C(psi)_run1.txt',skiprows=5,delimiter=',')

Corr1 = CorrMat[:,0]
C_err1 = CorrMat[:,1]

new_Corr = (np.array(Corr)*200.+Corr1*100)/300
GSS = (Corr-new_Corr)**2.0*200+(Corr1-new_Corr)**2.0*100
C_err = (np.array(C_err)*200+C_err1*100+GSS)/300

print np.mean(C_err)

fig = plt.figure()
plt.errorbar(psi,new_Corr,yerr = C_err)
plt.title('4-point correlator')
plt.xlabel('psi')
plt.ylabel('C')
fig.savefig('Corr_all.png')
plt.close(fig)