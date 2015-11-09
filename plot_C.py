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

##############################################################################
# Code
##############################################################################

CorrMat = np.loadtxt(os.getcwd()+'/output_data/C(psi).txt',skiprows=1,delimiter=',')

Corr = CorrMat[:,0]
C_err = CorrMat[:,1]
cos_psi = CorrMat[:,2]
psi = CorrMat[:,3]

fig = plt.figure()
plt.errorbar(psi,Corr,yerr = C_err)
plt.title('4-point correlator')
plt.xlabel('psi')
plt.ylabel('C')
fig.savefig('Corr.png')
plt.close(fig)