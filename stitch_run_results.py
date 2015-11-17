##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# stitch together correlators computed from differen simulation runs
#
#############################################################################


##############################################################################
# Imports
##############################################################################
import numpy as np
import os
import csv
import sys
##############################################################################
# Code
##############################################################################

runs = ['run1','run2','run3','run4']
file_list = ['/corr_'+run_name+'.csv' for run_name in runs]

data_list = []
q1_list =[]
q2_list =[]
phi_list=[]
total_frames = 0
all_data = []
for this_file in file_list:
    with open(os.getcwd()+'/computed_results/'+this_file,'r') as f:
        csvreader = csv.reader(f,delimiter=' ',quotechar='|')
        q1_list.append(next(csvreader))
        q2_list.append(next(csvreader))
        phi_list.append(next(csvreader))
        frames = next(csvreader)

        for row in csvreader:
            total_frames += 1
            all_data.append(row)
        
print 'Total frame number is %d. check if it is correct.' % total_frames
# check is all the parameters used for calculations are equal
for ii in range(len(file_list)-1):
    if ((q1_list[ii] == q1_list[ii+1]) and (q2_list[ii] == q2_list[ii+1])) and (phi_list[ii] == phi_list[ii+1]):
        pass
    else:
        print "Parameters in this list of files are not equal. Check the black sheep %s."%file_list[ii+1]
        sys.exit()
        
with open(os.getcwd()+'/computed_results/combined_corr_all.csv','w') as f:
    csvwriter = csv.writer(f, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    csvwriter.writerow(q1_list[0])
    csvwriter.writerow(q2_list[0])
    csvwriter.writerow(phi_list[0])
    for row in all_data:
        csvwriter.writerow(row)