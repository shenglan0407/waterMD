##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
# 
# stitch together correlators computed from the same run
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

run_name = 'run8'
qs = [0.23,0.292,0.306,0.322,0.34,0.36,0.38,0.405,0.433,0.465]
frame_starts = [1000,5000]
for this_q in qs:
    file_list = ['corr_'+run_name+'_'+str(this_q)+'q_30p.csv']
    for this_start in frame_starts:
        file_list.append('corr_'+run_name+'_'+str(this_q)+'q_30p_'+str(this_start)+'.csv')
#     file_list = ['q0.465_formfactor_1.csv','q0.465_formfactor_2.csv']
    data_list = []
    q1_list =[]
    q2_list =[]
    phi_list=[]
    frame_list = []
    all_data = []
    for this_file in file_list:
        with open(os.getcwd()+'/computed_results/'+this_file,'r') as f:
            csvreader = csv.reader(f,delimiter=' ',quotechar='|')
            q1_list.append(next(csvreader))
            q2_list.append(next(csvreader))
            phi_list.append(next(csvreader))
            frames = next(csvreader)

            computed_data = []
            for row in csvreader:
                computed_data.append(row)
            data_list.append(computed_data)
            n_frames = len(computed_data)
            frame_list.append([int(frames[0]),int(frames[n_frames-1])])
    frame_list= np.array(frame_list)

    for ii in range(len(file_list)-1):
        if ((q1_list[ii] == q1_list[ii+1]) and (q2_list[ii] == q2_list[ii+1])) and (phi_list[ii] == phi_list[ii+1]):
            pass
        else:
            print "Parameters in this list of files are not equal. Check the black sheep %s."%file_list[ii+1]
            sys.exit()

    order = np.argsort(frame_list[:,0])

    for ii in range(len(order)-1):
        if frame_list[order[ii]][1]>=frame_list[order[ii+1]][0]:
            begin = frame_list[order[ii]][0]
            end=frame_list[order[ii+1]][0]-begin
            all_data.extend(data_list[order[ii]][:end])
        else:
            all_data.extend(data_list[order[ii]])

    all_data.extend(data_list[order[-1]])

    with open(os.getcwd()+'/computed_results/combined_corr_'+run_name+'_'+str(this_q)+'q_30p.csv','w') as f:
        csvwriter = csv.writer(f, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow(q1_list[0])
        csvwriter.writerow(q2_list[0])
        csvwriter.writerow(phi_list[0])
        for row in all_data:
            csvwriter.writerow(row)
    
    print("For q_inverse = %.3g \
    total number of simulation frame is %d. check if this is correct."%(this_q,len(all_data)))
        
