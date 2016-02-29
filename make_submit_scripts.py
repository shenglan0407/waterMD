import numpy as np
import os

run_name  = 'thor_run7_1'
# output = 'thor_run7_1_autocorr'
stride  = 2
n_q=50
ind_pairs=[(ii,ii+5) for ii in np.arange(0,n_q,5)]

# q_invs_float = [0.292,0.299,0.306,
#                 0.314,0.322,0.331,0.34,
#                 0.35,0.36,0.37,0.38,0.393,
#                 0.405,0.419,0.433,0.449,0.465]
                
with open(os.getcwd()+'/barley_scripts/submitBatch.sh','w') as submit:
    for this_pair in ind_pairs:
        sub_file = 'job_'+str(this_pair[0])+'_'+str(this_pair[1])+'.sh'
        output = 'thor_run7_1_autocorr_'+str(this_pair[0])
        with open(os.getcwd()+'/barley_scripts/'+sub_file,'w') as this_file:
            this_file.write('#!/bin/bash\n')
            this_file.write('#$ -N corr_'+str(this_pair[0])+'\n')
            this_file.write('#$ -wd /srv/zfs01/user_data/shenglan/waterMD/\n')
            this_file.write('#$ -M shenglan@stanford.edu\n')
            this_file.write('#$ -m besa\n')
            this_file.write('#$ -j y\n')
            this_file.write('#$ -o /srv/zfs01/user_data/shenglan/waterMD/barley_scripts/out_'+str(this_pair[0])+'.log\n')
            this_file.write('python thor_compute_autocorr.py -i %s -o %s -s %d -e %d -p %d\n'%(run_name,output,this_pair[0],this_pair[1],stride))
            this_file.write('hostname\n')
        submit.write('qsub '+sub_file+'\n')
