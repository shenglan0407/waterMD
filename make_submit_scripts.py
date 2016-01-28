import numpy as py
import os

q_invs_float = [0.292,0.299,0.306,
                0.314,0.322,0.331,0.34,
                0.35,0.36,0.37,0.38,0.393,
                0.405,0.419,0.433,0.449,0.465]
                
with open(os.getcwd()+'/barley_scripts/submitBatch.sh','w') as submit:
    for this_q in q_invs_float:
        sub_file = 'job_q'+str(this_q)+'.sh'
        with open(os.getcwd()+'/barley_scripts/'+sub_file,'w') as this_file:
            this_file.write('#!/bin/bash\n')
            this_file.write('#$ -N corr_q'+str(this_q)+'\n')
            this_file.write('#$ -wd /srv/zfs01/user_data/shenglan/waterMD/\n')
            this_file.write('#$ -M shenglan@stanford.edu\n')
            this_file.write('#$ -m besa\n')
            this_file.write('#$ -j y\n')
            this_file.write('#$ -o /srv/zfs01/user_data/shenglan/waterMD/barley_scripts/out_q'+str(this_q)+'.log\n')
            this_file.write('python compute_correlator.py -i run7 -q '+str(this_q)+' -p 100 -s all\n')
            this_file.write('hostname\n')
        submit.write('qsub '+sub_file+'\n')
