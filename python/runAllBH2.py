#!/usr/bin/python3

import os, sys
from CorrItems import allBH2vars

samples = ['y0_1000GeV', 'zp_1000GeV', 'xdxdtt_100GeV']

scriptname = '_run_all_BH2.sh'
os.system('rm -f {}'.format(scriptname))
scriptfile = open(scriptname, 'w')


for sample in samples:
    listfile = 'lists/list_{}.txt'.format(sample)
    for var in allBH2vars:
        scriptfile.write('nice -n 10 ./python/XsectStackReplicas2D.py {} {} >& logs/log_bh2_{}_{}'.format(listfile, var, sample, var) + '\n')

scriptfile.close()
print('./scripts/jobStarter.sh --jobs 8 {}'.format(scriptname))
