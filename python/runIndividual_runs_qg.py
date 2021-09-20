#!/usr/bin/python3

import os, sys

datadir='/home/public/data/ttbar/data/qitek/DelphesOut/'

tags = [ 
    'gg2ttbarj_allhad_NLO_ATLAS',
    'qg2ttbarj_allhad_NLO_ATLAS',
    'qq2ttbarj_allhad_NLO_ATLAS',

    'pp_2tj_allhad_NLO_nocuts_14TeV_ATLAS',
    'pp_2tj_incl_NLO_nocuts_14TeV_ATLAS'
] 

print('You can run:')
for tag in tags:
    #print('*** Processing {}'.format(tag))
    scriptname = '_run_{}.sh'.format(tag)
    os.system('rm -f {}'.format(scriptname))
    scriptfile = open(scriptname, 'w')
    ftag = (tag + '').replace('_weighted','')
    for line in os.popen('cd {} ; ls | grep {}'.format(datadir,ftag)).readlines():
        item = line[:-1]  
        run = item  + ''
        run = run.replace('out_boosted_AtKt4_and_10_', '').replace('.root', '')
        sample = run + ''
        weight = ''
        if '_weighted' in tag:
            weight = 'weighted'
            sample = sample.replace('_weighted','')
        wtag = ''
        if len(weight) > 0:
            wtag = '_' + weight
        scriptfile.write('nice -n 10 ./bin/AnalyzeBoosted_x {} {} Trimmed {} >& log_{}{}.txt \n'.format(datadir, sample, weight, run,wtag))
    scriptfile.close()
    #print('done!')
    #print('Now you can run:')
    print('./scripts/jobStarter.sh --jobs 8 {}'.format(scriptname))

#print('DONE!')
