#!/usr/bin/python3
# jk III/2020

import os, sys

runsdir = 'runs/'
cmd = 'mkdir -p {} ; mv *run_*.root {}'.format(runsdir, runsdir)
print(cmd)
os.system(cmd)

argv = sys.argv
reqtag = ''
if len(argv) > 1:
    reqtag = argv[1]
    print('OK, will merge samples containing {:}'.format(reqtag))


tags = [

    'pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted',
    'pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_weighted',
    'pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted',

          ]

for tag in tags:
    # HACK!!!
    #if tag != 'zp_ttbarj_allhad_1000GeV_14TeV_ATLAS' and  tag != 'zp_ttbarj_allhad_1250GeV_14TeV_ATLAS' and tag != 'pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS':
    #if tag != 'pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS':
    #if tag != 'zp_ttbarj_allhad_1000GeV_NEW_14TeV_ATLAS':
    #    #'pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV_ATLAS':
    #if tag != 'pp_2tj_allhad_NLO_ptj1min100_ptj2min100_14TeV_ATLAS':
    #if not 'pp_2tj_allhad_NLO_ptj1min200' in tag:
    #    continue
    #if not 'xd_100' in tag:
    #    continue
    #if not '2tj' in tag:
    #    continue
    #if not 'ptj1min200_ptj2min60max200' in tag:
    #    continue

    if reqtag != '' and not reqtag in tag:
        continue
        
    print('***processing {}'.format(tag))
    files = []
    #addGrepCmd = ''
    #if 'xdxd' in tag:
    #    # early xdxd samples have bugged particle jets as xd were included in the ptcl jets clustering
    #    addGrepCmd = ' | egrep "run_3|run_4" | egrep -v "run_30|run_31"'
    addGrepCmd = ''
    if not 'weighted' in tag:
        addGrepCmd = ' | grep -v weighted'
    for line in os.popen('cd {} ; ls *{}.root | grep run_ {}'.format(runsdir, tag, addGrepCmd)).readlines():
        #print(line)
        files.append(line[:-1])
        #print(files)
    print('Files: {}'.format(len(files)))


    flists = [ [], [] ]
    iif = 0
    ihalf = 0
    for fname in files:
        if iif == int(len(files)/2.):
            ihalf += 1
        flists[ihalf].append(fname)
        iif += 1
    print(flists)
    if len(flists[0]) != len(flists[1]):
        print('ERROR, GOT NON EQUAL NUMBER OF FILES IN HALF0 AND HALF1 SAMPLES! {} vs {}'.format(len(flists[0]), len(flists[1])))
    
    ihalf=0
    for flist in flists:
        fline = ''
        for item in flist:
            print(item)
            fline += ' {}{}'.format(runsdir,item)
        cmd = 'hadd -f analyzed_histos_{}_half{}.root {}'.format(tag, ihalf, fline)
        print(cmd)
        os.system(cmd)
        ihalf += 1
    cmd = 'hadd -f analyzed_histos_{}_all.root analyzed_histos_{}_half0.root analyzed_histos_{}_half1.root'.format(tag, tag, tag)
    os.system(cmd)

# SPECIAL: 29.10.2020
# hadd the full weighted ttbar sample!
if reqtag == '' or '2tj' in reqtag or 'weigh' in reqtag:
    os.system('hadd -f analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted_all.root analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_weighted_all.root analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted_all.root')
    pass

print('DONE!')

#print('Calling ./scripts/renameNEW.sh')
#os.system('./scripts/renameNEW.sh')

