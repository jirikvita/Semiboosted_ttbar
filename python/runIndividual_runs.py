#!/usr/bin/python3

import os, sys

# DEFAULT:
#skipExistingFiles=False
# if you really know what you're doing!
# this does not re-run, only run additinal meanwhile-produced samples!!!
skipExistingFiles=False

datadir='/home/public/data/ttbar/data/qitek/DelphesOut/'

tags = [ #'2tj_allhad_NLO_ATLAS', # no cuts, equiv to next sample:)
         #'2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS', # not really pTheavy...
	 #'2t_allhad_14TeV_LO_ATLAS',
          #'pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS',
         #'pp_2tj_allhad_NLO_ptj1min100_ptj2min100_14TeV_ATLAS',
         #'pp_2tj_allhad_NLO_ptj1min200_ptj2min200_14TeV_ATLAS',
	 'pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS',
         'pp2y02tt_y0_1000GeV_width100GeV_NLO_14TeV_ATLAS',
         'pp2y02tt_y0_1000GeV_width300GeV_NLO_14TeV_ATLAS',
	 #'pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS',
	 #'pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS',
         'pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV_ATLAS',
	 'pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV_ATLAS',
         'pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad_14TeV_ATLAS',
    
    'zp_ttbarj_allhad_1500GeV_14TeV_ATLAS',
         #'zp_ttbarj_allhad_1250GeV_14TeV_ATLAS',
    #'zp_ttbarj_allhad_1000GeV_14TeV_ATLAS',
    'zp_ttbarj_allhad_1000GeV_NEW_14TeV_ATLAS',
    'zp_ttbarj_allhad_1250GeV_NEW_14TeV_ATLAS',
         'zp_ttbarj_allhad_800GeV_14TeV_ATLAS',
         'zp_ttbarj_allhad_900GeV_14TeV_ATLAS',
         'zp_ttbarj_allhad_750GeV_14TeV_ATLAS',
    #'zp_ttbarj_allhad_700GeV_14TeV_ATLAS',
         #'pp2jjjMatched_pTj1j2min_60GeV_14TeV_ATLAS',
         #'pp2jjjjMatched_pTj1j2min_60GeV_14TeV_ATLAS',
         #'pp2jjjj_pTj1min400GeV_pTj2min60GeV_14TeV_ATLAS',
         'pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV_ATLAS',
    
    'pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS',
    'pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS',
    'pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS',

    'pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted',
    'pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_weighted',
    'pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted',

    'pp_zp_ttbar_allhad_1000GeV_14TeV_ATLAS',
    'pp2wwbb_allhad_pTj1j2min_60GeV_14TeV_ATLAS',

    # 30.6.2021:
    'pp_4j_ptj1j2min60_ptj1j2max200_14TeV_ATLAS',
    'pp_4j_ptj1j2min200_14TeV_ATLAS',
    'pp_4j_ptj1min200_ptj2min60max200_14TeV_ATLAS',

    # 1.7.2021:
    'pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV_ATLAS',
    'pp_2b2j_LO_matched_ptj1j2min200_14TeV_ATLAS',
    'pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV_ATLAS',

    # 13.7.2021
    'pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted',
    'pp_2b2j_LO_matched_ptj1j2min200_14TeV_ATLAS_weighted',
    'pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted',

    
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
        doRun = True
        outfname = 'analyzed_histos_{}.root'.format(sample)
        outputExist = os.path.isfile('./{}'.format(outfname)) or os.path.isfile('./runs/{}'.format(outfname)) 
        if skipExistingFiles and outputExist:
            doRun = False
        if doRun:
            scriptfile.write('nice -n 10 ./bin/AnalyzeBoosted_x {} {} Trimmed {} >& log_{}{}.txt \n'.format(datadir, sample, weight, run, wtag))
    scriptfile.close()
    #print('done!')
    #print('Now you can run:')
    print('./scripts/jobStarter.sh --jobs 8 {}'.format(scriptname))

#print('DONE!')
