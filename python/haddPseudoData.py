#!/usr/bin/python3

# jk III/2020, coronavirus times
# obsolete, adding now done based on xsect weights in XsectStack.py

import os, sys

toMerge = { #'analyzed_histos_pseudo_data_tt_ttj_y0_y0_1000GeV_half0.root' : ['analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
            #                                                                 'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
            #                                                                 'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
            #                                                                 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root'],
            # 
            # 'analyzed_histos_pseudo_data_ttj_y0_y0_1000GeV_half0.root' : ['analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
            #                                                               'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
            #                                                               'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root'],
            # 
            'analyzed_histos_pseudo_data_tt_ttj_half0.root' : ['analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                               'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                               'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                               'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],
              
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_y0_y0_1000GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                    'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                    'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                    'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root'],
            
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_y0_y0_1000GeV_width100GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                                'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2y02tt_y0_1000GeV_width100GeV_NLO_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_y0_y0_1000GeV_width300GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                                'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2y02tt_y0_1000GeV_width300GeV_NLO_14TeV_ATLAS_half0.root'],
            
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_xdxd_y0_1000GeV_xd_10GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                               'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                               'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                               'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                               'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                               'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                               'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV_ATLAS_half0.root'],

            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_xdxd_y0_1000GeV_xd_100GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                                'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV_ATLAS_half0.root'],
            
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_xdxd_y0_1000GeV_xd_300GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                                'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                                'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad_14TeV_ATLAS_half0.root'],

            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_1500GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_1500GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_1250GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_1250GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_1000GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_900GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_900GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_800GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_800GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_750GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_750GeV_14TeV_ATLAS_half0.root'],
            'analyzed_histos_pseudo_data_Wj_Wbb_tt_ttj_Zprime_700GeV_half0.root' : ['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root',
                                                                                     'analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root',
                                                                                    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
                                                                                     'analyzed_histos_zp_ttbarj_allhad_700GeV_14TeV_ATLAS_half0.root'],
                                                                                     
            
            
            
}


for ihalf in range(0,2):
    print('*** Processing half {}'.format(ihalf))
    for starget in toMerge:

        # HACK!!!
        #if not 'width300' in starget:
        #    #if not 'Zp' in starget:
        #    continue
        
        target = starget + ''
        target = starget.replace('half0', 'half{}'.format(ihalf))
        files = ''
        for fname in toMerge[starget]:
            files += '{} '.format(fname)
        files = files.replace('half0', 'half{}'.format(ihalf))
        cmd = '  hadd -f {} {}'.format(target, files)
        print(cmd)
        os.system(cmd)
