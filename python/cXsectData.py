#!/usr/bin/python

# Jiri Kvita, 14.4.2020
# file to store xsect help and stack properties classes and xsect data


#########################################

class cXsectHold:
    def __init__(self, name, xsect):
        self.name = name
        self.xsect = xsect

#########################################
        
# format: sample, xsect in pb
# obtained as:
#
# cd /home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/
# /home/qitek/install/mg5analysis/boosted_gacr/scripts/GetXsectMergedAver.sh
#
# from files like
# /home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1min60_ptj2min60/Events/run_22/tag_1_merged_xsecs.txt
# containing info like
# Merging scale       Cross-section [pb]  MC uncertainty [pb] 
# 45                  1.649391e+02        2.02e-01            
# 67.5                1.625231e+02        2.00e-01            
# 90                  1.534332e+02        1.94e-01 
# and taking the middle 67.5 GeV scale xsection;)

xsdata = [

    # DM model
    cXsectHold('pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad',  0.002440),
    cXsectHold('pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad', 0.002292),
    cXsectHold('pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad', 0.001239),

    # pp -> y0 -> ttbar
    # width 10 GeV # some problematic runs removed in Events/:
    cXsectHold('pp2y02tt_y0_1000GeV_NLO',               0.304690), 
    # larger widths:
    cXsectHold('pp2y02tt_y0_1000GeV_width100GeV_NLO',   0.031113),
    cXsectHold('pp2y02tt_y0_1000GeV_width300GeV_NLO',   0.010017),
    
    # ttbar

    # OLD main sample
    # cXsectHold('pp_2tj_allhad_NLO_ptj1min60_ptj2min60',  308.582143),
    # semiboosted ttbar sample, 21.6.2020
    # cXsectHold('pp_2tj_allhad_NLO_ptj1min100_ptj2min100',   92.096275),
    # boosted ttbar sample, 21.6.2020
    # cXsectHold('pp_2tj_allhad_NLO_ptj1min200_ptj2min200',   6.933571),

    
    # Main sample (a) X/2020
    cXsectHold('pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200', 137.756388),  # Fixed, unmerged: 271.090909
    # Main sample (b) X/2020
    cXsectHold('pp_2tj_allhad_NLO_ptj1j2min200', 4.606223),                # Fixed, unmerged:   7.281100
    # Main sample (c) X/2020
    cXsectHold('pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200', 19.484641), # Fixed, unmerged:  29.936500


    # alternative samples:
    # cXsectHold('pp_2t_allhad_14TeV_LO',          246.081818), # REMOVED ONE ODD RUN! after Delphes
    # cXsectHold('pp_2tj_allhad_NLO',              531.642857),
    # obsolete: cXsectHold('pp_2tj_allhad_NLO_ptheavy50GeV', 531.317073),

    # Zprime:
    cXsectHold('zp_ttbarj_allhad_1500GeV', 0.000037),  # Fixed; unmerged: 0.000066
    cXsectHold('zp_ttbarj_allhad_1250GeV', 0.000065),  # Fixed; unmerged: 0.000113
    cXsectHold('zp_ttbarj_allhad_1000GeV', 0.000121),  # Fixed; unmerged: 0.000206
    cXsectHold('zp_ttbarj_allhad_900GeV',  0.000158),  # Fixed; unmerged: 0.000267
    cXsectHold('zp_ttbarj_allhad_800GeV',  0.000210),  # Fixed; unmerged: 0.000350
    cXsectHold('zp_ttbarj_allhad_750GeV',  0.000244),  # Fixed; unmerged: 0.000403
    cXsectHold('zp_ttbarj_allhad_700GeV',  0.000284),  # Fixed; unmerged: 0.000466

    
    # Wjets
    #cXsectHold('pp2Wbb_pTj1min50GeV_allhad', 50.633810),
    #cXsectHold('pp2Wj_pTj1min50GeV_allhad', 12162.380952),
    cXsectHold('pp2wbbjjMatched_allhad_pTj1j2min_60GeV', 82.589070), # fixed; unmerged: 143.895000

    cXsectHold('pp2wwbb_allhad_pTj1j2min_60GeV', 126.100000),


    

    # QCD:
    # watch out the by-hand factor for jjj and jjjj samples!;)
    #cXsectHold('pp2jjjMatched_pTj1j2min_60GeV',    15854000., ),
    #cXsectHold('pp2jjjjMatched_pTj1j2min_60GeV', 17305375.00000),
    #cXsectHold('pp2jjjj_pTj1min400GeV_pTj2min60GeV',  20377.300000),
    # BAD: cXsectHold('pp2jjjj_mjjmin400GeV_14TeV', 193550000.000000),

    # NEW QCD 30.6.2021
    #pp_4j_ptj1j2min200
    #136040.010000 pb / 4 = 34010.002500 pb
    #----------------------------------------------
    #pp_4j_ptj1j2min60_ptj1j2max200
    #44062290.000000 pb / 5 = 8812458.000000 pb
    #----------------------------------------------
    #pp_4j_ptj1min200_ptj2min60max200
    #144231.370000 pb / 4 = 36057.842500 pb

    # 30.6.2021
    # QCD light flavoured samples
    cXsectHold('pp_4j_ptj1j2min200',  39212.325000 ),
    cXsectHold('pp_4j_ptj1j2min60_ptj1j2max200',  6311169.000000 ),
    cXsectHold('pp_4j_ptj1min200_ptj2min60max200',  26378.370000 ), 

    # 1.7.2021
    # matched bb + bbj + bbjj:
    # final numbers: 7.9.2021
    cXsectHold('pp_2b2j_LO_matched_ptj1j2min200', 604.473870),
    cXsectHold('pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200',  21773.217852),
    cXsectHold('pp_2b2j_LO_matched_ptj1min200_ptj2min60max200',  2659.769536),

] # end
