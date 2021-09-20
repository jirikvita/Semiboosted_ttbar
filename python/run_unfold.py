#!/usr/bin/python3

# jk 25.2.2018, 7.3.2018, 3.4.2018, 11.6.2019, 4.6.2020
# 2020 TODO: add bg-subtracted spectra, acc and eff ingredients!


from Unfold import *

fnames = [

    ['analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root',
     'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half1.root',
     #[],
     3., 3., -1, False],
    
]

dirnames = ['2B0S/',
            '1B1S/',
            '0B2S/'
]
migranames = [ 
    'DiTopMass_Particle_Detector',
    'TopPt_Particle_Detector',
]


for fname in fnames:
    fname_data = fname[0]
    fname_mc = fname[1]
    #SF = fname[2]
    ValSFdown = fname[2]
    ValSFup = fname[3]
    ValAddSFLastFirst = fname[4]
    smear = fname[5]
    woption = 'recreate'
    rebin = -1
    #ValSFdown = SF # sf to restrict the range of truth to be sampled around reco
    #ValSFup = SF # sf to restrict the range of truth to be sampled around reco
    corrs = {}
    for dirname in dirnames:
        for migraname in migranames:
            print('Unfolding projections of %s from file %s' % (fname_data, migraname,))
            outtag = dirname + migraname
            outtag = outtag.replace('/','_')
            outfilename = fname_data.replace('analyzed_histos', 'unfolded_' + outtag)
            ftag = ''
            #corrs['acc'] =
            #corrs['eff'] = 
            Unfold(fname_data, fname_mc, dirname, migraname, outfilename, woption, smear, rebin, ValSFdown, ValSFup, ftag, ValAddSFLastFirst)
            #woption = 'update'
