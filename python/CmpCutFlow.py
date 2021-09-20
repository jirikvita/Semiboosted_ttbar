#!/usr/bin/python

from __future__ import print_function

# jk 16.9.2018

from ROOT import *
from Tools import MakeNiceLegendEntry

superdir = './'
fnames = [ superdir + 'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root',

    superdir + 'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_all.root',

           #superdir + 'analyzed_histos_zp_ttbarj_allhad_700GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_750GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_800GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_900GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_1250GeV_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_zp_ttbarj_allhad_1500GeV_14TeV_ATLAS_all.root',

    superdir + 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp2y02tt_y0_1000GeV_width100GeV_NLO_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp2y02tt_y0_1000GeV_width300GeV_NLO_14TeV_ATLAS_all.root',
    
    superdir + 'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV_ATLAS_all.root',
    superdir + 'analyzed_histos_pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad_14TeV_ATLAS_all.root',
    
    superdir + 'analyzed_histos_pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV_ATLAS_all.root',
    ]


# HACK!
#fnames  = [superdir + 'analyzed_histos_run_190_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_origFwk.root',
#           superdir + 'analyzed_histos_run_190_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS.root',
#]

lsts = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,2,1,2,1,2,1,2,1,2,1,2]
lcols = [kBlue,
         kBlue+1, kBlue+2, kBlue-2,
         kRed-1, kRed-2, kRed+1, kRed+2, kRed+3, kRed,
         kViolet, kViolet + 2, kViolet -2,
         kBlack, kGray, kGray+2,
         kGreen+2]

levels = [   
    ###!!!'Ptcl',
     'Det',
]

hname = 'CutFlow'

sameopt = ''
opt = 'hist'

inormbin = 2
iInitial = 2

canname = 'CutflowCmp'
can = TCanvas(canname, canname, 0, 0, 1600, 800)
#nx = len(levels)
#ny = len (fnames)
nx = 2 # 3
ny = 1 # 2

stuff = []

can.Divide(nx,ny)
for i in range(0, nx*ny):
    can.cd(i+1)
    #gPad.SetLogy()

rfiles = []
for fname in fnames:
    rfiles.append(TFile(fname, 'read'))

#idir = 1
legs = {}
for rfile,lst,lcol in zip(rfiles,lsts,lcols):


    #print('Working on file {}'.format(rfile.GetName()))
    for level in levels:

        if rfile == rfiles[0]:
            leg = TLegend(0.35, 0.6, 0.88, 0.88)
            leg.SetBorderSize(0)
            legs[level] = leg
        idir = levels.index(level)
        can.cd(idir+1)
        #idir = idir+1
        cname = '{}{}'.format(hname,level)
        #print('...working on {}'.format(cname))
        h = rfile.Get(cname)
        if 'Initial' in h.GetXaxis().GetBinLabel(iInitial):
            print('{} : {} initial events'.format(rfile.GetName(), h.GetBinContent(iInitial)))

        if inormbin > 0:
            val = h.GetBinContent(inormbin)
            if val > 0.:
                h.Scale(1./val)
        #print('Filename: {}'.format(rfile.GetName()))
        ftag = rfile.GetName().split('/')[-1].replace('histos','').replace('.root','')
        legs[level].AddEntry(h, MakeNiceLegendEntry(ftag), 'PL')
        if inormbin > 0:
            h.SetMinimum(0.)
            h.SetMaximum(2);
        else:
            h.SetMinimum(0.)
            h.SetMaximum(1.e7);

        h.SetStats(0)
        h.SetLineColor(lcol)
        h.SetLineStyle(lst)
        h.SetLineWidth(2)
        h.Draw(opt + sameopt)
        h.GetYaxis().SetMoreLogLabels()
    sameopt = 'same'

for rfile,lst,lcol in zip(rfiles,lsts,lcols):
    for level in levels:
        idir = levels.index(level)
        can.cd(idir+1)
        legs[level].Draw()
        gPad.Update()

can.Print('png/' + can.GetName() + '.png')
can.Print('pdf/' + can.GetName() + '.pdf')
gApplication.Run()
