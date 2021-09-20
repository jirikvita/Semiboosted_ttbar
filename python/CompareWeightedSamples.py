#!/usr/bin/python

# jk 11.1.2019, 13.7.2021

from __future__ import print_function

from collections import OrderedDict

import os, sys

import ROOT

from Tools import MakeNiceLegendEntry

divideByBinWidth = True

from Tools import DivideByBinWidth

stuff = []
hists = []

sels = ['Detector', 'Particle']

ROOT.gStyle.SetOptTitle(False)

dirname='NoCuts/'

nSamples = 3
iWeighted = -1
iUnweighted = 1

fnames_2tj = [
    'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_unweighted_all.root',
    'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root',
]

# 13.7.2021
fnames_2b2j = [
    'analyzed_histos_pp_2b2j_LO_matched_ALL_14TeV_ATLAS_unweighted_all.root',
    'analyzed_histos_pp_2b2j_LO_matched_ptj1j2min200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted_all.root',
    'analyzed_histos_pp_2b2j_LO_matched_ALL_14TeV_ATLAS_weighted_all.root',
    ]


### STEERING CHOICE!!!
fnames = fnames_2tj
#fnames = fnames_2b2j

lsts = [1, 1, 1, 1, 1]
lws =  [1, 2, 2, 2, 3]
cols = [ ROOT.kBlack, ROOT.kRed+2, ROOT.kRed+1, ROOT.kOrange, ROOT.kRed]
#dfst = 3001
fsts = [-1, 3325, 3352, 3350, -1]

pngdir='png/'
pdfdir='pdf/' 
os.system('mkdir -p {:}'.format(pngdir,))
os.system('mkdir -p {:}'.format(pdfdir,))

Cans = {}
Legs = {}
FileTagDict = OrderedDict()

hnames = [#'JetPt',
          #'Jet2Pt',
    #'Jet1Rapidity',
    'Jet1Pt',
    #'LJet1Rapidity',
    'LJet1Pt',
]

ymin = 1e-3
ysf = 20.

Cans = {}
Hists = {}
for sel in sels:
    Cans[sel] = []
    Hists[sel] = []
    Legs = {}

    sameopt = ''
    for fname,lst,lw in zip(fnames,lsts,lws):
        rfile = ROOT.TFile(fname, 'read')
        ftag=fname
        ftag=ftag.replace('analyzed_histos_', '').replace('.root', '')
        ltag = MakeNiceLegendEntry(ftag)
        FileTagDict[ltag] = rfile


    for hname in hnames:


        leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
        leg.SetBorderSize(0)
        Legs[sel] = leg
        
        canname = 'can_WeightCmp_{}_{}'.format(sel,hname)
        can = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
        ROOT.gPad.SetLogy(1)
        Cans[sel].append(can)

        sameopt = ''
        opt = 'hist'
        
        for ftag,col,lw,lst,fst in zip(FileTagDict,cols,lws,lsts,fsts):
            rfile = FileTagDict[ftag]
            fullhname = '{:}{:}{:}'.format(dirname,sel,hname)
            print('Trying to get {}'.format(fullhname))
            histo = rfile.Get(fullhname)
            Hists[sel].append(histo)
            histo.SetStats(0)
            histo.SetLineColor(col)        
            histo.SetLineWidth(lw)
            histo.SetLineStyle(lst)
            if fst >= 0:
                histo.SetFillStyle(fst)
                histo.SetFillColor(col)
            if divideByBinWidth:
                DivideByBinWidth(histo)

            histo.SetMinimum(ymin)
            histo.SetMaximum(ysf*histo.GetMaximum())
            histo.Draw(opt + sameopt)
            sameopt = 'same'
            leg.AddEntry(histo, ftag, 'LF')

        leg.Draw()
        ROOT.gPad.Update()
        tag = ''
        if '2b2j' in fnames[0]:
            tag = '_2b2j'
        ROOT.gPad.Print(pngdir + can.GetName() + tag + '.png')
        ROOT.gPad.Print(pdfdir + can.GetName() + tag + '.pdf')

        

ROOT.gApplication.Run()
