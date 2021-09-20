#!/usr/bin/python
# Tue 23 Jul 09:00:00 CEST 2021

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp, pi
import os, sys, getopt

from FitTools import *


###############################################
###############################################
###############################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    gTag=''
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
        
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    bgfilename = 'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_half0.root'
    bgfile = ROOT.TFile(bgfilename, 'read')
    sigfilename = 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root'
    sigfile = ROOT.TFile(sigfilename, 'read')

    bgfilenameAlt = 'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_half1.root'
    bgfileAlt = ROOT.TFile(bgfilenameAlt, 'read')
    sigfilenameAlt = 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half1.root'
    sigfileAlt = ROOT.TFile(sigfilenameAlt, 'read')
    
    draw = True

    Results = []

    hname ='scoreh0'
    scoreh0 = ROOT.TH1D(hname, hname + ';score0;replicas', 50, -10., 10)
    hname ='scoreh1'
    scoreh1 = ROOT.TH1D(hname, hname + ';score1;replicas', 50, -10., 10)
    
    #for irep in range(1,2):
    for irep in range(1,99):
        print('=== Processinfg relica {}'.format(irep))
        
        # KEY STEETING -- HISTOGRAMME CHOICE!
        
        # 1D play"
        sigSF = 0.155
        hname = '1B1S/replicas/DetectorDiTopMass_rep{}'.format(irep)
        
        # 2D play:
        #sigSF = 0.015
        #hname = '1B1S/replicas/DetectorHTjVsDetectorDiTopMass_rep{}'.format(irep)

        fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1, stuff = Analyze(hname, sigSF, bgfile, sigfile, bgfileAlt, sigfileAlt, draw)
        Results.append([fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1])
        draw = False
        scoreh0.Fill(logt0)
        scoreh1.Fill(logt1)

    canname = 'betterScores'
    if 'Vs' in hname:
        canname = canname + '2d'
    else:
        canname = canname + '1d'
    scan = ROOT.TCanvas(canname, canname, 600, 300, 1200, 600)
    scan.Divide(2,1)
    cans.append(scan)


    hists = [scoreh0, scoreh1]
    cols = [ROOT.kSpring, ROOT.kCyan, ROOT.kOrange, ROOT.kGray]
    for i in range(0, len(hists)):
        hist = hists[i]
        hist.SetFillColor(cols[i])
        scan.cd(i+1)
        hist.Draw('hist')
    
    scan.Print(scan.GetName() + '.pdf')
    scan.Print(scan.GetName() + '.png')


    ROOT.gApplication.Run()
        
    return

###################################
###################################
###################################

if __name__ == '__main__':
    # execute only if run as a script'
    main(sys.argv)
    
###################################
###################################
###################################

