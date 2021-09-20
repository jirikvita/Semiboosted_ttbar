#!/usr/bin/python
# Tue 20 Jul 08:28:08 CEST 2021

# references:
# https://root.cern/doc/master/rf707__kernelestimation_8py.html
# https://root.cern/doc/v610/rf402__datahandling_8C_source.html

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

##########################################
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

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    canname = 'can'
    can = ROOT.TCanvas(canname, canname)
    cans.append(can)

    # later open half0 and half1 files:
    filename = 'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root'
    bgfile = ROOT.TFile(filename, 'read')
    filename = 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_all.root'
    sigfile = ROOT.TFile(filename, 'read')

    # 1D play"

    hname = '1B1S/DetectorDiTopMass'
    hbg = bgfile.Get(hname)
    hsig = sigfile.Get(hname)

    # make data by a simple addition:
    hdata = hbg + hsig 
    
    stuff.append([hbg, hsig])
    xmin = hsig.GetXaxis().GetXmin()
    xmax = hsig.GetXaxis().GetXmax()
    x = ROOT.RooRealVar('x', 'x', xmin, xmax)
    data = ROOT.RooDataHist('data', 'data', ROOT.RooArgList(x), ROOT.RooFit.Import(hdata))
    bg = ROOT.RooDataHist('bg', 'bg', ROOT.RooArgList(x), ROOT.RooFit.Import(hbg))
    sig = ROOT.RooDataHist('sig', 'sig', ROOT.RooArgList(x), ROOT.RooFit.Import(hsig))
    
    # make kernel estimate smoothed pdf
    # An adaptive kernel estimation pdf on the same data without mirroring option
    # for comparison
    # bg:

    # JK 20.7.2021
    # seems this ain't work, RooFit supports only the unbinned RooDataSet as input to get the kernel smoothed PDF...
    kest_bg = ROOT.RooKeysPdf('kestbg1', 'kestbg1', x, bg, ROOT.RooKeysPdf.NoMirror)
    # sig:
    kest_sig = ROOT.RooKeysPdf('kestsig1', 'kestsig1', x, sig, ROOT.RooKeysPdf.NoMirror)
    
    frac1 = ROOT.RooRealVar('frac1', 'frac1', 0.2, 0., 1.);
    nEvt = hdata.GetEntries()

    n1 = ROOT.RooRealVar('n1', 'n1', nEvt*frac1.getVal(), 0., 1000000);
    n2 = ROOT.RooRealVar('n2', 'n2', nEvt*(1.-frac1.getVal()), 0., 1000000);
    
    # fit data to model:
    model = ROOT.RooAddPdf('model', 'kestbg1 + kestsig1', RooArgSet(kestbg1,kestsig1), RooArgSet(frac1))
    #model = ROOT.RooAddPdf model('model', 'gsigPdf1 + gsigPdf2', RooArgSet(sigPdf1, sigPdf2), RooArgSet(n1,n2))
    model.fitTo(data)
    model.plotOn(frame)
    
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

