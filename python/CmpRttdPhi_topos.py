#!/usr/bin/python
# Sun 23 May 13:03:25 CEST 2021

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

    filename = ''
    topos = ['2B0S/', '1B1S/', '0B2S/']
    basehnames = [ 'DetectorRttbarVsDetectorDiTopDeltaPhi' ]
    filename = 'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root'
    rfile = ROOT.TFile(filename, 'read')
    can2s = []
    cols = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2]
    for hname in  basehnames:
        canname = 'can_{}'.format(hname)
        can2 = ROOT.TCanvas(canname, canname, 0, 0, 1200, 400)
        can2.Divide(3,1)
        can2s.append(can2)
        opt = 'col' # col
        itopo = -1
        for col,topo in zip(cols,topos):
            itopo = itopo + 1
            can2.cd(itopo+1)
            fullhname = topo + hname
            print('Trying to get {}'.format(fullhname))
            h = rfile.Get(fullhname)
            h.SetMarkerColor(col)
            h.SetMarkerSize(0.2)
            #h.SetLineColor(col)
            #h.SetStats(0)
            h.Scale(1./h.GetEntries())
            #h.SetMaximum(0.01)
            #h.SetMinimum(1.e-5)
            h.Draw(opt)
            ROOT.gPad.Update()
            st = h.FindObject('stats')
            try:
                wx = 0.19
                wy = 0.13
                x0 = 0.12
                y0 = 0.12
                st.SetX1NDC(x0)
                st.SetY1NDC(y0)
                st.SetX2NDC(x0 + wx)
                st.SetY2NDC(y0 + wy)
            except:
                print('Failed to get the histo stat object!')
                
            ROOT.gPad.SetLogz(1)
            #opt = 'hist same'
        can2.Update()
    
    ROOT.gApplication.Run()
    return

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

