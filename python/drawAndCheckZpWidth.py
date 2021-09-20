#!/usr/bin/python
# Fri 23.3.2021

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []
rfiles = []

def GetHisto(fname, path, mass, eps = 0.025):
    rfile = ROOT.TFile(path + fname, 'read')
    rfiles.append(rfile)
    Delphes = rfile.Get("Delphes")
    hname = 'zp_mass_{}'.format(mass)
    histo = ROOT.TH1D(hname, hname, 100, mass - eps*mass/1250., mass + eps*mass/1250.)
    Delphes.Draw("GenZPrime.Mass >> {}".format(hname))
    return histo


    

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
    gr = ROOT.TGraphErrors()
    gr.SetName('zpw')
    gr.SetTitle(';m_{Z\'} [GeV];#Gamma_{Z\'} [GeV]')
    wparam = 0.0008252

    path = '/home/public/data/ttbar/data/qitek/DelphesOut/'
    fnames = [    'out_boosted_AtKt4_and_10_run_11_zp_ttbarj_allhad_750GeV_14TeV_ATLAS.root', 
                  'out_boosted_AtKt4_and_10_run_11_zp_ttbarj_allhad_800GeV_14TeV_ATLAS.root', 
                  'out_boosted_AtKt4_and_10_run_11_zp_ttbarj_allhad_900GeV_14TeV_ATLAS.root', 
                  'out_boosted_AtKt4_and_10_run_13_zp_ttbarj_allhad_1000GeV_NEW_14TeV_ATLAS.root', 
                  'out_boosted_AtKt4_and_10_run_13_zp_ttbarj_allhad_1250GeV_NEW_14TeV_ATLAS.root', 
                  'out_boosted_AtKt4_and_10_run_11_zp_ttbarj_allhad_1500GeV_14TeV_ATLAS.root'
                  ]


    massPoints = [750, 800, 900, 1000, 1250, 1500]

    Data = {}
    ip = 0
    histos = {}
    canname = 'zpMasses'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    can.Divide(3, 2)
    stuff.append(can)
    ican = 1
    for mass in massPoints:
        for fname in fnames:
            if '{}'.format(mass) in fname:
                can.cd(ican)
                ROOT.gPad.SetLogy()
                ican = ican + 1

                histo = GetHisto(fname, path, mass)
                histo.SetMarkerColor(1)
                histo.SetMarkerSize(0.75)
                histo.SetMarkerStyle(24)
                fitname = 'fit_{}'.format(mass)
                fun = ROOT.TF1(fitname, '[0]*[2] / ( (x^2 - [1]^2 )^2 + [1]^2*[2]^2 / 4)', histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax())
                fun.SetParameters(0.1*histo.GetEntries() * histo.GetRMS() , histo.GetMean(), 200*histo.GetRMS() / mass)
                fun.SetNpx(1000)
                fun.SetLineWidth(1)
                histo.Fit(fitname)
                #histo.Fit(fitname)
                #histo.Fit(fitname)
                ROOT.gStyle.SetOptFit(11111)
                histo.Draw('e1')
                fun.Draw('same')
                histo.Draw('e1 same')
                histos[mass] = histo
                # wobs = histo.GetMean()
                wobs = fun.GetParameter(2)
                werr = fun.GetParError(2)
                Data[mass] = wobs
                wcomp = wparam*mass*mass
                print('mass: {}, wobs: {}, wcopm={} ratio={}'.format(mass, wobs, wcomp,wobs/wcomp))
                gr.SetPoint(ip, mass, wobs)
                gr.SetPointError(ip, 0., werr)
                ip = ip + 1

    stuff.append(histos)


    canname = 'zpMasses_gr'
    gcan = ROOT.TCanvas(canname, canname, 100, 100, 800, 800)
    stuff.append(gcan)
    
    gr.SetMarkerSize(1)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(ROOT.kBlack)
    gr.Draw('AP')
    gr.Fit('pol2')
    stuff.append(gr)
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

