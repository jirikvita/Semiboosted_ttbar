#!/usr/bin/python
# Sat 29 May 10:16:43 CEST 2021

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from Tools import DivideByBinWidth

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

    cw = 800
    ch = 600

    ROOT.gStyle.SetOptTitle(0)
    
    scanname = 'shat'
    scan = ROOT.TCanvas(scanname, scanname, 0, 0, cw, ch)
    cans.append(scan)

    canname = 'flavour_fractions'
    fcan = ROOT.TCanvas(canname, canname, 150, 150, cw, ch)
    cans.append(fcan)

    filename = 'foo.root'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    rfile = ROOT.TFile(filename, 'read')

    dirname='pdfinfo/'
    hname = 'InitPartonsFracs'
    hInitPartonsFracs = rfile.Get(dirname + hname)
    hInitPartonsFracs.Scale(1./hInitPartonsFracs.Integral())

    canname = 'ffractions'
    gcan = ROOT.TCanvas(canname, canname, 250, 250, cw, ch)
    cans.append(gcan)
    gcan.cd()
    hInitPartonsFracs.SetMarkerSize(1.5)
    hInitPartonsFracs.SetLineColor(ROOT.kBlack)
    hInitPartonsFracs.SetMarkerStyle(20)
    hInitPartonsFracs.SetMinimum(0.)
    hInitPartonsFracs.SetMaximum(1.)
    hInitPartonsFracs.GetXaxis().SetLabelSize(0.08)
    hInitPartonsFracs.SetStats(0)
    hInitPartonsFracs.Draw('e1')
    gcan.Update()
 
    
    qgs = ['gg', 'qg', 'qq']
    dcols = {'gg' : ROOT.kPink, 'qq' : ROOT.kSpring, 'qg' : ROOT.kCyan+2} 
    hshat = {}
    ints = {}

    s1 = 0.
    s2 = 1.2 # TeV in range of shat
    
    htot = ROOT.TH1D('dum', 'dum', 10, s1, s2)
    for qg in qgs:
        basename = 'InitPartonsFracsShat_'
        h = rfile.Get(dirname + basename + qg)
        DivideByBinWidth(h)
        #h.Rebin(2)
        h.SetStats(0)
        h.GetXaxis().SetTitle('#hat{s} [TeV]')
        h.GetYaxis().SetTitle('normalized events')
        h.SetLineColor(dcols[qg])
        hshat[qg] = h
        ints[qg] = h.Integral(0, h.GetNbinsX(), 'width')
        if htot.GetEntries() <= 1e-4:
            htot = h.Clone('shat_tot')
        else:
            htot.Add(h)

    scan.cd()
    opt = ''
    ROOT.gPad.SetLogy(1)
    for qg in qgs:
        h = hshat[qg]
        h.SetLineWidth(2)
        h.GetXaxis().SetRangeUser(s1,s2)
        h.Draw('hist' + opt)
        opt = 'same'
    
            
    sname = 'shat_stack'
    hstack = ROOT.THStack(sname, sname)
    hfrac = {}
    delta = 0.07
    #fleg = ROOT.TLegend(0.80, 0.65, 0.88, 0.82)
    fleg = ROOT.TLegend(0.91, 0.50 - delta, 1., 0.50 + delta)
    fleg.SetBorderSize(0)
    lqg = []
    for qg in qgs:
        lqg.append(qg)
        ratio = hshat[qg].Clone(hshat[qg].GetName() + '_ratio')
        ratio.Divide(htot)
        ratio.SetFillColor(dcols[qg])
        ratio.SetFillStyle(1001)
        ratio.SetLineWidth(1)
        hfrac[qg] = ratio
        hstack.Add(hfrac[qg], 'hist')
    lqg.reverse()
    for qg in lqg:
        fleg.AddEntry(hfrac[qg], qg, 'F')

    fcan.cd()
    s2 = 1. # TeV!
    eps = 0.004
    h2 = ROOT.TH2D('dum2', 'dum2;#hat{s} [TeV];fraction', 100, s1, s2, 100, 0. - eps, 1. + eps)
    h2.SetStats(0)
    h2.Draw()
    hstack.Draw('same')
    #hstack.GetXaxis().SetRangeUser(s1,s2)
    #hstack.GetYaxis().SetRangeUser(0.,1.)
    fleg.Draw()
    ROOT.gPad.Update()

    scan.cd()
    fleg.Draw()
    
    
    pdgid = [21, -4, -3, -2, -1, 1, 2, 3, 4]
    hx = {}
    for pdg in pdgid:
        hx[pdg] = rfile.Get(dirname + 'h_xPdfHisto_' + str(pdg))
        DivideByBinWidth(hx[pdg])

    canname = 'pdfx'
    canx = ROOT.TCanvas(canname, canname, 300, 300, cw, ch)
    cans.append(canx)
    ROOT.gPad.SetLogy(1)
    cols = [ROOT.kRed, ROOT.kBlack, ROOT.kGreen + 2, ROOT.kBlue, ROOT.kViolet]
    opt = ''
    ymax = -1
    for pdg in hx:
        hmax = hx[pdg].GetMaximum()
        #hx[pdg].Rebin(10)
        hx[pdg].SetStats(0)
        hx[pdg].GetXaxis().SetTitle('x')
        hx[pdg].Scale(1./hx[pdg].Integral(0, hx[pdg].GetNbinsX(), 'width'))
        if ymax < hmax:
            ymax = 1.1*hmax
    xleg = ROOT.TLegend(0.75, 0.55, 0.88, 0.88)
    for pdg in hx:
        i = abs(pdg)
        if i > 20:
            i = 0
        h = hx[pdg]
        h.SetMaximum(ymax)
        h.SetLineColor(cols[i])
        h.SetLineWidth(2)
        if pdg < 0:
            h.SetLineStyle(2)
        h.Draw('hist' + opt)
        xleg.AddEntry(h, str(pdg), 'L')
        opt = 'same'
    xleg.Draw()


    for can in cans:
        can.Print('pdf/' + can.GetName() + '.pdf')
        can.Print('png/' + can.GetName() + '.png')
        
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

