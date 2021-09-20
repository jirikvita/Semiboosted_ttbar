#!/usr/bin/python
# Sun 26 Jul 16:07:35 CEST 2020

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

##########################################
def MakeLine(x1, y1 ,x2, y2, col, lst = 2, lw = 3):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(col)
    line.SetLineWidth(lw)
    line.SetLineStyle(lst)
    line.Draw()
    stuff.append(line)
    return line


##########################################
def PlotH2AndProfile(h, yBandWidth = 0.15, center = 1., delta = 0.05, lc = ROOT.kRed, mc = ROOT.kRed, ms = 20, lcol = ROOT.kBlack):
    hname = h.GetName()
    h.SetStats(0)
    meany = h.GetMean(2)

    if 'JetsPt' in hname:
        h.GetXaxis().SetTitle('p_{T} [GeV]')
        h.GetXaxis().SetMoreLogLabels()
    
    h.GetXaxis().SetTitleOffset(1.5)
    p = h.ProfileX(hname + '_profX')
    if 'R_E' in h.GetYaxis().GetTitle():
        h.GetYaxis().SetTitle(h.GetYaxis().GetTitle().replace('R_E','R_{E}'))
    h.Draw('colz')
    print('Setting y axis rantge to {},{}'.format(center - yBandWidth,  center + yBandWidth))
    h.GetYaxis().SetRangeUser(center - yBandWidth,  center + yBandWidth)

    
    #p.SetName(hname + '_profX')

    if ('JetsPt' in hname or 'JetsE' in hname ) and not 'JetsEta' in hname:
        # to be removed!
        # h.GetXaxis().SetRangeUser(25., h.GetXaxis().GetXmax())
        ROOT.gPad.SetLogx(1)
        #ROOT.gPad.SetLogz(1)
        h.GetXaxis().SetMoreLogLabels()
    p.SetMarkerColor(mc)
    p.SetLineColor(lc)
    p.SetMarkerStyle(ms)
    p.Draw('e1 hist same')
    #p.SetLineColor(mc)
    #p.DrawCopy('hist same')
    tex = ROOT.TLatex(0.35, 0.93, '<{}> = {:1.4f}'.format(h.GetYaxis().GetTitle().split()[0],meany))
    tex.SetNDC()
    #tex.SetTextColor(lc)
    tex.Draw()
    stuff.append(tex)

    line1 = MakeLine(h.GetXaxis().GetXmin(), center, h.GetXaxis().GetXmax(), center, lcol)
    line1.SetLineStyle(1)
    line1.Draw()
    line2 = MakeLine(h.GetXaxis().GetXmin(), center+delta, h.GetXaxis().GetXmax(), center+delta, lcol)  
    line3 = MakeLine(h.GetXaxis().GetXmin(), center-delta, h.GetXaxis().GetXmax(), center-delta, lcol)  

    res = [p, line1, line2, line3]
    stuff.append(res)
    
    return

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

    # https://root.cern.ch/doc/master/classTColor.html
    #ROOT.gStyle.SetPalette(ROOT.kDeepSea)
    #ROOT.gStyle.SetPalette(ROOT.kCool)
    #ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    #ROOT.gStyle.SetPalette(ROOT.kCopper)
    #ROOT.gStyle.SetPalette(ROOT.kBird) # std, nice
    #ROOT.gStyle.SetPalette(ROOT.kAvocado)
    #ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
    #ROOT.gStyle.SetPalette(ROOT.kLake)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)
    
    ROOT.gStyle.SetOptTitle(0)
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    if len(argv) < 2:
        print('Usage: {} filename.root'.format(argv[0]))
        #print('E.g.:')
        #print('{} data/analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_all_forJES.root'.format(argv[0]))
        exit(1)
    
    filename = argv[1]

    jets = ['Jets', 'LJets']

    quantity = 'Closure'
    yBandWidth = 0.25 # 0.15
    if not 'Closure' in filename:
        quantity = ''
        #yBandWidth = 0.30
    
    for jet in jets:
    
        rfile = ROOT.TFile(filename, 'read')
        stuff.append(rfile)
        nvars = 2
        hnames = [###'JESRPt' + quantity + jet + 'Eta',
                  ###'JESRPt' + quantity + jet + 'Pt',
                  ###'JESRPt' + quantity + jet + 'E',
                  'JESRE' + quantity + jet +  'Eta',
                  'JESRE' + quantity + jet +  'Pt',
                  ###'JESRE' + quantity + jet +  'E',
                  
                  # eta and phi diff closure:
                  #'JESDEta' + quantity + jet + 'Eta',
                  #'JESDEta' + quantity + jet + 'Pt',
                  #'JESDEta' + quantity + jet +  'E',
                  #'JESDPhi' + quantity + jet + 'Eta',
                  #'JESDPhi' + quantity + jet + 'Pt',
                  #'JESDPhi' + quantity + jet +  'E',
        ]
        ncorrs = len(hnames) / nvars
        canname = 'JES' + quantity + jet
        #can = ROOT.TCanvas(canname, canname, 0, 0, 800, 1200)
        #can.Divide(2,3)
        can = ROOT.TCanvas(canname, canname, 0, 0, 400*ncorrs, 400*nvars)
        can.Divide(ncorrs,nvars)
        cans.append(can)
        

        dirname = 'JES/'
        histos = []
        for hname in hnames:
            print('Getting {}'.format(dirname + hname))
            h = rfile.Get(dirname + hname)
            histos.append(h)
            print(h)
        stuff.append(histos)

        for h in histos:
            ii = histos.index(h)
            ican = ii+1
            # compact:
            if len(hnames) <= 6:
                ican = ncorrs*ii+1
                if ii > ncorrs:
                    ican = ncorrs*(ii-nvars)+ncorrs
            can.cd(ican)
            yband = 1.*yBandWidth
            center = 1.
            delta = 0.05
            if 'DEta' in h.GetName() or 'DPhi' in h.GetName():
                print('OK, zooming the y axis!')
                yband = 0.05
                center = 0.
                delta = 0.01
            res = PlotH2AndProfile(h, yband, center, delta)
            
        can.Update()
        can.Print('png/' + canname + '.png')
        can.Print('pdf/' + canname + '.pdf')
        
    ROOT.gApplication.Run()
    return cans

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

