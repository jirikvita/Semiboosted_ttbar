#!/usr/bin/python
# Tue 20 Jul 14:15:00 CEST 2021

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp, pi
import os, sys, getopt

cans = []
stuff = []

##########################################

class fitres:
    def __init__(self, result,error,fitstatus): #,fitchi2,fitprob):
        self.result = result
        self.error = error
        self.fitstatus = fitstatus
        #self.fitchi2 = fitchi2
        #self.fitprob = fitprob
##########################################



##########################################
def GetFractionFitterResult(hdata, samples):
    xmin = hdata.GetXaxis().GetXmin()
    xmax = hdata.GetXaxis().GetXmax()

    nsamples = len(samples)
    mc_array = ROOT.TObjArray(nsamples)
    for sample in samples:
        mc_array.Add(sample)

    
    fitter = ROOT.TFractionFitter(hdata, mc_array)
    # constrain fractions to be between 0 and something
    for i in range(0, nsamples):
        fitter.Constrain(i, 0.0, 1.1)

    # Fit! ;-)
    fitstatus = fitter.Fit()

    result = ROOT.Double()
    error = ROOT.Double()
    Result = []
    for i in range(0, nsamples):
      fitter.GetResult(i, result, error)
      Result.append([1.*result, 1.*error])
      print('sample {} fraction: {} error: {} signif.: {}'.format(i, result, error, result/error))
    
    #fitchi2 = fitter.GetChisquare()
    #fitprob = fitter.GetProb()

    return fitres(Result[0][0],Result[0][1],fitstatus) #,fitchi2,fitprob)

###############################################

def Analyze(hname, sigSF, bgfile, sigfile, bgfileAlt, sigfileAlt, draw = True):

    is2d = 'Vs' in hname
    
    hbg = bgfile.Get(hname)
    hsig = sigfile.Get(hname)
    hsig.Scale(sigSF)

    print('Integrals sig: {} bg: {}'.format(hsig.Integral(), hbg.Integral()))
    
    hbg.SetFillStyle(3345)
    hbg.SetFillColor(ROOT.kAzure+7)
    hbg.SetLineColor(hbg.GetFillColor())
    hsig.SetFillStyle(3354)
    hsig.SetFillColor(ROOT.kPink+10)
    hsig.SetLineColor(hsig.GetFillColor())
    
    stack = ROOT.THStack()
    stack.Add(hbg)
    stack.Add(hsig)
    
    hbgAlt = bgfileAlt.Get(hname)
    hsigAlt = sigfileAlt.Get(hname)
    hsigAlt.Scale(sigSF)
        

    # make data by a simple addition:
    # works only for 1D: hdata = hbg + hsig

    # so far the null hypothesis:
    hdata = hbgAlt.Clone(hbgAlt.GetName() + '_data')
    hdata.SetFillStyle(-1)
    hdata.SetLineColor(ROOT.kBlack)
    hdataNull = hdata.Clone(hdata.GetName() + 'copyNull')

    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetMarkerSize(1)
    hdata.SetMarkerStyle(20)
    stuff.append([hdata, hdataNull, hbg, hsig])
    
    samples = [hsig, hbg]
    # bg only fit:
    #samples = [hbg]
    fitresultNull = GetFractionFitterResult(hdata, samples)
    
    hdata.Add(hsigAlt)

    can = ROOT.TObject()
    leg = ROOT.TObject()
    opt = 'e1x0'

    legopt = 'F'
    if is2d:
        opt = 'hist'
        hbg.SetMarkerColor(hbg.GetFillColor())
        hsig.SetMarkerColor(hsig.GetFillColor())
        hdata.SetMarkerStyle(24)
        hbg.SetMarkerStyle(25)
        hsig.SetMarkerStyle(26)
        legopt = 'P'
    
    if draw:
        canname = 'FitFracNullAltPval'
        if is2d:
            canname = canname + '2d'
        else:
            canname = canname + '1d'
        can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 600)
        can.Divide(2,1)
        cans.append(can)
        can.cd(1)
        leg = ROOT.TLegend(0.55, 0.5, 0.88, 0.88)
        leg.SetBorderSize(0)
        hdata.Draw(opt)
        stack.Draw('hist same')
        hdata.Draw(opt + ' same')
        hdataNull.Draw('hist same')
        ROOT.gPad.RedrawAxis()
        leg.AddEntry(hdata, 'pseudodata', 'P')
        if not is2d:
            leg.AddEntry(hdataNull, 'Null pseudodata', 'L')
        leg.AddEntry(hsig, 'Signal', legopt)
        leg.AddEntry(hbg, 'Background', legopt)
        leg.Draw()
        #ROOT.gPad.Update()
        stuff.append(leg)
    
    fitresultSig = GetFractionFitterResult(hdata, samples)

    x1 = 0.
    x2 = 1.
    gform = '[0]*exp(-(x-[1])^2/(2*[2]^2))'
    
    pdfNull = ROOT.TF1('pdfNull', gform, x1, x2)
    sigmaNull = fitresultNull.error
    # note the factor of 2. as we center this gauss on 0.!
    pdfNull.SetParameters(2./sqrt(2*pi)/sigmaNull, 0., sigmaNull)
    pdfNull.SetLineColor(ROOT.kBlack)
    pdfNull.SetNpx(1000)
    
    pdfSig = ROOT.TF1('pdfSig', gform, x1, x2)
    mean = fitresultSig.result
    sigma = fitresultSig.error
    pdfSig.SetParameters(1./sqrt(2*pi)/sigma, mean, sigma)
    pdfSig.SetNpx(1000)

    p0  = 2*pdfNull.Integral(mean, 1.)
    t = 999.
    logt = 999.
    if p0 > 0:
        t = -log(p0)
        if t > 0.:
            logt = log(t)
    print('Null hypothesis probability: {}'.format(p0))
    print('t = -log(p0): {}'.format(t))
    print('log(t): {}'.format(logt))
    
    if draw:
        can.cd(2)
        pdfleg = ROOT.TLegend(0.5, 0.5, 0.88, 0.88)
        pdfleg.SetBorderSize(0)
        pdfNull.Draw('')
        pdfNull.GetXaxis().SetTitle('signal fraction')
        pdfSig.Draw('same')

        pdfleg.AddEntry(pdfNull, 'Null test stats p.d.f.', 'L')
        pdfleg.AddEntry(pdfSig, 'Sigal test stats p.d.f.', 'L')
        pdfleg.Draw()
        #ROOT.gPad.Update()
        stuff.append(pdfleg)
        
        can.Print(can.GetName() + '.pdf')
        can.Print(can.GetName() + '.png')

    return fitresultNull, fitresultSig, p0, t, logt, stuff
    
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

    hname ='scoreh'
    scoreh = ROOT.TH1D(hname, hname + ';score;replicas', 25, 0., 10)
    hname = 'nullsigmah'
    nullsigmah = ROOT.TH1D(hname, hname + ';null fit sigma;replicas', 25, 0., 0.25)
    hname = 'nullsigmah'
    sigmah = ROOT.TH1D(hname, hname + ';S+B fit sigma;replicas', 50, 0., 0.06)
    hname = 'nullsigmah'
    frach = ROOT.TH1D(hname, hname + ';S+B fitted signal fraction;replicas', 25, 0., 0.5)


    #for irep in range(1,2):
    for irep in range(1,99):
        print('=== Processinfg relica {}'.format(irep))
        
        # KEY STEETING -- HISTOGRAMME CHOICE!
        
        # 1D play"
        sigSF = 0.035
        hname = '1B1S/replicas/DetectorDiTopMass_rep{}'.format(irep)
        
        # 2D play:
        #sigSF = 0.015
        #hname = '1B1S/replicas/DetectorHTjVsDetectorDiTopMass_rep{}'.format(irep)

        fitresultNull, fitresultSig, p0, t, logt, stuff = Analyze(hname, sigSF, bgfile, sigfile, bgfileAlt, sigfileAlt, draw)
        Results.append([fitresultNull, fitresultSig, p0, t, logt])
        draw = False
        scoreh.Fill(logt)
        nullsigmah.Fill(fitresultNull.error)
        sigmah.Fill(fitresultSig.error)
        frach.Fill(fitresultSig.result)

    canname = 'scores'
    if 'Vs' in hname:
        canname = canname + '2d'
    else:
        canname = canname + '1d'
    scan = ROOT.TCanvas(canname, canname, 300, 0, 1200, 1200)
    scan.Divide(2,2)
    cans.append(scan)


    hists = [nullsigmah, sigmah, frach, scoreh]
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

