#!/usr/bin/python
# Tue 23 Jul 09:00:00 CEST 2021

# compute the fit p-value
# either from the TFraction fitter
# or, in case of checking data to the background onyl hypothesis,
# as just a chi2 probability comparison test


from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp, pi

from BumpSignifTools import kBHInfty, kBHZero, kBHSFcmp

cans = []
stuff = []



##########################################

class fitres:
    def __init__(self, result,error,fitstatus,fitchi2,fitprob, fittedData):
        self.result = result
        self.error = error
        self.fitstatus = fitstatus
        self.fitchi2 = fitchi2
        self.fitprob = fitprob
        self.fittedData = fittedData
##########################################



##########################################
def GetFractionFitterResult(hdata, samples):

    #hdata = hdataOrig.Clone(hdataOrig.GetName() + '_forfit')
    xmin = hdata.GetXaxis().GetXmin()
    xmax = hdata.GetXaxis().GetXmax()
    
    nsamples = len(samples)
    mc_array = ROOT.TObjArray(nsamples)
    for sample in samples:
        mc_array.Add(sample)
    
    fitstatus = -1
    fitchi2 = -1
    fitprob = -1

    Result = []
    
    if len(samples) > 1:

        fitter = ROOT.TFractionFitter(hdata, mc_array)
        # constrain fractions to be between 0 and something
        for i in range(0, nsamples):
            fitter.Constrain(i, 0.0, 1.1)
        
        # Fit! ;-)
        fitstatus = fitter.Fit()
        fitchi2 = fitter.GetChisquare()
        fitprob = fitter.GetProb()
        result = ROOT.Double()
        error = ROOT.Double()
        for i in range(0, nsamples):
            fitter.GetResult(i, result, error)
            Result.append([1.*result, 1.*error])
            print('sample {} fraction: {} error: {} signif.: {}'.format(i, result, error, result/error))
    else:
        # do a chi2 comparison between histogrammes:
        model = samples[0]
        intmodel = model.Integral()
        intdata = hdata.Integral()
        # adjust to over/underflow bins and 2D?
        frac = intmodel / intdata
        error = 0. # try error propagation?
        ### !!! hacked by jk 29.7.2021
        # hdata.Scale(frac)
        # https://root.cern.ch/doc/master/classTH1.html#a6c281eebc0c0a848e7a0d620425090a5
        fitchi2 = hdata.Chi2Test(model, 'WW CHI2')
        fitprob = hdata.Chi2Test(model, 'WW') # by default returns the pval
        Result.append([frac, error]) # attempt some error

    return fitres(Result[0][0],Result[0][1],fitstatus,fitchi2,fitprob, hdata)

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

    # full data:
    hdata = hbgAlt.Clone(hbgAlt.GetName() + '_data')
    hdata.Add(hsigAlt)
    hdata.SetFillStyle(-1)
    hdata.SetLineColor(ROOT.kBlack)

    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetMarkerSize(1)
    hdata.SetMarkerStyle(20)
    stuff.append([hdata, hbg, hsig])
    
    # bg only fit:
    samples = [hbg]
    fitresultNull = GetFractionFitterResult(hdata, samples)
    
    p0  = fitresultNull.fitprob
    t0 = 999.
    logt0 = 999.
    if p0 > 0:
        t0 = -log(p0)
        if t0 > 0.:
            logt0 = log(t0)
    print('Null hypothesis probability p0: {}'.format(p0))
    print('t0 = -log(p0): {}'.format(t0))
    print('log(t0): {}'.format(logt0))

    samples = [hsig, hbg]
    
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
        canname = 'betterFitFracNullAltPval'
        if is2d:
            canname = canname + '2d'
        else:
            canname = canname + '1d'
        can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 600)
        can.Divide(2,1)
        cans.append(can)
        can.cd(2)
        leg = ROOT.TLegend(0.55, 0.6, 0.88, 0.88)
        leg.SetBorderSize(0)
        hdata.Draw(opt)
        stack.Draw('hist same')
        hdata.Draw(opt + ' same')
        ROOT.gPad.RedrawAxis()
        leg.AddEntry(hdata, 'pseudodata', 'P')
        leg.AddEntry(hsig, 'Signal', legopt)
        leg.AddEntry(hbg, 'Background', legopt)
        leg.Draw()
        #ROOT.gPad.Update()
        stuff.append(leg)

        can.cd(1)
        leg2 = ROOT.TLegend(0.35, 0.65, 0.88, 0.88)
        leg2.SetBorderSize(0)
        fittedData = fitresultNull.fittedData
        ymax = max(fittedData.GetMaximum(), hbg.GetMaximum())
        fittedData.SetMaximum(1.5*ymax)
        hbg.SetMaximum(1.5*ymax)
        fittedData.Draw(opt)
        hbg.Draw('hist same')
        ROOT.gPad.RedrawAxis()
        leg2.AddEntry(hdata, 'pseudodata scaled to bg.', 'P')
        leg2.AddEntry(hbg, 'Background', legopt)
        leg2.Draw()
        #ROOT.gPad.Update()
        stuff.append(leg2)

        can.Print(can.GetName() + '.pdf')
        can.Print(can.GetName() + '.png')

        
        
    # S+B fit to full data:
    fitresultSig = GetFractionFitterResult(hdata, samples)
    p1  = fitresultSig.fitprob
    t1 = 999.
    logt1 = 999.
    if p1 > 0:
        t1 = -log(p1)
        if t1 > 0.:
            logt1 = log(t1)
    print('Alternative hypothesis probability p1: {}'.format(p1))
    print('t1 = -log(p1): {}'.format(t1))
    print('log(t1): {}'.format(logt1))
    
    return fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1, stuff
    

###############################################
# jk VII/2021
def AnalyzeCompatibilityH0H1(hdata, hbg, hsig):

    hname = hdata.GetName()
    is2d = 'Vs' in hname
    print('Integrals sig: {} bg: {}'.format(hsig.Integral(), hbg.Integral()))
    
    # bg only fit:
    samples = [hbg]
    fitresultNull = GetFractionFitterResult(hdata, samples)
    
    p0  = fitresultNull.fitprob
    t0 = 999.
    logt0 = 999.
    if p0 > 0:
        t0 = -log(p0)
        if t0 > 0.:
            logt0 = log(t0)
    #print('Null hypothesis probability p0: {}'.format(p0))
    #print('t0 = -log(p0): {}'.format(t0))
    #print('log(t0): {}'.format(logt0))

    samples = [hsig, hbg]
    
    # S+B fit to full data:
    fitresultSig = GetFractionFitterResult(hdata, samples)
    p1  = fitresultSig.fitprob
    t1 = 999.
    logt1 = 999.
    if p1 > 0:
        t1 = -log(p1)
        if t1 > 0.:
            logt1 = log(t1)
    print('Alternative hypothesis probability p1: {}'.format(p1))
    print('t1 = -log(p1): {}'.format(t1))
    print('log(t1): {}'.format(logt1))
    
    return fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1, stuff


##########################################
def AnalyzeDataBgCompatibility(hdata, hbg, x, y, title = '', tst = 0.025):

    hname = hdata.GetName()
    is2d = 'Vs' in hname
    
    # bg only fit:
    samples = [hbg]
    fitresultNull = GetFractionFitterResult(hdata, samples)
    
    p0  = fitresultNull.fitprob
    t0 = 999.
    logt0 = 999.
    if p0 > 0:
        t0 = -log(p0)
        if t0 > 0.:
            logt0 = log(t0)
    #print('Null hypothesis probability p0: {}'.format(p0))
    #print('t0 = -log(p0): {}'.format(t0))
    #print('log(t0): {}'.format(logt0))

    tex = ROOT.TLatex(x, y, '')
    if t0 > kBHInfty*kBHSFcmp:
        tex.Selogt0ext(x, y, title + '+#infty')
    elif t0 < kBHZero:
        tex.SetText(x, y, title + '-#infty')
    else:
        tex.SetText(x, y, title + '{:2.1f}'.format(logt0))
    tex.SetTextSize(tst)
    tex.SetNDC()
    
    return fitresultNull, p0, t0, logt0, stuff, tex
    
##########################################
##########################################
##########################################
