#!/usr/bin/python
# jk 30.11.2017, Feb/March 2018, 

from __future__ import print_function

from math import sqrt, pow

import ROOT

from Tools import DivideByBinWidth

gRand = ROOT.TRandom3()

Legs = []
stuff = []

###################################################################
def ComputeChi2WrtVal(ratio, Normalized, const = 1.):
    nx = ratio.GetXaxis().GetNbins()
    chi2 = 0.
    ndf = 0
    for i in range(1, nx+1):
        val = ratio.GetBinContent(i)
        err = ratio.GetBinError(i)        
        if err > 0:
            chi2 = chi2 + pow((const - val)/err,2)
            ndf = ndf + 1
    if Normalized:
        ndf = ndf - 1
    return chi2,ndf

###################################################################
def MyMultiply(h, g, power = 1, tag = '_clone'):
    f = h.Clone(h.GetName() + tag)
    nx = h.GetXaxis().GetNbins()
    nx1 = g.GetXaxis().GetNbins()
    if nx != nx1:
        print('ERROR! histograms of different binnings! {:} vs {:}'.format(nx,nx1))
    for i in range(1, nx+1):
        val0 = h.GetBinContent(i)
        err0 = h.GetBinError(i)
        val1 = g.GetBinContent(i)
        err1 = g.GetBinError(i)
        val = val0*val1
        err = 0.
        if power == -1:
           if val1 > 0:
              val = val0/val1
           else:
              val = 0.
              err = 0.
        if val0 > 0 and val1 > 0:
            err = val*sqrt( pow(err0/val0,2) + pow(err1/val1,2) )
        #print(val0,err0,val1,err1,val,err)
        f.SetBinContent(i, val)
        f.SetBinError(i, err)
    f.Scale(1)
    return f

              
###################################################################        
# warning, fragile! Works for uniform binning only! 28.2.2018
def RemoveFirstFewBins(h2, nfirst):
    origname = h2.GetName()
    nx = h2.GetXaxis().GetNbins()
    ny = h2.GetYaxis().GetNbins()
    xmax = h2.GetXaxis().GetXmax()
    xmin = h2.GetXaxis().GetBinLowEdge(1+nfirst)
    ymax = h2.GetYaxis().GetXmax()
    ymin = h2.GetYaxis().GetBinLowEdge(1+nfirst)
    cropped = ROOT.TH2D(origname + '_cropped', origname, nx - nfirst, xmin, xmax, ny - nfirst, ymin, ymax)

    for i in range(1, nx+1):
        if i <= nfirst: continue
        for j in range(1, ny+1):
            if j <= nfirst: continue
            val = h2.GetBinContent(i,j)
            err = h2.GetBinError(i,j)
            cropped.SetBinContent(i-nfirst,j-nfirst,val)
            cropped.SetBinError(i-nfirst,j-nfirst,err)

    cropped.Scale(1.)
    h2.SetName(origname + '_orig')
    cropped.SetName(origname)
    return cropped

###################################################################
def PrintBinContent(histo):
    nx = histo.GetXaxis().GetNbins()
    line=''
    prefix = ''
    for binx in range(0,nx+2):
        line = '{:}{:} {:4.3f}'.format(line, prefix, histo.GetBinContent(binx),)
        prefix=','
    print(line)


###################################################################
def PrintCan(can, dirname = ''):
    xdirname = dirname + ''
    xdirname = xdirname.replace('pdf','png')
    can.Print(xdirname + can.GetName() + '.png')
    xdirname = xdirname.replace('png','pdf')
    can.Print(xdirname + can.GetName() + '.pdf')

###################################################################
def DrawRatios(hists,labels,idiv, Normalize = True, canname = 'Ratios', xlabel = '', addTag = '', dirname = ''):
    if idiv < 0 or idiv >= len(hists):
        idiv = 0
    can = ROOT.TCanvas(canname, canname, 200, 200, 800, 800)
    stuff.append(can)
    can.cd()
    ratios = []
    span1 = 1.
    span2 = 1.
    tmp = ROOT.TH2D('Ratios' + canname , 'Ratios;' + xlabel + ';Ratio to {:}'.format(labels[idiv].lower()),
                    100, hists[0].GetXaxis().GetXmin(), hists[0].GetXaxis().GetXmax(),
                    1000, 1.-span1, 1.+span2)
    tmp.SetStats(0)
    tmp.Draw()
    stuff.append(tmp)


    if Normalize:
        for hist in hists:
            hist.Scale(1./hist.Integral())
    
    opt = 'e1same'
    yshift = 0.44
    xshift = 0.0
    leg = ROOT.TLegend(0.15+xshift, 0.29+yshift, 0.55+xshift, 0.45+yshift)
    leg.SetName('RatioLeg')
    leg.SetBorderSize(0)
    ih = -1
    for hist,label in zip(hists,labels):
        ih = ih+1
        if ih == idiv:
            continue
        ratio = hist.Clone(hist.GetName() + '_ratio')
        ratio.Divide(hists[idiv])
        ratios.append(ratio)

        #if label != 'Unfolded':
        #    opt = 'sameL'
        #    leg.AddEntry(ratio, label, 'L')
        #else:
        #    opt = 'PLsame'

        chi2,ndf = ComputeChi2WrtVal(ratios[-1], Normalize, 1.)
        print('chi2={:} ndf={:} chi2/ndf={:}'.format(chi2, ndf, chi2/ndf))
        label = label + '  #chi^{2}/ndf' + ' = {:1.2f}'.format(chi2/ndf)
        
        leg.AddEntry(ratio, label, 'PL')
        print('Drawing ratio {:} with option {:}'.format(label, opt))
        ratio.Draw(opt + 'PL')
        PrintBinContent(ratio)
        opt = 'e1same'

    ROOT.gPad.Update()
    stuff.append(ratios)
    
    latex2 = ROOT.TLatex(0.60, 0.78, addTag)
    latex2.SetNDC()
    latex2.SetTextSize(0.035)
    latex2.Draw()
    stuff.append(latex2)
    line = ROOT.TLine(tmp.GetXaxis().GetXmin(), 1., tmp.GetXaxis().GetXmax(), 1.)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)
    line.Draw()
    leg.Draw()
    can.Update()
    PrintCan(can, dirname)
    stuff.append(leg)
    stuff.append(line)
    return can,ratios,[leg,latex2,line,tmp]
    
def PlotHisto(hist, opt = 'hist', canname = ''):
    if canname == '':
        canname = 'Can_' + hist.GetName()
    can = ROOT.TCanvas(canname, canname, 100, 100, 800, 800)
    stuff.append(can)
    can.cd()
    hist.Draw(opt)
    PrintCan(can)
    return can

def SmearData(hist, smear):
    data = hist.Clone(hist.GetName() + '_smeared')
    for i in range(1, hist.GetXaxis().GetNbins()+1):
        val = gRand.Poisson(hist.GetBinContent(i))
        err = sqrt(val)
        if smear:
            data.SetBinContent(i, val)
            data.SetBinError(i, err)
        else:
            # ok, just a clone, actually;)
            data.SetBinContent(i, hist.GetBinContent(i))
            data.SetBinError(i, hist.GetBinError(i))
    data.Scale(1.)
    return data



def MakeUnfoldedHisto(reco4bins, h1s, tag = 'h_unfolded'):
    hname = tag
    nbins = len(h1s)
    print ('Length of the h1s: {:}'.format(len(h1s)) )
    hist = reco4bins.Clone(hname)
    hist.Reset()
    i = -1
    for h1 in h1s:
        i = i+1
        val = h1.GetMean()
        hist.SetBinContent(i+1, val)
        hist.SetBinError(i+1, h1.GetRMS())
    stuff.append(hist)    
    return hist

###################################################################
def MakeTH1Ds(trace, nbins = 600, tag = 'trace'):
    h1s = []
    gmax = -999
    gmin = 1e19
    # find global min and max of the unfolded bin content:
    for line in trace:
        xmin = min(line)
        xmax = max(line)
        gmin = min(xmin, gmin)
        gmax = max(xmax, gmax)
    i = -1
    ### HACK!!!
    gmax = -999
    gmin = 1e19
    
    for line in trace:
        ### HACK!!!
        gxmin = min(line)
        gxmax = max(line)
        i = i+1
        hname = tag + '_{:}'.format(i)
        h1 = ROOT.TH1D(hname, hname, nbins, gmin, gmax)
        for val in line:
            h1.Fill(val)
        h1s.append(h1)
    stuff.append(h1s)
    return h1s

###################################################################
def PlotPosterior(psts, canname = 'PosteriorsDiag'):
    # only used in unfolding code, just after unfolding finishes
    can = ROOT.TCanvas(canname, canname, 100, 100, 800, 800)
    stuff.append(can)
    nx = int(sqrt(1.*len(psts)))
    ny = nx
    while nx*ny < len(psts):
        ny = ny+1
    can.Divide(nx,ny)
    i = -1
    for pst in psts:
        i = i+1
        can.cd(i+1)
        pst.Draw('hist')
    PrintCan(can)
    return can

###################################################################
def PlotIngredients(h_data, h_reco, h_ptcl, h_migra, canname = 'Ingredients', cantitle = 'Ingredients', level = 'Particle',
                    dirname = '', ms = 1.5):
    can = ROOT.TCanvas(canname, cantitle, 0, 0, 1200, 600)
    stuff.append(can)
    can.Divide(2,1)

    can.cd(1)
    ROOT.gPad.SetLogy(1)
   
    h_ptcl.SetLineColor(ROOT.kRed)
    h_ptcl.SetMarkerColor(ROOT.kRed)
    h_ptcl.SetMarkerStyle(20)
    h_ptcl.SetMarkerSize(ms)
    h_ptcl.Draw('e1')
    h_ptcl.GetYaxis().SetMoreLogLabels()

    # can.cd(2)
    h_reco.SetLineColor(ROOT.kBlue)
    h_reco.SetMarkerColor(ROOT.kBlue)
    h_reco.SetMarkerStyle(25)
    h_reco.SetMarkerSize(ms)
    h_reco.Draw('e1 same')
    
    #can.cd(3)
    h_data.SetMarkerSize(1.)
    h_data.SetMarkerStyle(26)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerSize(ms)
    h_data.Draw('e1same')
    
    leg = ROOT.TLegend(0.5, 0.5, 0.85, 0.85)
    leg.AddEntry(h_ptcl, level, 'PL')
    leg.AddEntry(h_reco, 'Detector', 'PL')
    leg.AddEntry(h_data, 'Pseudo data', 'PL')
    leg.SetBorderSize(0)
    Legs.append(leg)

    #h_response.Draw('colz')
    can.cd(2)
    h_migra.SetMinimum(0.)
    h_migra.SetMaximum(1.)
    h_migra.Draw('colz text')

    PrintCan(can, dirname)
    can.Update()
    return can,leg

###################################################################
def MakeProjectionXFromHisto(h2, tag = '_projX'):
    h2.ClearUnderflowAndOverflow()
    h2.GetXaxis().SetRange(1, h2.GetXaxis().GetNbins() )
    h2.GetYaxis().SetRange(1, h2.GetYaxis().GetNbins() )
    hist = h2.ProjectionX(h2.GetName() + tag)
    return hist

###################################################################
def MakeProjectionYFromHisto(h2, tag = '_projY'):
    h2.ClearUnderflowAndOverflow()
    h2.GetXaxis().SetRange(1, h2.GetXaxis().GetNbins() )
    h2.GetYaxis().SetRange(1, h2.GetYaxis().GetNbins() )
    hist = h2.ProjectionY(h2.GetName() + tag)
    return hist

###################################################################
def MakeListFromHisto(hist):
    vals = []
    for i in range(1, hist.GetXaxis().GetNbins()+1):
        val = hist.GetBinContent(i)
        vals.append(val)
    return vals

###################################################################
def NormalizeResponse(h2, tag = '_migra'):
    migra = h2.Clone(h2.GetName() + tag)
    for i in range(1, h2.GetXaxis().GetNbins()+1):
        sum = 0.
        for j in range(1, h2.GetYaxis().GetNbins()+1):
            val = h2.GetBinContent(i,j)
            sum = sum + val
        if sum > 0.:
            for j in range(1, h2.GetYaxis().GetNbins()+1):
                migra.SetBinContent(i,j,migra.GetBinContent(i,j) / sum)
    return migra

###################################################################
def MakeListResponse(h2):
    vals = []
    for i in range(1, h2.GetXaxis().GetNbins()+1):
        column = []
        for j in range(1, h2.GetYaxis().GetNbins()+1):
            val = h2.GetBinContent(i,j)
            column.append(val)
        vals.append(column)
    return vals
