#!/usr/bin/python
# jk 2020

from __future__ import print_function

import ROOT

from Tools import *

from math import sqrt, pow, log

k2B0S = '2B0S'
k1B1S = '1B1S'
k0B2S = '0B2S'
kAnySel = 'AnySel'
kNoCuts = 'NoCuts'

# power of 2 to scale the signal in addition:
addSignalSFsPower = {#'2B0S' : 0.0625, '1B1S' : 0.125, '0B2S' : 1.,
    # to be read from list config files!
                k2B0S : 0., k1B1S : 0., k0B2S : 0,
    # zero means no scaling, not read from file:
                kAnySel : 0., kNoCuts : 0.}
    
##########################################
def ReplaceInStringToEmpty(ss, tags):
    rs = ss + ''
    for tag in tags:
        rs = rs.replace(tag, '')
    return rs

##########################################
def ComputeChi2AndKS(hdata, htot, x=0.13, y=0.73):
    # compute chi2 between data and MC:
    chi2 = ROOT.Double(0.)
    # ndf = hdata.GetNbinsX()
    # chi2 = hdata.Chi2Test(htot, "UU")
    ndf,chi2 = GetChi2(hdata, htot) # Tools
    ks = hdata.KolmogorovTest(htot)
    #if normalize == 'normalize':
    #    ndf = ndf - 1
    ctex = ROOT.TLatex(x, y, '')
    if ndf > 0:
        chi2ndf = chi2 / ndf
        ctex.SetText(x, y, '#chi^{2}/ndf=' + '{:2.2f}'.format(chi2ndf) + ' KS={:1.2f}'.format(ks))
        ctex.SetTextSize(0.05) # 0.055
        ctex.SetNDC()
    return ndf,chi2,ctex

##########################################
def DrawNice2DRatio(htot, hdata, gzratioMin, gzratioMax, stuff, iplot = 0, opt = 'box',  ratioTag = 'Pseudo-data / Prediction'):
            
        ratio = hdata.Clone(hdata.GetName() + '_ratio2d{}'.format(iplot))
        ratio.Divide(htot)
        # Print2DHisto(ratio)

        ratio.GetZaxis().SetRangeUser(gzratioMin, gzratioMax)
        ratio.Draw(opt)
        return ratio

#########################################################
def CopyStyle(h1, h2, copyFillAtts = True):
    h2.SetLineColor(h1.GetLineColor())
    h2.SetLineStyle(h1.GetLineStyle())
    h2.SetMarkerColor(h1.GetMarkerColor())
    h2.SetMarkerSize(h1.GetMarkerSize())
    h2.SetMarkerStyle(h1.GetMarkerStyle())
    if copyFillAtts:
        h2.SetFillStyle(h1.GetMarkerStyle())
        h2.SetFillColor(h1.GetMarkerStyle())



#########################################################
def MakeRatioHisto( data, prediction, newname = 'clone' ):
    ratio = data.Clone(newname)
    ratio.Reset()
    CopyStyle(data, ratio, False)
    nbins = data.GetNbinsX()
    #print 'nbins: %i' % (nbins,)
    i = 0
    for n in range( nbins ):
        #x_mc = ROOT.Double()
        y_mc = ROOT.Double()
        #x_data = ROOT.Double()
        y_data = ROOT.Double()
        #x_mc = prediction.GetBinCenter( n+1 )
        y_mc = prediction.GetBinContent( n+1 )   
        #print x_mc, y_mc
        if y_mc == 0.:
            continue
        x_data = data.GetBinCenter( n+1 )
        y_data = data.GetBinContent( n+1 )
        y_err = data.GetBinError( n+1 )
        ratio.SetBinContent( n+1, y_data/y_mc )
        ratio.SetBinError( n+1, y_err / y_mc)
        
        i += 1
    ratio.Scale(1.)
    return ratio


    
##########################################

def DrawNiceRatioWithBand(htot, hdata, hxmin, hxmax, gyratioMin, gyratioMax, stuff, iplot = 0, ratioTag = 'Pseudo-data / Prediction'):
            
        band = MakeOneWithErrors(htot)
        band.SetFillColor(ROOT.kYellow)
        #print('band:')
        #PrintBinContent(band)

        # double counting of visual errors on data and perdiction band!!!
        #ratio = hdata.Clone(hdata.GetName() + '_ratio{}'.format(iplot))
        #ratio.Divide(htot)
        ratio = MakeRatioHisto(hdata, htot, hdata.GetName() + '_ratio{}'.format(iplot))
        
        #print('ratio:')
        #PrintBinContent(ratio)
        ratioScaleHisto = ROOT.TH2D(ratio.GetName() + '_tmp', ratio.GetName() + '_tmp' + ';;' + ratioTag,
                                    ratio.GetNbinsX(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax(),
                                    100, gyratioMin, gyratioMax)
        ratioScaleHisto.SetStats(0)
        ratioScaleHisto.GetYaxis().SetTitleOffset(0.45)
        ratioScaleHisto.GetXaxis().SetLabelSize(0.085)
        ratioScaleHisto.GetYaxis().SetLabelSize(0.085)
        ratioScaleHisto.GetXaxis().SetTitle(MakePrettyTitle(ratio.GetXaxis().GetTitle()))
        ratioScaleHisto.GetXaxis().SetTitleSize(0.095)
        ratioScaleHisto.GetYaxis().SetTitleSize(0.095)
        ratioScaleHisto.Draw()
        ratioScaleHisto.GetXaxis().SetRangeUser(hxmin, hxmax)

        band.Draw('same e2')
        line = ROOT.TLine(hxmin, 1., hxmax, 1.)
        line.SetLineColor(ROOT.kRed)
        line.Draw()
        ratio.Draw('e1 X0 same')
        arrows = DrawArrowForPointsOutsideYAxisRange(ratio, ratioScaleHisto, hxmin, hxmax)
        stuff.append(arrows)

        stuff.append(line)
        stuff.append(ratioScaleHisto)
        
        
        return ratio, band, ratioScaleHisto


##########################################

def DivideByErrorBars(h, diviser):
    for i in xrange(1,h.GetXaxis().GetNbins()+1):
        hval = h.GetBinContent(i)
        herr = h.GetBinError(i)
        bgerr = diviser.GetBinError(i)
        if bgerr > 0:
            err = sqrt( pow(herr,2) + pow(bgerr,2) )
            h.SetBinContent(i, hval/err)
            h.SetBinError(i, 0.)
        else:
            h.SetBinContent(i, 0.)
            h.SetBinError(i, 0.)            
            print('ERROR: in DivideByErrorBars of histo {} hbg error is {}!'.format(h.GetName(), bgerr))
        h.Scale(1.)

##########################################        

def DrawSignificance(hbg, hdata, hxmin, hxmax, yMin, yMax, stuff, iplot = 0, ratioTag = 'Signal signif.'):
            
        signifh = hdata.Clone(hdata.GetName() + '_signif{}'.format(iplot))
        signifh.Add(hbg, -1.)
        DivideByErrorBars(signifh, hbg)
        #print('signifh:')
        #PrintBinContent(signifh)
        ratioScaleHisto = ROOT.TH2D(signifh.GetName() + '_tmp', signifh.GetName() + '_tmp' + ';;' + ratioTag,
                                    signifh.GetNbinsX(), signifh.GetXaxis().GetXmin(), signifh.GetXaxis().GetXmax(),
                                    100, yMin, yMax)
        ratioScaleHisto.SetStats(0)
        ratioScaleHisto.GetYaxis().SetTitleOffset(0.45)
        ratioScaleHisto.GetXaxis().SetLabelSize(0.085)
        ratioScaleHisto.GetYaxis().SetLabelSize(0.085)
        ratioScaleHisto.GetXaxis().SetTitle(MakePrettyTitle(signifh.GetXaxis().GetTitle()))
        ratioScaleHisto.GetXaxis().SetTitleSize(0.15)
        ratioScaleHisto.GetYaxis().SetTitleSize(0.095)
        ratioScaleHisto.Draw()
        ratioScaleHisto.GetXaxis().SetRangeUser(hxmin, hxmax)

        line0 = ROOT.TLine(hxmin, 0., hxmax, 0.)
        line0.SetLineColor(ROOT.kRed)
        line0.Draw()
        line1p = ROOT.TLine(hxmin, 1., hxmax, 1.)
        line1p.SetLineColor(ROOT.kRed)
        line1p.SetLineStyle(2)
        line1p.Draw()
        line1n = ROOT.TLine(hxmin, -1., hxmax, -1.)
        line1n.SetLineColor(ROOT.kRed)
        line1n.SetLineStyle(2)
        line1n.Draw()

        signifh.SetFillColor(ROOT.kRed)
        signifh.SetFillStyle(1111)
        signifh.Draw('hist same X0')
        
        stuff.append([line0, line1p, line1n])
        stuff.append(ratioScaleHisto)
        
        return signifh, ratioScaleHisto

    
##########################################
def ScaleHistAndRebin(hist, hname, w, rebin = True):
    if rebin and IsUniformlyBinned(hist):
        if not 'N' in hname:
            if 'DiTopM' in hname or 'Tau' in hname:
                    hist.Rebin(5)
            else:
                if not ('Delta' in hname or 'CosTheta' in hname or 'Yboost' in hname or 'Chittbar' in hname):
                    #hist.Rebin(4) ## default was 2!
                    #else:
                    if 'Pt' in hname:
                            #if 'LJet' in hname and 'Mass' in hname:
                            hist.Rebin(10)
                    elif 'Mass' in hname or 'Pout' in hname:
                            #if 'LJet' in hname and 'Mass' in hname:
                            hist.Rebin(5)
                    else:
                            hist.Rebin(20)
            
    #norm = ROOT.Double(hist.Integral())
    #if norm > 0.:
    #hist.Scale(w / norm)
    hist.SetStats(0)
    hist.Scale(w)
