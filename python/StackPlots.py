#!/usr/bin/python


# jk 3.3.2017, 30.3.2017, 7.4.2017, 6.9.2017, 2018, 2020
# Example running: ./python/StackPlots.py ./python/list.txt "_mytag" notnorm batch

import ROOT
import sys, os
from Tools import *
from stackPlotItems import *
from math import *

from array import array

from ROOT import std,Double,AddressOf

stuff = []

###!!!
#batch='runTheApp'
batch='batch'

###!!!
bw=0
#bw=1

#normalize='normalize'
### default
normalize='NO'

tag=''
doratio='ratio'
ratioTag = 'Pseudo data / Prediction'
idata = 0

#idataAlt = -1 # 1 or -1
#gopt = ['e1 X0', 'e1', 'F', 'F', 'F', 'F']
#gropt = ['e1 X0', 'e1', 'F', 'F', 'F', 'F']
#glegopt = ['P', 'P', 'F', 'F', 'F', 'F' ]

idataAlt = -1
gopt = ['e1 X0', 'F', 'F', 'F', 'F', 'F', 'F', 'F']
gropt = ['e1 X0', 'F', 'F', 'F', 'F', 'F', 'F', 'F']
glegopt = ['P', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F' ]


gcol = [
    #standard:
    ROOT.kBlack,
    ###ROOT.kRed,
    ROOT.kGreen+1,
    ROOT.kGreen+2,
    ROOT.kAzure+2,
    ROOT.kAzure,
    ROOT.kAzure+5,
    ROOT.kAzure+6,
    ROOT.kRed+1,
    ROOT.kSpring,
    ROOT.kSpring-4,
    #ROOT.kGreen+2,
    ROOT.kMagenta+2,
    ROOT.kCyan+1,
    ROOT.kBlack,
    ROOT.kPink+2,
    ROOT.kGray,
    ROOT.kYellow+2,
    ROOT.kOrange+2,
    
    ROOT.kBlue, #ROOT.kBlue+2, ROOT.kBlue+2,
    ROOT.kGreen+2,
    #ROOT.kGreen+2,
    ROOT.kRed,
    ROOT.kMagenta+2,
    ROOT.kCyan+1,
    ROOT.kBlack,
    ROOT.kPink+2,
    ROOT.kGray,
    ROOT.kYellow+2,
    ROOT.kOrange+2,


]

gwidth = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2,2,2]
gls = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2]
#gls = [1, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1]
gmsz = [1., 0., 0., 0., 0., 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
gmst = [20, 25, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 20, 21]
gyratioMin = 0.6
gyratioMax = 1.4
gySFlin = 1.9
gySFlog = 500.
gLegAddText = ['', '', '', '', '', '', '', '', '', '', '' ]


# Tools:
SampleWeightsParticle = [ 1., 1., 1., 1., 1., 1., 1., 1., 1.]
SampleWeightsDetector = [ 1., 1., 1., 1., 1., 1., 1., 1., 1.]


ROOT.gStyle.SetPalette(1)
# for precision of numbers in 2D when using the TEXT option:
ROOT.gStyle.SetPaintTextFormat("1.2f")

pngdir='python/stack_png/'
pdfdir='python/stack_pdf/' 
os.system('mkdir -p %s' % (pngdir,))
os.system('mkdir -p %s' % (pdfdir,))


# users, add your specific settings if needed:
if os.getenv('USER') == 'qitek':
    # JK specific settings
    pass
if os.getenv('USER') == 'pebaron':
    # add your specific settings
    pass

#os.system('notify-send "Running %s"' % (sys.argv[0],))
    
# get command line arguments:
cantag='cmp_'
if len(sys.argv) > 1:
    flist = sys.argv[1]
    cantag = flist
    cantag = cantag.replace('python/', '')
    cantag = cantag.replace('/', '_')
    cantag = cantag.replace('.txt', '')
    cantag = cantag.replace('list_', '')
    cantag = cantag.replace('._', '')
    print 'OK, using user-defined histograms root file list %s' % (flist,)
    # read and parse the list:
    fnames = []
    listfile = open(flist, 'r')
    print '  ...reading lines'
    for line in listfile.readlines():
        ffname = line[:-1]
        if len(ffname) > 5 and ffname[0] != '#' :
            fnames.append(ffname) # remove end line
    print '  Read:'
    for fname in fnames:
        print '   %s' % (fname,)
    
if len(sys.argv) > 2:
    tag = sys.argv[2]
    print 'OK, using user-defined histograms tag for output pngs %s' % (tag,)

if len(sys.argv) > 3:
    normalize = sys.argv[3]
    print 'OK, using user-defined normalize %s' % (normalize,)

if len(sys.argv) > 4:
    batch = sys.argv[4]
    print 'OK, using user-defined batch mode %s' % (batch,)


if batch == 'batch':
    ROOT.gROOT.SetBatch(1)
 
print '*** Settings:'
print 'tag=%s, normalize=%s, batch=%s' % (tag, normalize, batch,)
print 'files:'
print fnames
print


cans = []
Hists = []
Objs = []
rfiles = []
Legs = []
tags = {}
Pads = []

cw = 800
ch = 800


ROOT.gStyle.SetOptTitle(0)

Chi2Hists = {}
for key in ChiKeys:
    name = 'chi2_{}'.format(key)
    title = name + ';#chi^{2}'
    Chi2Hists[key] = ROOT.TH1D(name, title, 20, 0, 10)

for fname in fnames:
    rfile = ROOT.TFile(fname, 'read')
    rfiles.append(rfile)
    tag = fname.replace('analyzed_histos_', '').replace('histos_', '').replace('.root','')
    tags[fname] = tag

ToPlot = []
Pads = []
Stacks = []
for names in names1d:
    for name in names:
        print 'Processing 1D histo %s' % (name,)
        canname = cantag
        canname = canname + '_' + name
        canname = canname.replace('/', '_')
        can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
        cans.append(can)
        
        # prepare TPads:
        pad1,pad2,pad_inset = MakePadsStack(can, 'centre', 0.40, 0., 0., 0.)
        Pads.append([can,pad1,pad2,pad_inset])

        print '  Processing histo %s' % (name,)
        leg = ROOT.TLegend(0.42, 0.57, 0.88, 0.42+0.49 )
        legh = MakeLegHeader(name)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        Legs.append(leg)
        
        stack = ROOT.THStack(name + '_stack', name + '_stack')
        hists = []        
        for rfile in rfiles:
            ifile = rfiles.index(rfile)
            print '    Processing file %s' % (rfile.GetName())
            hname = name
            print '       getting %s' % (hname,)
            hist = rfile.Get(hname)
            hist.SetStats(0)
            hist.SetLineColor(gcol[ifile])
            hist.SetLineWidth(gwidth[ifile])
            hist.SetLineStyle(gls[ifile])
            if ifile == idata or ifile == idataAlt:
                hist.SetMarkerSize(gmsz[ifile])
                hist.SetMarkerColor(gcol[ifile])
                hist.SetMarkerStyle(gmst[ifile])
                print '     ...this should be data...'
            else:
                hist.SetFillColor(gcol[ifile])
                hist.SetFillStyle(1001)
                print '     ...this should be MC...'
            
            fracttxt = ' ({:.0f})'.format(hist.GetEntries())
            sampleWeight = 1.
            if 'Particle' in hist.GetName():
                 sampleWeight = SampleWeightsParticle[ifile]
            else:
                 sampleWeight = SampleWeightsDetector[ifile]
            if normalize == 'norm' or normalize == 'normalize' and ifile != idata and ifile != idataAlt:
                fracttxt = ' f={:1.2f}'.format(sampleWeight)    
            leg.AddEntry(hist, MakeNiceLegendEntry(tags[rfile.GetName()] + gLegAddText[ifile] + fracttxt, 0), glegopt[ifile])
            if not 'N' in hname:
                if 'DiTopM' in hname:
                    hist.Rebin(4)
                else:
                    if not ('Delta' in hname or 'CosTheta' in hname):
                        #hist.Rebin(4) ## default was 2!
                        #else:
                        hist.Rebin(2)
            hists.append(hist)  
            if normalize == 'norm' or normalize == 'normalize' and ifile != idata and ifile != idataAlt:
                # fragile code, means idata must be 0;)
                norm = Double(hist.Integral())
                if norm > 0.:
                    hist.Scale(sampleWeight*hists[idata].Integral() / norm)
            if ifile != idata and ifile != idataAlt:
                stack.Add(hist)
        # loop over files
        #if idataAlt >= 0:
        ToPlot.append( [ stack, hists[idata], [leg, legh], Pads[-1], hists[idataAlt], hname ] )
        #else:
        #    ToPlot.append( [ stack, hists[idata], leg, Pads[-1] ] )
        Hists.append(hists)
        Stacks.append(stack)

### now draw the stack and data ###
tmpScale = []
lines = []
ratios = []
for toplot in ToPlot:
    stack = toplot[0]
    hdata = toplot[1]
    hdataAlt = toplot[4]
    leg = toplot[2][0]
    legh = toplot[2][1]
    pads = toplot[3]
    can = pads[0]
    pad = pads[1]
    fullhname = toplot[5]
    
    pad.cd()
    hdata.SetMaximum(gySFlin*hdata.GetMaximum())
    hdata.SetMinimum(0.)
    if 'DiTopM' in hdata.GetName():
        hdata.GetXaxis().SetRangeUser(300, 2000.)
    hdata.Draw(gopt[idata])
    stack.Draw('samehist')
    hdata.Draw(gopt[idata] + 'same')
    if idataAlt >= 0:
        hdataAlt.Draw(gopt[idataAlt] + 'same')

    leg.Draw()
    
    ltex = ROOT.TLatex(0.13, 0.86, 'pp #sqrt{s} = 14 TeV')
    ltex.SetTextSize(0.055)
    ltex.SetNDC()
    ltex.Draw()
    Ltex.append(ltex)
    gentxt = 'MadGraph5'
    if 'Detector' in hdata.GetName():
        gentxt += ' + Delphes'
    ltex = ROOT.TLatex(0.13, 0.9395, gentxt)
    ltex.SetTextSize(0.07)
    ltex.SetNDC()
    ltex.Draw()
    Ltex.append(ltex)
    Legs.append(leg)

    # the total stack!
    last = stack.GetHists().At(stack.GetNhists()-1)
    htot = last.Clone(last.GetName() + '_tot')
    for ih in range(0, stack.GetNhists()-1):
        htot.Add(stack.GetHists().At(ih))

    # compute chi2 between data and MC:
    chi2 = ROOT.Double(0.)
    # ndf = hdata.GetNbinsX()
    # chi2 = hdata.Chi2Test(htot, "UU")
    ndf,chi2 = GetChi2(hdata, htot) # Tools
    ks = hdata.KolmogorovTest(htot)
    if normalize == 'normalize':
        ndf = ndf - 1
    if ndf > 0:
        chi2ndf = chi2 / ndf
        ctex = ROOT.TLatex(0.13, 0.73, '#chi^{2}/ndf=' + '{:2.2f}'.format(chi2ndf) + ' KS={:1.2f}'.format(ks))
        ctex.SetTextSize(0.055)
        ctex.SetNDC()
        ctex.Draw()
        Ltex.append(ctex)
        for key in ChiKeys:
            if key in fullhname:
                print('OK, filling chi2 hist with {}'.format(chi2/ndf))
                Chi2Hists[key].Fill(chi2 / ndf)


    # channel name into plot
    chtex = ROOT.TLatex(0.13, 0.79, legh)
    chtex.SetTextSize(0.055)
    chtex.SetNDC()
    chtex.Draw()
    Ltex.append(chtex)
    
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    #ROOT.gPad.GetFrame().Draw()

    
    rpad = pads[2]
    rpad.cd()
    #print 'hdata:'
    #PrintBinContent(hdata)
    #print 'Stacked histos: {:}'.format(stack.GetNhists())
   
    
  
    #print 'htot:'
    #PrintBinContent(htot)
    band = MakeOneWithErrors(htot)
    band.SetFillColor(ROOT.kYellow)
    #print 'band:'
    #PrintBinContent(band)
    # todo: TDir is removed from name, mem leak in replacement...:
    ratio = hdata.Clone(hdata.GetName() + '_ratio')
    #print 'ratio:'
    #PrintBinContent(ratio)
    ratio.Divide(htot)
    ratios.append(band)
    ratios.append(ratio)

    if idataAlt >= 0:
        ratioAlt = hdataAlt.Clone(hdataAlt.GetName() + '_ratioAlt')
        ratioAlt.Divide(htot)
        ratios.append(ratioAlt)
    
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
    if 'DiTopM' in htot.GetName():
        ratioScaleHisto.GetXaxis().SetRangeUser(300, 2000.)

    band.Draw('same e2')
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
    line.SetLineColor(ROOT.kRed)
    line.Draw()
    ratio.Draw('e1 X0 same')
    arrows = DrawArrowForPointsOutsideYAxisRange(ratio, ratioScaleHisto)
    stuff.append(arrows)
    if idataAlt >= 0:
        ratioAlt.Draw('e1 same')

    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    #ROOT.gPad.GetFrame().Draw()

    pad.cd()
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    
    lines.append(line)
    tmpScale.append(ratioScaleHisto)

    # ptag = tag
    # better, shorter yet sufficient name;-)
    ptag = ''
    can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
    can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')
    hdata.SetMaximum(gySFlog*hdata.GetMaximum())
    hdata.SetMinimum(0.1)
    pad.SetLogy(1)
    pad.Update()
    ROOT.gPad.RedrawAxis()
    can.Print(pngdir + can.GetName() + ptag + '_logy' + '.png')
    can.Print(pdfdir + can.GetName() + ptag + '_logy' + '.pdf')


canname = 'Chi2Hists'
can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
cans.append(can)
opt = 'hist'
chileg = ROOT.TLegend(0.5, 0.5, 0.84, 0.84)
chileg.SetBorderSize(0)
chicols = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+3, ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+3]
i = 0
ymax  =-1
for key in Chi2Hists:
    chi2h = Chi2Hists[key]
    val = chi2h.GetMaximum()
    if val > ymax:
        ymax = val
for key in Chi2Hists:
    chi2h = Chi2Hists[key]
    chi2h.SetLineColor(chicols[i])
    chi2h.SetLineWidth(2)
    chi2h.SetStats(0)
    chi2h.SetMaximum(ymax*1.1)
    chi2h.Draw(opt)
    chileg.AddEntry(chi2h, '{:10} #mu={:2.2f}'.format(key, chi2h.GetMean()) , 'L')
    opt = 'hist same'
    i += 1
chileg.Draw()
can.Print(pngdir + can.GetName() + '.png')
can.Print(pdfdir + can.GetName() + '.pdf')

os.system('mkdir {}/{}/'.format(pngdir,cantag))
os.system('mkdir {}/{}/'.format(pdfdir,cantag))
os.system('mv {}*.png {}/{}/'.format(pngdir,pngdir,cantag))
os.system('mv {}*.pdf {}/{}/'.format(pdfdir,pdfdir,cantag))


os.system('notify-send "DONE! %s"' % (sys.argv[0],))
    
if batch != 'batch':
    ROOT.gApplication.Run()

