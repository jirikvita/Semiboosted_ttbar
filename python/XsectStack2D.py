#!/usr/bin/python

from __future__ import print_function

# jk 14.4.2020 (covid-19 homeoffice), based on code of StackPlots.py
# Example running: ./python/XsectStack.py ./python/list.txt "_mytag" notnorm batch
# 10.9.2020


import ROOT

from xSectTools import *
from BumpSignifTools import *
from FitTools import *
from CorrItems import *
from cTopPlot import *

from LhoodFitTools import *

from math import *
import sys, os
from array import array


##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    stuff = []

    debug = 0
    
    #batch='runTheApp'
    batch='batch'
    
    normalize='normalize'
    # scale to this lumi!
    lumi = 10000 # 1000000 # pb! 
    csf = 1.

    doDivideByBinArea = False
    printStats = False
    plotRatios = False
    plotProjections = False
    
    tag=''
    idata = 0

    gyratioMin = 0.2
    gyratioMax = 1.8
    gySFlin = 1.9
    gySFneg = 0.25
    gySFlog = 500.
    gyLogMin = 1.e-5

    # for 2D
    ROOT.gStyle.SetPalette(1)
    # for precision of numbers in 2D when using the TEXT option:
    ROOT.gStyle.SetPaintTextFormat("1.2f")
    ROOT.gStyle.SetPadTopMargin(0.25)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    
    pngdir='python/stack2d_png/'
    pdfdir='python/stack2d_pdf/' 
    os.system('mkdir -p {}'.format(pngdir,))
    os.system('mkdir -p {}'.format(pdfdir,))

    # users, add your specific settings if needed:
    if os.getenv('USER') == 'qitek':
        # JK specific settings
        pass

    #os.system('notify-send "Running {}"'.format(argv[0],))
    
    # get command line arguments:
    cantag='cmp_'
    signalFileName = ''
    sigLegTitle = ''
    SignalSF = ''
    sigTag = '' # later a short signal tag for TLaTeX
    
    if len(argv) > 1:
        flist = argv[1]
        cantag = flist
        cantag = cantag.replace('/', '_') # the order is important!
        cantag = ReplaceInStringToEmpty(cantag, ['python_', '.txt', 'lists_', 'list_', '._'])
        print('Reading signal file details from config {}'.format(flist,))
        # read and parse
        print('  ...reading lines')
        lines = []
        cfgfile = open(flist, 'r')
        for line in cfgfile.readlines():
            sline = line[:-1]
            lines.append(sline) # remove end line
        print('  Read:')
        for line in lines:
            print('   {}'.format(line,))
        if len(lines) < 3:
            print('Config file should contain signal file name, legend title and signal SF, i.e. at last 3 lines!')
            return 
        signalFileName = lines[0]
        sigLegTitle = lines[1]
        
        signalSFsStr = lines[2].split()
        SignalSF = int(signalSFsStr[0])
        # this order must be respected in file!
        addSignalSFsPower[k2B0S] = float(signalSFsStr[1])
        addSignalSFsPower[k1B1S] = float(signalSFsStr[2])
        addSignalSFsPower[k0B2S] = float(signalSFsStr[3])
        
        if len(lines) > 3:
            sigTag = lines[3]
    else:
        print('Usage: {} config.txt'.format(argv[0]))
        return
    # NEW! 14.4.2020:
    print('*** Initializing samples and cross-section weigths... ***')
    stackItems = MakeStackItems(signalFileName,
                                sigLegTitle,
                                lumi, csf, SignalSF) # lumi, csf ... common SF, ssf = signal SF

    xsdict = MakeXsextDic()
    #print(xsdict)
    fnames = stackItems.keys()
    fnames.pop(fnames.index(gStackName))
    # make sure signal is the last in list!
    isig = fnames.index(signalFileName)
    if isig != len(fnames)-1:
        fnames[-1],fnames[isig] = fnames[isig],fnames[-1]
    
    if batch == 'batch':
        ROOT.gROOT.SetBatch(1)

    ROOT.gStyle.SetPalette(ROOT.kSolar)
    #ROOT.gStyle.SetPalette(ROOT.kCherry)
    #ROOT.gStyle.InvertPalette()
    print('*** Settings:')
    print('tag={}, normalize={}, batch={}'.format(tag, normalize, batch,))
    print('files:')
    print(fnames)
    print('')

    cans = []
    Hists = []
    Objs = []
    rfiles = []
    rfilesStatIndep = []
    Legs = []
    tags = {}
    Pads = []

    cw = 800
    ch = 800
    if plotRatios or plotProjections:
        ch = 2*cw
    lx1 = 0.49
    ly1 = 0.57
    lx2 = 0.88
    ly2 = 0.42+0.49

    ROOT.gStyle.SetOptTitle(0)

    #Chi2Hists = {}
    #for key in ChiKeys:
    #    name = 'chi2_{}'.format(key)
    #    title = name + ';#chi^{2}'
    #    Chi2Hists[key] = ROOT.TH1D(name, title, 50, 0, 25)

    for fname in fnames:
        rfile = ROOT.TFile(fname, 'read')
        rfiles.append(rfile)
        tag = fname.replace('analyzed_histos_', '').replace('histos_', '').replace('.root','')
        tags[fname] = tag

        fnameStatIndep = fname + ''
        if 'half0' in fname:
            fnameStatIndep = fnameStatIndep.replace('half0','half1')
        elif 'half1' in fname:
            fnameStatIndep = fnameStatIndep.replace('half1','half0')
        else:
            print('ERROR finding half0 or half1 in filename! Need this to define the complementary sample!')
            return
        rfileStatIndep = ROOT.TFile(fnameStatIndep, 'read')
        rfilesStatIndep.append(rfileStatIndep)

    ToPlot = []
    Pads = []
    sPads = []
    Stacks = []

      
    print('*** Processing 2D histograms to stack;-) ***')
    Topos = [ '0B2S', '1B1S', '2B0S']
    
    # HACK!!! for quick studies
    names = [ #'TopPtVsDiTopMass',
        #'CosThetaStarVsDiTopPout',
        'HTjVsDiTopMass',
    ]
    #for xname in names:
    for xname in hnamesDict:
        for topo in Topos:
            name = topo + '/' + ('Detector' + xname).replace('Vs', 'VsDetector')

            # HACK!!!
            #if not ('DiTopMass' in name or 'DiTopPt' in name):
            #    continue

            
            print('Processing 1D histo {}'.format(name,))
            canname = cantag
            canname = canname + '_' + name
            canname = canname.replace('/', '_')
            can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
            cans.append(can)

            # prepare TPads:

            bcanratio = 0.55
            PadSeparation = 0.0
            UpperPadBottomMargin = 0.1
            LowerPadTopMargin = 0.1
            if not (plotRatios or plotProjections):
                bcanratio = 0.
                PadSeparation = 0.0
                UpperPadBottomMargin = 0.2
                LowerPadTopMargin = 0.

                
            pad1,pad2,pad_inset = MakePadsStack(can, 'centre', bcanratio,  PadSeparation, UpperPadBottomMargin, LowerPadTopMargin)
            Pads.append([can,pad1,pad2,pad_inset])
            
            # prepare TPads, for bg-subtraction, i.e. signal-only
            scanname = canname + '_bgSub'
            scan = ROOT.TCanvas(scanname, scanname, 200, 200, cw, ch)
            cans.append(scan)
            spad1,spad2,spad_inset = MakePadsStack(scan, 'centre', bcanratio, PadSeparation, UpperPadBottomMargin, LowerPadTopMargin)
            sPads.append([scan,spad1,spad2,spad_inset])

            print('  Processing histo {}'.format(name,) )
            leg = ROOT.TLegend(lx1, ly1, lx2, ly2 )
            legitems = []
            legh = MakeLegHeader(name)
            topotag = MakeTopoTag(name)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            Legs.append(leg)

            stack = ROOT.THStack(name + '_stack', name + '_stack')
            hists = []
            stackStatIndep = ROOT.THStack(name + '_stack_StatIndep', name + '_stack_stack_StatIndep')
            for rfile,rfileStatIndep in zip(rfiles, rfilesStatIndep):
                ifile = rfiles.index(rfile)
                fname = rfile.GetName()
                isSignal = rfile.GetName() == signalFileName
                print('    Processing file {}'.format(rfile.GetName()))
                hname = name
                print('       getting {}'.format(hname,))
                hist = rfile.Get(hname)
                # HACK jk 21.4.2021
                if hist.GetXaxis().GetTitle() == "Mass [GeV]":
                    hist.GetXaxis().SetTitle("DiTop m [GeV]")
                if hist.GetYaxis().GetTitle() == "Mass [GeV]":
                    hist.GetYaxis().SetTitle("DiTop m [GeV]")
                if doDivideByBinArea:
                    DivideByBinArea(hist)
                hist.SetLineColor(stackItems[fname].lcol)
                hist.SetLineWidth(stackItems[fname].lw)
                hist.SetLineStyle(stackItems[fname].lst)
                hist.SetFillColor(stackItems[fname].fcol)
                hist.SetFillStyle(stackItems[fname].fst)

                # now get the other version of the histogram from the complementary sample:
                histStatIndep = rfileStatIndep.Get(hname)
                if doDivideByBinArea:
                    DivideByBinArea(histStatIndep)
                
                sampleWeight = stackItems[fname].weight
                if isSignal:
                    print('       USING ADDITIONAL SAMPLE WEIGHT 2^{} based on topology {}'.format(addSignalSFsPower[topotag], topotag))
                    sampleWeight = sampleWeight*pow(2,addSignalSFsPower[topotag])

                fracttxt = ''
                if fabs(stackItems[fname].sf - 1.) > 1.e-4:
                    powtag = '*2^{' + str(addSignalSFsPower[topotag]) + '}'
                    if abs(addSignalSFsPower[topotag]) < 1.e-5:
                        powtag = ''
                    fracttxt = fracttxt + ' (#times{}'.format(stackItems[fname].sf) + powtag + ')'
                legitems.append(cLegItem(hist, stackItems[fname].legtag + fracttxt, stackItems[fname].lopt))

                # HACK: do not use this sample except 0B2S topology due to stat insuff. and fluctuations in 1B1S abd 2B0S!
                # jk 19.7.2021, here from 8.7.2021
                if 'pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200' in fname and not '0B2S' in hname:
                    continue
                
                if debug > 0:
                        print('* Scaling histos by weight {}'.format(sampleWeight))
                ScaleHistAndRebin(hist, hname, sampleWeight, False)
                ScaleHistAndRebin(histStatIndep, hname, sampleWeight, False)
                if debug > 1:
                        print('    stack0 adding nbins: {} N={} mu={} sigma={} min={} max={}'.format(hist.GetXaxis().GetNbins(), hist.GetEntries(),
                                                                                                     hist.GetMean(), hist.GetRMS(),
                                                                                                     hist.GetMinimum(), hist.GetMaximum()) )

                        print('    stack1 adding nbins: {} N={} mu={} sigma={} min={} max={}'.format(histStatIndep.GetXaxis().GetNbins(), histStatIndep.GetEntries(),
                                                                                                     histStatIndep.GetMean(), histStatIndep.GetRMS(),
                                                                                                     histStatIndep.GetMinimum(), histStatIndep.GetMaximum()) )
                stack.Add(hist)
                stackStatIndep.Add(histStatIndep)
                hists.append([hist, histStatIndep])  
                    
            # loop over files
            print('Pads:')
            print(Pads[-1])
            ToPlot.append( cTopPlot(stack, stackStatIndep, leg, legh, legitems, Pads[-1], sPads[-1], hname) )
            Hists.append(hists)
            Stacks.append(stack)

    ### now draw the stack and data! ;-) ###
    ### TODO: subtract BG!
    
    ratios = []
    iplot = -1
    print('*** Plotting stacked histograms;) ***')
    for toplot in ToPlot:
        iplot = iplot+1
        stack = toplot.stack
        hdata = toplot.hdata
        hdata.SetMarkerStyle(stackItems[gStackName].mst)
        hdata.SetMarkerSize(stackItems[gStackName].msz)
        hdata.SetMarkerColor(stackItems[gStackName].mcol)
        hdata.SetLineColor(stackItems[gStackName].lcol)
        leg = toplot.leg
        legh = toplot.legh

        pads = toplot.pads
        can = pads[0]  # canvas
        pad = pads[1]  # upper pad
        rpad = pads[2] # ratio pad
        print(pads)
        
        # signal only, background-subtracted:
        spads = toplot.spads
        scan = spads[0]  # canvas
        spad = spads[1]  # upper
        srpad = spads[2] # ratio
        
        fullhname = toplot.fullhname

        pad.cd()
        hdata.SetMaximum(gySFlin*hdata.GetMaximum())
        #hdata.SetMinimum(0.1)
        
        hxmin = hdata.GetXaxis().GetXmin()
        hxmax = hdata.GetXaxis().GetXmax()

        if DivideByBinArea:
            hdata.GetZaxis().SetTitle('Expected events / #Delta^{2}_{ij}')
        else:
            hdata.GetZaxis().SetTitle('Expected events')
            
        # the total stack!
        last = stack.GetHists().At(stack.GetNhists()-1)
        htot = last.Clone(last.GetName() + '_tot')
        hbg = last.Clone(last.GetName() + '_bg')
        hbg.Reset()
        for ih in range(0, stack.GetNhists()-1):
            htot.Add(stack.GetHists().At(ih))
            hbg.Add(stack.GetHists().At(ih))
        hsig = toplot.hsig
        dopt = 'box' #stackItems[gStackName].dopt
        #hdata.Draw(dopt)
        #stack.Draw(dopt + 'same')
        #hdata.Draw(dopt + 'same')

        hsig.SetLineColor(ROOT.kRed)
        hsig.SetFillColor(ROOT.kRed)
        htot.SetFillColor(ROOT.kBlue)
        htot.SetFillStyle(1111)
        htot.SetLineColor(ROOT.kBlack)
        
        hsig.Scale(1.)
        htot.Scale(1.)
        
        pad.cd()
        #htot.SetLineColor(ROOT.kGreen+2)
        #htot.Draw('box')
        FixXaxisTitle(htot, 'x')
        FixXaxisTitle(htot, 'y')
        htot.Draw(dopt)
        hsig.Draw(dopt + ' same')

        ctex_corr = ROOT.TLatex(0.02, 0.02,'Correlations:')
        ctex_corr.SetNDC()
        ctex_corr.SetTextSize(0.040)
        
        corr_bg = hbg.GetCorrelationFactor()
        ctex_bg = ROOT.TLatex(0.32, 0.02,'bg. #rho=' + '{:1.2f}'.format(corr_bg))
        ctex_bg.SetNDC()
        ctex_bg.SetTextSize(0.040)
        ctex_bg.SetTextColor(htot.GetFillColor())
        
        corr_sig = hsig.GetCorrelationFactor()
        ctex_sig = ROOT.TLatex(0.54, 0.02,'sig. #rho=' + '{:1.2f}'.format(corr_sig))
        ctex_sig.SetNDC()
        ctex_sig.SetTextSize(0.040)
        ctex_sig.SetTextColor(hsig.GetFillColor())

        ctex_corr.Draw()
        ctex_bg.Draw()
        ctex_sig.Draw()
            

        stuff.append([ctex_sig, ctex_bg])
        
        #PrintBinContent2D(hbg)
        #PrintBinContent2D(hsig)
        ROOT.gPad.Update()
        
        nb = hdata.GetNbinsX()
        stats = ''
        #if printStats:
        #    stats =  ' N={:.0f} I={:.0f}'.format(hdata.GetEntries(),hdata.Integral(0, nb+1))
        leg.AddEntry(hdata, stackItems[gStackName].legtag + stats, stackItems[gStackName].lopt)
        
        for ileg in range(len(legitems)-1, -1, -1):
            legitem = toplot.legitems[ileg]
            stats = ''
            #if printStats:
            #    stats = ' N={:.0f} I={:.0f}'.format(legitem.hist.GetEntries(), legitem.hist.Integral(0, nb+1))
            leg.AddEntry(legitem.hist, legitem.legtag + stats, legitem.lopt)
        ###leg.Draw()

        ltex = ROOT.TLatex(0.02, 0.02, 'pp #sqrt{s} = 14 TeV  ' + 'L = {:.0f} ab'.format(lumi/1e6) + '{}^{-1}')
        ltex.SetTextSize(0.045)
        ltex.SetNDC()
        ###ltex.Draw()
        Ltex.append(ltex)

        gentxt = 'MadGraph5'
        if 'Detector' in hdata.GetName():
            gentxt += ' + Delphes'
        gltex = ROOT.TLatex(0.16, 0.9395, gentxt)
        gltex.SetTextSize(0.045)
        gltex.SetNDC()
        gltex.Draw()
        Ltex.append(gltex)
        Legs.append(leg)

        # channel name into plot
        chtex = ROOT.TLatex(0.60, 0.9395, legh)
        chtex.SetTextSize(0.045)
        chtex.SetNDC()
        chtex.Draw()
        Ltex.append(chtex)

        ndf,chi2,ctex = ComputeChi2AndKS(hdata, htot, 0.13, 0.73)

        # BumpHunter2D:
        ytxt = 0.09 # 0.08, 0.02

        bhcol = ROOT.kGreen + 1 #2
        bhhcol = ROOT.kGreen + 1 #2

        bhcolx =  ROOT.kOrange + 7 #ROOT.kYellow + 2
        bhhcolx = ROOT.kOrange + 7 #ROOT.kYellow + 2

        bhcoly = ROOT.kMagenta  # + 2
        bhhcoly = ROOT.kMagenta # + 2

        bhtx = 0.32
        bhresult,bhtex = FindBestBumpHunter2DArea(hdata, hbg, bhtx, ytxt, 'xy:   ', 0.025)
        used2dBH = bhresult.bins
        usedLinesBH = DrawUsedBinsGrid(hdata, used2dBH, bhresult.i, bhresult.j, bhcol, bhhcol, 1)

        # JK 29.7.2021
        # compute the full 2D range chi2 probability and its log of data compatibility with the null (H0) or alternatve (H1) hypothesis
        #fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1, stuffbg = AnalyzeCompatibilityH0H1(hdata, hbg, hsig)
        # compatibility of data shape with bg only:
        fitresultNull, p0, t0, logt0, stuffbg, bgtext2d = AnalyzeDataBgCompatibility(hdata, hbg, bhtx, 0.06, 'xy: ')
        print('2D Null hypothesis probability p0: {}'.format(p0))
        print('2D t0 = -log(p0): {}'.format(t0))
        print('2D log(t0): {}'.format(logt0))

        
        # Sinificance:
        # ytxt = 0.08
        #scol = ROOT.kMagenta + 2
        #shcol = ROOT.kMagenta + 2
        #significance, stex = ComputeExcessSignificance2D(hdata, hbg, 0.02, 0.08)
        #significance, stex, n2d, used2d, i0, j0 = ComputeExcessSignificance2DIteratively(hdata, hbg, 0.02, ytxt, 'Signif. xy: ' )
        #usedLines = DrawUsedBinsGrid(hdata, used2d, i0, j0, scol, shcol, 2)
        
         #1D's:
        hdatax = hdata.ProjectionX(hdata.GetName() + '_projX')
        hdatay = hdata.ProjectionY(hdata.GetName() + '_projY')
        hbgx = hbg.ProjectionX(hbg.GetName() + '_projX')
        hbgy = hbg.ProjectionY(hbg.GetName() + '_projY')
        ## old, full range: significance1Dx, stex1Dx = ComputeExcessSignificance1D(hdatax, hbgx, 0.02, 0.02, '1D signif. x: ')
        ## old, full range: significance1Dy, stex1Dy = ComputeExcessSignificance1D(hdatay, hbgy, 0.55, 0.02, '1D signif. y: ')
        #significance1Dx, stex1Dx, nx, usedx, i0x = ComputeExcessSignificance1DIteratively(hdatax, hbgx, 0.345, ytxt, 'x: ')
        #significance1Dy, stex1Dy, ny, usedy, i0y = ComputeExcessSignificance1DIteratively(hdatay, hbgy, 0.500, ytxt, 'y: ')

        # BH 1D:
        #bhresultx,bhtexx = FindBestBumpHunter1DInterval(hdatax, hbgx, 0.54, ytxt, 'x: ')
        #bhresulty,bhtexy = FindBestBumpHunter1DInterval(hdatay, hbgy, 0.74, ytxt, 'y: ')
        bhresultx,bhtexx = FindBestBumpHunter1DInterval(hdatax, hbgx, 0.54, ytxt, 'x: ', 0.025)
        bhresulty,bhtexy = FindBestBumpHunter1DInterval(hdatay, hbgy, 0.74, ytxt, 'y: ', 0.025)
        
        usedLinesBHx = DrawUsedBinsLines(hdatax, hdatay, bhresultx.bins, bhresultx.i, True,  0.08, 0.005, bhcolx, bhhcolx, 1, 3, 3)
        usedLinesBHy = DrawUsedBinsLines(hdatay, hdatax, bhresulty.bins, bhresulty.i, False, 0.10, 0.005, bhcoly, bhhcoly, 1, 3, 3)
        #usedLinesx   = DrawUsedBinsLines(hdatax, hdatay, usedx, i0x,  True,  0.10, 0.005, scol, shcol, 1)
        #usedLinesy   = DrawUsedBinsLines(hdatay, hdatax, usedy, i0y,  False, 0.12, 0.005, scol, shcol, 1)

        # compatibility of data shape with bg only 1Dx:
        fitresultNullx, p0x, t0x, logt0x, stuffbgx, bgtext1dx = AnalyzeDataBgCompatibility(hdatax, hbgx, 0.54, 0.06, 'x: ')
        print('1Dx Null hypothesis probability p0: {}'.format(p0x))
        print('1Dx t0 = -log(p0): {}'.format(t0x))
        print('1Dx log(t0): {}'.format(logt0x))
        # compatibility of data shape with bg only 1Dy:
        fitresultNully, p0y, t0y, logt0y, stuffbgy, bgtext1dy = AnalyzeDataBgCompatibility(hdatay, hbgy, 0.74, 0.06, 'y: ')
        print('1Dy Null hypothesis probability p0: {}'.format(p0y))
        print('1Dy t0 = -log(p0): {}'.format(t0y))
        print('1Dy log(t0): {}'.format(logt0y))

        bgtext2d.Draw()
        bgtext1dx.Draw()
        bgtext1dy.Draw()
        
        #stex1Dx.Draw()
        #stex1Dy.Draw()
        #if significance > significance1Dy and significance > significance1Dx:
        #    stex.SetTextColor(hsig.GetFillColor())
        #stex.Draw()

        bhtex.SetTextColor(bhhcol)
        ctex_bh = ROOT.TLatex(0.02, ytxt,'BH log(t)')
        if bhresult.t > bhresultx.t and bhresult.t > bhresulty.t:
            ctex_bh.SetTextColor(bhhcol)
        elif bhresultx.t > bhresulty.t:
            ctex_bh.SetTextColor(bhhcolx)
        else:
            ctex_bh.SetTextColor(bhhcoly)

        ctex_bh.SetNDC()
        ctex_bh.SetTextSize(0.040)
        ctex_bh.Draw()

        bhtex.Draw()

        bhtexx.Draw()
        bhtexx.SetTextColor(bhhcolx)
        bhtexy.Draw()
        bhtexy.SetTextColor(bhhcoly)
        stuff.append([bhtex, bhtexx, bhtexy])
        
        """
        if ndf > 0:
            ctex.Draw()
            Ltex.append(ctex)
            ctex.Draw()
            for key in ChiKeys:
                if key in fullhname:
                    print('filling chi2 hist with {}'.format(chi2/ndf) )
                    Chi2Hists[key].Fill(chi2 / ndf)
        """ 

        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        if plotRatios or plotProjections:
            rpad.cd()
            if plotRatios:
                ratio = DrawNice2DRatio(htot, hdata, gyratioMin, gyratioMax, stuff, 10000 + iplot, 'colz')
            else:
                hdatax.SetLineColor(ROOT.kRed)
                hdatax.SetMarkerColor(ROOT.kRed)
                hdatax.SetStats(0)
                hdatax.Draw('hbar')
                hbgx.SetLineColor(ROOT.kBlue)
                hbgx.SetFillColor(ROOT.kBlue)
                hbgx.Draw('hbar same')
                
            """
            ratio,band,tmp = DrawNiceRatioWithBand(htot, hdata, hxmin, hxmax, gyratioMin, gyratioMax, stuff, 10000 + iplot)
            ratios.append(band)
            ratios.append(ratio)
            ROOT.gPad.Update()
            ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        pad.cd()
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis()
        """

        ptag = ''
        can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
        can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')

        hdata.SetMaximum(gySFlog*hdata.GetMaximum())
        hdata.SetMinimum(gyLogMin*lumi)
        #ROOT.gPad.SetLogz(1)
        #pad.Update()
        #ROOT.gPad.RedrawAxis()
        #can.Print(pngdir + can.GetName() + ptag + '_logy' + '.png')
        #can.Print(pdfdir + can.GetName() + ptag + '_logy' + '.pdf')

        ##########################################
        #    and now subtract the background!    #
        ##########################################
        
        spad.cd()
        hdataBgSub = hdata.Clone(hdata.GetName() + '_bgSub')
        hdataBgSub.Add(toplot.hbg, -1.)
        hdataBgSub.SetMinimum(-gySFneg*hdataBgSub.GetMaximum())
        hdataBgSub.SetMaximum(gySFlin*hdataBgSub.GetMaximum())
        stuff.append([hsig, hdataBgSub])
        
        dopt = 'box' # stackItems[gStackName].dopt
        hdataBgSub.Draw(dopt)
        #hsig.Draw(stackItems[signalFileName].dopt + 'same') # why this ain't work?
        hsig.Draw(dopt + ' same')
        hdataBgSub.Draw(dopt + 'same')
        #sarrowsUp = DrawArrowForPointsOutsideYAxisRange(hdataBgSub, upHelpHisto, hxmin, hxmax)
        #stuff.append(sarrowsUp)
        ltex.Draw()
        gltex.Draw()
        chtex.Draw()

        # and now the ratio:
        if plotRatios:
            srpad.cd()
            sratio = DrawNice2DRatio(hsig, hdataBgSub, gyratioMin, gyratioMax, stuff, iplot, 'colz')
            ROOT.gPad.Update()

        sndf,schi2,sctex = ComputeChi2AndKS(hsig, hdataBgSub, 0.13, 0.73)
        if sndf > 0:
            sctex.Draw()
            Ltex.append(sctex)
            sctex.Draw()

        sleg = ROOT.TLegend(lx1, ly1 + (ly2-ly1)/2., lx2, ly2)
        sleg.SetBorderSize(0)
        sleg.AddEntry(hdataBgSub, stackItems[gStackName].legtag, stackItems[gStackName].lopt)
        # fragile, but again, relying on the fact that signal is expected the last sample in stack
        sleg.AddEntry(hsig, legitems[-1].legtag, stackItems[signalFileName].lopt)
        #sleg.Draw()
        
        scan.Print(pngdir + scan.GetName() + ptag + '_liny' + '.png')
        scan.Print(pdfdir + scan.GetName() + ptag + '_liny' + '.pdf')

    """
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

    os.system('mkdir -p {}/{}/'.format(pngdir,cantag))
    os.system('mkdir -p {}/{}/'.format(pdfdir,cantag))
    os.system('mv {}*.png {}/{}/'.format(pngdir,pngdir,cantag))
    os.system('mv {}*.pdf {}/{}/'.format(pdfdir,pdfdir,cantag))
    """

    print('DONE!')
    os.system('notify-send "DONE! {}"'.format(argv[0],))

    if batch != 'batch':
        ROOT.gApplication.Run()

    # kill oneself:
    os.system('killall XsectStack2D.py')

    return 

#########################################
#########################################
#########################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
#########################################
#########################################
#########################################



