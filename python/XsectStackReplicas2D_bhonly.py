#!/usr/bin/python

from __future__ import print_function

# jk 14.4.2020 (covid-19 homeoffice), based on code of StackPlots.py
# Example running: ./python/XsectStack.py ./python/list.txt "_mytag" notnorm batch
# 10.9.2020, 7.11.2020


import ROOT

from collections import OrderedDict

from xSectTools import *
from BumpSignifTools import *

from CorrItems import *

from cTopPlot import *

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
    #ROOT.gStyle.SetPadTopMargin(0.25)
    #ROOT.gStyle.SetPadLeftMargin(0.15)
    #ROOT.gStyle.SetPadRightMargin(0.05)
    #ROOT.gStyle.SetPadBottomMargin(0.15)
    
    pngdir='python/stack2d_png/'
    pdfdir='python/stack2d_pdf/' 
    os.system('mkdir -p {}'.format(pngdir,))
    os.system('mkdir -p {}'.format(pdfdir,))
    os.system('mkdir -p tex')
    
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
    
    if len(argv) > 2:
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
        print('Usage: {} config.txt variable'.format(argv[0]))
        return
    toSkipOrReqStr = ''
    if len(argv) > 2:
        toSkipOrReqStr = argv[2]
        
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

    Topos = [ '0B2S',
              '1B1S',
              '2B0S' ]
    v1Dx,v1Dy,v2D = '1Dx','1Dy','2D'
    BHversions = [v1Dx,v1Dy,v2D]
    vcols = {v1Dx : ROOT.kOrange + 7, v1Dy : ROOT.kMagenta, v2D : ROOT.kGreen + 1}
    #vcols = {v1Dx : ROOT.kYellow + 2, v1Dy : ROOT.kMagenta + 2, v2D : ROOT.kGreen + 2}
    #vcols = {v1Dx : ROOT.kGreen+2, v1Dy : ROOT.kBlue, v2D : ROOT.kRed}
    vfills = {v1Dx : 3345, v1Dy : 3354, v2D : 3305}
    
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

    nReplicas = 100 # !!! 100
    printAnyway = False

    SignifHists = OrderedDict()
    rhoVals = OrderedDict()
    rhoValsAver = OrderedDict()
    
    print('*** Processing 2D histograms to stack;-) ***')

    names = [ 'CosThetaStarVsDiTopPout',
              #'RttbarVsDiTopMass',
              #'DiTopMassVsDiTopPout',
              #'YboostVsChittbar',
    ]
    texTag = '_all'

    
    # HACK!!! for quick studies
    #for xname in names:
    ### DEFAULT!
    print('Working with toSkipOrReqStr: {}'.format(toSkipOrReqStr))
    for xname in hnamesDict:
        
        print('PROCESSING {}'.format(xname))

        """
        if toSkipOrReqStr != '':
            if not toSkipOrReqStr in xname:
                continue
            else:
                texTag = '_' + toSkipOrReqStr
        else:
            if 'DiTopMass' in xname or 'TopPt' in xname or 'Delta' in xname:
                continue
            else:
                texTag = '_noDiTopMass_noTopPt_noDelta'
        """
        if xname != toSkipOrReqStr:
            continue
        else:
            texTag = '_' + xname
        


        sname = 'BHscore_{}'.format(xname)
        sname = sname.replace('/','_')
        barehname = xname.split('/')[-1]
        title = sname + ';BH score;Replicas'
        # to change naming to max of BH score = log(t)!
        signifMin = 2
        signifMax = 8 #!!!8 # 7
        nsignifBins = 6*(signifMax - signifMin) #6
        # hack:
        #nsignifBins = 100*(signifMax - signifMin) #6

        for topo in Topos:
            try:
                print(SignifHists[topo])
            except:
                SignifHists[topo] = OrderedDict()
                rhoVals[topo] = OrderedDict()
                rhoValsAver[topo] = OrderedDict()

            # create a triplet of BH score histograms: 1Dx, 1Dy, 2D and correlations for averaging
            rhoVals[topo][barehname] = []
            SignifHists[topo][barehname] = OrderedDict()
            # TODO: rename the name and title to respective x and y variables!
            print('*** Creating histos for {} {}'.format(topo,barehname))
            SignifHists[topo][barehname][v1Dx] = ROOT.TH1D(sname + '_' + topo + '_' + v1Dx, title, nsignifBins, signifMin, signifMax)
            SignifHists[topo][barehname][v1Dy] = ROOT.TH1D(sname + '_' + topo + '_' + v1Dy, title, nsignifBins, signifMin, signifMax)
            SignifHists[topo][barehname][v2D]  = ROOT.TH1D(sname + '_' + topo + '_' + v2D,  title, nsignifBins, signifMin, signifMax)
            for vkey in SignifHists[topo][barehname]:
                hh = SignifHists[topo][barehname][vkey]
                hh.SetStats(0)
                #hh.SetLineWidth(1)
                #if 'VsDiTopMass' in hh.GetName() and vkey == v1Dx:
                hh.SetLineWidth(3)
                hh.SetLineColor(vcols[vkey])
                hh.SetFillColor(vcols[vkey])
                hh.SetFillStyle(vfills[vkey])
                
                
            print('=== current SignifHists ===')
            print(SignifHists)
            
        #for irep in range(35, 45):
        for irep in range(0, nReplicas):

            for topo in Topos:
                name = topo + '/replicas/' + ('Detector' + xname).replace('Vs', 'VsDetector') + '_rep{}'.format(irep)
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
                    if doDivideByBinArea:
                        DivideByBinArea(hist)
                    hist.SetLineColor(stackItems[fname].lcol)
                    hist.SetLineWidth(stackItems[fname].lw)
                    hist.SetLineStyle(stackItems[fname].lst)
                    hist.SetFillColor(stackItems[fname].fcol)
                    hist.SetFillStyle(stackItems[fname].fst)

                    print('    ...got histo {} of integral {} from file {}'.format(hist.GetName(), hist.Integral(), rfile.GetName()))
                    
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
                
                #print('Pads:')
                #print(Pads[-1])
                #print('Histos in stack: {}, histo0 entries: {}'.format(stack.GetNhists(), stack.GetHists().At(0).GetEntries()))
                ToPlot.append( cTopPlot(stack, stackStatIndep, leg, legh, legitems, Pads[-1], sPads, hname, barehname) )
                print('...checking hdata integral: {}'.format(ToPlot[-1].hdata.Integral()))
                Hists.append(hists)
                Stacks.append(stack)

    ### now draw the stack and data! ;-) ###
    ### TODO: subtract BG!

    print('=== Significance histos to be filled ===')
    print(SignifHists)
    
    ratios = []
    iplot = -1
    print('*** Plotting stacked histograms;) ***')
    print('Plot items to go through: {}'.format(len(ToPlot)))
    #print(ToPlot)
    for toplot in ToPlot:
        print(' --- processing {} {} ---'.format(toplot.basename, toplot.fullhname))
        iplot = iplot+1
        stack = toplot.stack
        hdata = toplot.hdata
        print('hdata {} integral: {}'.format(toplot.fullhname, hdata.Integral()))
        if hdata.Integral() < 1.e-6:
            print('ERROR, got data integral of ZERO, skipping this case...')
            continue
        #hdata.SetMarkerStyle(stackItems[gStackName].mst)
        #hdata.SetMarkerSize(stackItems[gStackName].msz)
        #hdata.SetMarkerColor(stackItems[gStackName].mcol)
        hdata.SetLineColor(stackItems[gStackName].lcol)
        leg = toplot.leg
        legh = toplot.legh

        pads = toplot.pads
        can = pads[0]  # canvas
        pad = pads[1]  # upper pad
        rpad = pads[2] # ratio pad
        #print(pads)
        
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
        print('Last histo integral: {}'.format(last.Integral()))
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
        htot.Draw(dopt)
        hsig.Draw(dopt + ' same')

        corr_bg = hbg.GetCorrelationFactor()
        ctex_bg = ROOT.TLatex(0.65, 0.08,'#rho=' + '{:1.2f}'.format(corr_bg))
        ctex_bg.SetNDC()
        ctex_bg.SetTextSize(0.045)
        ctex_bg.SetTextColor(htot.GetFillColor())
        
        corr_sig = hsig.GetCorrelationFactor()
        ctex_sig = ROOT.TLatex(0.81, 0.08,'#rho=' + '{:1.2f}'.format(corr_sig))
        ctex_sig.SetNDC()
        ctex_sig.SetTextSize(0.045)
        ctex_sig.SetTextColor(hsig.GetFillColor())

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
        ytxt = 0.08 # 0.02
        bhcol = vcols[v2D]
        bhhcol = vcols[v2D]
        print('Computing 2D BH test from data,bg of integrals {},{}'.format(hdata.Integral(), hbg.Integral()))
        bhresult,bhtex = FindBestBumpHunter2DArea(hdata, hbg, 0.02, ytxt, 'BH log(t) xy: ')
        used2dBH = bhresult.bins
        usedLinesBH = DrawUsedBinsGrid(hdata, used2dBH, bhresult.i, bhresult.j, bhcol, bhhcol, 1)
        
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
        bhresultx,bhtexx = FindBestBumpHunter1DInterval(hdatax, hbgx, 0.345, ytxt, 'x: ')
        bhresulty,bhtexy = FindBestBumpHunter1DInterval(hdatay, hbgy, 0.500, ytxt, 'y: ')
        usedLinesBHx = DrawUsedBinsLines(hdatax, hdatay, bhresultx.bins, bhresultx.i, True,  0.08, 0.005, bhcol, bhhcol, 1, 3, 3)
        usedLinesBHy = DrawUsedBinsLines(hdatay, hdatax, bhresulty.bins, bhresulty.i, False, 0.10, 0.005, bhcol, bhhcol, 1, 3, 3)
        #usedLinesx   = DrawUsedBinsLines(hdatax, hdatay, usedx, i0x,  True,  0.10, 0.005, scol, shcol, 1)
        #usedLinesy   = DrawUsedBinsLines(hdatay, hdatax, usedy, i0y,  False, 0.12, 0.005, scol, shcol, 1)
        
        #stex1Dx.Draw()
        #stex1Dy.Draw()
        #if significance > significance1Dy and significance > significance1Dx:
        #    stex.SetTextColor(hsig.GetFillColor())
        #stex.Draw()

        barehname = toplot.basename
        topo = toplot.fullhname.split('/')[0]
        print('   ...filling {} {}'.format(topo,barehname))
        SignifHists[topo][barehname]['1Dx'].Fill(ComputeSafeLog(bhresultx.t))
        SignifHists[topo][barehname]['1Dy'].Fill(ComputeSafeLog(bhresulty.t))
        SignifHists[topo][barehname]['2D'].Fill(ComputeSafeLog(bhresult.t))
        # data correlation factor
        rhoVals[topo][barehname].append(hdata.GetCorrelationFactor())

        if bhresult.t > bhresultx.t and bhresult.t > bhresulty.t:
            bhtex.SetTextColor(hsig.GetFillColor())
        bhtex.Draw()
        bhtexx.Draw()
        bhtexy.Draw()
        stuff.append([bhtex, bhtexx, bhtexy])
        
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

        if printAnyway:
            ptag = ''
            can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
            can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')
            #hdata.SetMaximum(gySFlog*hdata.GetMaximum())
            #hdata.SetMinimum(gyLogMin*lumi)


    # now plot the distributions of the BH log(t) over spectra and topologies:)
    cans = {}
    legs = {}
    topots = {}
    baseopt = 'hist'
    sameopt = 'same'
    
    print('Plotting the BHscore distributions')

    #sigfiletag = sigTag.replace(' ','').replace(',','_').replace('#','').replace('{','').replace('}','').replace('=','').replace('m_','').replace('y_0','y0').replace("Z'","zp")
    
    for topo in Topos:
        texFileName = 'tex/table_{}_BH_{}{}.tex'.format(cantag,topo,texTag) # sigfiletag,
        texTableFile = open(texFileName, 'w')

        for barehname in SignifHists[topo]:
            VsIndex = barehname.find('Vs')
            vtags = OrderedDict()
            vtags[v1Dx] = ivarsLabelsDict[barehname[VsIndex+2:]]
            vtags[v1Dy] = ivarsLabelsDict[barehname[:VsIndex]]
            vtags[v2D]  = v2D # + ' ' + barehname
            vtagsTeX = OrderedDict()
            keyx = barehname[VsIndex+2:]
            keyy = barehname[:VsIndex]
            print('keyx,keyy: {},{}. Labels {},{}'.format(keyx, keyy, ivarsLabelsDictTeX[keyx], ivarsLabelsDictTeX[keyy]))
            vtagsTeX[v1Dx] = '{}'.format(ivarsLabelsDictTeX[keyx])
            vtagsTeX[v1Dy] = '{}'.format(ivarsLabelsDictTeX[keyy])
            vtagsTeX[v2D]  = '{} vs. {}'.format(vtagsTeX[v1Dy], vtagsTeX[v1Dx])
            print('...processing {} {}'.format(topo, barehname))
            opt = baseopt + sameopt
            ckey = barehname + '_' + topo
            try:
                cans[ckey].cd()
            except:
                canname = 'BHcmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cans[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                cans[ckey].cd()
                #legs[ckey] = ROOT.TLegend(0.60, 0.60, 0.89, 0.89)
                legs[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legs[ckey].SetNColumns(3)
                #legs[ckey].SetHeader(topo, 'C')
                legs[ckey].SetBorderSize(0)
                rho = sum(rhoVals[topo][barehname]) / len(rhoVals[topo][barehname])
                rhoValsAver[topo][barehname] = 1.*rho
                topot = ROOT.TLatex(0.30, 0.92, '{}'.format(topo) + '         #bar{#rho}=' + '{:+1.2f}'.format(rho))
                topot.SetNDC()
                topot.SetTextSize(0.06)
                topots[ckey] = topot
                opt = baseopt
            maxScore = -1000
            for version in BHversions:
                hh = SignifHists[topo][barehname][version]
                score = hh.GetMean()
                if score > maxScore:
                    maxScore = score
            for version in BHversions:
                hh = SignifHists[topo][barehname][version]
                #hh.SetMaximum(2*hh.GetMaximum())
                hh.SetMaximum(0.60*nReplicas)
                hh.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hh.GetName(), hh.GetEntries(), hh.Integral(0, hh.GetXaxis().GetNbins()+1)))
                opt = baseopt + sameopt
                score = hh.GetMean()
                legs[ckey].AddEntry(hh, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], score, hh.GetRMS()), 'F')
                if version == BHversions[0]:
                    texTableFile.write( ' {} & '.format(topo))
                font = ''
                if abs(score - maxScore) < 1.e-5:
                    font = '\\mathbf'
                texTableFile.write( ' {} & $'.format(vtagsTeX[version]) + font + '{' + '{:1.2f}'.format(score) + '} \\pm ' + '{:1.2f}$ '.format(hh.GetRMS()))
                if version == BHversions[-1]:
                    texTableFile.write( r' \\')
                else:
                    texTableFile.write( r' &')
            texTableFile.write(' % rho % ' + ' {:1.3f}  '.format(rhoValsAver[topo][barehname]))
            texTableFile.write( '\n')
        # hline after each topology:
        texTableFile.write( r' \hline' + '\n')
        texTableFile.close()
        print('Written {}'.format(texFileName))
    for ckey in cans:
        ptag = ''
        cans[ckey].cd()
        legs[ckey].Draw()
        topots[ckey].Draw()
        #cans[ckey].Print(pngdir + 'model_' + sigfiletag + '_' + cans[ckey].GetName() + ptag + '_liny' + '.png')
        #cans[ckey].Print(pdfdir + 'model_' + sigfiletag + '_' + cans[ckey].GetName() + ptag + '_liny' + '.pdf')
        cans[ckey].Print(pngdir + cans[ckey].GetName() + ptag + '_liny' + '.png')
        cans[ckey].Print(pdfdir + cans[ckey].GetName() + ptag + '_liny' + '.pdf')
            
    print('DONE!')
    os.system('notify-send "DONE! {}"'.format(argv[0],))

    if batch != 'batch':
        ROOT.gApplication.Run()

    # kill oneself:
    # os.system('killall XsectStackReplicas2D.py')

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



