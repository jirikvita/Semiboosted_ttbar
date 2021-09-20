#!/usr/bin/python

from __future__ import print_function

# jk 14.4.2020 (covid-19 homeoffice), based on code of StackPlots.py
# Example running: ./python/XsectStack.py ./python/list.txt "_mytag" notnorm batch

import ROOT

from xSectTools import *
from BumpSignifTools import *
from stackPlotItems import *
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
    lumi = 1000000 # pb! 
    csf = 1.

    doDivideByBinWidth = True
    printStats = False
    
    tag=''
    idata = 0

    gyratioMin = 0.65
    gyratioMax = 1.35

    gysignifMin = -3.
    gysignifMax = 9.

    gySFlin = 1.9
    gySFneg = 0.25
    gySFlog = 500.
    gyLogMin = 1.e-5

    # for 2D, obsolete here
    # ROOT.gStyle.SetPalette(1)
    # for precision of numbers in 2D when using the TEXT option:
    # ROOT.gStyle.SetPaintTextFormat("1.2f")

    pngdir='python/stack_png/'
    pdfdir='python/stack_pdf/' 
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
    
    cw = 900
    ch = 900
    lx1 = 0.49
    ly1 = 0.57
    lx2 = 0.88
    ly2 = 0.42+0.49

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadRightMargin(0.02)

    Chi2Hists = {}
    for key in ChiKeys:
        name = 'chi2_{}'.format(key)
        title = name + ';#chi^{2}'
        Chi2Hists[key] = ROOT.TH1D(name, title, 50, 0, 25)

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

    names1d = []
    names1d = names1d + MakeGlobalPlotItems()
    names1d = names1d + MakeRelPlotItems()
    names1d = names1d + MakeSpecialPlotItems()
    #names1d = names1d + MakeNPlotItems()
    names1d = names1d + MakeMassPlotItems()
    names1d = names1d + MakeMainKinemPlotItems()
    ###names1d = names1d + MakeCutPlotItems()
    #names1d = names1d + MakeFirstCutPlotItems()

    yieldsPrinted = {}
    yieldsTeXfile = {}
    for topo in ChiKeys:
        yieldsPrinted[topo] = False
        os.system('mkdir -p yields')
        yieldsTeXfile[topo] = open('yields/yields_{}_{}.tex'.format(cantag,topo), 'w')

    print('*** Processing histograms to stack;-) ***')
    for names in names1d:
        for name in names:
            # if name == names[0]:
            # HACK!!!
            #if not ('DiTopMass' in name or 'DiTopPt' in name):
            #    #if not ('DiTopMass' in name or 'DiTopPt' in name):
            #    continue
            print('Processing 1D histo {}'.format(name,))
            canname = cantag
            canname = canname + '_' + name + '_RS' # ratio and significance;) jk 15.9.2020
            canname = canname.replace('/', '_')
            can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
            cans.append(can)

            # prepare TPads:
            
            pads,pad_inset = MakeMultiSubPads(can, [0.70, 0.12, 0.18])
            Pads.append([can,pads,pad_inset])
            
            # prepare TPads, for bg-subtraction, i.e. signal-only
            scanname = canname + '_bgSub'
            scan = ROOT.TCanvas(scanname, scanname, 200, 200, cw, ch)
            cans.append(scan)
            spad1,spad2,spad_inset = MakePadsStack(scan, 'centre', 0.40, 0., 0., 0.)
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
                isSignal = rfile.GetName() == signalFileName
                ifile = rfiles.index(rfile)
                fname = rfile.GetName()
                print('    Processing file {}'.format(rfile.GetName()))
                hname = name
                print('       getting {}'.format(hname,))
                hist = rfile.Get(hname)
                if doDivideByBinWidth:
                    DivideByBinWidth(hist)
                hist.SetLineColor(stackItems[fname].lcol)
                hist.SetLineWidth(stackItems[fname].lw)
                hist.SetLineStyle(stackItems[fname].lst)
                hist.SetFillColor(stackItems[fname].fcol)
                hist.SetFillStyle(stackItems[fname].fst)

                # now get the other version of the histogram from the complementary sample:
                histStatIndep = rfileStatIndep.Get(hname)
                if doDivideByBinWidth:
                    DivideByBinWidth(histStatIndep)

                print('...Sample {} entries histo {} half0/half1 = {}/{} = {}'.format(fname, hname, hist.GetEntries(), histStatIndep.GetEntries(), hist.GetEntries() / histStatIndep.GetEntries() ) )
                    
                # note the additional SF here and for legend!
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

                # TO MOVE / adjust, add not to dra bool? 
                legitems.append(cLegItem(hist, stackItems[fname].legtag + fracttxt, stackItems[fname].lopt))
                # HACK: do not use this sample except 0B2S topology due to stat insuff. and fluctuations in 1B1S abd 2B0S!
                # jk 19.7.2021
                if 'pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200' in fname and not '0B2S' in hname:
                    continue

                if debug > 0:
                        print('* Scaling histos by weight {}'.format(sampleWeight))
                ScaleHistAndRebin(hist, hname, sampleWeight)
                ScaleHistAndRebin(histStatIndep, hname, sampleWeight)
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
    
    ratios = []
    signifratios = []
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
        pad = pads[1][0]  # upper pad
        rpad = pads[1][1] # ratio pad
        signifpad = pads[1][2] # significance pad
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

        #if 'DiTopM' in hdata.GetName():
        #    hxmin = 300
        #    hxmax = 2000
        #    hdata.GetXaxis().SetRangeUser(hxmin, hxmax)
        
        if doDivideByBinWidth:
            hdata.GetYaxis().SetTitle('Expected events / #Delta')
        else:
            hdata.GetYaxis().SetTitle('Expected events')
        # DiTop replace?
        FixXaxisTitle(hdata)
        hdata.Draw(stackItems[gStackName].dopt)
        stack.Draw('samehist')
        hdata.Draw(stackItems[gStackName].dopt + 'same')

        nb = hdata.GetNbinsX()
        stats = ''
        # option for integration:
        iopt = ''
        if doDivideByBinWidth:
            iopt = 'width'
        if printStats:
            stats =  ' N={:.0f} I={:.0f}'.format(hdata.GetEntries(),hdata.Integral(0, nb+1, iopt))
        leg.AddEntry(hdata, stackItems[gStackName].legtag + stats, stackItems[gStackName].lopt)

        topotag = ''
        for topo in ChiKeys:
            if topo in legh:
                topotag = topo + ''
        hhname = hdata.GetName()
        printYieldsNow = 'DiTopMass' in hhname and not 'Particle' in hhname and not 'Geo' in hhname
        if  printYieldsNow and not yieldsPrinted[topotag]:
            yieldsTeXfile[topotag].write('%%% {} %%%\n'.format(hhname))
            yieldsTeXfile[topotag].write(' \\begin{tabular}{l|rrr} \n')
            yieldsTeXfile[topotag].write(' \\hline \\hline \n')
            yieldsTeXfile[topotag].write(' {} topology & Events yield & Weighted yeild & Stat. unc.\\\\ \\hline \n'.format(topotag))

        sumEntries = 0.
        sumYields = 0.
        sumYieldsErr2 = 0.
        iopt = ''
        if doDivideByBinWidth:
            iopt = 'width'
        for ileg in range(len(legitems)-1, -1, -1):
            legitem = toplot.legitems[ileg]
            stats = ''
            if printStats:
                stats = ' N={:.0f} I={:.0f}'.format(legitem.hist.GetEntries(), legitem.hist.Integral(0, nb+1, iopt))
            # HACK JK VII/2021
            if not ('bbjj' in legitem.legtag and 'p_{T}^{j1,j2} #in (60,200)' in legitem.legtag and not '0B2S' in topotag):
                leg.AddEntry(legitem.hist, legitem.legtag + stats, legitem.lopt)
                if printYieldsNow and not yieldsPrinted[topotag]:
                    textag = legitem.legtag.replace('#','\\').replace('{T}','\\mathrm{T}').replace('GeV','\\,\\mathrm{GeV}').replace('TeV','\\,\\mathrm{TeV}') # .replace(',','\\,,\\quad')
                    # sample yield, weighted:
                    syield = legitem.hist.Integral(0, nb+1, iopt)
                    sentries = legitem.hist.GetEntries()
                    sweight = syield / sentries
                    print('    ...yields: adding sample {} entries: {} yield: {} weight: {}'.format(legitem.legtag, sentries, syield, sweight))
                    sentriesErr = sqrt(sentries)
                    srelErr = sentriesErr / sentries 
                    syieldErr = syield * srelErr
                    yieldsTeXfile[topotag].write(' ${}$ & ${:.0f} \\pm {:.0f}$ & ${:.0f} \\pm {:.0f}$ &  {:2.1f}\\% \\\\ \n'.format(textag,
                                                                                                                                    myround(sentries), myround(sentriesErr),
                                                                                                                                    myround(syield), myround(syieldErr),
                                                                                                                                    srelErr*100))
                    sumEntries = sumEntries + sentries
                    sumYields = sumYields + syield
                    sumYieldsErr2 = sumYieldsErr2 + syield*sweight
        if printYieldsNow and not yieldsPrinted[topotag]:
            yieldsTeXfile[topotag].write(' \\hline \n')
            dyield = hdata.Integral(0, nb+1, iopt)
            dentries = hdata.GetEntries()
            dentriesErr = sqrt(dentries)
            #drelErr = dentriesErr / dentries 
            #yieldsTeXfile[topotag].write(' pseudodata & ${:.0f} \\pm {:.0f}$ & ${:.0f} \\pm {:.0f}$ &  {:2.1f}\\%  \\\\ \n '.format(myround(entries),myround(entriesErr),myround(dyield),myround(intErr),relErr*100))
            yieldsTeXfile[topotag].write(' pseudodata & ${:.0f}$ & ${:.0f}$  &   \\\\ \n '.format(int(dentries),myround(dyield)))
            
            entriesErr = sqrt(sumEntries)
            yieldsErr = sqrt(sumYieldsErr2)
            relErr = yieldsErr / sumYields
            yieldsTeXfile[topotag].write(' prediction & ${:.0f} \\pm {:.0f}$ & ${:.0f} \\pm {:.0f}$  &  {:2.1f}\\% \\\\ \n '.format(int(sumEntries), myround(entriesErr),
                                                                                                                                    myround(sumYields), myround(yieldsErr),
                                                                                                                                    relErr*100))
            
            yieldsTeXfile[topotag].write(' pseudodata / prediction & {:1.2f} & {:1.2f} & \\\\ \n '.format(dentries / (1.*sumEntries),
                                                                                                          dyield / (1.*sumYields) ) )
            yieldsTeXfile[topotag].write(' \\hline \\hline \n')
            yieldsTeXfile[topotag].write(' \\end{tabular} \n')
            yieldsPrinted[topotag] = True
        if printYieldsNow:
            yieldsTeXfile[topotag].close()

        leg.Draw()

        ltex = ROOT.TLatex(0.13, 0.86, 'pp #sqrt{s} = 14 TeV  ' + 'L = {:.0f} ab'.format(lumi/1e6) + '{}^{-1}')
        ltex.SetTextSize(0.05)
        ltex.SetNDC()
        ltex.Draw()
        Ltex.append(ltex)

        gentxt = 'MadGraph5'
        if 'Detector' in hdata.GetName():
            gentxt += ' + Delphes'
        gltex = ROOT.TLatex(0.15, 0.9395, gentxt)
        gltex.SetTextSize(0.07)
        gltex.SetNDC()
        gltex.Draw()
        Ltex.append(gltex)
        Legs.append(leg)

        # channel name into plot
        chtex = ROOT.TLatex(0.13, 0.79, legh)
        chtex.SetTextSize(0.055)
        chtex.SetNDC()
        chtex.Draw()
        Ltex.append(chtex)

        # the total stack!
        last = stack.GetHists().At(stack.GetNhists()-1)
        htot = last.Clone(last.GetName() + '_tot')
        hbg = last.Clone(last.GetName() + '_bg')
        hbg.Reset()
        for ih in range(0, stack.GetNhists()-1):
            htot.Add(stack.GetHists().At(ih))
            hbg.Add(stack.GetHists().At(ih))

        ndf,chi2,ctex = ComputeChi2AndKS(hdata, htot, 0.13, 0.73)
        significance, stex = ComputeExcessSignificance1D(hdata, hbg)
        stex.Draw()
        if ndf > 0:
            ctex.Draw()
            Ltex.append(ctex)
            ctex.Draw()
            for key in ChiKeys:
                if key in fullhname:
                    print('filling chi2 hist with {}'.format(chi2/ndf) )
                    Chi2Hists[key].Fill(chi2 / ndf)

        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        rpad.cd()
        ratio,band,hratio = DrawNiceRatioWithBand(htot, hdata, hxmin, hxmax, gyratioMin, gyratioMax, stuff, 10000 + iplot, 'ratio')
        hratio.GetYaxis().SetTitleSize(0.235)
        hratio.GetYaxis().SetTitleOffset(0.17)
        ratios.append(band)
        ratios.append(ratio)
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        signifpad.cd()
        signifratio,hsratio = DrawSignificance(hbg, hdata, hxmin, hxmax, gysignifMin, gysignifMax, stuff, 100000 + iplot, 'signif.')
        hsratio.GetYaxis().SetTitleSize(0.16)
        hsratio.GetYaxis().SetTitleOffset(0.22)
        signifratios.append(signifratio)
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()
        
        
        pad.cd()
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis()

        ptag = ''
        can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
        can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')

        hdata.SetMaximum(gySFlog*hdata.GetMaximum())
        hdata.SetMinimum(gyLogMin*lumi)
        ROOT.gPad.SetLogy(1)
        pad.Update()
        ROOT.gPad.RedrawAxis()
        can.Print(pngdir + can.GetName() + ptag + '_logy' + '.png')
        can.Print(pdfdir + can.GetName() + ptag + '_logy' + '.pdf')

        ####################################
        # and now subtract the background!
        spad.cd()
        hdataBgSub = hdata.Clone(hdata.GetName() + '_bgSub')
        hsig = toplot.hsig
        hdataBgSub.Add(toplot.hbg, -1.)
        hdataBgSub.SetMinimum(-gySFneg*hdataBgSub.GetMaximum())
        hdataBgSub.SetMaximum(gySFlin*hdataBgSub.GetMaximum())

        stuff.append([hsig, hdataBgSub])
        
        # for out-of-range arrows:
        upHelpHisto = ROOT.TH2D(hdataBgSub.GetName() + '_upHelp', '',
                                hdataBgSub.GetXaxis().GetNbins(), hxmin, hxmax,
                                100, hdataBgSub.GetMinimum(), hdataBgSub.GetMaximum())
        upHelpHisto.SetStats(0)
        upHelpHisto.Draw()
        hdataBgSub.Draw(stackItems[gStackName].dopt + 'same')
        #hsig.Draw(stackItems[signalFileName].dopt + 'same') # why this ain't work?
        hsig.Draw('hist same')
        hdataBgSub.Draw(stackItems[gStackName].dopt + 'same')
        sarrowsUp = DrawArrowForPointsOutsideYAxisRange(hdataBgSub, upHelpHisto, hxmin, hxmax)
        stuff.append(sarrowsUp)
        ltex.Draw()
        gltex.Draw()
        chtex.Draw()

        # and now the ratio:
        srpad.cd()
        print('  ...making subtracted ratios...')
        sratio,sband,htmp = DrawNiceRatioWithBand(hsig, hdataBgSub, hxmin, hxmax, gyratioMin, gyratioMax, stuff, iplot, 'Pseudo (Data - Bg) / Sig')
        ratios.append(sband)
        ratios.append(sratio)
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()
        spad.cd()
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis()

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
        sleg.Draw()
        
        scan.Print(pngdir + scan.GetName() + ptag + '_liny' + '.png')
        scan.Print(pdfdir + scan.GetName() + ptag + '_liny' + '.pdf')

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

    print('DONE!')
    os.system('notify-send "DONE! {}"'.format(argv[0],))

    if batch != 'batch':
        ROOT.gApplication.Run()

        
    # kill oneself:
    os.system('killall XsectStackRatioSignif.py')

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



