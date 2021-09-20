#!/usr/bin/python

from __future__ import print_function

# jk 14.4.2020 (covid-19 homeoffice), based on code of StackPlots.py
# Example running: ./python/XsectStack.py ./python/list.txt "_mytag" notnorm batch
# rev. 3.9.2020

import ROOT

from xSectTools import *
from cTopPlot import *
from BumpSignifTools import *

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

    gyratioMin = 0.6
    gyratioMax = 1.4
    gySFlin = 1.9
    gySFneg = 0.25
    gySFlog = 500.
    gyLogMin = 1.e-5

    # for 2D, obsolete here
    # ROOT.gStyle.SetPalette(1)
    # for precision of numbers in 2D when using the TEXT option:
    # ROOT.gStyle.SetPaintTextFormat("1.2f")

    pngdir='python/signif_png/'
    pdfdir='python/signif_pdf/' 
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

    cw = 800
    ch = 800
    lx1 = 0.49
    ly1 = 0.57
    lx2 = 0.88
    ly2 = 0.42+0.49

    ROOT.gStyle.SetOptTitle(0)


    Topos = [ '2B0S',
              '1B1S',
              '0B2S',
    #'AnySel'
    ]


    
    Chi2Hists = {}
    for key in Topos:
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

    

    # TO FINISH!
    suff = '_denser'
    names1d = []
    for topo in Topos:
        names1d.append(topo + '/DetectorDiTopMass' + suff)
        names1d.append(topo + '/DetectorDiTopPt' + suff)
        names1d.append(topo + '/DetectorDiTopCosThetaStar' + suff)
        names1d.append(topo + '/DetectorDiTopPout' + suff)
        names1d.append(topo + '/DetectorTopPt' + suff)
            #topo + '/ParticleDiTopMass' + suff,


    nReplicas = 100 # !!!
    printAnyway = False
    
    print('*** Processing histograms to stack;-) ***')
    SignifHists = {}

    for names in names1d:
        sname = 'signif_{}'.format(names)
        sname = sname.replace('/','_')
        barehname = names.split('/')[-1]
        title = sname + ';Signal significance;Replicas'
        signifMax = 2048+1
        signifMin = 1
        nsignifBins = 256
        #if not 'xd' in flist:
        #    signifMax = 150
            
        for key in Topos:
            if key in names:
                try:
                    print(SignifHists[key])
                except:
                    SignifHists[key] = {}
                SignifHists[key][barehname] = ROOT.TH1D(sname, title, nsignifBins, signifMin, signifMax)
            
        for irep in range(0, nReplicas):
            name = names
            if 'denser' in name:
                name = name.replace('_denser',  '_rep{}_denser'.format(irep))
            else:
                name = names + '_rep{}'.format(irep)
            name = name.replace('/','/replicas/')
            
            print('Processing 1D histo {}'.format(name,))
            canname = cantag
            canname = canname + '_' + name
            canname = canname.replace('/', '_')
            can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
            cans.append(can)

            # prepare TPads:
            pad1,pad2,pad_inset = MakePadsStack(can, 'centre', 0.40, 0., 0., 0.)
            Pads.append([can,pad1,pad2,pad_inset])
            
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
                legitems.append(cLegItem(hist, stackItems[fname].legtag + fracttxt, stackItems[fname].lopt))

                # HACK: do not use this sample except 0B2S topology due to stat insuff. and fluctuations in 1B1S abd 2B0S!
                # jk 19.7.2021, here from 8.7.2021
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
            ToPlot.append( cTopPlot(stack, stackStatIndep, leg, legh, legitems, Pads[-1], sPads[-1], hname, barehname) )
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
        basehname = toplot.basename

        pad.cd()
        hdata.SetMaximum(gySFlin*hdata.GetMaximum())
        #hdata.SetMinimum(0.1)
        
        hxmin = hdata.GetXaxis().GetXmin()
        hxmax = hdata.GetXaxis().GetXmax()

        #if 'DiTopM' in hdata.GetName():
        #    hxmin = 300
        #    hxmax = 2000
        #    hdata.GetXaxis().SetRangeUser(hxmin, hxmax)
        
        if DivideByBinWidth:
            hdata.GetYaxis().SetTitle('Expected events / #Delta')
        else:
            hdata.GetYaxis().SetTitle('Expected events')
        hdata.Draw(stackItems[gStackName].dopt)
        stack.Draw('samehist')
        hdata.Draw(stackItems[gStackName].dopt + 'same')

        nb = hdata.GetNbinsX()
        stats = ''
        if printStats:
            stats =  ' N={:.0f} I={:.0f}'.format(hdata.GetEntries(),hdata.Integral(0, nb+1))
        leg.AddEntry(hdata, stackItems[gStackName].legtag + stats, stackItems[gStackName].lopt)
        
        for ileg in range(len(legitems)-1, -1, -1):
            legitem = toplot.legitems[ileg]
            stats = ''
            if printStats:
                stats = ' N={:.0f} I={:.0f}'.format(legitem.hist.GetEntries(), legitem.hist.Integral(0, nb+1))
            leg.AddEntry(legitem.hist, legitem.legtag + stats, legitem.lopt)
        leg.Draw()

        ltex = ROOT.TLatex(0.13, 0.86, 'pp #sqrt{s} = 14 TeV  ' + 'L = {:.0f} ab'.format(lumi/1e6) + '{}^{-1}')
        ltex.SetTextSize(0.055)
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
            for topo in Topos:
                if topo in fullhname:
                    print('filling chi2 hist with {}'.format(chi2/ndf) )
                    Chi2Hists[topo].Fill(chi2 / ndf)
                    SignifHists[topo][basehname].Fill(significance)
                    

        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        rpad.cd()
        ratio,band = DrawNiceRatioWithBand(htot, hdata, hxmin, hxmax, gyratioMin, gyratioMax, stuff, 10000 + iplot)
        ratios.append(band)
        ratios.append(ratio)
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()        #ROOT.gPad.GetFrame().Draw()

        pad.cd()
        ROOT.gPad.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis()

        ptag = ''
        if printAnyway:
            can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
            can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')

        hdata.SetMaximum(gySFlog*hdata.GetMaximum())
        hdata.SetMinimum(gyLogMin*lumi)
        ROOT.gPad.SetLogy(1)
        pad.Update()
        ROOT.gPad.RedrawAxis()
        if printAnyway:
            can.Print(pngdir + can.GetName() + ptag + '_logy' + '.png')
            can.Print(pdfdir + can.GetName() + ptag + '_logy' + '.pdf')

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
        sratio,sband = DrawNiceRatioWithBand(hsig, hdataBgSub, hxmin, hxmax, gyratioMin, gyratioMax, stuff, iplot, 'Pseudo (Data - Bg) / Sig')
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
        if printAnyway:
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
    for topo in Chi2Hists:
        chi2h = Chi2Hists[topo]
        val = chi2h.GetMaximum()
        if val > ymax:
            ymax = val
    for topo in Chi2Hists:
        chi2h = Chi2Hists[topo]
        chi2h.SetLineColor(chicols[i])
        chi2h.SetLineWidth(2)
        chi2h.SetStats(0)
        chi2h.SetMaximum(ymax*1.1)
        chi2h.Draw(opt)
        chileg.AddEntry(chi2h, '{:10} #mu={:2.2f}'.format(topo, chi2h.GetMean()) , 'L')
        opt = 'hist same'
        i += 1
    chileg.Draw()
    can.Print(pngdir + can.GetName() + '.png')
    can.Print(pdfdir + can.GetName() + '.pdf')



    hkeys = [] # histogram names as keys, regardless the topology
    for topo in SignifHists:
        for hkey in SignifHists[topo]:
            if not hkey in hkeys:
                hkeys.append(hkey)
    print(hkeys)
    for hkey in hkeys:

        canname = 'SignifHists_{}'.format(hkey)
        canname = canname.replace('/','_')
        can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
        cans.append(can)
        opt = 'hist'
        signifleg = ROOT.TLegend(0.45, 0.55, 0.88, 0.88)
        signifleg.SetBorderSize(0)
        legh = hkey + ''
        legh = legh.replace('Detector','Detector ').replace('Particle',' Particle ').replace('_denser','')
        signifleg.SetHeader(legh)
        stuff.append(signifleg)
        signifcols = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+3, ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+3]
        ymax  = -1
        for key in SignifHists:
            for hkey2 in SignifHists[key]:
                if hkey != hkey2:
                    continue
                signifh = SignifHists[key][hkey2]
                val = signifh.GetMaximum()
                if val > ymax:
                    ymax = val
        i = 0
        for key in SignifHists:
            for hkey2 in SignifHists[key]:
                if hkey != hkey2:
                    continue
                signifh = SignifHists[key][hkey]
                signifh.SetLineColor(chicols[i])
                signifh.SetLineWidth(2)
                signifh.SetStats(0)
                signifh.SetMaximum(ymax*1.9)
                signifh.Draw(opt)
                ROOT.gPad.SetLogx(1)
                ltag = key + ''
                #ltag = ltag.replace('/',' ')
                signifleg.AddEntry(signifh, '{:10} #mu={:2.2f}'.format(ltag, signifh.GetMean()) , 'L')
                opt = 'hist same'
                i += 1
        signifleg.Draw()
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
    os.system('killall signifXsectStackReplicas1D.py')

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



