#!/usr/bin/python

# jk 3.3.2017, 30.3.2017, 7.4.2017, 6.9.2017, 7--16.1.2019, May 2019, May 2020
# Example running: ./python/CmpPlotHistos_GenDet.py ./python/list_cmp.txt "_allcmp" norm batch
# compares side-by-side particle (generator) and detector (delphes) level histograms in root files as provided in list

# TODO: option to divide by the bin area in case of TH2!;-)
# TODO: convert TProfileY to a transposed TGraphErrors, and similar for TProjectionY
# same-plot the profileY graph into the 2D canvas
# TRY also plotting projections into vertical and horizontal side canvases?


#############################################################

from __future__ import print_function

import ROOT
import sys, os, getopt
# or? import argparse ?;-)

from Tools import *
from math import *

from array import array

from ROOT import std,Double,AddressOf


#############################################################
#############################################################
#############################################################

    
def main(argv):

    from PlotItems import *

    stuff = []
    
    #############################################################

    gDrawOpt = 'hist ' # 'C hist'
    gAxisOpt = ''

    addRatioTitle = False
    # draw diagonal in migration matrices:
    drawDiag = False

    #PlotAlsoProfile = False
    PlotAlsoProfile = True

    #PlotAlsoProjection = False
    PlotAlsoProjection = True

    divideByBinWidth = True # TODO!!! so far effective only on projections!
    intInLeg = False
    resoInLeg = False
    #############################################################

    bw=0
    bwtag = ''

    doratio='ratio'
    ratioIndex = 0 # file index to do ratios w.r.t.
    ratioTag = 'ratio to t#bar{t}   '

    # some more nasty global variables;)
    gcol = [
        #standard:
        ROOT.kRed,
        ROOT.kBlue,
        #ROOT.kBlue,
        #ROOT.kBlue+2, ROOT.kBlue+2,
        #ROOT.kGreen+2,
        ROOT.kGreen+2,
        #ROOT.kMagenta+2,
        ROOT.kMagenta+2,
        #ROOT.kCyan+1,
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

    LevelLabel = ['particle', 'detector']

    gwidth = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    #gls = [1, 1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,]
    #gls = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,]
    # DEFAULT:
    gls = [1, 2, 2, 3, 4, 5, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,]
    #gls = [1, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1]
    gmsz = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    gmst = [20, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    gyratioMin = 0. # 0.55
    gyratioMax = 2. # 1.45
    gySFlin = 1.9
    gySFlog = 100.




    #############################################################
    # command-line controllable:
    gBatch        = True
    gNormalize    = False
    gNormalizeMigra = False
    gStopAfter1D  = False
    gPngTag = ''
    
    ROOT.gStyle.SetPalette(1)
    # for precision of numbers in 2D when using the TEXT option:
    ROOT.gStyle.SetPaintTextFormat("1.2f")
    
    if bw == 1:
        gcol = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        bwtag = '_bw'
        ROOT.gStyle.SetPalette(52)


    pngdir='python/png/'
    pdfdir='python/pdf/' 
    os.system('mkdir -p {:}'.format(pngdir,))
    os.system('mkdir -p {:}'.format(pdfdir,))


    # users, add your specific settings if needed:
    if os.getenv('USER') == 'qitek':
        # JK specific settings
        pass

    # get command line arguments:
    cantag='cmp_'

    print(argv)
    if len(argv) > 1:
        flist = argv[1]
        cantag = flist
        cantag = cantag.replace('python/', '')
        cantag = cantag.replace('lists/', '')
        cantag = cantag.replace('/', '_')
        cantag = cantag.replace('.txt', '')
        cantag = cantag.replace('list_', '')
        print('OK, using user-defined histograms root file list {:}'.format(flist,) )
        # read and parse the list:
        fnames = []
        listfile = open(flist, 'r')
        print('  ...reading lines' )
        for line in listfile.readlines():
            ffname = line[:-1]
            if len(ffname) > 5 and ffname[0] != '#' :
                fnames.append(ffname) # remove end line
        print('  Read:' )
        for fname in fnames:
            print('   {:}'.format(fname,) )
    else:
        print('ERROR: you MUST provide a filename containing a list of ROOT files to process as the first argument!')
        sys.exit(1)
            
            
    
    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    print(argv[2:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbnmst:', ['help','batch','norm', 'mnorm', 'stop', 'tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -n --norm -s --stop -tSomeTag -m --mnorm --tag="MyCoolPngs"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -b --batch -n --norm -s --stop -tSomeTag -m --mnorm --tag="MyCoolPngs"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-n", "--norm"):
            gNormalize = True
        elif opt in ("-m", "--mnorm"):
            gNormalizeMigra = True
        elif opt in ("-s", "--stop"):
            gStopAfter1D = True
        elif opt in ("-t", "--tag"):
            gPngTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gPngTag,) )

    print('Done.')

    if gNormalizeMigra:
        print('Asked for normalizing migrations, NOT plotting profile nor projection!')
        PlotAlsoProjection = False
        PlotAlsoProfile = False

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    print('*** Settings:')
    print('tag={:}, normalize={:}, batch={:}, StopAfter1D={:}'.format(gPngTag, gNormalize, gBatch,gStopAfter1D))
    print('files:')
    print(fnames)
    print('')

    if len(fnames) < 1:
        print('ERROR: Got no files in config, nothing to draw, exiting!')
        exit(1)
    
    if len(fnames) < 2:
        print('Got only 1 sample, assuming user wants only the 2D migrations, removing 1D plotting!')
        # names1d = []
        names1d = names1d_Glob[:1]
        
    
    cans = []
    Hists = []
    Objs = []
    rfiles = []
    Legs = []
    tags = {}
    Pads = []
    txts = []
    
    cw = 1200
    ch = 750

    ROOT.gStyle.SetOptTitle(0)

    for fname in fnames:
        rfile = ROOT.TFile(fname, 'read')
        rfiles.append(rfile)
        gPngTag = fname.replace('analyzed_histos_', '').replace('.root','')
        tags[fname] = gPngTag

    for name in names1d:
        iname = names1d.index(name)
        DualHists = []
        opt = {}
        print('Processing 1D histos {:}'.format(name,))
        nlevels = len(name) # is usually 2: generator and detector levels;-)
        canname = cantag
        for i in range(0, nlevels):
            canname = canname + '_' + name[i]
        canname = canname.replace('/', '_').replace('._','')
        can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
        cans.append(can)
        DualPads = []
        side='left'

        # prepare TPads:
        for i in range(0, nlevels):
            pad1,pad2,pad_inset = MakePads(can, side, 0.40, 0., 0., 0.)
            DualPads.append([pad1,pad2,pad_inset])
            side='right'
        Pads.append(DualPads)

        # old: can.Divide(nlevels,1)
        legs = []

        # find also global maximum over particle and detector levels!
        maxy = -999
        for i in range(0, nlevels):
            hists = []
            opt[name[i]] = ''
            print('  Processing histo {:}'.format(name[i],))
            leg = ROOT.TLegend(0.12, 0.15+0.44+0.05, 0.88, 0.40+0.44+0.05)
            topotag = name[i].split('/')[0]
            leg.SetHeader('Selection: {}'.format(topotag))
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            legs.append(leg)
            for rfile in rfiles:
                ifile = rfiles.index(rfile)
                print('    Processing file {:}'.format(rfile.GetName()))
                hname = name[i]
                
                # HACK!!!
                #if not 'DiTopMass' in hname:
                #    continue
                
                print('       getting {:}'.format(hname,))
                hist = rfile.Get(hname)
                hist.SetStats(0)
                hist.SetLineColor(gcol[ifile])
                #hist.SetMarkerSize(gmsz[ifile])
                #hist.SetMarkerColor(gcol[ifile])
                #hist.SetMarkerStyle(gmst[ifile])
                hist.SetLineWidth(gwidth[ifile])
                hist.SetLineStyle(gls[ifile])
                inttag = ""
                if intInLeg:
                    inttag = ' I={:.0f}'.format(hist.Integral())
                elif resoInLeg:
                    inttag = ' #sigma/#mu={:1.3f}'.format(hist.GetRMS()/hist.GetMean())
                leg.AddEntry(hist, MakeNiceLegendEntry(tags[rfile.GetName()], i) + inttag, 'L')
                print('LEG ENTRY {}'.format( MakeNiceLegendEntry(tags[rfile.GetName()], i) + inttag))
                ### REBINNING HERE!!!
                if not hname.find('JetN') >= 0 and IsUniformlyBinned(hist):
                    #print('OK, DOING rebinning {:}'.format(hname))
                    if 'TopM' in hname or 'WM' in hname:
                        hist.Rebin(4)
                    else:
                        hist.Rebin(2)
                else:
                    pass #print('OK, NOT rebinning {:}'.format(hname))

                ylabel = 'events'
                if gNormalize:
                    norm = Double(hist.Integral())
                    if norm > 0.:
                        hist.Scale(1./norm)
                    else:
                        print('ERROR! Histogram {:} from file {:} has ZERO integral!'.format(hist.GetName(), rfile.GetName(), ))
                if divideByBinWidth:
                    if gNormalize:
                        ylabel = 'fraction of ' + ylabel
                    DivideByBinWidth(hist)
                    hist.GetYaxis().SetTitle(ylabel + ' / #Delta')
                else:
                    hist.GetYaxis().SetTitle(ylabel)
                hist.GetYaxis().SetTitleSize(0.05)
                hist.GetYaxis().SetTitleOffset(1.)
                if hist.GetMaximum() > maxy:
                    maxy = hist.GetMaximum()
                
                # old:
                #can.cd(i+1)
                #print(Pads)
                Pads[-1][i][0].cd()
                #ROOT.gPad.SetLogy()
                ROOT.gPad.SetTicks(1,1)
                ###if opt[name[i]].find('same') < 0:
                ###    hist.SetMaximum(gySFlin*hist.GetMaximum())
                FixXaxisTitle(hist)
                hist.Draw(gDrawOpt + gAxisOpt + opt[name[i]])
                #hist.GetYaxis().SetMoreLogLabels()
                opt[name[i]] = ' same'
                hists.append(hist)
            # end of loop over files
            DualHists.append(hists)
        # end of loop over levels
        Hists.append(DualHists)
        for dh in DualHists:
            for h in dh:
                h.SetMaximum(maxy*gySFlin)
        for i in range(0, nlevels):
            # old: can.cd(i+1)
            Pads[-1][i][0].cd()
                
            ROOT.gPad.SetTicks(1,1)
            legs[i].Draw()
            delphestag = ''
            if 'detector' in LevelLabel[i]:
                delphestag = 'Delphes '
            ltex = ROOT.TLatex(0.10, 0.95, '{}{:} level'.format(delphestag, LevelLabel[i]))
            ltex.SetTextSize(0.06)
            ltex.SetNDC()
            ltex.Draw()
            Ltex.append(ltex)

            etex = ROOT.TLatex(0.6, 0.95, 'pp #sqrt{s} = 14 TeV')
            etex.SetTextSize(0.055)
            etex.SetNDC()
            etex.Draw()
            Ltex.append(etex)
            
        Legs.append(legs)


    Ratios = []
    for name in names1d:
        DualRatios = []
        iname = names1d.index(name)
        ###print('* Hists[iname]: ', Hists[iname])
        #print(len(Hists[iname]))
        for i in range(0, nlevels):
            ###print('** Hists[iname][i]: ', Hists[iname][i])
            #print(len(Hists[iname][i]))
            ratios = []
            if doratio == 'ratio' and ratioIndex < len(Hists[iname][i]) and ratioIndex >= 0:
                #print('OK, doing the ratios...')
                for ifile in range(0,len(Hists[iname][i])):
                    ### LET'S DO also the ratio to the denominatror, to have the error band;-)
                    if ratioIndex != ifile:
                        ratio = Hists[iname][i][ifile].Clone(Hists[iname][i][ifile].GetName() + '_ratio_{:}_{:}_{:}'.format(iname, i, ifile) )
                        ratio.Divide(Hists[iname][i][ratioIndex])
                        ratio.Scale(1.)
                        ratios.append(ratio)
                    else:
                        ratio = MakeOneWithErrors(Hists[iname][i][ifile])
                        ratios.append(ratio)
                else:
                    pass #print('SKIPPING the ratio!')
                DualRatios.append(ratios)
        Ratios.append(DualRatios)


    # now plot the ratios:
    print('*** Ratios to plot')
    #print(Ratios)
    #print(len(Ratios))
    tmpScale = []
    lines = []
    k = -1
    for iname in range(0,len(Ratios)):
        k = k+1
        opt = {}
        #print(Ratios[iname])
        #print(len(Ratios[iname]))
        for i in range(0, len(Ratios[iname])):
            opt[name[i]] = ''
            iratio = -1
            ratioScaleHisto = ROOT.TObject()
            for ratio in Ratios[iname][i]:
                iratio = iratio+1
                Pads[iname][i][1].cd()
                ROOT.gPad.SetTicks(1,1)
                if opt[name[i]] == '':
                    #print("OK, DRAWING THE RATIO DUMMY HISTO!")
                    ratioScaleHisto = ROOT.TH2D(ratio.GetName()+'_tmp_{}_{}_{}'.format(i,iratio,k), ratio.GetName()+'_tmp' + ';;' + ratioTag,
                                                ratio.GetNbinsX(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax(),
                                                100, gyratioMin, gyratioMax)
                    ratioScaleHisto.SetStats(0)
                    stuff.append(ratioScaleHisto)
                    ratioScaleHisto.GetYaxis().SetTitleOffset(0.55)
                    if addRatioTitle:
                        ratioScaleHisto.GetYaxis().SetTitle('Ratio to {:}'.format(tags[rfiles[ratioIndex].GetName()]))
                    ratioScaleHisto.GetXaxis().SetLabelSize(0.085)
                    ratioScaleHisto.GetYaxis().SetLabelSize(0.085)
                    ratioScaleHisto.GetXaxis().SetTitle(MakePrettyTitle(ratio.GetXaxis().GetTitle()))
                    ratioScaleHisto.GetXaxis().SetTitleSize(0.095)
                    ratioScaleHisto.GetYaxis().SetTitleOffset(0.7)
                    ratioScaleHisto.GetYaxis().SetTitleSize(0.08) # 0.095
                    ratioScaleHisto.Draw("")
                    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
                    line.SetLineColor(gcol[ratioIndex])
                    line.Draw()
                    lines.append(line)
                    tmpScale.append(ratioScaleHisto)
                    opt[name[i]] = 'same'
                if iratio == ratioIndex:
                    ratio.SetFillStyle(3001)
                    ratio.SetFillColor(ROOT.kGray)
                    ratio.Draw( gAxisOpt + 'e2' + opt[name[i]]) # error band!
                else:
                    ratio.Draw(gDrawOpt + gAxisOpt + opt[name[i]])
                Pads[iname][i][1].Update()
                #ratioScaleHisto.GetYaxis().SetRangeUser(0.11+gyratioMin, 0.9*gyratioMax)
                ratioScaleHisto.GetYaxis().SetRangeUser(0.05+gyratioMin, 0.96*gyratioMax)
                ROOT.gPad.Update()


    ican = 0
    for can in cans:

        normtag = ''
        if gNormalize:
            normtag = '_norm'

        can.Update()
        can.Print(pngdir + can.GetName() + '_liny' + bwtag + normtag + '.png')
        can.Print(pdfdir + can.GetName() + '_liny' + bwtag + normtag + '.pdf')

        extraSFleft = 1.
        extraSFright = 1.
        Hists[ican][0][0].SetMaximum(extraSFleft*gySFlog*Hists[ican][0][0].GetMaximum())
        Pads[ican][0][0].SetLogy(1)
        Hists[ican][1][0].SetMaximum(extraSFright*gySFlog*Hists[ican][1][0].GetMaximum())
        Pads[ican][1][0].SetLogy(1)
        can.Print(pngdir + can.GetName() + '_logy' + bwtag + normtag + '.png')
        can.Print(pdfdir + can.GetName() + '_logy' + bwtag + normtag + '.pdf')
        ican = ican + 1

    ################################################################    
    # stop here if you do not want to plot all the migrations;)
    if gStopAfter1D:
        print('DONE!')
        if not gBatch:
            ROOT.gApplication.Run()
        return 0

    ##########
    #   2D   #
    ##########

    cw = 800
    ch = 800
    if len(rfiles) > 4 :
        ch = 1200
    print('Drawing migrations...')

    for name in names2d:
        n = len(names2d)
        normtag = ''
        if gNormalizeMigra:
            normtag = 'norm_'
        canname = 'migra_' + normtag + cantag
        for i in range(0, nlevels):
            canname = canname + '_' + name ## TODO???
        canname = canname.replace('/', '_')
        can = ROOT.TCanvas(canname, canname, 10, 10, cw, ch)

        Red    = [ 1.00, 0.43, 0.05]
        Green  = [ 1.00, 0.78, 0.59]
        Blue   = [ 1.00, 0.69, 0.53]

        Red2    = [ 1.00, 0.23, 0.05]
        Green2  = [ 1.00, 0.98, 0.99]
        Blue2   = [ 1.00, 0.99, 0.93]

        Number = len (Red)
        Length = [ 0.00, 0.25, 1.00 ]
        nb=100

        if bw != 1 and gNormalizeMigra:
            tcol = ROOT.TColor()
            if canname.find('Parton_Particle') >= 0:
                tcol.CreateGradientColorTable(Number, array( 'd', Length),
                                              array( 'd', Blue),
                                              array( 'd', Red),
                                              array( 'd', Red),
                                              nb)
            elif canname.find('Parton_Detector') >= 0:
                tcol.CreateGradientColorTable(Number, array( 'd', Length),
                                              array( 'd', Red),
                                              array( 'd', Blue),
                                              array( 'd', Green),
                                              nb)
            elif canname.find('Particle_Detector') >= 0:
                tcol.CreateGradientColorTable(Number, array( 'd', Length),
                                              array( 'd', Red),
                                              array( 'd', Green),
                                              array( 'd', Blue),
                                              nb)
            elif canname.find('ParticleMatch_DetectorMatch') >= 0:
                tcol.CreateGradientColorTable(Number, array( 'd', Length),
                                              array( 'd', Red2),
                                              array( 'd', Green2),
                                              array( 'd', Blue2),
                                              nb)



        nn = int(sqrt(len(rfiles)))
        nx = nn
        ny = nx
        while nx*ny < len(rfiles):
            #print('    {:} {:}'.format(nx,ny))
            if nx < ny:
                nx = nx+1
            else:
                ny = ny+1
        #print('Dividing Canvas by {:} {:}'.format(nx,ny))
        if (len(rfiles) < 3):
            can.Divide(len(rfiles),1)
        else:
            can.Divide(nx,ny)
        cans.append(can)

        hname = name
        print('Working on {:}'.format(hname,))
        for rfile in rfiles:
            ifile = rfiles.index(rfile)
            ftag = rfile.GetName()
            ftag = ftag.replace('analyzed_histos_','').replace('_',' ').replace('.root','').replace('ATLAS','Delphes ATLAS').replace('CMS','Delphes CMS')
            ftag = ftag.replace('ptj1min60 ptj2min60', 'p_{T}^{j1,j2} > 60 GeV').replace('allhad','').replace('all','')
            hist = rfile.Get(hname)
            hists.append(hist)
            hist.SetStats(0)
            ifile = rfiles.index(rfile)
            can.cd(ifile+1)
            if not gNormalizeMigra:
                ROOT.gPad.SetLogz(1)
            corr = hist.GetCorrelationFactor()
            if gNormalizeMigra:
                NormalizeByColumns(hist)
                hist.SetMinimum(0.)
                hist.SetMaximum(1.)
                hist.GetXaxis().SetTitle(MakePrettyTitle(hist.GetXaxis().GetTitle()))
                hist.GetYaxis().SetTitle(MakePrettyTitle(hist.GetYaxis().GetTitle()))
            addopt = ''    
            if 'SelectionMigra' in hist.GetName() or gNormalizeMigra:
                addopt = ' text'
            if 'SelectionMigra' in hist.GetName():
                hist.SetMarkerSize(2)
            hist.Draw("colz" + gAxisOpt + addopt)
            if 'SelectionMigra' in hist.GetName():
                textX = ROOT.TLatex(0.55, 0.02, 'Particle level')
                textX.SetNDC()
                textX.Draw()
                txts.append(textX)
                textY = ROOT.TLatex(0.03, 0.55, 'Detector level')
                textY.SetTextAngle(90)
                textY.SetNDC()
                textY.Draw()
                txts.append(textY)
            hist.GetYaxis().SetTitleOffset(1.3)
            if PlotAlsoProfile:
                prof = hist.ProfileX(hist.GetName() + ftag + '_prof_{:}'.format(len(Objs)) )
                Objs.append(prof)
                profcol = ROOT.kBlack # ROOT.kWhite
                prof.SetStats(0)
                prof.SetMarkerColor(profcol)
                prof.SetLineColor(profcol)
                prof.SetMarkerSize(1)
                prof.SetMarkerStyle(2)
                #prof.Draw('e1same')
                prof.Draw('Psame' + gAxisOpt)
                Objs.append(prof)
            if PlotAlsoProjection:
                proj = hist.ProjectionX(hist.GetName() + ftag + '_proj_{:}'.format(len(Objs)) )
                if divideByBinWidth:
                    DivideByBinWidth(proj)
                if proj.GetMaximum() > 0.:
                    proj.Scale(hist.GetYaxis().GetXmax()*0.90/proj.GetMaximum())
                Objs.append(proj)
                projcol = ROOT.kBlack # ROOT.kWhite
                proj.SetStats(0)
                #proj.SetMarkerColor(projcol)
                proj.SetLineColor(projcol)
                #proj.SetMarkerSize(1)
                proj.SetMarkerStyle(4)
                #proj.Draw('e1same' + gAxisOpt)
                #proj.Draw('Psame' + gAxisOpt)
                proj.Draw('hist same' + gAxisOpt)
                Objs.append(proj)



            #hist.Draw("colz text")
            xx = 0.14
            yy = 0.95 # 852
            text = ROOT.TLatex(xx, yy, '{:}'.format(ftag,))
            text.SetNDC()
            text.SetTextSize(0.035)
            ### jk 10.3.2021 text.Draw()
            Objs.append(text)
            xx = 0.65
            yy = 0.162
            rtext = ROOT.TLatex(xx, yy, '#rho = {:1.2f}'.format(corr,))
            rtext.SetNDC()
            rtext.SetTextSize(0.06)
            rtext.Draw()
            Objs.append(rtext)
            if drawDiag:
                line = ROOT.TLine(hist.GetXaxis().GetXmin(), hist.GetYaxis().GetXmin(),hist.GetXaxis().GetXmax(), hist.GetYaxis().GetXmax() )
                if bw == 1:
                    line.SetLineColor(ROOT.kGray+2)
                else:
                    line.SetLineColor(ROOT.kRed)
                line.SetLineWidth(2)
                line.Draw()
                Objs.append(line)
            ROOT.gPad.RedrawAxis()
            ROOT.gPad.Update()

        can.Print(pngdir + canname + bwtag +  '.png')
        can.Print(pdfdir + canname + bwtag + '.pdf')


    print('DONE!')
    if not gBatch:
        ROOT.gApplication.Run()


#############################################################
#############################################################
#############################################################

if __name__ == "__main__":
    main(sys.argv)

#############################################################
#############################################################
#############################################################
