#!/usr/bin/python

# jk 3.3.2017, 30.3.2017, 7.4.2017, 6.9.2017, 7--16.1.2019, May 2019, May 2020, Sept 2020
# Example running: ./python/compareReplicas.py rootfile.root

#############################################################

from __future__ import print_function

import ROOT
import sys, os, getopt
# or? import argparse ?;-)

from Tools import *
from replicaItems import *
from math import *

from array import array

from ROOT import std,Double,AddressOf

stuff = []

### !!! ;-)
nReplicas = 100


#############################################################
#############################################################
#############################################################

    
def main(argv):


    #############################################################

    gDrawOpt = ' PLC ' # 'C hist' # PLC PMC AMC PFC
    gAxisOpt = ''

    addRatioTitle = False
    # draw diagonal in migration matrices:
    drawDiag = False

    #PlotAlsoProfile = False
    PlotAlsoProfile = True

    #PlotAlsoProjection = False
    PlotAlsoProjection = True

    divideByBinWidth = True # TODO: so far effective only on projections!
    intInLeg = False
    #############################################################

    bw=0
    bwtag = ''

    doratio='ratio'
    ratioIndex = 0 # file index to do ratios w.r.t.
    ratioTag = 'ratio to replica 0   '

    # some more nasty global variables;)

    indices = range(0, 5) + range(-9,-1)
    baseCols = [
        ROOT.kRed,
        ROOT.kBlue, 
        ROOT.kGreen,
        ROOT.kMagenta,
        ROOT.kCyan,
        ROOT.kPink,
        ROOT.kYellow,
        ROOT.kSpring]
    gcol = [    ]
    for bc in baseCols:
        gcol = gcol + [ bc + ind for ind in indices]

    LevelLabel = ['detector', 'detector']

    gyratioMin = 0. # 0.55
    gyratioMax = 2. # 1.45
    gySFlin = 1.9
    gySFlog = 100.

    #############################################################
    # command-line controllable:
    gBatch        = False
    gNormalize    = False
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
    cantag='replicas'

    print(argv)
    fname = ''
    if len(argv) > 1:
        fname = argv[1]
        print('OK, using user-defined histograms root file {:}'.format(fname,) )
    else:
        print('ERROR: you MUST provide a ROOT filename!')
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
        elif opt in ("-t", "--tag"):
            gPngTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gPngTag,) )

    print('Done.')

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    print('*** Settings:')
    print('tag={:}, normalize={:}, batch={:}'.format(gPngTag, gNormalize, gBatch,))
    print('files:')
    print(fname)
    print('')

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

    rfile = ROOT.TFile(fname, 'read')
    gPngTag = fname.replace('analyzed_histos_', '').replace('.root','')
    stuff.append(rfile)
    
    for name in names1d:
        iname = names1d.index(name)
        DualHists = []
        opt = {}
        print('Processing 1D histos {:}'.format(name,))
        nlevels = len(name) # is usually 2: generator and detector levels;-)
        canname = cantag
        for i in range(0, nlevels):
            canname = canname + '_' + name[i]
        canname = canname.replace('/', '_')
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
            leg = ROOT.TLegend(0.12, 0.15+0.44+0.20, 0.88, 0.40+0.44+0.07)
            topotag = name[i].split('/')[0]
            leg.SetHeader('Selection: {}'.format(topotag))
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            legs.append(leg)
            # TODO:
            # make this a while loop, trying to get replicas,
            # break at exception
            for irep in range(0,nReplicas):
                hname = name[i]
                if 'denser' in hname:
                    hname = hname.replace('_denser',  '_rep{}_denser'.format(irep))
                else:
                    hname = name[i] + '_rep{}'.format(irep)
                print('       getting {:}'.format(hname,))
                hist = rfile.Get(hname)
                if divideByBinWidth:
                    DivideByBinWidth(hist)
                hist.SetStats(0)
                #hist.SetLineColor(gcol[irep])
                #hist.SetLineWidth(1)
                #hist.SetLineStyle(gls[ifile])
                inttag = ""
                if intInLeg:
                    inttag = ' I={:.0f}'.format(hist.Integral())

                if gNormalize:
                    norm = Double(hist.Integral())
                    if norm > 0.:
                        hist.Scale(1./norm)
                    else:
                        print('ERROR! Histogram {:} from file {:} has ZERO integral!'.format(hist.GetName(), rfile.GetName(), ))

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
                hist.Draw(opt[name[i]] + gDrawOpt + gAxisOpt)
                hist.GetYaxis().SetMoreLogLabels()
                opt[name[i]] = ' same '
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
                    ratioScaleHisto.GetYaxis().SetTitleOffset(0.55)
                    if addRatioTitle:
                        ratioScaleHisto.GetYaxis().SetTitle('Ratio to {:}'.format(tags[rfiles[ratioIndex].GetName()]))
                    ratioScaleHisto.GetXaxis().SetLabelSize(0.085)
                    ratioScaleHisto.GetYaxis().SetLabelSize(0.085)
                    ratioScaleHisto.GetXaxis().SetTitle(MakePrettyTitle(ratio.GetXaxis().GetTitle()))
                    ratioScaleHisto.GetXaxis().SetTitleSize(0.095)
                    ratioScaleHisto.GetYaxis().SetTitleSize(0.095)
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


    ican = 0
    for can in cans:

        normtag = ''
        if gNormalize:
            normtag = '_norm'

        Pads[ican][0][0].Update()
        Pads[ican][1][0].Update()
        can.Print(pngdir + can.GetName() + '_liny' + bwtag + normtag + '.png')
        can.Print(pdfdir + can.GetName() + '_liny' + bwtag + normtag + '.pdf')

        Hists[ican][0][0].SetMaximum(gySFlog*Hists[ican][0][0].GetMaximum())
        Pads[ican][0][0].SetLogy(1)
        Pads[ican][0][0].Update()
        Hists[ican][1][0].SetMaximum(gySFlog*Hists[ican][1][0].GetMaximum())
        Pads[ican][1][0].SetLogy(1)
        Pads[ican][1][0].Update()
        can.Print(pngdir + can.GetName() + '_logy' + bwtag + normtag + '.png')
        can.Print(pdfdir + can.GetName() + '_logy' + bwtag + normtag + '.pdf')
        ican = ican + 1

    ################################################################    

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
