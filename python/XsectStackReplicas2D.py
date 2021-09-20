#!/usr/bin/python

from __future__ import print_function

# jk 14.4.2020 (covid-19 homeoffice), based on code of StackPlots.py
# Example running: ./python/XsectStack.py ./python/list.txt "_mytag" notnorm batch
# 10.9.2020, 7.11.2020
# 12.8.2021, 26.8.2021

import ROOT

from collections import OrderedDict

from xSectTools import *
from BumpSignifTools import *
from FitTools import *
from LhoodFitTools import *

from CorrItems import *

from cTopPlot import *

from math import *
import sys, os
from array import array

##########################################
class resultsClass:
    # to store histograms of BumpHunter scores, background compatibility scores
    # but also of the fitted signal strength mu and the related zero mu compatibility score
    def __init__(self):
        # BH and bg. hypothese check scores
        self.BHscoreHists = OrderedDict()
        self.BGscoreHists = OrderedDict()

        self.muscoreHists = OrderedDict()
        self.muHists = OrderedDict()

        self.bhmuscoreHists = OrderedDict()
        self.bhmuHists = OrderedDict()


        # originally were histos for data, will make these histos for axes;-)
        self.bhmuscorevsmuscoreHists2d = OrderedDict()
        self.BHvsBGHists2d = OrderedDict()
        self.muscorevsBGHists2d = OrderedDict()
        self.muscorevsBHHists2d = OrderedDict()
        self.bhmuscorevsBGHists2d = OrderedDict()
        self.bhmuscorevsBHHists2d = OrderedDict()

        # make these also TGraphs?
        # ok, it's a graph, not a tgraph2d but I will keep the 2d tag in name;-)
        self.bhmuscorevsmuscoreGr2d = OrderedDict()
        self.BHvsBGGr2d = OrderedDict()
        self.muscorevsBGGr2d = OrderedDict()
        self.muscorevsBHGr2d = OrderedDict()
        self.bhmuscorevsBGGr2d = OrderedDict()
        self.bhmuscorevsBHGr2d = OrderedDict()


    def MakeAlsoGraphDictsFromHisto2dDicts(self):
        # make these lists:
        dictsToGo = [[  self.bhmuscorevsmuscoreHists2d , self.bhmuscorevsmuscoreGr2d],
                     [  self.BHvsBGHists2d             , self.BHvsBGGr2d],
                     [  self.muscorevsBGHists2d        , self.muscorevsBGGr2d],
                     [  self.muscorevsBHHists2d        , self.muscorevsBHGr2d],
                     [  self.bhmuscorevsBGHists2d      , self.bhmuscorevsBGGr2d],
                     [  self.bhmuscorevsBHHists2d      , self.bhmuscorevsBHGr2d],
        ]
        for dicts in dictsToGo:
            hdict = dicts[0]
            grdict = dicts[1]
            for key in hdict:
                grdict[key] = ROOT.TGraph()
                grdict[key].SetName('gr_' + hdict[key].GetName())
                grdict[key].SetLineColor(hdict[key].GetMarkerColor())
                grdict[key].SetMarkerStyle(4)
                grdict[key].SetMarkerSize(1)
        return
                
    
##########################################
def SetGraphPoint(gr, x, y):
    gr.SetPoint(gr.GetN(), x, y)
    return

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

    cansBH = []
    cansBG = []
    cansmu = []
    cansmuscore = []
    
    Hists = []
    Objs = []
    rfiles = []
    rfilesStatIndep = []
    Legs = []
    Pads = []
    tags = {}


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

    nReplicas = 100 # HACKS !!! DEFAULT: 100
    printAnyway = False

    ScoreHists = OrderedDict()
    rhoVals = OrderedDict()
    rhoValsAver = OrderedDict()
    
    print('*** Processing 2D histograms to stack;-) ***')

    names = [ #'CosThetaStarVsDiTopPout',
              #'RttbarVsDiTopMass',
              'DiTopMassVsDiTopPout',
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
        scoreMin = 2
        scoreMax = 8 #!!!8 # 7
        nscoreBins = int(6*(scoreMax - scoreMin)) #6
        # hack:
        #nscoreBins = 100*(scoreMax - scoreMin) #6

        bgscoreMin = -4
        bgscoreMax = 4 #!!!8 # 7
        nbgscoreBins = int(6*(bgscoreMax - bgscoreMin)) #6
        bgsname = sname.replace('BHscore', 'BGscore')
        bgtitle = bgsname + ';BG score;Replicas'

        muMin = -0.5
        muMax = 5.
        nmuBins = 90 # 6*(muMax - muMin) #6        
        musname = sname.replace('BHscore', 'mu')
        mutitle = musname + ';fitted #mu;Replicas'

        musMin = -4
        musMax = 8
        nmusBins = 6*(musMax - musMin) #6        
        mussname = sname.replace('BHscore', 'mu score')
        mustitle = musname + ';#mu score;Replicas'

        # BH assisted mu score = bhmu and bhmuscore
        bhmuMin = -0.5
        bhmuMax = 5.
        nbhmuBins = 101 # 6*(bhmuMax - bhmuMin) #6        
        bhmusname = sname.replace('BHscore', 'bhmu')
        bhmutitle = bhmusname + ';fitted #mu_{BH};Replicas'

        bhmusMin = -4
        bhmusMax = 8
        nbhmusBins = 6*(bhmusMax - bhmusMin) #6        
        bhmussname = sname.replace('BHscore', 'bhmu score')
        bhmustitle = bhmusname + ';#mu_{BH} score;Replicas'

        for topo in Topos:
            try:
                print(ScoreHists[topo])
            except:
                ScoreHists[topo] = OrderedDict()
                rhoVals[topo] = OrderedDict()
                rhoValsAver[topo] = OrderedDict()

            # create a triplet of BH score histograms: 1Dx, 1Dy, 2D and correlations for averaging
            rhoVals[topo][barehname] = []
            ScoreHists[topo][barehname] = resultsClass()
            # TODO: rename the name and title to respective x and y variables!
            print('*** Creating histos for {} {}'.format(topo,barehname))

            ScoreHists[topo][barehname].BHscoreHists[v1Dx] = ROOT.TH1D(sname + '_' + topo + '_' + v1Dx, title, nscoreBins, scoreMin, scoreMax)
            ScoreHists[topo][barehname].BHscoreHists[v1Dy] = ROOT.TH1D(sname + '_' + topo + '_' + v1Dy, title, nscoreBins, scoreMin, scoreMax)
            ScoreHists[topo][barehname].BHscoreHists[v2D]  = ROOT.TH1D(sname + '_' + topo + '_' + v2D,  title, nscoreBins, scoreMin, scoreMax)

            ScoreHists[topo][barehname].BGscoreHists[v1Dx] = ROOT.TH1D(bgsname + '_' + topo + '_' + v1Dx, bgtitle, nbgscoreBins, bgscoreMin, bgscoreMax)
            ScoreHists[topo][barehname].BGscoreHists[v1Dy] = ROOT.TH1D(bgsname + '_' + topo + '_' + v1Dy, bgtitle, nbgscoreBins, bgscoreMin, bgscoreMax)
            ScoreHists[topo][barehname].BGscoreHists[v2D]  = ROOT.TH1D(bgsname + '_' + topo + '_' + v2D,  bgtitle, nbgscoreBins, bgscoreMin, bgscoreMax)

            # lhood mu
            # mu histos
            ScoreHists[topo][barehname].muHists[v1Dx] = ROOT.TH1D(musname + '_' + topo + '_' + v1Dx, mutitle, nmuBins, muMin, muMax)
            ScoreHists[topo][barehname].muHists[v1Dy] = ROOT.TH1D(musname + '_' + topo + '_' + v1Dy, mutitle, nmuBins, muMin, muMax)
            ScoreHists[topo][barehname].muHists[v2D]  = ROOT.TH1D(musname + '_' + topo + '_' + v2D,  mutitle, nmuBins, muMin, muMax)

            ScoreHists[topo][barehname].bhmuHists[v1Dx] = ROOT.TH1D(bhmusname + '_' + topo + '_' + v1Dx, bhmutitle, nbhmuBins, bhmuMin, bhmuMax)
            ScoreHists[topo][barehname].bhmuHists[v1Dy] = ROOT.TH1D(bhmusname + '_' + topo + '_' + v1Dy, bhmutitle, nbhmuBins, bhmuMin, bhmuMax)
            ScoreHists[topo][barehname].bhmuHists[v2D]  = ROOT.TH1D(bhmusname + '_' + topo + '_' + v2D,  bhmutitle, nbhmuBins, bhmuMin, bhmuMax)

            # lhood mu pval score histos:
            ScoreHists[topo][barehname].muscoreHists[v1Dx] = ROOT.TH1D(mussname + '_' + topo + '_' + v1Dx, mustitle, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].muscoreHists[v1Dy] = ROOT.TH1D(mussname + '_' + topo + '_' + v1Dy, mustitle, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].muscoreHists[v2D]  = ROOT.TH1D(mussname + '_' + topo + '_' + v2D,  mustitle, nmusBins, musMin, musMax)

            # lhood bhmu pval score histos:
            ScoreHists[topo][barehname].bhmuscoreHists[v1Dx] = ROOT.TH1D(bhmussname + '_' + topo + '_' + v1Dx, bhmustitle, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscoreHists[v1Dy] = ROOT.TH1D(bhmussname + '_' + topo + '_' + v1Dy, bhmustitle, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscoreHists[v2D]  = ROOT.TH1D(bhmussname + '_' + topo + '_' + v2D,  bhmustitle, nbhmusBins, bhmusMin, bhmusMax)


            # 2D histos:
            # TO CHECK!!!
            
            ScoreHists[topo][barehname].muscorevsBHHists2d[v1Dx]  = ROOT.TH2D(mussname + '_vs_' + sname + '_' + topo + '_' + v1Dx,
                                                                              ';BH score;#mu score;',
                                                                              nscoreBins, scoreMin, scoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].muscorevsBGHists2d[v1Dx]  = ROOT.TH2D(mussname + '_vs_' + bgsname + '_' + topo + '_' + v1Dx,
                                                                              ';BG score;#mu score;',
                                                                              nbgscoreBins, bgscoreMin, bgscoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].bhmuscorevsBHHists2d[v1Dx]  = ROOT.TH2D(bhmussname + '_vs_' + sname + '_' + topo + '_' + v1Dx,
                                                                              ';BH score;#mu_{BH} score;',
                                                                              nscoreBins, scoreMin, scoreMax, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscorevsBGHists2d[v1Dx]  = ROOT.TH2D(bhmussname + '_vs_' + bgsname + '_' + topo + '_' + v1Dx,
                                                                              ';BG score;#mu_{BH} score;',
                                                                              nbgscoreBins, bgscoreMin, bgscoreMax, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].BHvsBGHists2d[v1Dx]  = ROOT.TH2D(sname + '_vs_' + bgsname + '_' + topo + '_' + v1Dx,
                                                                         ';BG score;BH score;',
                                                                         nbgscoreBins, bgscoreMin, bgscoreMax, nscoreBins, scoreMin, scoreMax)

            ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d[v1Dx]  = ROOT.TH2D(musname + '_vs_' + bhmusname + '_' + topo + '_' + v1Dx,
                                                                                     ';#mu score;#mu_{BH} score',
                                                                                     nmusBins, musMin, musMax, nbhmusBins, bhmusMin, bhmusMax)
            
            ScoreHists[topo][barehname].muscorevsBHHists2d[v1Dy]  = ROOT.TH2D(mussname + '_vs_' + sname + '_' + topo + '_' + v1Dy,
                                                                              ';BH score;#mu score;',
                                                                              nscoreBins, scoreMin, scoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].muscorevsBGHists2d[v1Dy]  = ROOT.TH2D(mussname + '_vs_' + bgsname + '_' + topo + '_' + v1Dy,
                                                                              ';BG score;#mu score;',
                                                                              nbgscoreBins, bgscoreMin, bgscoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].bhmuscorevsBHHists2d[v1Dy]  = ROOT.TH2D(bhmussname + '_vs_' + sname + '_' + topo + '_' + v1Dy,
                                                                                ';BH score;#mu_{BH} score;',
                                                                              nscoreBins, scoreMin, scoreMax, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscorevsBGHists2d[v1Dy]  = ROOT.TH2D(bhmussname + '_vs_' + bgsname + '_' + topo + '_' + v1Dy,
                                                                              ';BG score;#mu_{BH} score;',
                                                                              nbgscoreBins, bgscoreMin, bgscoreMax, nbhmusBins, bhmusMin, bhmusMax)



            
            ScoreHists[topo][barehname].BHvsBGHists2d[v1Dy]  = ROOT.TH2D(sname + '_vs_' + bgsname + '_' + topo + '_' + v1Dy,
                                                                         ';BG score;BH score;',
                                                                         nbgscoreBins, bgscoreMin, bgscoreMax, nscoreBins, scoreMin, scoreMax)
            ScoreHists[topo][barehname].BHvsBGHists2d[v2D]  = ROOT.TH2D(sname + '_vs_' + bgsname + '_' + topo + '_' + v2D,
                                                                        ';BG score;BH score;',
                                                                        nbgscoreBins, bgscoreMin, bgscoreMax, nscoreBins, scoreMin, scoreMax)
            ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d[v1Dy]  = ROOT.TH2D(musname + '_vs_' + bhmusname + '_' + topo + '_' + v1Dy,
                                                                                     ';#mu score;#mu_{BH} score',
                                                                                     nmusBins, musMin, musMax, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d[v2D]  = ROOT.TH2D(musname + '_vs_' + bhmusname + '_' + topo + '_' + v2D,
                                                                                     ';#mu score;#mu_{BH} score',
                                                                                     nmusBins, musMin, musMax, nbhmusBins, bhmusMin, bhmusMax)



            ScoreHists[topo][barehname].muscorevsBHHists2d[v2D]  = ROOT.TH2D(mussname + '_vs_' + sname + '_' + topo + '_' + v2D,
                                                                             ';BH score;#mu score;',
                                                                             nscoreBins, scoreMin, scoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].muscorevsBGHists2d[v2D]  = ROOT.TH2D(mussname + '_vs_' + bgsname + '_' + topo + '_' + v2D,
                                                                             ';BG score;#mu score;',
                                                                             nbgscoreBins, bgscoreMin, bgscoreMax, nmusBins, musMin, musMax)
            ScoreHists[topo][barehname].bhmuscorevsBHHists2d[v2D]  = ROOT.TH2D(bhmussname + '_vs_' + sname + '_' + topo + '_' + v2D,
                                                                             ';BH score;#mu_{BH} score;',
                                                                             nscoreBins, scoreMin, scoreMax, nbhmusBins, bhmusMin, bhmusMax)
            ScoreHists[topo][barehname].bhmuscorevsBGHists2d[v2D]  = ROOT.TH2D(bhmussname + '_vs_' + bgsname + '_' + topo + '_' + v2D,
                                                                             ';BG score;#mu_{BH} score;',
                                                                             nbgscoreBins, bgscoreMin, bgscoreMax, nbhmusBins, bhmusMin, bhmusMax)

            ScoreHists[topo][barehname].MakeAlsoGraphDictsFromHisto2dDicts()
            
            
            # set colors:
            dictsToGo = [ScoreHists[topo][barehname].BHscoreHists,
                         ScoreHists[topo][barehname].BGscoreHists,
                         ScoreHists[topo][barehname].muHists,
                         ScoreHists[topo][barehname].muscoreHists,
                         ScoreHists[topo][barehname].bhmuHists,
                         ScoreHists[topo][barehname].bhmuscoreHists
            ]
            for scorehist in dictsToGo:
                for vkey in scorehist:
                    hh = scorehist[vkey]
                    hh.SetStats(0)
                    #hh.SetLineWidth(1)
                    #if 'VsDiTopMass' in hh.GetName() and vkey == v1Dx:
                    hh.SetLineWidth(3)
                    hh.SetLineColor(vcols[vkey])
                    hh.SetFillColor(vcols[vkey])
                    hh.SetFillStyle(vfills[vkey])
            
            dictsToGo2d = [ScoreHists[topo][barehname].muscorevsBHHists2d,
                           ScoreHists[topo][barehname].muscorevsBGHists2d,
                           ScoreHists[topo][barehname].bhmuscorevsBHHists2d,
                           ScoreHists[topo][barehname].bhmuscorevsBGHists2d,
                           ScoreHists[topo][barehname].BHvsBGHists2d,
                           ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d,
            ]
            for scorehist in dictsToGo2d:
                for vkey in scorehist:
                    h2 = scorehist[vkey]
                    h2.SetStats(0)
                    h2.SetMarkerColor(vcols[vkey])
                    h2.SetMarkerSize(0.7)
                    h2.SetMarkerStyle(4)
            grsToGo2d = [ScoreHists[topo][barehname].muscorevsBHGr2d,
                           ScoreHists[topo][barehname].muscorevsBGGr2d,
                           ScoreHists[topo][barehname].bhmuscorevsBHGr2d,
                           ScoreHists[topo][barehname].bhmuscorevsBGGr2d,
                           ScoreHists[topo][barehname].BHvsBGGr2d,
                           ScoreHists[topo][barehname].bhmuscorevsmuscoreGr2d,
            ]
            for Gr in grsToGo2d:
                for vkey in Gr:
                    gr = Gr[vkey]
                    gr.SetMarkerColor(vcols[vkey])
                    gr.SetMarkerSize(0.7)
                    gr.SetMarkerStyle(4)

        ###############################
        #     LOOP OVER REPLICAS!     #
        ###############################

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
                        if abs(addSignalSFsPower[topotag]) < kEpsilon:
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
                
                #print('Pads:')
                #print(Pads[-1])
                #print('Histos in stack: {}, histo0 entries: {}'.format(stack.GetNhists(), stack.GetHists().At(0).GetEntries()))
                ToPlot.append( cTopPlot(stack, stackStatIndep, leg, legh, legitems, Pads[-1], sPads, hname, barehname) )
                print('...checking hdata integral: {}'.format(ToPlot[-1].hdata.Integral()))
                Hists.append(hists)
                Stacks.append(stack)

    ### now draw the stack and data! ;-) ###
    ### TODO: subtract BG!

    #print('=== Scoreicance histos to be filled ===')
    #print(ScoreHists)
    
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

        ### I. BUMP HUNTER
        
        # BumpHunter2D:
        ytxt = 0.08 # 0.02
        bhcol = vcols[v2D]
        bhhcol = vcols[v2D]
        print('Computing 2D BH test from data,bg of integrals {},{}'.format(hdata.Integral(), hbg.Integral()))
        bhresult,bhtex = FindBestBumpHunter2DArea(hdata, hbg, 0.02, ytxt, 'BH log(t) xy: ')
        used2dBH = bhresult.bins
        usedLinesBH = DrawUsedBinsGrid(hdata, used2dBH, bhresult.i, bhresult.j, bhcol, bhhcol, 1)

        # JK 29.7.2021, 11.8.2021
        # compute the full 2D range chi2 probability and its log of data compatibility with the null (H0) or alternatve (H1) hypothesis
        #fitresultNull, fitresultSig, p0, t0, logt0, p1, t1, logt1, stuffbg = AnalyzeCompatibilityH0H1(hdata, hbg, hsig)
        # compatibility of data shape with bg only:
        bhtx = 0.32
        fitresultNull, p0, t0, logt0, stuffbg, bgtext2d = AnalyzeDataBgCompatibility(hdata, hbg, bhtx, 0.06, 'xy: ')
        print('2D Null hypothesis probability p0: {}'.format(p0))
        print('2D t0 = -log(p0): {}'.format(t0))
        print('2D log(t0): {}'.format(logt0))

        # REVIVE???
        # Sinificance:
        # ytxt = 0.08
        #scol = ROOT.kMagenta + 2
        #shcol = ROOT.kMagenta + 2
        #scoreicance, stex = ComputeExcessScoreicance2D(hdata, hbg, 0.02, 0.08)
        #scoreicance, stex, n2d, used2d, i0, j0 = ComputeExcessScoreicance2DIteratively(hdata, hbg, 0.02, ytxt, 'Score. xy: ' )
        #usedLines = DrawUsedBinsGrid(hdata, used2d, i0, j0, scol, shcol, 2)
        
         #1D's:
        hdatax = hdata.ProjectionX(hdata.GetName() + '_projX')
        hdatay = hdata.ProjectionY(hdata.GetName() + '_projY')
        hbgx = hbg.ProjectionX(hbg.GetName() + '_projX')
        hbgy = hbg.ProjectionY(hbg.GetName() + '_projY')
        hsigx = hsig.ProjectionX(hsig.GetName() + '_projX')
        hsigy = hsig.ProjectionY(hsig.GetName() + '_projY')
        
        
        ## old, full range: scoreicance1Dx, stex1Dx = ComputeExcessScoreicance1D(hdatax, hbgx, 0.02, 0.02, '1D score. x: ')
        ## old, full range: scoreicance1Dy, stex1Dy = ComputeExcessScoreicance1D(hdatay, hbgy, 0.55, 0.02, '1D score. y: ')
        #scoreicance1Dx, stex1Dx, nx, usedx, i0x = ComputeExcessScoreicance1DIteratively(hdatax, hbgx, 0.345, ytxt, 'x: ')
        #scoreicance1Dy, stex1Dy, ny, usedy, i0y = ComputeExcessScoreicance1DIteratively(hdatay, hbgy, 0.500, ytxt, 'y: ')

        # BH 1D:
        bhresultx,bhtexx = FindBestBumpHunter1DInterval(hdatax, hbgx, 0.345, ytxt, 'x: ')
        bhresulty,bhtexy = FindBestBumpHunter1DInterval(hdatay, hbgy, 0.500, ytxt, 'y: ')
        usedLinesBHx = DrawUsedBinsLines(hdatax, hdatay, bhresultx.bins, bhresultx.i, True,  0.08, 0.005, bhcol, bhhcol, 1, 3, 3)
        usedLinesBHy = DrawUsedBinsLines(hdatay, hdatax, bhresulty.bins, bhresulty.i, False, 0.10, 0.005, bhcol, bhhcol, 1, 3, 3)
        #usedLinesx   = DrawUsedBinsLines(hdatax, hdatay, usedx, i0x,  True,  0.10, 0.005, scol, shcol, 1)
        #usedLinesy   = DrawUsedBinsLines(hdatay, hdatax, usedy, i0y,  False, 0.12, 0.005, scol, shcol, 1)
        
        #stex1Dx.Draw()
        #stex1Dy.Draw()
        #if scoreicance > scoreicance1Dy and scoreicance > scoreicance1Dx:
        #    stex.SetTextColor(hsig.GetFillColor())
        #stex.Draw()

        # II. compatibility of data shape with bg only 1Dx:
        fitresultNullx, p0x, t0x, logt0x, stuffbgx, bgtext1dx = AnalyzeDataBgCompatibility(hdatax, hbgx, 0.54, 0.06, 'x: ')
        print('1Dx Null hypothesis probability p0: {}'.format(p0x))
        print('1Dx t0 = -log(p0): {}'.format(t0x))
        print('1Dx log(t0): {}'.format(logt0x))
        # compatibility of data shape with bg only 1Dy:
        fitresultNully, p0y, t0y, logt0y, stuffbgy, bgtext1dy = AnalyzeDataBgCompatibility(hdatay, hbgy, 0.74, 0.06, 'y: ')
        print('1Dy Null hypothesis probability p0: {}'.format(p0y))
        print('1Dy t0 = -log(p0): {}'.format(t0y))
        print('1Dy log(t0): {}'.format(logt0y))
        
        # III. Lhood and mu score by null mu probability:
        # 2d:
        print('Computing the lhood fit mu...')
        mu,muerr = fitSignalStrength(hsig, hbg, hdata, False)
        pmu0, logtmu0 = ComputeZeroCompatibility(mu, muerr)
        # 1dx:
        mux,muerrx = fitSignalStrength(hsigx, hbgx, hdatax, True)
        pmu0x, logtmu0x = ComputeZeroCompatibility(mux, muerrx)
        # 1dy:
        muy,muerry = fitSignalStrength(hsigy, hbgy, hdatay, True)
        pmu0y, logtmu0y = ComputeZeroCompatibility(muy, muerry)
        

        # IV?
        # Lhood and mu score by null mu probability:
        # 2d:
        print('Computing the lhood fit bhmu...')
        bhmu,bhmuerr = fitSignalStrength(hsig, hbg, hdata, False, bhresult.bins)
        pbhmu0, logtbhmu0 = ComputeZeroCompatibility(bhmu, bhmuerr)
        print('Likelihood 2d fit result:   mu={} +/- {}, p0={} log(-log(p0))={}'.format(mu, muerr, pmu0, logtmu0))
        print('Likelihood 2d fit result: bhmu={} +/- {}, p0={} log(-log(p0))={}'.format(bhmu, bhmuerr, pbhmu0, logtbhmu0))
        # 1dx:
        bhmux,bhmuerrx = fitSignalStrength(hsigx, hbgx, hdatax, True, bhresultx.bins)
        pbhmu0x, logtbhmu0x = ComputeZeroCompatibility(bhmux, bhmuerrx)
        print('Likelihood 1Dx fit result:   mu={} +/- {}, p0={} log(-log(p0))={}'.format(mux, muerrx, pmu0x, logtmu0x))
        print('Likelihood 1Dx fit result: bhmu={} +/- {}, p0={} log(-log(p0))={}'.format(bhmux, bhmuerrx, pbhmu0x, logtbhmu0x))
        # 1dy:
        bhmuy,bhmuerry = fitSignalStrength(hsigy, hbgy, hdatay, True, bhresulty.bins)
        pbhmu0y, logtbhmu0y = ComputeZeroCompatibility(bhmuy, bhmuerry)
        print('Likelihood 1Dy fit result:   mu={} +/- {}, p0={} log(-log(p0))={}'.format(muy, muerry, pmu0y, logtmu0y))
        print('Likelihood 1Dy fit result: bhmu={} +/- {}, p0={} log(-log(p0))={}'.format(bhmuy, bhmuerry, pbhmu0y, logtbhmu0y))
        
        
        barehname = toplot.basename
        topo = toplot.fullhname.split('/')[0]
        print('   ...filling {} {}'.format(topo,barehname))

        bhlogtx = ComputeSafeLog(bhresultx.t)
        bhlogty = ComputeSafeLog(bhresulty.t)
        bhlogt = ComputeSafeLog(bhresult.t)
    
        bglogtx = ComputeSafeLog(t0x)
        bglogty = ComputeSafeLog(t0y)
        bglogt = ComputeSafeLog(t0)

        # TODO: use also the ComputeSafeLog for mu score ?!?!
        
        ScoreHists[topo][barehname].BHscoreHists['1Dx'].Fill(bhlogtx)
        ScoreHists[topo][barehname].BHscoreHists['1Dy'].Fill(bhlogty)
        ScoreHists[topo][barehname].BHscoreHists['2D'].Fill(bhlogt)

        ScoreHists[topo][barehname].BGscoreHists['1Dx'].Fill(bglogtx)
        ScoreHists[topo][barehname].BGscoreHists['1Dy'].Fill(bglogty)
        ScoreHists[topo][barehname].BGscoreHists['2D'].Fill(bglogt)

        ScoreHists[topo][barehname].muHists['1Dx'].Fill(mux)
        ScoreHists[topo][barehname].muHists['1Dy'].Fill(muy)
        ScoreHists[topo][barehname].muHists['2D'].Fill(mu)

        ScoreHists[topo][barehname].muscoreHists['1Dx'].Fill(logtmu0x)
        ScoreHists[topo][barehname].muscoreHists['1Dy'].Fill(logtmu0y)
        ScoreHists[topo][barehname].muscoreHists['2D'].Fill(logtmu0)

        ScoreHists[topo][barehname].bhmuHists['1Dx'].Fill(bhmux)
        ScoreHists[topo][barehname].bhmuHists['1Dy'].Fill(bhmuy)
        ScoreHists[topo][barehname].bhmuHists['2D'].Fill(bhmu)

        ScoreHists[topo][barehname].bhmuscoreHists['1Dx'].Fill(logtbhmu0x)
        ScoreHists[topo][barehname].bhmuscoreHists['1Dy'].Fill(logtbhmu0y)
        ScoreHists[topo][barehname].bhmuscoreHists['2D'].Fill(logtbhmu0)

        # 2D:
        """
        ScoreHists[topo][barehname].muscorevsBHHists2d['1Dx'].Fill(bhlogtx, logtmu0x)
        ScoreHists[topo][barehname].muscorevsBHHists2d['1Dy'].Fill(bhlogty, logtmu0y)
        ScoreHists[topo][barehname].muscorevsBHHists2d['2D'].Fill(bhlogt, logtmu0)
        
        ScoreHists[topo][barehname].muscorevsBGHists2d['1Dx'].Fill(bglogtx, logtmu0x)
        ScoreHists[topo][barehname].muscorevsBGHists2d['1Dy'].Fill(bglogty, logtmu0y)
        ScoreHists[topo][barehname].muscorevsBGHists2d['2D'].Fill(bglogt, logtmu0)

        ScoreHists[topo][barehname].bhmuscorevsBHHists2d['1Dx'].Fill(bhlogtx, logtbhmu0x)
        ScoreHists[topo][barehname].bhmuscorevsBHHists2d['1Dy'].Fill(bhlogty, logtbhmu0y)
        ScoreHists[topo][barehname].bhmuscorevsBHHists2d['2D'].Fill(bhlogt, logtbhmu0)
        
        ScoreHists[topo][barehname].bhmuscorevsBGHists2d['1Dx'].Fill(bglogtx, logtbhmu0x)
        ScoreHists[topo][barehname].bhmuscorevsBGHists2d['1Dy'].Fill(bglogty, logtbhmu0y)
        ScoreHists[topo][barehname].bhmuscorevsBGHists2d['2D'].Fill(bglogt, logtbhmu0)

        ScoreHists[topo][barehname].BHvsBGHists2d['1Dx'].Fill(bglogtx, bhlogtx)
        ScoreHists[topo][barehname].BHvsBGHists2d['1Dy'].Fill(bglogty, bhlogty)
        ScoreHists[topo][barehname].BHvsBGHists2d['2D'].Fill(bglogt, bhlogt)

        ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d['1Dx'].Fill(logtmu0x, logtbhmu0x)
        ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d['1Dy'].Fill(logtmu0y, logtbhmu0y)
        ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d['2D'].Fill(logtmu0, logtbhmu0)
        """
        
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBHGr2d['1Dx'], bhlogtx, logtmu0x)
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBHGr2d['1Dy'], bhlogty, logtmu0y)
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBHGr2d['2D'], bhlogt, logtmu0)
        
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBGGr2d['1Dx'], bglogtx, logtmu0x)
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBGGr2d['1Dy'], bglogty, logtmu0y)
        SetGraphPoint(ScoreHists[topo][barehname].muscorevsBGGr2d['2D'], bglogt, logtmu0)

        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBHGr2d['1Dx'], bhlogtx, logtbhmu0x)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBHGr2d['1Dy'], bhlogty, logtbhmu0y)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBHGr2d['2D'], bhlogt, logtbhmu0)
        
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBGGr2d['1Dx'], bglogtx, logtbhmu0x)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBGGr2d['1Dy'], bglogty, logtbhmu0y)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsBGGr2d['2D'], bglogt, logtbhmu0)

        SetGraphPoint(ScoreHists[topo][barehname].BHvsBGGr2d['1Dx'], bglogtx, bhlogtx)
        SetGraphPoint(ScoreHists[topo][barehname].BHvsBGGr2d['1Dy'], bglogty, bhlogty)
        SetGraphPoint(ScoreHists[topo][barehname].BHvsBGGr2d['2D'], bglogt, bhlogt)

        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsmuscoreGr2d['1Dx'], logtmu0x, logtbhmu0x)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsmuscoreGr2d['1Dy'], logtmu0y, logtbhmu0y)
        SetGraphPoint(ScoreHists[topo][barehname].bhmuscorevsmuscoreGr2d['2D'], logtmu0, logtbhmu0)

        
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
    cansBH = {}
    legsBH = {}
    cansBG = {}
    legsBG = {}
    cansmu = {}
    legsmu = {}
    cansmuscore = {}
    legsmuscore = {}

    cansbhmu = {}
    legsbhmu = {}
    cansbhmuscore = {}
    legsbhmuscore = {}

    cans2d = {}
    
    topots = {}
    baseopt = 'hist'
    sameopt = 'same'
    
    print('Plotting the score distributions...')

    #sigfiletag = sigTag.replace(' ','').replace(',','_').replace('#','').replace('{','').replace('}','').replace('=','').replace('m_','').replace('y_0','y0').replace("Z'","zp")

    rholeg_BHvsBG = OrderedDict()
    rholeg_bhmuscorevsmuscore = OrderedDict()
    rholeg_muscorevsBH = OrderedDict()
    rholeg_muscorevsBG = OrderedDict()
    rholeg_bhmuscorevsBH = OrderedDict()
    rholeg_bhmuscorevsBG = OrderedDict()

    stuff.append([rholeg_muscorevsBG,
                  rholeg_muscorevsBH,
                  rholeg_bhmuscorevsBG,
                  rholeg_bhmuscorevsBH,
                  rholeg_BHvsBG,
                  rholeg_bhmuscorevsmuscore,
    ])

    
    for topo in Topos:
        texFileName = 'tex/table_{}_BH_{}{}.tex'.format(cantag,topo,texTag) # sigfiletag,
        texTableFile = open(texFileName, 'w')
        # legends for the scores 2D plots correlations:
        
        for barehname in ScoreHists[topo]:
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
                cansBH[ckey].cd()
            except:

                for version in BHversions:
                    rholeg_BHvsBG[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_BHvsBG[ckey].SetBorderSize(0)
                    rholeg_BHvsBG[ckey].SetHeader(ckey.split('_')[0])

                    rholeg_bhmuscorevsmuscore[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_bhmuscorevsmuscore[ckey].SetBorderSize(0)
                    rholeg_bhmuscorevsmuscore[ckey].SetHeader(ckey.split('_')[0])

                    
                    rholeg_muscorevsBH[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_muscorevsBH[ckey].SetBorderSize(0)
                    rholeg_muscorevsBH[ckey].SetHeader(ckey.split('_')[0])
                    rholeg_muscorevsBG[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_muscorevsBG[ckey].SetBorderSize(0)
                    rholeg_muscorevsBG[ckey].SetHeader(ckey.split('_')[0])

                    rholeg_bhmuscorevsBH[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_bhmuscorevsBH[ckey].SetBorderSize(0)
                    rholeg_bhmuscorevsBH[ckey].SetHeader(ckey.split('_')[0])
                    rholeg_bhmuscorevsBG[ckey] = ROOT.TLegend(0.14, 0.75, 0.45, 0.88)
                    rholeg_bhmuscorevsBG[ckey].SetBorderSize(0)
                    rholeg_bhmuscorevsBG[ckey].SetHeader(ckey.split('_')[0])

                canname = 'BHcmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansBH[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsBH[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsBH[ckey].SetNColumns(3)
                legsBH[ckey].SetBorderSize(0)
                
                canname = 'BGcmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansBG[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsBG[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsBG[ckey].SetNColumns(3)
                legsBG[ckey].SetBorderSize(0)

                canname = 'mucmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansmu[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsmu[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsmu[ckey].SetNColumns(3)
                legsmu[ckey].SetBorderSize(0)

                canname = 'muscorecmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansmuscore[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsmuscore[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsmuscore[ckey].SetNColumns(3)
                legsmuscore[ckey].SetBorderSize(0)

                canname = 'bhmucmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansbhmu[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsbhmu[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsbhmu[ckey].SetNColumns(3)
                legsbhmu[ckey].SetBorderSize(0)

                canname = 'bhmuscorecmp_model_{}_{}_{}'.format(cantag, barehname, topo)
                cansbhmuscore[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 800, 800)
                legsbhmuscore[ckey] = ROOT.TLegend(0.12, 0.80, 0.89, 0.89)
                legsbhmuscore[ckey].SetNColumns(3)
                legsbhmuscore[ckey].SetBorderSize(0)

                
                canname = 'cmp2d_model_{}_{}_{}'.format(cantag, barehname, topo)
                cans2d[ckey] = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
                cans2d[ckey].Divide(3,2)
                
                rho = sum(rhoVals[topo][barehname]) / len(rhoVals[topo][barehname])
                rhoValsAver[topo][barehname] = 1.*rho
                topot = ROOT.TLatex(0.30, 0.92, '{}'.format(topo) + '         #bar{#rho}=' + '{:+1.2f}'.format(rho))
                topot.SetNDC()
                topot.SetTextSize(0.06)
                topots[ckey] = topot
                opt = baseopt
                
            maxBHScore = -1000
            maxBGScore = -1000
            maxmuScore = -1000
            maxbhmuScore = -1000
            for version in BHversions:
                hh = ScoreHists[topo][barehname].BHscoreHists[version]
                score = hh.GetMean()
                if score > maxBHScore:
                    maxBHScore = score
                hh = ScoreHists[topo][barehname].BGscoreHists[version]
                score = hh.GetMean()
                if score > maxBGScore:
                    maxBGScore = score
                hh = ScoreHists[topo][barehname].muscoreHists[version]
                score = hh.GetMean()
                if score > maxmuScore:
                    maxmuScore = score
                hh = ScoreHists[topo][barehname].bhmuscoreHists[version]
                score = hh.GetMean()
                if score > maxbhmuScore:
                    maxbhmuScore = score

         
            
            for version in BHversions:
                
                hhBH = ScoreHists[topo][barehname].BHscoreHists[version]
                hhBG = ScoreHists[topo][barehname].BGscoreHists[version]
                hhmu = ScoreHists[topo][barehname].muHists[version]
                hhmuscore = ScoreHists[topo][barehname].muscoreHists[version]
                hhbhmu = ScoreHists[topo][barehname].bhmuHists[version]
                hhbhmuscore = ScoreHists[topo][barehname].bhmuscoreHists[version]

                # 2d: histos, now just for axis
                hhmuscorevsBH = ScoreHists[topo][barehname].muscorevsBHHists2d[version]
                hhmuscorevsBG = ScoreHists[topo][barehname].muscorevsBGHists2d[version] 
                hhbhmuscorevsBH = ScoreHists[topo][barehname].bhmuscorevsBHHists2d[version]
                hhbhmuscorevsBG = ScoreHists[topo][barehname].bhmuscorevsBGHists2d[version] 
                hhBHvsBG = ScoreHists[topo][barehname].BHvsBGHists2d[version]
                hhbhmuscorevsmuscore = ScoreHists[topo][barehname].bhmuscorevsmuscoreHists2d[version] 

                # grs: for data storage:)
                grmuscorevsBH = ScoreHists[topo][barehname].muscorevsBHGr2d[version]
                grmuscorevsBG = ScoreHists[topo][barehname].muscorevsBGGr2d[version] 
                grbhmuscorevsBH = ScoreHists[topo][barehname].bhmuscorevsBHGr2d[version]
                grbhmuscorevsBG = ScoreHists[topo][barehname].bhmuscorevsBGGr2d[version] 
                grBHvsBG = ScoreHists[topo][barehname].BHvsBGGr2d[version]
                grbhmuscorevsmuscore = ScoreHists[topo][barehname].bhmuscorevsmuscoreGr2d[version] 
                
                hhs = [hhBH, hhBG, hhmuscore, hhbhmuscore]
                for hh in hhs:
                    #hh.SetMaximum(2*hh.GetMaximum())
                    hh.SetMaximum(0.60*nReplicas)
                hhs = [hhmu, hhbhmu]
                for hh in hhs:
                    #hh.SetMaximum(2*hh.GetMaximum())
                    hh.SetMaximum(0.5*nReplicas)
                    
                cansBH[ckey].cd()
                hhBH.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhBH.GetName(), hhBH.GetEntries(), hhBH.Integral(0, hhBH.GetXaxis().GetNbins()+1)))
                
                cansBG[ckey].cd()
                hhBG.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhBG.GetName(), hhBG.GetEntries(), hhBG.Integral(0, hhBG.GetXaxis().GetNbins()+1)))

                cansmu[ckey].cd()
                hhmu.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhmu.GetName(), hhmu.GetEntries(), hhmu.Integral(0, hhmu.GetXaxis().GetNbins()+1)))

                cansmuscore[ckey].cd()
                hhmuscore.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhmuscore.GetName(), hhmuscore.GetEntries(), hhmuscore.Integral(0, hhmuscore.GetXaxis().GetNbins()+1)))


                cansbhmu[ckey].cd()
                hhbhmu.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhbhmu.GetName(), hhbhmu.GetEntries(), hhbhmu.Integral(0, hhbhmu.GetXaxis().GetNbins()+1)))

                cansbhmuscore[ckey].cd()
                hhbhmuscore.Draw(opt)
                print('Drawing {} with N={} I={}'.format(hhbhmuscore.GetName(), hhbhmuscore.GetEntries(), hhbhmuscore.Integral(0, hhbhmuscore.GetXaxis().GetNbins()+1)))

                
                cans2d[ckey].cd(1)
                hhBHvsBG.Draw(opt)
                grBHvsBG.Draw('P')
                rholeg_BHvsBG[ckey].AddEntry(grBHvsBG, '{} #rho={:1.2f}'.format(version, grBHvsBG.GetCorrelationFactor()), 'P')
                
                cans2d[ckey].cd(2)
                hhmuscorevsBH.Draw(opt)
                grmuscorevsBH.Draw('P')
                rholeg_muscorevsBH[ckey].AddEntry(grmuscorevsBH, '{} #rho={:1.2f}'.format(version, grmuscorevsBH.GetCorrelationFactor()), 'P')
                
                cans2d[ckey].cd(3)
                hhmuscorevsBG.Draw(opt)
                grmuscorevsBG.Draw('P')
                rholeg_muscorevsBG[ckey].AddEntry(grmuscorevsBG, '{} #rho={:1.2f}'.format(version, grmuscorevsBG.GetCorrelationFactor()), 'P')
                
                cans2d[ckey].cd(4)
                hhbhmuscorevsBH.Draw(opt)
                grbhmuscorevsBH.Draw('P')
                rholeg_bhmuscorevsBH[ckey].AddEntry(grbhmuscorevsBH, '{} #rho={:1.2f}'.format(version, grbhmuscorevsBH.GetCorrelationFactor()), 'P')
                
                cans2d[ckey].cd(5)
                hhbhmuscorevsBG.Draw(opt)
                grbhmuscorevsBG.Draw('P')
                rholeg_bhmuscorevsBG[ckey].AddEntry(grbhmuscorevsBG, '{} #rho={:1.2f}'.format(version, grbhmuscorevsBG.GetCorrelationFactor()), 'P')

                cans2d[ckey].cd(6)
                hhbhmuscorevsmuscore.Draw(opt)
                grbhmuscorevsmuscore.Draw('P')
                rholeg_bhmuscorevsmuscore[ckey].AddEntry(grbhmuscorevsmuscore, '{} #rho={:1.2f}'.format(version, grbhmuscorevsmuscore.GetCorrelationFactor()), 'P')
                
                opt = baseopt + sameopt
                
                BHscore = hhBH.GetMean()
                legsBH[ckey].AddEntry(hhBH, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], BHscore, hhBH.GetRMS()), 'F')
                
                BGscore = hhBG.GetMean()
                legsBG[ckey].AddEntry(hhBG, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], BGscore, hhBG.GetRMS()), 'F')

                mu = hhmu.GetMean()
                legsmu[ckey].AddEntry(hhmu, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], mu, hhmu.GetRMS()), 'F')
                
                muscore = hhmuscore.GetMean()
                legsmuscore[ckey].AddEntry(hhmuscore, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], muscore, hhmuscore.GetRMS()), 'F')

                bhmu = hhbhmu.GetMean()
                legsbhmu[ckey].AddEntry(hhbhmu, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], bhmu, hhbhmu.GetRMS()), 'F')
                
                bhmuscore = hhbhmuscore.GetMean()
                legsbhmuscore[ckey].AddEntry(hhbhmuscore, '{}: {:1.1f} #pm {:1.1f} '.format(vtags[version], bhmuscore, hhbhmuscore.GetRMS()), 'F')


                # LaTeX printing:
                
                if version == BHversions[0]:
                    texTableFile.write( ' {} & '.format(topo))
                font = ''
                if abs(BHscore - maxBHScore) < kEpsilon:
                    font = '\\mathbf'
                texTableFile.write( ' {} & $'.format(vtagsTeX[version]) + font + '{' + '{:1.2f}'.format(BHscore) + '} \\pm ' + '{:1.2f}$ '.format(hhBH.GetRMS()))
                font = ''
                if abs(BGscore - maxBGScore) < kEpsilon:
                    font = '\\mathbf'
                texTableFile.write( ' & ' + font + '{' + '{:1.2f}'.format(BGscore) + '} \\pm ' + '{:1.2f}$ '.format(hhBG.GetRMS()))
                font = ''
                if abs(muscore - maxmuScore) < kEpsilon:
                    font = '\\mathbf'
                texTableFile.write( ' & ' + font + '{' + '{:1.2f}'.format(muscore) + '} \\pm ' + '{:1.2f}$ '.format(hhmuscore.GetRMS()))
                font = ''
                if abs(bhmuscore - maxbhmuScore) < kEpsilon:
                    font = '\\mathbf'
                texTableFile.write( ' & ' + font + '{' + '{:1.2f}'.format(bhmuscore) + '} \\pm ' + '{:1.2f}$ '.format(hhbhmuscore.GetRMS()))

                # and finally also mu and bhmu;-) 2.9.2021
                texTableFile.write( ' & ' + font + '{' + '{:1.2f}'.format(mu) + '} \\pm ' + '{:1.2f}$ '.format(hhmu.GetRMS()))
                texTableFile.write( ' & ' + font + '{' + '{:1.2f}'.format(bhmu) + '} \\pm ' + '{:1.2f}$ '.format(hhbhmu.GetRMS()))
                
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
    for ckey in cansBH:
        ptag = ''
        cansToPrint = [cansBH[ckey], cansBG[ckey], cansmu[ckey], cansmuscore[ckey], cansbhmu[ckey], cansbhmuscore[ckey]]
        legsToPrint = [legsBH[ckey], legsBG[ckey], legsmu[ckey], legsmuscore[ckey], legsbhmu[ckey], legsbhmuscore[ckey]]
        for can,leg in zip(cansToPrint, legsToPrint):
            can.cd()
            leg.Draw()
            topots[ckey].Draw()
            can.Print(pngdir + can.GetName() + ptag + '_liny' + '.png')
            can.Print(pdfdir + can.GetName() + ptag + '_liny' + '.pdf')

        for ic in range(0,3*2):
            cans2d[ckey].cd(ic+1)
            topots[ckey].Draw()
        cans2d[ckey].cd(1)
        rholeg_BHvsBG[ckey].Draw()
        cans2d[ckey].cd(2)
        rholeg_muscorevsBH[ckey].Draw()
        cans2d[ckey].cd(3)
        rholeg_muscorevsBG[ckey].Draw()
        cans2d[ckey].cd(4)
        rholeg_bhmuscorevsBH[ckey].Draw()
        cans2d[ckey].cd(5)
        rholeg_bhmuscorevsBG[ckey].Draw()
        cans2d[ckey].cd(6)
        rholeg_bhmuscorevsmuscore[ckey].Draw()

        
        cans2d[ckey].Print(pngdir + cans2d[ckey].GetName() + ptag + '_liny' + '.png')
        cans2d[ckey].Print(pdfdir + cans2d[ckey].GetName() + ptag + '_liny' + '.pdf')

            
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



