#!/usr/bin/python

# Jiri Kvita, 14.4.2020
# file to store xsect help and stack properties classes and xsect data

from __future__ import print_function

from collections import OrderedDict

from cXsectData import *

import ROOT

gStackName = 'STACK'
# 1.7.2020
# IMPORTANT: cut flow bin to take the Pythia weighted events from, for the cross-section normalization
kWeightedSumEventsCutFlowBin = 2
kCutFlowHname = 'CutFlowDet'

#########################################

def FindXsectData(name):
    for xsd in xsdata:
        if xsd.name == name:
            return xsd
    print('ERROR getting xsection data for {}!'.format(name))
    return 0

#########################################
# better this way:
def MakeXsextDic():
    xsdict = {}
    for xsd in xsdata:
        xsdict[xsd.name] = xsd
    return xsdict


#########################################
class cLegItem:
    def __init__(self, hist, legtag, lopt):
        self.hist = hist
        self.legtag = legtag
        self.lopt = lopt


#########################################

class cTopPlot:
    def __init__(self, stack, stackStatIndep, leg, legh, legitems, pads, spads, fullhname, basename = ''):

        # for bg_sig stack:
        self.stack = stack
        # for data:
        self.stackStatIndep = stackStatIndep
        
        # we expect the last stack item to be the signal, all preceeding to add up to background
        self.hsig = stackStatIndep.GetHists().At(stackStatIndep.GetNhists()-1)
        self.hdata = self.hsig.Clone(self.hsig.GetName() + '_tot')
        for ih in range(0, stackStatIndep.GetNhists()-1):
            self.hdata.Add(stackStatIndep.GetHists().At(ih))
        # but we want the signal model to be stat indepent, so let's switch it:
        self.hsig = stack.GetHists().At(stack.GetNhists()-1)
        self.hbg = self.hsig.Clone(self.hsig.GetName() + '_bg')
        self.hbg.Reset()
        for ih in range(0, stack.GetNhists()-1):
            self.hbg.Add(stack.GetHists().At(ih))
        self.leg = leg
        self.legh = legh
        self.legitems = legitems
        # stack pads:
        self.pads = pads
        # bg subtracted = signal only pads
        self.spads = spads
        self.fullhname = fullhname
        self.basename = basename
        return

#########################################

class cSampleProp:
    # dopt = draw option, lopt = legend option nev = [run1evt, run2evt, ...]
    def __init__(self, mcol, mst, msz, lcol, lst, lw, fcol, fst,
                 xs, filetagForNevts, lumi, sf, legtag, dopt, lopt):
        self.mcol = mcol
        self.mst = mst
        self.msz = msz
        self.lcol = lcol
        self.lst = lst
        self.lw = lw
        self.fcol = fcol
        self.fst = fst
        self.xs = xs
        #self.nev = nev
        #if len(nev) > 0:

        self.nev = 1.
        fname = 'analyzed_histos_{}_14TeV_ATLAS_all.root'.format(filetagForNevts)
        try:
            tmpfile = ROOT.TFile(fname, 'read')
            cutflowh = tmpfile.Get(kCutFlowHname)
            self.nev = cutflowh.GetBinContent(kWeightedSumEventsCutFlowBin) / 2. # here the division by 2 assumes ~stat equivalent half0 and half1 files!
            print('cSampleProp: got the total weighted nEvets from cutflow histo for {} as {}'.format(fname,self.nev))
        except:
            if gStackName in fname:
                print('OK, not computing xsect weights as here we define only colors and styles for the total stack;-)')
            else:
                print('ERROR getting the number of weighted events from file {}'.format(fname))
            
        self.lumi = lumi
        self.sf = sf
        self.legtag = legtag
        self.dopt = dopt
        self.lopt = lopt
        print('Weighting ingredients for {}: xs={} sumN={} lumi={} sf={}'.format(fname, self.xs, self.nev, self.lumi, self.sf)) # legtag
        self.weight = self.xs / self.nev * self.lumi * self.sf
        print('Initialized process {:30} with weight {}'.format(legtag, self.weight))
        return

#########################################
    
def MakeStackItems(signalFile = 'analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root',
                   legtag = 'y_{0} #rightarrow t#bar{t}, m_{y0} = 1000 GeV, #Gamma_{y0} = 10 GeV',
                   lumi = 1., csf = 1., ssf = 1): # csf ... common SF, ssf = signal SF


    stackItems = OrderedDict() # jk 18.6.2020
    stag = signalFile + ''
    stag = stag.replace('analyzed_histos_', '').replace('_14TeV_ATLAS_half0.root', '')
    xsdict = MakeXsextDic()


    # QCD bbjj samples:
    # 1.7.2021
    stackItems['analyzed_histos_pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGray+2, 1001,
                                                                                                                   xsdict['pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200'].xsect,
                                                                                                                   'pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200',
                                                                                                                   lumi, csf, 'bbjj, p_{T}^{j1,j2} #in (60,200) GeV',
                                                                                                                   'F', 'F')
    
    stackItems['analyzed_histos_pp_2b2j_LO_matched_ptj1j2min200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGray+1, 1001,
                                                                                                       xsdict['pp_2b2j_LO_matched_ptj1j2min200'].xsect,
                                                                                                       'pp_2b2j_LO_matched_ptj1j2min200',
                                                                                                       lumi, csf, 'bbjj, p_{T}^{j1,j2} #geq 200 GeV',
                                                                                                       'F', 'F')
    
    stackItems['analyzed_histos_pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGray, 1001,
                                                                                                                     xsdict['pp_2b2j_LO_matched_ptj1min200_ptj2min60max200'].xsect,
                                                                                                                     'pp_2b2j_LO_matched_ptj1min200_ptj2min60max200',
                                                                                                                     lumi, csf, 'bbjj, p_{T}^{j1} #geq 200 GeV p_{T}^{j2} #in (60,200) GeV',
                                                                                                                     'F', 'F')



    
     # OLD Wbb and Wj samples
    """
    stackItems['analyzed_histos_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGreen+2, 1001,
                                                                                                  xsdict['pp2Wbb_pTj1min50GeV_allhad_14TeV'].xsect,
                                                                                                  'pp2Wbb_pTj1min50GeV_allhad_14TeV',
                                                                                                  lumi, csf, 'Wbb p_{T}^{j1} > 50 GeV', 'F', 'F')
    stackItems['analyzed_histos_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGreen+1, 1001,
                                                                                                 xsdict['pp2Wj_pTj1min50GeV_allhad_14TeV'].xsect,
                                                                                                 'pp2Wj_pTj1min50GeV_allhad_14TeV',
                                                                                                 lumi, csf, 'Wj p_{T}^{j1} > 50 GeV', 'F', 'F')
    """

    # Wbb
    stackItems['analyzed_histos_pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGreen+3, 1001,
                                                                                                xsdict['pp2wbbjjMatched_allhad_pTj1j2min_60GeV'].xsect,
                                                                                                'pp2wbbjjMatched_allhad_pTj1j2min_60GeV',
                                                                                                lumi, csf, 'Wbbjj, p_{T}^{j1,j2} > 60 GeV', 'F', 'F')



    # WWbb
    stackItems['analyzed_histos_pp2wwbb_allhad_pTj1j2min_60GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kGreen+2, 1001,
                                                                                                      xsdict['pp2wwbb_allhad_pTj1j2min_60GeV'].xsect,
                                                                                                      'pp2wwbb_allhad_pTj1j2min_60GeV',
                                                                                                      lumi, csf, 'WWbb, p_{T}^{j1,j2} > 60 GeV', 'F', 'F')


    
    # QCD samples:
    """
    stackItems['analyzed_histos_pp2jjjMatched_pTj1j2min_60GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kMagenta, 1001,
                                                                                                     xsdict['pp2jjjMatched_pTj1j2min_60GeV'].xsect,
                                                                                                     'pp2jjjMatched_pTj1j2min_60GeV',
                                                                                                     lumi, csf, 'Matched 3j p_{T}^{j1,j2} > 60 GeV', 'F', 'F')
    """

    # best to use with more stats...
    """
    
    stackItems['analyzed_histos_pp2jjjjMatched_pTj1j2min_60GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kMagenta+1, 1001,
                                                                                                      xsdict['pp2jjjjMatched_pTj1j2min_60GeV'].xsect,
                                                                                                      'pp2jjjjMatched_pTj1j2min_60GeV',
                                                                                                      lumi, csf, 'Matched 4j p_{T}^{j1,j2} > 60 GeV', 'F', 'F')  
    
    stackItems['analyzed_histos_pp2jjjj_pTj1min400GeV_pTj2min60GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kMagenta+2, 1001,
                                                                                                          xsdict['pp2jjjj_pTj1min400GeV_pTj2min60GeV'].xsect,
                                                                                                          'pp2jjjj_pTj1min400GeV_pTj2min60GeV',
                                                                                                          lumi, csf, 'Matched 4j p_{T}^{j1(j2)} > 400 (60) GeV', 'F', 'F')


    """

    # BAD: stackItems['analyzed_histos_pp2jjjj_mjjmin400GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kMagenta, 1001,
    #                                                                                        xsdict['pp2jjjj_mjjmin400GeV_14TeV'].xsect,
    #                                                                                        'pp2jjjj_mjjmin400GeV_14TeV',
    #                                                                                        lumi, csf, 'Matched 4j m_{jj} > 400', 'F', 'F')

    # OLD tt:
    
    """
    OLD tt samples:
    stackItems['analyzed_histos_2t_allhad_14TeV_LO_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+6, 1001,
                                                                                    xsdict['pp_2t_allhad_14TeV_LO'].xsect,
                                                                                    'pp_2t_allhad_14TeV_LO',
                                                                                    lumi, csf, 't#bar{t} LO', 'F', 'F')
    stackItems['analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+4, 1001,
                                                                                                   xsdict['pp_2tj_allhad_NLO_ptheavy50GeV'].xsect,
                                                                                                   'pp_2tj_allhad_NLO_ptheavy50GeV',
                                                                                                   lumi, csf, 't#bar{t} NLO II', 'F', 'F')
    stackItems['analyzed_histos_2tj_allhad_NLO_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+2, 1001,
                                                                                xsdict['pp_2tj_allhad_NLO'].xsect,
                                                                                'pp_2tj_allhad_NLO',
                                                                                lumi, csf, 't#bar{t} NLO I', 'F', 'F')
    """



 
    # QCD light flavoured samples, 30.6.2021
    #stackItems['analyzed_histos_pp_4j_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kOrange, 1001,
    #                                                                                                    xsdict['pp_4j_ptj1j2min60_ptj1j2max200'].xsect,
    #                                                                                                    'pp_4j_ptj1j2min60_ptj1j2max200',
    #                                                                                                    lumi, csf, 'QCD, 60 GeV #leq p_{T}^{j1,j2} #leq 200 GeV',
    #                                                                                                    'F', 'F')
    #stackItems['analyzed_histos_pp_4j_ptj1j2min200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kOrange-3, 1001,
    #                                                                                                    xsdict['pp_4j_ptj1j2min200'].xsect,
    #                                                                                                    'pp_4j_ptj1j2min200',
    #                                                                                                    lumi, csf, 'QCD, p_{T}^{j1,j2} #geq 200 GeV',
    #                                                                                                    'F', 'F')
    #stackItems['analyzed_histos_pp_4j_ptj1min200_ptj2min60max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kOrange+7, 1001,
    #                                                                                                    xsdict['pp_4j_ptj1min200_ptj2min60max200'].xsect,
    #                                                                                                    'pp_4j_ptj1min200_ptj2min60max200',
    #                                                                                                    lumi, csf, 'QCD, p_{T}^{j1} #geq 200 GeV p_{T}^{j2} #in (60,200) GeV',
    #                                                                                                    'F', 'F')
    
   

    # FOR ADDING MORE:
    #stackItems['analyzed_histos__14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kOrange+7, 1001,
    #                                                                    xsdict[''].xsect,
    #                                                                    '',
    #                                                                    lumi, csf, '',
    #                                                                    'F', 'F')



    
    # ttbar base semiboosted sample
    #stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure, 1001,
    #                                                                                                         xsdict['pp_2tj_allhad_NLO_ptj1min60_ptj2min60'].xsect,
    #                                                                                                         'pp_2tj_allhad_NLO_ptj1min60_ptj2min60',
    #                                                                                                         lumi, csf, 't#bar{t}, p_{T}^{j1,j2} #geq 60 GeV',
    #                                                                                                         'F', 'F')
    # semiboosted ttbar sample, 21.6.2020
    #stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1min100_ptj2min100_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+5, 1001,
    #                                                                                                           xsdict['pp_2tj_allhad_NLO_ptj1min100_ptj2min100'].xsect,
    #                                                                                                           'pp_2tj_allhad_NLO_ptj1min100_ptj2min100',
    #                                                                                                           lumi, csf, 't#bar{t}, p_{T}^{j1,j2} #geq 100 GeV',
    #                                                                                                           'F', 'F')
    
                                                                                                             
    # boosted ttbar sample, 21.6.2020
    #stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+10, 1001,
    #                                                                                                           xsdict['pp_2tj_allhad_NLO_ptj1min200_ptj2min200'].xsect,
    #                                                                                                           'pp_2tj_allhad_NLO_ptj1min200_ptj2min200',
    #                                                                                                           lumi, csf, 't#bar{t}, p_{T}^{j1,j2} #geq 200 GeV',
    #                                                                                                           'F', 'F')


    # main ttbar sample (a), 23.10.2020
    stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure-6, 1001,
                                                                                                                  xsdict['pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200'].xsect,
                                                                                                                  'pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200',
                                                                                                                  lumi, csf, 't#bar{t}, p_{T}^{j1,j2} #in (60, 200) GeV',
                                                                                                                  'F', 'F')
    # main ttbar sample (b), 23.10.2020
    stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure-1, 1001,
                                                                                                      xsdict['pp_2tj_allhad_NLO_ptj1j2min200'].xsect,
                                                                                                      'pp_2tj_allhad_NLO_ptj1j2min200',
                                                                                                      lumi, csf, 't#bar{t}, p_{T}^{j1,j2} #geq 200 GeV',
                                                                                                      'F', 'F')
    # main ttbar sample (c), 23.10.2020
    stackItems['analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_half0.root'] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kAzure+7, 1001,
                                                                                                                    xsdict['pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200'].xsect,
                                                                                                                    'pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200',
                                                                                                                    lumi, csf, 't#bar{t}, p_{T}^{j1} #geq 200 GeV p_{T}^{j2} #in (60,200) GeV',
                                                                                                                    'F', 'F')


 
    

    
    # signal file plotting properties
    stackItems[signalFile] = cSampleProp(0, 0, 0, 0, 0, 0, ROOT.kRed+1, 1001,
                                         xsdict[stag].xsect,
                                         stag,
                                         lumi, ssf, legtag, 'F', 'F')

    # stack plotting properties
    stackItems[gStackName] = cSampleProp(ROOT.kBlack, 20, 1, ROOT.kBlack, 1, 1, 0, 0,
                                         1.,
                                         gStackName,
                                         1., 1., 'Pseudo-data',
                                         'e1 X0', 'P')

    return stackItems

#########################################
