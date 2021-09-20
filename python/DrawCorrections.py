#!/usr/bin/python

# jiri kvita, 2015, 2019, 4.6.2020

from __future__ import print_function

import ROOT
from Tools import *

import os, sys
from array import array

kEff = 0
kAcc = 2

_cans = []
_files = []
_corrs = []
_Corrs = []
_legs = []

Fst = [3305, 3335, 3353]
# TDirectory paths
Paths = [ '0B2S',
          '1B1S',
          '2B0S',
]


# Col = [, ROOT.kGreen+2, ROOT.kAzure+2 ]
Cols = { kEff : [ROOT.kOrange-3, ROOT.kRed,  ROOT.kRed+2],
         kAcc : [ROOT.kAzure+9,  ROOT.kBlue, ROOT.kMagenta+3,]
}


ObjNames = {
             'LJet1' : 't1',
             'LJet2' : 't2',
             'Dijet' : 't#bar{t}',
             'difference' : ''}
TitleNames = { 'pt' : [  'p_{T}', '[GeV]' ], 
               'm' : [  'm', '[GeV]' ], 
               'absrap' : [  '|y|', '' ], 
               'rapidity' : [  'y', '' ],
               'Pout' : [  'p_{out}', '[GeV]' ],
               'absPout' : [  '|p_{out}|', '[GeV]' ],
               'z_ttbar' : [  'z_{#hat{t}_{l}#hat{t}_{h}}', '' ],
               'Yboost' : [  'y_{boost}', '' ],
               'Chi_ttbar' : [  '#chi_{#hat{t}_{l}#hat{t}_{h}}', '' ],
               'dPhi_ttbar' : [  '#Delta#phi_{#hat{t}_{l}#hat{t}_{h}}', '' ],
               'Salam_ttbar' : [  'S_{#hat{t}_{l}#hat{t}_{h}}', '' ],
               'HT_ttbar' : [  'H_{T}^{t#bar{t}}', '[GeV]' ],
               'HT_pseudo' : [  'H_{T}^{pseudo}', '[GeV]' ],
               'R_lb' : [  '[p_{T}^{j1} + p_{T}^{j2}] / [ p_{T}^{b,lep} + p_{T}^{b,had}]', '' ],
               'R_Wb_had' : [  'p_{T}^{W,had} / p_{T}^{b,had}', '' ],
               'R_Wb_lep' : [  'p_{T}^{W,lep} / p_{T}^{b,lep}', '' ],
               'R_Wt_had' : [  'p_{T}^{W,had} / p_{T}^{t,had}', '' ],
               'R_Wt_lep' : [  'p_{T}^{W,lep} / p_{T}^{t,lep}', '' ],
               }
CorrNames = {
    #'eff' : 'Efficiency #varepsilon', 
    ##'match' : 'Matching correction f_{match}', 
    #'acc' : 'Acceptance correction f_{acc}'
    'eff' : 'Efficiency', 
    #'match' : 'Matching correction f_{match}', 
    'acc' : 'Acceptance correction'
}


#################
def ZeroErrorBars(corr):
    for i in range(0,corr.GetN()):
        corr.SetPointError(i,0., 0., 0., 0.)


#################
def CheckAcc(acc,name):
    #vals = array('d', [0.,0.])
    x = Double(0.) # ROOT.Double
    eff = Double(0.)
    for i in range(0,acc.GetN()):
        #acc.GetPoint(i,vals[0],vals[1])
        acc.GetPoint(i,x,eff)
        if eff > 1.:
            print('  ERROR: acceptance={:4.2f} in bin {:} of {:}!'.format(eff,i,name) )
    return

def GetTag(objname, varname):
    tag = ''
    #if objname.find('top') >= 0 and varname != "absrap": tag = '_0'
    #if objname.find('tt') >= 0 and varname == "pt": tag = '_5'
    return tag



#################
def GetCorrection(rfile, pfile, objname, varname, icorr, basepath, migrabasepath = 'migrations'):
    # here access the particle spectrum
    h_part = pfile.Get(basepath + '/Particle' + objname + varname)

    # get the migration matrix:
    matrixPath =  basepath + '/' + migrabasepath + '/' + objname + varname + '_Particle_Detector'
    print('Trying to get {:}'.format(matrixPath) )
    h_matrix = rfile.Get( matrixPath )
    print('Trying to get {:}'.format(basepath + '/Detector' + objname + varname,) )
    h_reco = rfile.Get(basepath + '/Detector' + objname + varname)

    # to fix in future with unified bins:
    # old:
    #h_recopart_r = h_matrix.ProjectionY( "particle_recoandparticle", 1, h_matrix.GetNbinsY() )
    #h_recopart_p = h_matrix.ProjectionX( "particle_recoandparticle", 1, h_matrix.GetNbinsX() )
    h_matrix.ClearUnderflowAndOverflow()
    h_matrix.GetXaxis().SetRange(1, h_matrix.GetXaxis().GetNbins() )
    h_matrix.GetYaxis().SetRange(1, h_matrix.GetYaxis().GetNbins() )
    h_recopart_r = h_matrix.ProjectionY( "reco_recoandparticleX" )#, 1, h_matrix.GetNbinsX() )
    h_recopart_p = h_matrix.ProjectionX( "particle_recoandparticle") #, 1, h_matrix.GetNbinsX() )

    # wrong solution with bad over/over flow bins treatment:
    #h_recopart_r = rfile.Get(basepath + '/DetectorAndParticle_det' + objname + varname)
    #h_recopart_p = rfile.Get(basepath + '/DetectorAndParticle_ptcl' + objname + varname)

    print('h_part:')
    PrintBinContent(h_part)
    print('h_reco:')
    PrintBinContent(h_reco)
    print('h_recopart_p:')
    PrintBinContent(h_recopart_p)
    print('h_recopart_r:')
    PrintBinContent(h_recopart_r)

    print('  Making eff...')
    #eff = MakeRatio( h_recopart_p, h_part)
    eff = h_recopart_p.Clone(h_recopart_p.GetName() + '_eff')
    eff.Divide(h_part)

    print('  Making acc...')
    #acc = MakeRatio( h_recopart_r, h_reco)
    acc = h_recopart_r.Clone(h_recopart_r.GetName() + '_acc')
    acc.Divide(h_reco)
    
    #CheckAcc(acc,'{:} {:}'.format(h_recopart_r.GetName(),h_recopart_r.GetTitle()) )

    #print('  Making match...'}
    #print('    RMS check: {:} {:}'.format(h_match_r.GetRMS(), h_rp.GetRMS())}
    #match = MakeRatio( h_match_r,  h_rp)

    if icorr == kEff: return eff#,h_part#,h_match_p
    if icorr == kAcc: return acc#,h_rp,h_reco
    #if icorr == 1: return match,h_match_r,h_rp
    return

#################
def DrawCorrection(ll, rfile, pfile, objname, varname, icorr, tdpaths):
    print('DrawCorrection icor requested: {:}'.format(icorr,))

    tag = ''
    if icorr == kEff: tag = 'eff'
    if icorr == kAcc: tag = 'acc'
    #if icorr == 1: tag = 'match'

    canname = '{:}_{:}_{:}_{:}'.format(tag,objname,varname,ll)
    can = ROOT.TCanvas(canname, canname, 1, 1, 800, 800)
    _cans.append(can)
    can.cd()

    corrs = []
    
    opt = 'e3 x0'
    count = 0
    xmin = -1
    xmax = -1
    #for rfile,pfile in zip(rfiles,pfiles):
    leg = ROOT.TLegend(0.70, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    _legs.append(leg)
    ipath = -1
    for basepath in tdpaths:
        ipath = ipath + 1
        print('*** Processing {:} ***'.format(basepath))
        corr = GetCorrection(rfile, pfile, objname, varname, icorr, basepath)
        corrs.append(corr)
        print('Got correction:')
        print(corr)
        col = Cols[icorr][ipath]
        corr.SetLineColor(col)
        corr.SetLineStyle(1)
        corr.SetFillStyle(Fst[ipath])
        corr.SetFillColor(col)
        corr.SetLineWidth(2)
        #if count > 0:
            #corr.SetLineWidth(3)
            #corr.SetLineStyle(count)
            #corr.SetMarkerStyle(0)
            #corr.SetMarkerSize(0)
            #corr.SetMarkerStyle(0)
            #corr.SetMarkerColor(0)
            #ZeroErrorBars(corr) # worked only for graphs, now working with hists;-)
        # else
        corr.SetMarkerSize(1.5)
        corr.SetMarkerStyle(20 + count)
        corr.SetMarkerColor(col)


        if count == 0:
            leg.AddEntry(corr, basepath, 'P') # PLF
        else:
            leg.AddEntry(corr, basepath, 'P') # PLF
        
        if count == 0:
            xmin = corr.GetXaxis().GetXmin()
            xmax = corr.GetXaxis().GetXmax() 

        xtitle = corr.GetXaxis().GetTitle().replace('DetectorAndParticle_det','').replace('DetectorAndParticle_ptcl','')
        ytitle = CorrNames[tag]
        title = CorrNames[tag] + ';' + xtitle# + ';' + ytitle

        SetStyle(corr, xtitle, ytitle)
        corr.SetStats(0)
        corr.SetMinimum(0)
        corr.SetMaximum(1.1)
        corr.Draw(opt)
        corr.GetYaxis().SetTitleOffset(1)
        corr.GetXaxis().SetTitleOffset(1)
    
        _corrs.append(corr)
        opt = 'e3 x0 same'
        count = count+1
        
    leg.Draw()
    yy = 0.865
    yyoff = 0.80-0.75
    #myText(0.19, yy-yyoff, kBlack, "Allhad");

    can.Print('png/' + canname + '.png')
    can.Print('pdf/' + canname + '.pdf')
    #can.Print('C/' + canname + '.C' )
    return corrs

####################################################
####################################################
####################################################

ROOT.gStyle.SetOptTitle(0)


ptag=''
ftag=''

ppath='analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root'
#analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root'

#ppath='analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_all.root'
#ppath='analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'
#ppath='analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half1.root'
rpath=ppath

os.system('mkdir C png pdf')
ROOT.gROOT.SetBatch(0)


rfile = ROOT.TFile(rpath, 'read')
_files.append(rfile)
print('Opened file {:}'.format(rfile.GetName(),))

pfile = ROOT.TFile(ppath, 'read')
_files.append(pfile)
print('Opened file {:}'.format(pfile.GetName(),))



Obj = [    #'Top1',
           #'Top2',
           'Top',
           'DiTop',
]
Vars = ['Pt',
        'Mass', 
        'Rapidity',
        'Pout',
        'CosThetaStar',
        'DeltaPhi',
        ##'Yboost',
        ##'Chittbar',
        ##'Delta'
]


ll = 'allhad'

for obj in Obj:
    for var in Vars:
        if obj != 'DiTop' and (var != 'Pt' and var != 'Rapidity'):
            continue
        effs = DrawCorrection(ll, rfile, pfile, obj, var, kEff, Paths)
        #DrawCorrection(ll, rfile, pfile, obj, var, 1, Paths)
        accs = DrawCorrection(ll, rfile, pfile, obj, var, kAcc, Paths)
        _Corrs.append([effs,accs])
        pass

# not mainted for the moment:
#SpecObj = [] #['difference' ]
#SpecVars = [
#           'HTj', 'HTJ',
#           #'R_lb', 
#           #'R_Wb_lep', 
#           #'R_Wb_had',
#           #'R_Wt_lep', 
#           #'R_Wt_had' 
#]
#for obj in SpecObj:
#    for var in SpecVars:
#        DrawCorrection(ll, rfiles, pfiles, obj, var, 0, tdpath)
#        DrawCorrection(ll, rfiles, pfiles, obj, var, 1, tdpath)
#        #DrawCorrection(ll, rfiles, pfiles, obj, var, 2, tdpath)
#        pass


ROOT.gApplication.Run()

