#!/usr/bin/python

# jk 11.1.2019

from __future__ import print_function

import os, sys

import ROOT

from Tools import MakeNiceLegendEntry

divideByBinWidth = True

from Tools import DivideByBinWidth

P = ['Top', 'W']
ptags = ['t', 'W']


stuff = []
hists = []

sels = ['Detector', 'Particle']

ROOT.gStyle.SetOptTitle(False)

dirname='TagHistos/'
#fname='analyzed_histos_2tj_allhad_NLO_ATLAS.root'
#fname='analyzed_histos_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS.root'
###fname='analyzed_histos_2t_allhad_14TeV_LO_ATLAS.root'
#fname='analyzed_histos_4t_allhad_14TeV_ATLAS.root'
#fname='analyzed_histos_4t_incl_13TeV_ATLAS.root'
#fname='analyzed_histos_2tj_ljets_NLO_14TeV_ATLAS.root'
#fname='analyzed_histos_2tj_ljets_NLO_14TeV_ptheavy50GeV_ATLAS.root'
#fname='analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS.root'
# fname='analyzed_histos_zp_ttbarj_allhad_allM_14TeV_ATLAS.root'

fnames = [#'analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_all.root',
          #'analyzed_histos_zp_ttbarj_allhad_1250GeV_14TeV_ATLAS_all.root',
          #'analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS_all.root'

    #'analyzed_histos_runs_162_191_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS.root',
    #'analyzed_histos_pp_2tj_allhad_NLO_mergedUnweighted_14TeV_ATLAS_all.root',

    'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root',
    ###'analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS_all.root',
    ###'analyzed_histos_zp_ttbarj_allhad_1250GeV_14TeV_ATLAS_all.root',
    
    #'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV_ATLAS_weighted_all.root',
    #'analyzed_histos_pp_2tj_allhad_NLO_ptj1j2min200_14TeV_ATLAS_weighted_all.root',
    #'analyzed_histos_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV_ATLAS_weighted_all.root',

    # QCD bg:
    'analyzed_histos_pp_2b2j_LO_matched_ALL_14TeV_ATLAS_weighted_all.root'
    
          ]

lsts = [1, 2, 3, 4]
lws = [2, 3, 3, 3]

pngdir='png/'
pdfdir='pdf/' 
os.system('mkdir -p {:}'.format(pngdir,))
os.system('mkdir -p {:}'.format(pdfdir,))

Cans = {}

for sel in sels:
    cans = {}
    for ptag,p in zip(ptags,P):
        canname1 = 'can_tageff_{:}_{:}_spect'.format(p,sel)
        can1 = ROOT.TCanvas(canname1, canname1, 0, 0, 800, 800)
        canname2 = 'can_tageff_{:}_{:}_eff'.format(p,sel)
        can2 = ROOT.TCanvas(canname2, canname2, 100, 1000, 800, 800)
        cans[p] = [can1,can2]
    Cans[sel] = cans

print(Cans)
Legs = {}
Rfiles = []

for sel in sels:
    legs = {}
    for ptag,p in zip(ptags,P):
        leg1 = ROOT.TLegend(0.47, 0.65, 0.88, 0.88)
        leg1.SetBorderSize(0)
        ###leg1.SetHeader('{} {}'.format(sel,p))
        leg2 = ROOT.TLegend(0.17, 0.65, 0.88, 0.88)
        leg2.SetBorderSize(0)
        ###leg2.SetHeader(leg1.GetHeader())
        legs[p] = [leg1,leg2]

        sameopt = ''
        for fname,lst,lw in zip(fnames,lsts,lws):
            print('*** processing {} ***'.format(fname))
            rfile = ROOT.TFile(fname, 'read')
            Rfiles.append(rfile)
            ftag=fname
            ftag=ftag.replace('analyzed_histos_', '').replace('.root', '')
            ltag = MakeNiceLegendEntry(ftag)
        
            h_TagMatched = rfile.Get('{:}{:}Tag{:}truthMatchedLjetsPt{:}'.format(dirname,p,ptag,sel))
            h_TagNonMatched = rfile.Get('{:}{:}Tagnon{:}MatchedLjetsPt{:}'.format(dirname,p,ptag,sel))
            h_notTagMatched = rfile.Get('{:}not{:}Tag{:}truthMatchedLjetsPt{:}'.format(dirname,p,ptag,sel))
            h_notTagNonMatched = rfile.Get('{:}not{:}Tagnon{:}MatchedLjetsPt{:}'.format(dirname,p,ptag,sel))

            h_TagMatched.SetStats(0)
            h_TagMatched.GetXaxis().SetTitle('Large-R jet p_{T} [GeV]')
            h_TagMatched.GetXaxis().SetTitleOffset(1.2)
            h_TagMatched.SetLineColor(ROOT.kRed)        
            h_TagMatched.SetLineWidth(lw)
            h_TagMatched.SetLineStyle(lst)

            h_TagNonMatched.SetStats(0)
            h_TagNonMatched.SetLineColor(ROOT.kBlue)        
            h_TagNonMatched.SetLineWidth(lw)
            h_TagNonMatched.SetLineStyle(lst)

            if divideByBinWidth:
                DivideByBinWidth(h_TagMatched)
                DivideByBinWidth(h_TagNonMatched)
                DivideByBinWidth(h_notTagMatched)
                DivideByBinWidth(h_notTagNonMatched)


            h_sumMatched = h_TagMatched.Clone(h_TagMatched.GetName() + '_all')
            h_sumMatched.Add(h_notTagMatched)
            h_sumMatched.GetYaxis().SetTitle('dN/dp_{T}')
            h_sumMatched.GetYaxis().SetTitleOffset(1.2)
            h_sumMatched.SetStats(0)
            h_sumMatched.SetLineColor(ROOT.kGreen+2)        
            h_sumMatched.SetLineWidth(lw)
            h_sumMatched.SetLineStyle(lst)


            h_sumNonMatched = h_TagNonMatched.Clone(h_TagNonMatched.GetName() + '_all')
            h_sumNonMatched.Add(h_notTagNonMatched)
            h_sumNonMatched.GetYaxis().SetTitle('dN/dp_{T}')
            h_sumNonMatched.GetYaxis().SetTitleOffset(1.2)
            h_sumNonMatched.SetStats(0)
            h_sumNonMatched.SetLineColor(ROOT.kBlack)        
            h_sumNonMatched.SetLineWidth(lw)
            h_sumNonMatched.SetLineStyle(lst)

            h_eff_real = h_TagMatched.Clone(h_TagMatched.GetName() + '_eff_real')
            h_eff_real.Divide(h_sumMatched)

            h_eff_fake = h_TagNonMatched.Clone(h_TagNonMatched.GetName() + '_eff_fake')
            h_eff_fake.Divide(h_sumNonMatched)

            print(sel,p)
            Cans[sel][p][0].cd()
            ROOT.gPad.SetLogy(1)
            h_sumMatched.SetMaximum(h_sumMatched.GetMaximum()*100)
            h_sumMatched.SetMinimum(1.)
            h_sumMatched.Draw('hist' + sameopt)
            h_sumNonMatched.Draw('histsame')
            h_TagMatched.Draw('histsame')
            h_TagNonMatched.Draw('histsame')
            leg1.AddEntry(h_sumMatched, 'Sum {:} matched {}'.format(p,ltag), 'L')
            leg1.AddEntry(h_sumNonMatched, 'Sum {:} non-matched {}'.format(p,ltag), 'L')
            leg1.AddEntry(h_TagMatched, '{:} tagged, matched {}'.format(p,ltag), 'L')
            leg1.AddEntry(h_TagNonMatched, '{:} tagged, non-matched {}'.format(p,ltag), 'L')
            leg1.Draw()
            stuff.append(leg1)
            ROOT.gPad.Update()
            ROOT.gPad.Print(pngdir + Cans[sel][p][0].GetName() + '.png')
            ROOT.gPad.Print(pdfdir + Cans[sel][p][0].GetName() + '.pdf')
            
            
            Cans[sel][p][1].cd()


            h_eff_real.SetMinimum(0.)
            ytitle = '#epsilon_{' + ptag + '}'
            h_eff_real.GetYaxis().SetTitle(ytitle)
            h_eff_real.SetMaximum(1.1)
            if not '2b2j' in fname:
                leg2.AddEntry(h_eff_real, 'Real eff. {}'.format(ltag), 'L')
                h_eff_real.Draw('hist' + sameopt)

            if '2b2j' in fname:
                h_eff_fake.Draw('hist same')
                leg2.AddEntry(h_eff_fake, 'Fake eff. {}'.format(ltag), 'L')


            leg2.Draw()
            stuff.append(leg2)
            ROOT.gPad.Update()
            ROOT.gPad.Print(pngdir + Cans[sel][p][1].GetName() + '.png')
            ROOT.gPad.Print(pdfdir + Cans[sel][p][1].GetName() + '.pdf')

            hists.append([h_eff_fake, h_eff_real, h_TagNonMatched, h_TagMatched, h_sumMatched, h_sumNonMatched])

            #ROOT.gPad.SetGridy(1)
            sameopt = 'same'

    Legs[sel] = legs

ROOT.gApplication.Run()
