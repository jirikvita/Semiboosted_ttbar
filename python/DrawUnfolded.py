#!/usr/bin/python

from __future__ import print_function

import os

# jk 2.2.2018, 26.2.2018, 30th March 2018, June 2020

localstuff = []
fits = []

from ReadTools import *

from matplotlib import pyplot as plt

###################################################################

def DrawUnfolded(filenames, migraname = '', truthlevel = 'Particle', xlabel = '', rebin = -1,
                 Normalize = False, addTag = '',divideByBinWidth = True,):

    ROOT.gStyle.SetOptTitle(0)

    filename = filenames[0]

    print('Some getting...')
    unffile = ROOT.TFile(filename, 'read')
    h_data = unffile.Get(migraname + '_projY_smeared')
    h_ptcl = unffile.Get(migraname + '_projX')
    h_reco = unffile.Get(migraname + '_projY')
    h_unfolded = unffile.Get('h_unfolded')
    h_unfolded_gfit = h_unfolded.Clone('h_unfolded_gfit')
    h_unfolded_gfit.Reset()
    h_migra = unffile.Get(migraname + '_migra')
    # response matrix: BoostedTopPt_ptcl_det

    #print(unffile)
    #print(h_data)
    #print(h_ptcl)
    #print(h_reco)
    #print(h_unfolded)
    #print(h_migra)

    h_data.SetStats(0)
    h_ptcl.SetStats(0)
    h_reco.SetStats(0)
    h_unfolded.SetStats(0)
    h_migra.SetStats(0)

    if divideByBinWidth:
        DivideByBinWidth(h_ptcl)
        DivideByBinWidth(h_reco)
        DivideByBinWidth(h_data)
        DivideByBinWidth(h_unfolded)
    
    # add the unfolded data to the first comprison plot
    canname = filename.replace('unfolded_', 'Ingredients').replace('.root', '')
    cantitle = canname
    if '/' in canname:
        cantags = canname.split('/')
        canname = cantags[-1]
        cantitle = ''
        for cantag in cantags[:-2]:
            cantitle = cantitle + '_' + cantag
            
    origcanname = canname
    
    can,leg = PlotIngredients(h_data, h_reco, h_ptcl, h_migra, canname, cantitle, truthlevel, 'pdf/unfolded/')

    localstuff.append(leg)
    can.cd(1)
    unfcol = ROOT.kBlack # ROOT.kGreen+2
    h_unfolded.SetLineColor(unfcol)
    h_unfolded.SetLineStyle(1)
    h_unfolded.SetLineWidth(1)
    h_unfolded.SetMarkerColor(unfcol)
    h_unfolded.SetMarkerStyle(20)
    h_unfolded.SetMarkerSize(2)
    h_unfolded.SetStats(0)
    #h_unfolded_copy = h_unfolded.DrawCopy('hist same')
    #localstuff.append(h_unfolded_copy)
    #Legs[-1].AddEntry(h_unfolded, 'Unfolded', 'PL')
    Legs[-1].Draw()

    #can_unf = PlotHisto(h_unfolded, 'hist', 'Unfolded')
    #can_unf.SetLogy(1)

    #Hists = [h_ptcl, h_unfolded, h_reco, h_data]
    #Labels = [truthlevel, 'Unfolded', 'Detector', 'Pseudo data']
    Hists = [h_ptcl, h_unfolded,  h_data]
    Labels = [truthlevel, 'Unfolded', 'Detector level']
    #Hists = [h_ptcl, h_unfolded]
    #Labels = ['Particle', 'Unfolded']

    idiv = 0
    canname = canname.replace('Ingredients', 'Closure')
    can_div, ratios, widgets = DrawRatios(Hists, Labels, idiv, Normalize, canname, xlabel, addTag, 'pdf/unfolded/')

    if Normalize:
        h_reco.Scale(1./h_reco.Integral())
        
    #ratios[0].DrawCopy('PLsame')
    
    #canname = 'UnfoRatio'
    #can_ur = ROOT.TCanvas(canname, canname, 500, 200, 800, 800)
    #can_ur.cd()
    #ratios[0].DrawCopy('PL')  

    #Hists = [h_ptcl, h_reco, h_data]
    #Labels = ['Particle', 'Reco', 'Data']
    #idiv = 0
    #can_div = DrawRatios(Hists, Labels, idiv)

    #PrintCan(can)
    #PrintCan(can_unf)
    #PrintCan(can_div)

    canname = 'Posteriors_' + origcanname
    postcan = ROOT.TCanvas(canname, canname, 0, 0, 1000,1000)

    nn = h_data.GetNbinsX()
    nx = int(sqrt(1.*nn))
    ny = nx
    while nx*ny < nn:
        ny = ny+1
    postcan.Divide(nx,ny)
    posts = []
    for i in range(0, h_data.GetNbinsX()):
        post = unffile.Get('trace_{}'.format(i))
        postcan.cd(i+1)
        gfit = ROOT.TF1('fit_{}'.format(i), "[0]*exp(-(x-[1])^2 / (2*[2]^2))", post.GetXaxis().GetXmin(),  post.GetXaxis().GetXmax())
        gfit.SetParameters(post.GetMaximum()/2, post.GetMean(), post.GetRMS())
        post.Draw('hist')
        post.Fit(gfit)
        gfit.Draw('same')
        fits.append(gfit)
        posts.append(post)
        h_unfolded_gfit.SetBinContent(i+1, gfit.GetParameter(1))
        h_unfolded_gfit.SetBinError(i+1, gfit.GetParameter(2))
    localstuff.append([posts, postcan])
    h_unfolded_gfit.Scale(1.)

    postcan.Print('png/unfolded/' + canname + '.png')
    postcan.Print('pdf/unfolded/' + canname + '.pdf')

    # using gauss fit:
    #GHists = [h_ptcl, h_unfolded_gfit,  h_data]
    #can_div, ratios, widgets = DrawRatios(GHists, Labels, idiv, Normalize, canname, xlabel, addTag)
    
    
    # TO REVISE!!!
    if len(filenames) > 1:
        # OK, we assume root files to make corrections from
        # are also passed and a full closure including the acc and eff corrs
        # is supposed to be performed

        # so far, this is still unfolding the projection of migration matrix
        # onto the reco axis, i.e. after the acc correction, so we do not correct for acc
        # acc =

        truthfile_forCorrs = ROOT.TFile(filenames[1], 'read')
        truthfile_forClosure = ROOT.TFile(filenames[2], 'read')
        
        varname = migraname.split('_')[0]
        vartruth = varname + '_' + truthlevel
        truthhname = 'AllCuts/'
        if truthlevel == 'Parton':
            vartruth = truthlevel + '_' + varname
            truthhname = 'NoCuts/'
        truthhname = truthhname + vartruth
        # full phase space parton or particle level histo:
        print('Getting {:} from {:} '.format(truthhname, truthfile_forCorrs.GetName()))
        # truth spectrum at full phase space:
        h_truth_forCorrs = truthfile_forCorrs.Get(truthhname)

        truthhname_det = 'AllCuts/'
        if truthlevel == 'Parton':
            # parton spectrum after detector cuts:
            vartruth = truthlevel + '_det_' + varname
        truthhname_det = truthhname_det + vartruth
        print('Getting {:} from {:} '.format(truthhname_det, truthfile_forCorrs.GetName()))
        h_ptcl_forCorrs = truthfile_forCorrs.Get(truthhname_det)

        if rebin > 1:
            h_truth_forCorrs.Rebin(rebin)
            h_ptcl_forCorrs.Rebin(rebin)

        print('Making eff...')
        effcorr = MyMultiply(h_truth_forCorrs, h_ptcl_forCorrs, -1, '_EffCorr')
        eff = MyMultiply(h_ptcl_forCorrs, h_truth_forCorrs, -1, '_Eff')
        canname = origcanname.replace('Ingredients', 'eff_')
        can_eff = ROOT.TCanvas(canname, canname)
        eff.SetMinimum(0.)
        eff.SetMaximum(0.2)
        eff.Draw()
        PrintCan(can_eff)
        localstuff.append(can_eff)
        localstuff.append(effcorr)
        localstuff.append(eff)
        
        h_truth_forClosure = truthfile_forClosure.Get(truthhname)
        if rebin > 1:
            h_truth_forClosure.Rebin(rebin)
        print ('Making unfolded...')
        h_unfolded_corr = MyMultiply(h_unfolded, effcorr, 1, '_EffCorr')
        
        Hists = [h_truth_forClosure, h_unfolded_corr] #,  h_data]
        Labels = [truthlevel, 'Unfolded/eff', 'Detector level']
        #Hists = [h_ptcl, h_unfolded]
        #Labels = ['Particle', 'Unfolded']
        
        idiv = 0
        canname = origcanname.replace('Ingredients', 'ClosureInclEff')
        can_div_eff, ratios_eff = DrawRatios(Hists, Labels, idiv, Normalize, canname, xlabel, addTag)

        # PrintCan(can_div)

        localstuff.append(h_unfolded_corr)
        localstuff.append(h_truth_forClosure)
        localstuff.append(effcorr)

        localstuff.append(can_div_eff)
        localstuff.append(ratios_eff)
        
    localstuff.append(ratios)
    localstuff.append(unffile)
    localstuff.append(can)
    localstuff.append(can_div)
    localstuff.append(Hists)
    localstuff.append(fits)
    #ROOT.gApplication.Run()

    return localstuff
