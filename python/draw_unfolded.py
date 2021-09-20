#!/usr/bin/python
# jk 26.2.2018, 30.3.2018, 2.12.2019, 11.6.2020

from __future__ import print_function

from DrawUnfolded import *
from Tools import MakePrettyTitle

stuff = []

rebin = -1 # !!!

Ingredients = [
    ['unfolded_0B2S_TopPt_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],
   ['unfolded_0B2S_DiTopMass_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],

   ['unfolded_1B1S_TopPt_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],
    ['unfolded_1B1S_DiTopMass_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],

    ['unfolded_2B0S_DiTopMass_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],
    ['unfolded_2B0S_TopPt_Particle_Detector_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_half0.root'],
              ]

ROOT.gStyle.SetPaintTextFormat("1.2f")
ROOT.gStyle.SetPalette(1)

os.system('mkdir -p png/unfolded/')
os.system('mkdir -p pdf/unfolded/')

for ingredients in Ingredients:
    fname = ingredients[0]
    tokens = fname.split('_')
    hname = '{}_{}_{}'.format(tokens[2],tokens[3],tokens[4])
    print(hname)
    level = 'Particle'

    print('Processing {}'.format(fname) )
    print('Drawing closure for {} from file {}'.format(fname, hname,))
    xlabel = MakePrettyTitle(hname.split('_')[0])
    xlabel = xlabel.replace('Rapidity', ' rapidity')
    unfstuff = DrawUnfolded(ingredients, hname, level, xlabel, rebin)
    stuff.append(unfstuff)


ROOT.gApplication.Run()
