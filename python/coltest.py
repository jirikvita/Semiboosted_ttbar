#!/usr/bin/python

from __future__ import print_function

import ROOT

x1 = -2
x2 = 2
funname = 'fun'
fun = ROOT.TF1(funname, '[0]*exp( -(x-[1])^2 / (2*[2]^2))', x1, x2)


Pars = [ [1., 0, 0.5],
         [1., 0., 0.5],
         [1., 0., 0.5],
]

stuff = []

hs = []
opt = ''
N = 1000
#cols = [ ROOT.kViolet-2,
#         ROOT.kGray + 2,
#         ROOT.kAzure + 7]
cols = [          ROOT.kOrange + 7,
                  ROOT.kMagenta,
         ROOT.kGreen + 1]
fst = [ 3345, 3354, 3305]
ih = -1
nb = 10
maxy = -999
for pars in Pars:
    print(pars)
    ih = ih + 1
    hname = 'h{}'.format(ih)
    for i in range(0, len(pars)):
        print(i, pars[i])
        fun.SetParameter(i, pars[i])
    h = ROOT.TH1D(hname, hname, nb, x1, x2)
    h.SetStats(0)
    hs.append(h)
    h.FillRandom(funname, int(N*(1 + N/100*ih)) )
    h.Scale(1.)
    val = h.GetMaximum()
    if maxy < val:
        maxy = val
    h.SetFillStyle(fst[ih])
    h.SetFillColor(cols[ih])
    h.SetLineWidth(2)
    h.SetLineColor(cols[ih])
    h.Draw('hist' + opt)
    opt = 'same'
    ROOT.gPad.Update()
    
for h in hs:
    h.SetMaximum(1.25*maxy)
    
cand = ROOT.TCanvas()
eff = ROOT.TEfficiency(hs[0], hs[1])
cand.cd()
eff.Draw('')
stuff.append([eff, cand])

ROOT.gApplication.Run()
