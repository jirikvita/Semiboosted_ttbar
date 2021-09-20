#!/usr/bin/python

# jk 22-23.5.2019, Frascati, 3.7.2020 Olomouc

from __future__ import print_function

# IMPORTANT char** from python!;-)
# https://root-forum.cern.ch/t/creating-a-tpie-graph-from-pyroot/7306

import os, sys

import ROOT
from array import array

from math import sqrt, pow

#######################################################

def MakeTag(fname):
    tag = fname
    tag = tag.replace('analyzed_histos_', '').replace('_14TeV_ATLAS_all.root', '').replace('_', ' ').replace('zp ttbarj allhad',"Z'").replace('ptheavy','p_{T} > ' ).replace('2tj allhad', 't#bar{t}')
    return tag

#######################################################

def GetZprimeMass(sample):
    print(sample)
    if 't#bar{t}' in sample:
        return 2*172.5
    else:
        return float(sample.split()[1].replace('GeV', ''))

#######################################################
#######################################################
#######################################################

#ROOT.gStyle.SetOptTitle(0)

cans = []
stuff = []
hists = []

levels = [ 'Particle', 'Detector' ]
topos = ['0B2S', '1B1S', '2B0S']

fnames = [ 'analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root',
           #'analyzed_histos_zp_ttbarj_allhad_700GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_750GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_800GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_900GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_1250GeV_14TeV_ATLAS_all.root',
           'analyzed_histos_zp_ttbarj_allhad_1500GeV_14TeV_ATLAS_all.root',
           ]

# the main observable to get the efficiency from!
hnamebase = 'DiTopMass'           

SampleRates = {}
for fname in fnames:
    tag = MakeTag(fname)
    print('Opening file {}'.format(fname))
    rfile = ROOT.TFile(fname, 'read')
    Rates = {}
    for level in levels:
        rates = {}
        for topo in topos:
            hname = topo + '/' + level + hnamebase
            print('   ...getting {}'.format(hname))
            h = rfile.Get(hname)
            rates[topo] = h.Integral(0, h.GetNbinsX()+1)
        Rates[level] = rates
    SampleRates[tag] = Rates

print('Sample rates: ', SampleRates)

canname = 'TopologyFractions'
can = ROOT.TCanvas(canname, canname, 0, 0, 1828, 367)
can.Divide(len(fnames), len(levels))
cols = array('i', [ ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kViolet+1])
marks = { levels[0]: [20, 21, 22], levels[1] : [24, 25, 26] }
gSize = 2

gFracs = {}
gLines = {}

for level in levels:
    gfracs = []
    glines = []
    for topo in topos:
        i = topos.index(topo)
        gr = ROOT.TGraphErrors()
        gr.SetName('gfrac_{}_{}'.format(level,topo))
        gr.SetLineColor(cols[i])
        gr.SetMarkerColor(cols[i])
        gr.SetMarkerSize(1.5)
        gr.SetMarkerStyle(marks[level][i])
        #gfracs[topo] = gr
        gfracs.append(gr)
        glines.append(0)
    gFracs[level] = gfracs
    gLines[level] = glines

zpmax = 1600.
zpmin = 600.

# now make the total sum and compute the 2B0S 1B1S and 0B2S fractions;-)
# TODO!
ican = 0
for sample in SampleRates:
    for level in SampleRates[sample]:
        valsum = 0.
        vals = {}
        #print('--- {} {} ---'.format(sample, level))
        for topo in SampleRates[sample][level]:
            val = SampleRates[sample][level][topo]
            #print(val)
            valsum = valsum + val
            vals[topo] = val
        fracs = {}
        labels = []
        ff = []
        for topo in vals:
            frac = vals[topo]/valsum
            fracs[topo] = frac
            labels.append(array('c', topo + '\0'))
            ff.append(frac)
            gr = gFracs[level][topos.index(topo)]
            #gr = gFracs[topo]
            ip = gr.GetN()
            if 't#bar{t}' in sample:
                # draw a line;)
                xline = ROOT.TLine(zpmin, frac, zpmax, frac)
                gLines[level][topos.index(topo)] = xline
                continue
            mass = GetZprimeMass(sample)
            print('Setting point {} {} {} {}'.format(ip, mass, frac, sqrt(vals[topo]) / valsum))
            gr.SetPoint(ip, mass, frac)
            gr.SetPointError(ip, 0, sqrt(vals[topo])/valsum) # sqrt(frac*(1-frac)*valsum))
            
        #print(fracs)
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('{}, {} level'.format(sample, level))
        piename = 'pie_{}_{}'.format(sample, level).replace(' ', '')
        pietitle = sample + ' ' + level
        # https://root.cern.ch/doc/master/piechart_8C.html
        pie = ROOT.TPie(piename, pietitle, len(fracs), array('d', ff), cols)
        stuff.append(pie)
        pie.SetLabels(array( 'l', map( lambda x: x.buffer_info()[0], labels ) ) )
        for topo in fracs:
            print('{}  &  {:1.3f} \\\\'.format(topo, fracs[topo]))
        pie.SetEntryRadiusOffset(2,.05)
        #pie.SetEntryLineColor(2,2)
        #pie.SetEntryLineWidth(2,5)
        #pie.SetEntryLineStyle(2,2)
        #pie.SetEntryFillStyle(1,3030)
        pie.SetY(.32);
        pie.SetCircle(.5,.45,.3)
        # TODO: verbose labels?
        pie.SetLabelsOffset(.03)
        #pie.SetLabelFormat("#splitline{%val (%perc)}{%txt}")
        pie.SetLabelFormat("%txt %val")
        #pie.SetTextSize(0.05)
        
        utag = sample.split()[-1]
        tindex = -1
        for fname in fnames:
            if utag in fname:
                tindex = fnames.index(fname)
                break
        lf = len(fnames)
        jcan = tindex + levels.index(level)*lf + 1
        print('jcan etc: ', tag, tindex, ican, ican % lf, levels.index(level), jcan)

        can.cd( jcan  )
        ican = ican + 1
        pie.Draw('t nol <') # '3d '
        #pie.Draw("nol <")
        #pieleg = pie.MakeLegend()
        #pieleg.SetY1(.56)
        #pieleg.SetY2(.86)




gcanname = 'TopologyFractionsGraphs'
gcan = ROOT.TCanvas(gcanname, gcanname, 0, 0, 1000, 1000)
tmpname = 'htmp'
nbx = 100
nby = nbx
fmin = 0.
fmax = 1.
tmph = ROOT.TH2D(tmpname, tmpname + ";m_{Z'} [GeV];topology fraction", nbx, zpmin, zpmax, nby, fmin, fmax)
tmph.SetStats(0)
tmph.SetTitle('')
tmph.Draw()
opt = 'P'
funs = []
leg = ROOT.TLegend(0.12, 0.63, 0.35, 0.89)
leg.SetBorderSize(0)
for level in gFracs:
    ifrac = -1
    for gfrac in gFracs[level]:
        ifrac = ifrac + 1
        #gFracs[gfrac].Draw(opt)
        gfrac.Draw(opt)
        legtag = gfrac.GetName().replace('gfrac_','').replace('_',' ')
        leg.AddEntry(gfrac, legtag, 'PL')
        grname = gfrac.GetName()
        funname = 'fit_{}'.format(grname,)
        fun = ROOT.TF1(funname, '[0] + [1]*exp([2]*x)', zpmin, zpmax)
        fun.SetParameters(0., 1., -0.001)
        if '2B0S' in grname:
            fun.SetParameters(0., -1., -0.001)
        fun.SetLineColor(gfrac.GetLineColor())
        fun.SetLineWidth(gSize)
        if level == levels[1]:
            fun.SetLineStyle(2)
            gfrac.SetLineStyle(2)
        xline = gLines[level][ifrac]
        xline.SetLineColor(gfrac.GetLineColor())
        xline.SetLineWidth(gSize)
        xline.SetLineStyle(gfrac.GetLineStyle())
        xline.Draw()
        if level == levels[1]:
            xx = zpmax + 0.005*(zpmax - zpmin)
            topotag = legtag.replace(levels[0],'').replace(levels[1],'')
            txt = ROOT.TLatex(xx, xline.GetY1(), 't#bar{t} ' + topotag)
            txt.SetTextSize(0.029)
            txt.SetTextColor(xline.GetLineColor())
            txt.Draw()
            stuff.append(txt)
        fun.SetLineColor(gfrac.GetLineColor())
        gfrac.Fit(funname)
        funs.append(fun)
leg.Draw()
stuff.append(leg)

pngdir = 'png/'
pdfdir = 'pdf/'
can.Print(pngdir + canname + '.png')
can.Print(pdfdir + canname + '.pdf')
gcan.Print(pngdir + gcanname + '.png')
gcan.Print(pdfdir + gcanname + '.pdf')

ROOT.gApplication.Run()
