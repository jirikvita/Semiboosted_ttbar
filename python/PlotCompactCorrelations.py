#!/usr/bin/python
# Thu 30 Jul 21:23:32 CEST 2020, 5.8.2020

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from CorrItems import *

cans = []
stuff = []


topos = ['0B2S',
         '1B1S',
         '2B0S'
]

##########################################

def TransposeH2(h2, tag = ''):
    name = h2.GetName() + '_trans' + tag
    title = h2.GetTitle() + '_trans'
    ny = h2.GetNbinsY()
    nx = h2.GetNbinsX()
    print('Transposing {} mux={} muy={} sx={} sy={}'.format(h2.GetName(), h2.GetMean(1), h2.GetMean(2), h2.GetRMS(1), h2.GetRMS(2)) )
    trans = ROOT.TH2D(name, title,
                      ny, h2.GetYaxis().GetXbins().GetArray(),
                      nx, h2.GetXaxis().GetXbins().GetArray() )
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            val = h2.GetBinContent(i,j)
            err = h2.GetBinError(i,j)
            trans.SetBinContent(j,i, val)
            trans.SetBinError(j,i, err)
    trans.GetXaxis().SetTitle(h2.GetYaxis().GetTitle())
    trans.GetYaxis().SetTitle(h2.GetXaxis().GetTitle())
    trans.Scale(1.)
    return trans

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    gTag=''
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    marg = 0. # 0.12
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTopMargin(marg)
    ROOT.gStyle.SetPadBottomMargin(marg)
    ROOT.gStyle.SetPadLeftMargin(marg)
    ROOT.gStyle.SetPadRightMargin(marg)

    ROOT.gStyle.SetPaintTextFormat("1.2f")

    # https://root.cern.ch/doc/master/classTColor.html
    #ROOT.gStyle.SetPalette(ROOT.kDeepSea)
    #ROOT.gStyle.SetPalette(ROOT.kCool)
    #ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    #ROOT.gStyle.SetPalette(ROOT.kCopper)
    #ROOT.gStyle.SetPalette(ROOT.kBird) # std, nice
    #ROOT.gStyle.SetPalette(ROOT.kAvocado)
    #ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
    #ROOT.gStyle.SetPalette(ROOT.kLake)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    filename = 'foo.root'
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    rfile = ROOT.TFile(filename, 'read')

    levels = ['Detector' ] #, 'Particle']
    H2s = []

    corrhs = []
    for level in levels:
        h2s = []
        n = 0
        for var in ivarsLabelsDict:
            n = n+1
        n = n - 5 ### HACK!!! for variables removed from correlation studies! jk 6.3.2021
        for topo in topos:
            canname = 'corrs_{}_{}'.format(level,topo)
            can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 1200)
            mmarg = 0. # 0.01
            can.Divide(n,n, mmarg, mmarg)
            cans.append(can)
            cname = 'corrh2_{}_{}'.format(level, topo)
            m = n
            corrh = ROOT.TH2D(cname, cname, m,0,m, m,0,m)
            corrhs.append(corrh)
            for basehname in hnamesDict:
                lhname = basehname + ''
                
                hname = '{}/{}{}'.format(topo, level, lhname)
                hname = hname.replace('Vs', 'Vs' + level)
                h2 = rfile.Get(hname)
                try:
                    h2.SetStats(0)
                except:
                    print('ERROR getting {} from file {}'.format(hname, rfile.GetName()))
                    continue
                #h2.Scale(1./h2->Integral())
                h2s.append(h2)
                ii = hnamesDict[basehname]
                i,j = ii[0], ii[1]
                if ii[2]:
                    i,j = j,i
                ind = j*n + i + 1 # was (j-1)*n + i!
                print('{} {},{} ind={}'.format(basehname, i, j, ind))
                can.cd(ind)
                rho = 0.
                if ii[2]:
                    trans = TransposeH2(h2, '_' + level + '_' + topo)
                    trans.SetStats(0)
                    h2s.append(trans)
                    trans.Draw('colz')
                    rho = trans.GetCorrelationFactor()
                else:
                    h2.Draw('colz')
                    rho = h2.GetCorrelationFactor()
                corrh.SetBinContent(i+1,corrh.GetNbinsY() - j,rho)
                tex = ROOT.TLatex(0.17, 0.7, '{:1.2f}'.format(rho))
                tex.SetNDC()
                tex.SetTextSize(0.32)
                tex.SetTextColor(ROOT.kRed)
                tex.Draw()
                stuff.append(tex)
                stuff.append(h2s)

            # print variable label!
            for var in ivarsDict:
                ivar = ivarsDict[var]
                ind = (ivar)*n + ivar + 1  # was: (ivar-1)*n + ivar 
                can.cd(ind)
                tex = ROOT.TLatex(0.50 - 0.002*len(ivarsLabelsDict[var]), 0.45, '{}'.format(ivarsLabelsDict[var]))
                tex.SetNDC()
                tst = 0.20 #0.32
                tex.SetTextSize(tst)
                tex.SetTextColor(ROOT.kBlack)
                if ind != n*n:
                    tex.Draw()
                    stuff.append(tex)

                ind = (ivar)*n + ivar+1 # was (ivar-1)
                can.cd(ind)
                tex2 = ROOT.TLatex(0.30 - 0.002*len(ivarsLabelsDict[var]), 0.77, '{}'.format(ivarsLabelsDict[var]))
                tex2.SetNDC()
                tex2.SetTextSize(tst)
                tex2.SetTextColor(ROOT.kBlack)
                if ind > 1:
                    tex2.Draw()
                    stuff.append(tex2)
                corrh.GetXaxis().SetBinLabel(ivar+1, ivarsLabelsDict[var])
                corrh.GetYaxis().SetBinLabel(corrh.GetNbinsY() - ivar, ivarsLabelsDict[var])
                corrh.GetXaxis().ChangeLabel(ivar+1, 90)
                corrh.GetXaxis().SetTitleOffset(-2.)
                corrh.LabelsOption('v')
                pass
                
                
            H2s.append(h2s)

        corrcans = []
        itopo = -1
        for corrh in corrhs:
            itopo = itopo + 1
            canname = 'corrsh_{}_{}'.format(level,topos[itopo])
            corrcan = ROOT.TCanvas(canname, canname, 0, 0, 1100, 900)
            corrcans.append(corrcan)
            corrcan.cd()
            #ROOT.gPad.SetPad(0.1, 0.1, 0.9, 0.9)
            ROOT.gPad.SetBorderSize(0)
            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.15)
            ROOT.gPad.SetTopMargin(0.05)
            ROOT.gPad.SetBottomMargin(0.15)
            ROOT.gPad.SetFrameBorderMode(0)
            ROOT.gPad.SetGridx(1)
            ROOT.gPad.SetGridy(1)
            corrh.Scale(1.)
            corrh.SetMinimum(-1.)
            corrh.SetMaximum(1.)
            corrh.SetStats(0)
                
            corrh.Draw("colz text")
            corrcan.Print('png/' + corrcan.GetName() + '.png')
            corrcan.Print('pdf/' + corrcan.GetName() + '.pdf')

        for can in cans:
            can.Print('png/' + can.GetName() + '.png')
            can.Print('pdf/' + can.GetName() + '.pdf')


    stuff.append(rfile)
    ROOT.gApplication.Run()
    return

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

