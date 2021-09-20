#!/usr/bin/python
# Fri 19 Mar 11:42:07 CET 2021

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

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

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    canname = 'can'
    can = ROOT.TCanvas(canname, canname)
    cans.append(can)
    gr = ROOT.TGraphErrors()
    gr.SetName('zpw')
    gr.SetTitle(';m_{Z\'} [GeV];#Gamma_{Z\'} [GeV]')
    wparam = 0.0008252
    Data = {1500. : 0.0029,
            1250. : 0.00246,
            1000. : 0.001946,
             900. : 0.001795,
             800. : 0.001637,
             750. : 0.001611 
    }
    ip = 0
    for mass in Data:
        wobs = Data[mass]
        wcomp = wparam*mass*mass
        print('mass: {}, wobs: {}, wcopm={} ratio={}'.format(mass, wobs, wcomp,wobs/wcomp))
        gr.SetPoint(ip, mass, wobs)
        # just some dummy error for ~stable fit
        gr.SetPointError(ip, 0., wobs*0.01)
        ip = ip + 1

    gr.SetMarkerSize(1)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(ROOT.kBlack)
    gr.Draw('AP')
    gr.Fit('pol2')
    stuff.append(gr)
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

