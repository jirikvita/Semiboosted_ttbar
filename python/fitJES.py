#!/usr/bin/python
# Fri 24 Jul 16:13:30 CEST 2020

from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

kDiff = 0
kRelDiff = 1
kRatio = 2

##########################################
def PrintBinContent2D(ratio):
    for i in range(1,ratio.GetXaxis().GetNbins()+1):
        for j in range(1,ratio.GetYaxis().GetNbins()+1):
            val = ratio.GetBinContent(i,j)
            print(val)
    return

##########################################
def Make2DResiduals(h2, fun, dtype = kRelDiff, tag = '_reldif'):
    #print(h2)
    nx = h2.GetXaxis().GetNbins()
    ny = h2.GetYaxis().GetNbins()
    name = h2.GetName() + tag
    # no support for non-uniform binning yet...:(
    #res = ROOT.TH2D(name,'',
    #                nx,h2.GetXaxis().GetXbins().GetArray(),
    #                ny,h2.GetYaxis().GetXbins().GetArray())
    res = ROOT.TH2D(name,'residuals',
                    nx,h2.GetXaxis().GetXmin(),h2.GetXaxis().GetXmax(),
                    ny,h2.GetYaxis().GetXmin(),h2.GetYaxis().GetXmax() )
    res.GetXaxis().SetTitle(h2.GetXaxis().GetTitle())
    res.GetYaxis().SetTitle(h2.GetYaxis().GetTitle())
    res.GetZaxis().SetTitle(h2.GetZaxis().GetTitle())
    zmin = -0.5
    zmax = 0.5
    if dtype == kDiff:
        res.SetTitle('Difference')
    elif dtype == kRatio:
        res.SetTitle('Ratio')
        zmin = 0.8
        zmax = 1.2
    elif dtype == kRelDiff:
        res.SetTitle('Relative difference')
    else:
        print('ERROR: nonrecognized residual option {}, allowed are {},{},{}'.format(dtype,kDiff,kRatio,kRatio))
        return ROOT.TObject()
    
    # for x and y binned residuals:
    nz = 400
    resx = ROOT.TH2D(name + '_rx','residualsx',nx,h2.GetXaxis().GetXmin(),h2.GetXaxis().GetXmax(), nz, zmin, zmax)
    resx.GetXaxis().SetTitle(res.GetXaxis().GetTitle())
    resx.GetYaxis().SetTitle(res.GetZaxis().GetTitle())
    resy = ROOT.TH2D(name + '_ry','residualsy',ny,h2.GetYaxis().GetXmin(),h2.GetYaxis().GetXmax(), nz, zmin, zmax)
    resy.GetXaxis().SetTitle(res.GetYaxis().GetTitle())
    resy.GetYaxis().SetTitle(res.GetZaxis().GetTitle())
    resd = ROOT.TH1D(name + '_resall','residuals', nz, zmin, zmax)
    resd.GetXaxis().SetTitle(res.GetZaxis().GetTitle())
        
    # residuals distribution
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            val = h2.GetBinContent(i,j)
            if val > 0.:
                x = h2.GetXaxis().GetBinCenter(i)
                y = h2.GetYaxis().GetBinCenter(j)
                fitval = fun.Eval(x,y)
                err = h2.GetBinError(i,j)
                diff = val - fitval
                doset = True
                if dtype == kRatio or dtype == kRelDiff:
                    doset = False
                    if fitval > 0.:
                        doset = True
                        err = err / fitval
                        if dtype == kRatio:
                            diff = val / fitval
                            # print(diff)
                        elif dtype == kRelDiff:
                            diff = diff / fitval
                    else:
                        diff = 0.
                        err = 0.
                if doset:
                    #print('   ....setting {}'.format(diff))
                    res.SetBinContent(i,j,diff)
                    res.SetBinError(i,j,err)
                    resx.Fill(h2.GetXaxis().GetBinCenter(i), diff)
                    resy.Fill(h2.GetYaxis().GetBinCenter(i), diff)
                    resd.Fill(diff)
    res.Scale(1.)
    print(res)
    return res,resx,resy,resd


##########################################
def PrintChi2(fun, x = 0.15, y = 0.15):
    chi2 = fun.GetChisquare()
    ndf = fun.GetNDF()
    ss = ''
    if ndf > 0:
        ss = '#chi^{2}/ndf' + ' = {:3.1f}/{} = {:3.1f}'.format(chi2,ndf, chi2/ndf)
    else:
        ss = '#chi^{2}/ndf' + ' = {:3.1f}/{}'.format(chi2,ndf)
    print(ss)
    tex = ROOT.TLatex(x,y,ss)
    tex.SetNDC()
    tex.Draw()
    return chi2,ndf,tex

##########################################
def PrintProperties(h3):
    nx = h3.GetXaxis().GetNbins()
    ny = h3.GetYaxis().GetNbins()
    nz = h3.GetZaxis().GetNbins()
    print('The 3D histogram {} has nx*ny*nz bins: {}*{}*{} = {} total, nx*ny={}'.format(h3.GetName(),nx,ny,nz, nx*ny*nz,nx*ny))
    return
##########################################
##########################################
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


    ############################################
    #     JES fit to small and large jets!     #
    ############################################
    
    ##################################################################
    lstuffL = FitJES(argv, 'JESRPtLJetsEtaE', 'recreate')
    # alt.: lstuffL = FitJES(argv, 'JESRELJetsEtaE', 'recreate')
    ##################################################################
    lstuffS = FitJES(argv, 'JESRPtJetsEtaE', 'update')
    # alt.: lstuffS = FitJES(argv, 'JESREJetsEtaE', 'update')
    ##################################################################

    
    return [lstuffS,lstuffL]
    #return lstuffL
    
##########################################
def FitJES(argv, hname = 'JESRPtLJetsEtaE', fopt = 'recreate', dirname = "JES/"):

    lstuff = []
    # ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    ROOT.gStyle.SetOptTitle(0)
    if len(argv) < 2:
        print('Usage: {} filename.root'.format(argv[0]))
        print('E.g.:')
        print('{} data/analyzed_histos_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV_ATLAS_all_forJES.root'.format(argv[0]))
        exit(1)


    isLjets = 'LJets' in hname

    print('---------------- Welcome to JES fit! ;-) ----------------')

    cantag = '_smallJets'
    rbx = 4
    rby = 4
    if 'LJets' in hname:
        rbx = 4
        rby = 4
        cantag = '_largeJets'
        
    
    canname = 'JESfit' + cantag
    can = ROOT.TCanvas(canname, canname, 0, 0, 800, 1200)
    can.Divide(2, 3)
    cans.append(can)

    canname = 'tmpjesfit' + cantag
    cantmp = ROOT.TCanvas(canname, canname, 0, 0, 1200, 1200)
    cantmp.Divide(2, 2)
    cans.append(cantmp)
    
    filename = argv[1]
    print('Opening {}'.format(filename))
    rfile = ROOT.TFile(filename, 'read')
    stuff.append(rfile)
    
    ##########################################
    # 2D spectrum, rebin
    h3 = rfile.Get(dirname + hname)
    h3.SetMarkerColor(ROOT.kBlue)
    h3.SetStats(0)
    PrintProperties(h3)
    stuff.append(h3)
    
    h3.Rebin3D(rbx, rby, 1)
    PrintProperties(h3)
    
    can.cd(1)
    h3copy = h3.DrawCopy()
    stuff.append(h3)
    
    ph2 = h3.Project3DProfile('xy')
    ph2.SetStats(0)
    # it's swapped, for some reason...;-)
    ph2.GetXaxis().SetTitle(h3.GetYaxis().GetTitle())
    ph2.GetYaxis().SetTitle(h3.GetXaxis().GetTitle())
    ph2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())
    ph2.GetZaxis().SetRangeUser(0.8, 1.2)

    can.cd(3)
    ph2copy = ph2.DrawCopy('lego2')
    stuff.append(ph2copy)
    
    can.cd(5)
    ph2.DrawCopy('colz')

    x1 = ph2.GetXaxis().GetXmin()
    x2 = ph2.GetXaxis().GetXmax()
    y1 = ph2.GetYaxis().GetXmin()
    y2 = ph2.GetYaxis().GetXmax()

    fname = 'fit1dx'
    funx = ROOT.TObject()
    if True: #not isLjets:
        fform = '[2] + [0]*exp([1]*x)'
        funx = ROOT.TF1(fname, fform, x1, x2)
        funx.SetParameters(-0.09,-0.005, 0.8)
    else:
        #fform = '[5] + [0]*log([2]*(x-[1]) + [3]*(x-[1])^2 + [4]*(x-[1])^3)'
        #funx = ROOT.TF1(fname, fform, x1, x2)
        #funx.SetParameters(0.05, 20, 0.01, -0.00001, 0.000001, 0.8)

        fform = '[2] + [0]*log([1]*x)'
        funx = ROOT.TF1(fname, fform, x1, x2)
        funx.SetParameters(0.045, 1.e-3, 0.9)

        #fform = '[5] + [0]*log([2]*(x-[1]) + [3]*(x-[1])^2 + [4]*(x-[1])^3)'
        #funx = ROOT.TF1(fname, fform, x1, x2)
        #funx.SetParameters(0.05, 20, 0.01, -0.00001, 0.000001, 0.8)


        
    #can.cd(7)
    #funx.DrawCopy()

    fname = 'fit1dy'
    fform = '[3] + [0]*exp(-( abs(x)-[1])^2 / (2*[2]^2) )'
    funy = ROOT.TF1(fname, fform, y1, y2)
    funy.SetParameters( 0.05,1.7,0.5,  0.9)
    #can.cd(8)
    #funy.DrawCopy()
    
    cantmp.cd(1)
    p1x = h3.ProjectionX()
    p1x.GetXaxis().SetTitle(h3.GetXaxis().GetTitle())
    p1x.SetMinimum(0.)
    p1x.Draw('e1 hist')

    cantmp.cd(2)
    p1y = h3.ProjectionY()
    p1y.GetYaxis().SetTitle(h3.GetYaxis().GetTitle())
    p1y.SetMinimum(0.)
    p1y.Draw('e1 hist')

    ######################
    # 1D fits
    ######################
    print('----------------  1D fits  ----------------')
    cantmp.cd(3)
    ph1y = ph2.ProfileY(ph2.GetName() + '_ProfileY')
    ph1y.GetXaxis().SetTitle(h3.GetXaxis().GetTitle())
    print('...fitY')
    ph1y.Fit('fit1dy')
    ph1y.Draw('e1 hist')
    funy.Draw('same')
    chi2y,ndfy,texy = PrintChi2(funy)

    ipad = 4
    cantmp.cd(ipad)
    ph1x = ph2.ProfileX(ph2.GetName() + '_ProfileX')
    ph1x.GetXaxis().SetTitle(h3.GetYaxis().GetTitle())
    #ph1x.SetName(ph2.GetName() + '_ProfileX')
    print('...fitX')
    ph1x.Fit('fit1dx')
    ph1x.DrawCopy('e1 hist')
    funx.Draw('same')
    chi2x,ndfx,texx = PrintChi2(funx)
    
    #can.cd(6)
    #p1z = h3.Project3D('z')
    #p1z.Draw('e1 hist')

    ######################
    # 2D fit
    ######################
    print('----------------   2D fit  ----------------')
    fname = 'fit2d_' + hname 
    # x...E or pT
    # y ... eta

    fun = ROOT.TObject()
    if not isLjets:
        fform = '([9] + ([0]+[1]*abs(y))*exp(([2]+[3]*abs(y))*x)) * ([8] + ([4] + [5]*x)*exp(-(abs(y)-[6])^2 / (2*[7]^2) ) )'
        fun = ROOT.TF2(fname, fform, x1, x2, y1, y2)
        fun.SetParameters(-0.15,0.0,-0.004,0.0,     0.05,-0.0001,2.15,0.8,        0.8,0.8)
    else:
        fform = '([9] + ([0]+[1]*abs(y))*log(([2]+[3]*abs(y))*x)) * ([8] + ([4] + [5]*x)*exp(-(abs(y)-[6])^2 / (2*[7]^2) ) )'
        fun = ROOT.TF2(fname, fform, x1, x2, y1, y2)
        fun.SetParameters(0.05,0.0,0.001,0.0,     0.05,-0.0001,2.15,0.8,        0.8,0.8)
    fun.GetXaxis().SetTitle(h3.GetYaxis().GetTitle())
    fun.GetYaxis().SetTitle(h3.GetXaxis().GetTitle())
    ### IMPORTANT, fit is not sensitive to this parameter neither in exp nor in log parameterizations!
    fun.FixParameter(3, 0.)
    # main 2D fit!
    ph2.Fit(fun)
    ph2.Fit(fun)
    
    # print results:
    parnames = [r'\mathcal{A}_{E}^{0}', r'\mathcal{A}_{E}^{1}', r'\beta_0', r'\beta_1',
                r'\mathcal{D}^{0}', r'\mathcal{D}^{1}', r'\eta_0', r'\sigma_\eta',
                r'\mathcal{C}_E', r'\mathcal{C}_\eta'
                ]
    for i in range(0,fun.GetNpar()):
        #fun.SetParName(i,parnames[i])
        print('${}$ & {:1.5f} & {:1.5f} \\\\'.format(parnames[i], fun.GetParameter(i), fun.GetParError(i)))
    
    # does not work for z axis?
    fun.SetRange(x1, y1, 0.8, x2, y2, 1.2)

    can.cd(6)
    funcopy = fun.DrawCopy('colz')
    chi2,ndf,tex = PrintChi2(fun)
    stuff.append(funcopy)
    
    can.cd(4)
    fun.DrawCopy('lego2')

    
    #can.cd(2)
    #ph2 = h3.Project3DProfile('xy')
    #ph2.Draw('colz')
    
    ######################
    # fit residuals
    ######################
    print('---------------- residuals --------------')
    canname = 'JESfit2dresiduals' + cantag
    canres = ROOT.TCanvas(canname, canname, 0, 0, 1200, 1200)
    canres.Divide(2,2)
    cans.append(canres)

    reldif,rdx,rdy,rr = Make2DResiduals(ph2, fun, kRelDiff, 'reldif')
    #diff,dx,dy,rrd = Make2DResiduals(ph2, fun, kDiff, 'diff')
    ratio,rx,ry,rrr = Make2DResiduals(ph2, fun, kRatio, 'ratio')

    texs = [texx, texy, tex]
    ratios = [ratio, reldif] #, diff]
    for h in ratios:
        h.SetStats(0)
    #ratios.append([rrd,rr,rrr])
    
    #zz = 0.5
    #reldif.GetZaxis().SetRangeUser(-zz,zz)
    #diff.GetZaxis().SetRangeUser(-zz,zz)
    ratio.GetZaxis().SetRangeUser(0.9, 1.1)
    #canres.cd(1)
    #reldif.Draw('colz')
    #canres.cd(2)
    #diff.Draw('colz')
    canres.cd(1)
    ratio.Draw('colz')

    rs = [rx, ry]
    for rrs in rs:
        rrs.SetStats(0)
        rrs.GetYaxis().SetRangeUser(0.95, 1.05)
    
    canres.cd(3)
    rx.Draw('colz')
    prx = rx.ProfileX()
    prx.Draw('same')
    canres.cd(4)
    ry.Draw('colz')
    pry = ry.ProfileX()
    pry.Draw('same')
    canres.cd(2)
    ROOT.gPad.SetLogy(1)
    rrr.Draw('e1hist')
    rs.append([rr,rrr])

    # again...the troubling plot...
    cantmp.cd(ipad)
    #ph1x = ph2.ProfileX(ph2.GetName() + '_ProfileX')
    #ph1x.SetName(ph2.GetName() + '_ProfileX')
    #ph1x.Fit('fit1dx')
    ph1x.Draw('e1 hist')
    funx.Draw('same')
    #chi2x,ndfx,texx = PrintChi2(funx)
    
    ######################
    # print and save fit result
    ######################

    print('---------------- printing ----------------')
    can.Update()
    canres.Update()
    can.Print('png/' + can.GetName() + '.png')
    canres.Print('png/' + canres.GetName() + '.png')
    can.Print('pdf/' + can.GetName() + '.pdf')
    canres.Print('pdf/' + canres.GetName() + '.pdf')
    
    outfile = ROOT.TFile('JESfits.root', fopt)
    outfile.cd()
    fun.Write()
    # ph2.Write()
    outfile.Write()

    lstuff.append([cans,h3, ph2, p1x, p1y, ph1x, ph1y, fun, funx, funy])
    lstuff.append([ratios, rs, outfile])
    lstuff.append(texs)

    stuff.append(lstuff)
    
    #print(stuff)
    #print(ph1x)
    
    outfile.Close()
    return lstuff

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    llstuff = []
    local = main(sys.argv)
    stuff.append(local)
    llstuff.append(local)
    ROOT.gApplication.Run()

    
    
###################################
###################################
###################################

