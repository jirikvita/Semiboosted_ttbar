#!/usr/bin/python

# jk IX--XII 2020

import ROOT

from math import sqrt, pow, log

##########################################
#    https://arxiv.org/abs/1101.0390
##########################################
# 23.9.2020
# d ... data, should be counts
# b ... expected background, can be real

kBHInfty = 1.e100
kBHZero = 1.e-100

#kBHBigEpsilonInfty = 1.e289
kBHSFcmp = 0.9

##########################################
def ComputeSafeLog(tt, tmin = -10, tmax = 50):
    logt = 0.
    print('computing safe log of tt={}'.format(tt))
    if tt > kBHInfty*kBHSFcmp:
        logt = tmax
    elif tt < kBHZero:
        logt = tmin
    else:
        logt = log(tt)
    print('                 log(tt)={}'.format(logt))
    return logt

##########################################
def GetGammaP(d,b):
    P = 0.
    if d >= b:
        P = ROOT.TMath.Gamma(d,b)
    else:
        P = 1. - ROOT.TMath.Gamma(d + 1.,b)
    return P

##########################################
def GetBumpHuntPval(d,b,debug = 0):
    pval = 1.
    if d > b and b > 0:
        pval = GetGammaP(d,b)
    if debug:
        print('GetBumpHuntPval d={:4.1f}, b={:4.1f} pval={}'.format(d,b,pval))
    return pval
        
##########################################
def GetBumpHuntTestStat(d,b,debug = 0):
    # this is a natural logarithm;-)
    pval = GetBumpHuntPval(d,b,debug)
    if pval > 0:
        return -log(pval)
    else:
        return kBHInfty

##########################################
# classes to store BumpHunter results
##########################################
class cBumpRes1D:
    def __init__(self, i, bins, p, t):
        self.bins = bins
        self.i = i
        self.p = p
        self.t = t

##########################################
class cBumpRes2D:
    def __init__(self, i, j, bins, p, t):
        self.bins = bins
        self.i = i
        self.j = j
        self.p = p
        self.t = t

##########################################
def MarkAsUsedWindow1D(hdata, usedbins, i, iwwx):
    for ii in range(0,iwwx):
        if i + ii >= 0 and i + ii < len(usedbins) and hdata.GetBinContent(i+ii+1) > 0.: # jk 28.12.2020
            usedbins[i+ii] = 1
    return usedbins
##########################################
def MarkAsUsedWindow2D(hdata, usedbins, i, j, iwwx, iwwy):
    for ii in range(0,iwwx):
        for jj in range(0,iwwy):
            if i + ii >= 0 and i + ii < len(usedbins):
                if j + jj >= 0 and j + jj < len(usedbins[i+ii]) and hdata.GetBinContent(i+ii+1, j+jj+1) > 0.: # jk 28.12.2020
                    usedbins[i+ii][j+jj] = 1
    return usedbins
##########################################
def CountUsedBins1D(usedbins):
    nb = 0
    for bin in usedbins:
        if bin != 0:
            nb = nb + 1
    return nb

##########################################
def CountUsedBins2D(usedbins):
    nb = 0
    for bins in usedbins:
        for bin in bins:
            if bin != 0:
                nb = nb + 1
    return nb

##########################################
def UpdateIndex(extr, ii, studyMin):
    if studyMin:
        if ii < extr:
            extr = ii
    else:
        if ii > extr:
            extr = ii
    return extr

##########################################
def GetWidthXY(usedbins):
    # gives the max x and y span, does not check for a consistent set!
    imin = 999
    imax = -1
    jmin = 999
    jmax = -1
    for i in range(0, len(usedbins)):
        for j in range(0, len(usedbins[i])):
            if usedbins[i][j] != 0:
                imin = UpdateIndex(imin, i, studyMin = True)
                imax = UpdateIndex(imax, i, studyMin = False )
                jmin = UpdateIndex(jmin, j, studyMin = True)
                jmax = UpdateIndex(jmax, j, studyMin = False)

    wx = imax - imin + 1
    wy = jmax - jmin + 1
    return wx,wy
    
##########################################
def CheckAreaPropertiesOK1D(usedbins, nmax):
    nb = CountUsedBins1D(usedbins)
    return nb <= nmax

##########################################
def CheckAreaPropertiesOK2D(usedbins, nmax, wxmax, wymax):
    nb = CountUsedBins2D(usedbins)
    wx,wy = GetWidthXY(usedbins)
    return nb <= nmax and wx <= wxmax and wy <= wymax


##########################################
def MakeEmptyUsedBins1D(hdata):
    usedbins = []
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        usedbins.append(0)
    return usedbins

##########################################
def MakeEmptyUsedBins2D(hdata):
    usedbins = []
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        bins = []
        for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
            bins.append(0)
        usedbins.append(bins)
    return usedbins


##########################################
def MakeBestBumpHunterIntervalAroundWindow(hdata, hbg, i, iwwx, nmax, nmin, title = 'p-val: ', tst = 0.045):
    # find best BH interval around the initial bin i
    usedbins = MakeEmptyUsedBins1D(hdata)
    usedbins = MarkAsUsedWindow1D(hdata, usedbins, i, iwwx)
    # keep adding neighbouring bins as long as some criteria;)
    dsum = 0
    bsum = 0
    # 26.11.2020: need to add dsum and bsum based on the initial window!
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        if usedbins[ibin-1] == 0:
            if CheckAreaPropertiesOK1D(usedbins, nmax):
                data = hdata.GetBinContent(ibin)
                bg = hbg.GetBinContent(ibin)
                dsum = dsum + data
                bsum = bsum + bg
    # now try adding more bins:
    changed  = True
    while changed:
       changed = False
       for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
           if usedbins[ibin-1] == 0 and BinNeighboursUsedBins1D(ibin-1, usedbins):
               if CheckAreaPropertiesOK1D(usedbins, nmax):
                   data = hdata.GetBinContent(ibin)
                   bg = hbg.GetBinContent(ibin)
                   #edata = hdata.GetBinError(ibin)
                   #ebg = hbg.GetBinError(ibin)
                   #etotsq = pow(edata, 2) + pow(ebg, 2)
                   #if etotsq > 0. and data-bg > 0:
                   if data > 0 and GetBumpHuntTestStat(dsum + data, bsum + bg) > GetBumpHuntTestStat(dsum,bsum):
                       #signif = (data - bg) / sqrt(etotsq)
                       #if signif > thr:
                       usedbins[ibin-1] = 1
                       dsum = dsum + data
                       bsum = bsum + bg
                       changed = True

    p = GetBumpHuntPval(dsum,bsum)
    t = GetBumpHuntTestStat(dsum,bsum)
    result = cBumpRes1D(i, usedbins, p, t)
    return result

##########################################
def MakeBestBumpHunterAreaAroundWindow(hdata, hbg, i, j, iwwx, iwwy, nmax, nmin, wxmax, wymax):
    usedbins = MakeEmptyUsedBins2D(hdata)
    usedbins = MarkAsUsedWindow2D(hdata, usedbins, i, j, iwwx, iwwy)
    # keep adding neighbouring bins as long as some criteria;)
    dsum = 0
    bsum = 0

    # 26.11.2020: need to add dsum and bsum based on the initial window!
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
            if usedbins[ibin-1][jbin-1] == 0:
                if CheckAreaPropertiesOK2D(usedbins, nmax, wxmax, wymax):
                    dsum = dsum + hdata.GetBinContent(ibin,jbin)
                    bsum = bsum + hbg.GetBinContent(ibin,jbin)
    
    changed  = True
    while changed:
       changed = False
       for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
           for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
               if usedbins[ibin-1][jbin-1] == 0 and BinNeighboursUsedBins2D(ibin-1,jbin-1, usedbins):
                   if CheckAreaPropertiesOK2D(usedbins, nmax, wxmax, wymax):
                       data = hdata.GetBinContent(ibin,jbin)
                       bg = hbg.GetBinContent(ibin,jbin)
                       #edata = hdata.GetBinError(ibin,jbin)
                       #ebg = hbg.GetBinError(ibin,jbin)
                       #etotsq = pow(edata, 2) + pow(ebg, 2)
                       #if etotsq > 0. and data-bg > 0:
                       if data > 0 and GetBumpHuntTestStat(dsum + data, bsum + bg) > GetBumpHuntTestStat(dsum,bsum):
                           #signif = (data - bg) / sqrt(etotsq)
                           #if signif > thr:
                           usedbins[ibin-1][jbin-1] = 1
                           dsum = dsum + data
                           bsum = bsum + bg
                           changed = True

    p = GetBumpHuntPval(dsum,bsum)
    t = GetBumpHuntTestStat(dsum,bsum)
    result = cBumpRes2D(i, j, usedbins, p, t)
    return result

##########################################
def FindBestResultIndex(results):
    # find the best result:
    tmax = -999.
    ibest = -1
    ir = 0
    for result in results:
        if result.t > tmax:
            ibest = ir
            tmax = result.t
        ir = ir+1
    return ibest

##########################################
def FindBestBumpHunter1DInterval(hdata, hbg, x, y, title = 'log(t): ', tst = 0.045, showp = False):
    nx = hdata.GetXaxis().GetNbins()

    # limit parameters of the search area:
    nmax = nx / 2
    nmin = 3
    # initial window widths
    iwwx = 2
    # shift initial window
    results = []
    for i in range(1,nx+1-iwwx):
        result = MakeBestBumpHunterIntervalAroundWindow(hdata, hbg, i, iwwx, nmax, nmin)
        results.append(result)
    ibest = FindBestResultIndex(results)
    if ibest >= 0:
        bhtex = ROOT.TLatex(x, y, '')
        if showp:
            bhtex.SetText(x, y, title + '{:1.2f}'.format(results[ibest].p))
        else:
            tt = results[ibest].t
            if tt > kBHInfty*kBHSFcmp:
                bhtex.SetText(x, y, title + '+#infty')
            elif tt < kBHZero:
                bhtex.SetText(x, y, title + '-#infty')
            else:
                bhtex.SetText(x, y, title + '{:2.1f}'.format(log(tt)))
        bhtex.SetTextSize(tst)
        bhtex.SetNDC()
        print('BH 1D pval={} t={}'.format(results[ibest].p, results[ibest].t) )
        return results[ibest],bhtex
    else:
        return cBumpRes1D(-1, [], 1., 0.),ROOT.TObject()


##########################################
def GetNnonEmptyBisn(hdata):
    nx = hdata.GetXaxis().GetNbins()
    ny = hdata.GetYaxis().GetNbins()
    nnon = 0
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            if hdata.GetBinContent(i,j) > 0:
                nnon = nnon + 1
    return nnon
            
    
##########################################
def FindBestBumpHunter2DArea(hdata, hbg, x, y, title = 'log(t): ', tst = 0.045, showp = False, nNonEmpty = -1 ):
    nx = hdata.GetXaxis().GetNbins()
    ny = hdata.GetYaxis().GetNbins()

    # limit parameters of the search area:
    if nNonEmpty < 0:
        nNonEmpty = GetNnonEmptyBisn(hdata)
    nmax = min(nx*ny/4,nNonEmpty/4)
    nmin = min(4,nx*ny)
    wxmax = nx / 2
    wymax = ny / 2

    # initial window widths
    iwwx = 2
    iwwy = 2
    # shift initial window
    results = []
    for i in range(1,nx+1-iwwx):
        for j in range(1,ny+1-iwwy):
            result = MakeBestBumpHunterAreaAroundWindow(hdata, hbg, i,j, iwwx, iwwy, nmax, nmin, wxmax, wymax)
            results.append(result)
    ibest = FindBestResultIndex(results)
    if ibest >= 0:
        bhtex = ROOT.TLatex(x, y, '')
        if showp:
            bhtex.SetText(x, y, title + '{:1.2f}'.format(results[ibest].p))
        else:
            tt = results[ibest].t
            if tt > kBHInfty*kBHSFcmp:
                bhtex.SetText(x, y, title + '+#infty')
            elif tt < kBHZero:
                bhtex.SetText(x, y, title + '-#infty')
            else:
                #print('    t={:}'.format(results[ibest].t))
                bhscore = log(results[ibest].t)
                bhtex.SetText(x, y, title + '{:2.1f}'.format(bhscore))
        bhtex.SetTextSize(tst)
        bhtex.SetNDC()
        print('BH 2D pval={} t={}'.format(results[ibest].p, results[ibest].t) )
        return results[ibest],bhtex
    else:
        return cBumpRes2D(-1, -1, [], 1., 0.),ROOT.TObject()


    
##########################################
#             Significance 1D            #
##########################################
def FindMaxSignifBins1D(hdata, hbg):
    im = -1
    maxsignif = -999
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        tot = hdata.GetBinContent(ibin)
        bg = hbg.GetBinContent(ibin)
        edata = hdata.GetBinError(ibin)
        ebg = hbg.GetBinError(ibin)
        etot = sqrt( pow(edata, 2) + pow(ebg, 2) )
        if etot > 0.:
            signif = (tot - bg) / etot
            if signif > maxsignif:
                maxsignif = 1.*signif
                im = 1*ibin
    return im


##########################################
def BinNeighboursUsedBins1D(ibin, usedbins):
    for i in range(ibin-1,ibin+2):
        if i < 0 or i >= len(usedbins):
            continue
        # skip the tested bin itself
        if ibin == i:
            continue
        if usedbins[i] == 1:
            return True
    return False
            
##########################################
def ComputeExcessSignificance1DIteratively(hdata, hbg, x=0.13, y=0.66, title = '1D sig. signif.: ', thr = 2., tst = 0.045):
    sarea = 0.
    serrtotsq = 0.
    i0 = FindMaxSignifBins1D(hdata, hbg)
    usedbins = []
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        usedbins.append(0)
    usedbins[i0-1] = 1
    changed = True
    while changed:
        changed = False
        for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
            if usedbins[ibin-1] == 0 and BinNeighboursUsedBins1D(ibin-1, usedbins):
                tot = hdata.GetBinContent(ibin)
                bg = hbg.GetBinContent(ibin)
                edata = hdata.GetBinError(ibin)
                ebg = hbg.GetBinError(ibin)
                etotsq = pow(edata, 2) + pow(ebg, 2)
                if etotsq > 0. and tot-bg > 0:
                    signif = (tot - bg) / sqrt(etotsq)
                    if signif > thr:
                        usedbins[ibin-1] = 1
                        sarea = sarea + tot - bg
                        serrtotsq = serrtotsq + etotsq
                        changed = True

    # full significance
    nused = 0
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        if usedbins[ibin-1]:
            nused = nused + 1
    fsignif = 0.
    if serrtotsq > 0.:
        fsignif = sarea / sqrt(serrtotsq)

    print('{} used bins for significance: {}. Signif={}'.format(hdata.GetName(), nused, fsignif))
    stex = ROOT.TLatex(x, y, '')
    stex.SetText(x, y, title + '{:2.1f}'.format(fsignif))
    stex.SetTextSize(tst)
    stex.SetNDC()
    return fsignif, stex, nused, usedbins, i0


##########################################
def ComputeExcessSignificance1D(hdata, hbg, x=0.13, y=0.66, title = 'Sig. signif.: ', sumSignif = True):
    sarea = 0.
    errtotsq = 0.
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        tot = hdata.GetBinContent(ibin)
        bg = hbg.GetBinContent(ibin)
        edata = hdata.GetBinError(ibin)
        ebg = hbg.GetBinError(ibin)
        if sumSignif:
            if ebg > 0.:
                sarea = sarea + (tot-bg) / ebg
        else:
            sarea = sarea + tot - bg
            errtotsq = errtotsq + pow(edata, 2) + pow(ebg, 2)
    signif = 0.
    if sumSignif:
        signif = sarea
    else:
        if errtotsq > 0.:
            signif = sarea / sqrt(errtotsq)

    stex = ROOT.TLatex(x, y, '')
    stex.SetText(x, y, title + '{:2.1f}'.format(signif))
    stex.SetTextSize(0.05) # 0.055
    stex.SetNDC()
    return signif, stex




##########################################
def DrawUsedBinsLines(hdatax,hdatay,usedbins,i0, horiz, fflow = 0.07, fdy = 0.005, col = ROOT.kGreen, hcol = ROOT.kGreen+2, lst = 1, lw = 3, hlw = 4):
    xlines = []
    iline = -1
    lm = -1
    for i in range(0,len(usedbins)):
        if usedbins[i] == 1:
            iline = iline + 1
            x0 = hdatax.GetXaxis().GetBinLowEdge(i+1)
            wx = hdatax.GetXaxis().GetBinWidth(i+1)
            myxlines = []
            ymin = hdatay.GetXaxis().GetBinLowEdge(1)
            ymax = hdatay.GetXaxis().GetBinLowEdge(hdatay.GetXaxis().GetNbins()) + hdatay.GetXaxis().GetBinWidth(hdatay.GetXaxis().GetNbins())
            y0 =  ymin - fflow * (ymax - ymin)
            dy =  fdy * (ymax - ymin)
            if horiz:
                myxlines.append(ROOT.TLine(x0,y0,x0+wx,y0))
                # perp. edges:
                myxlines.append(ROOT.TLine(x0,y0-dy,x0,y0+dy))
                myxlines.append(ROOT.TLine(x0+wx,y0-dy,x0+wx,y0+dy))
            else:
                myxlines.append(ROOT.TLine(y0,x0,y0,x0+wx))
                # perp. edges:
                myxlines.append(ROOT.TLine(y0-dy,x0,y0+dy,x0))
                myxlines.append(ROOT.TLine(y0-dy,x0+wx,y0+dy,x0+wx))


            for xline in myxlines:
                xline.SetLineColor(hcol)
                xline.SetLineStyle(lst)
                if i+1 == i0:
                    xline.SetLineWidth(hlw)
                    xline.SetLineColor(col)
                    lm = iline
                else:
                    xline.SetLineWidth(lw)
                xline.Draw()
            xlines.append(myxlines)
    # redraw the box with the maximal significance:
    if lm >= 0 and lm < len(xlines):
        for xline in xlines[lm]:
            xline.Draw()
    return xlines


##########################################
#             Significance 2D            #
##########################################

def FindMaxSignifBins2D(hdata, hbg):
    im,jm = -1,-1
    maxsignif = -999
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
            tot = hdata.GetBinContent(ibin,jbin)
            bg = hbg.GetBinContent(ibin,jbin)
            edata = hdata.GetBinError(ibin,jbin)
            ebg = hbg.GetBinError(ibin,jbin)
            etot = sqrt( pow(edata, 2) + pow(ebg, 2) )
            if etot > 0.:
                signif = (tot - bg) / etot
                if signif > maxsignif:
                    maxsignif = 1.*signif
                    im = 1*ibin
                    jm = 1*jbin
    return im,jm


##########################################
def BinNeighboursUsedBins2D(ibin,jbin, usedbins):
    imin,imax = ibin-1,ibin+2
    for i in range(imin,imax):
        if i < 0 or i >= len(usedbins):
            continue
        jmin,jmax = jbin-1,jbin+2
        for j in range(jmin,jmax):
            if j < 0 or j >= len(usedbins[i]):
                continue
            # skip the tested bin itself
            if ibin == i and jbin == j:
                continue
            # skip diagonal neighbouring
            if abs(i - ibin) > 0 and abs(j - jbin) > 0:
                continue
            if usedbins[i][j] == 1:
                return True
    return False
            


##########################################
def ComputeExcessSignificance2DIteratively(hdata, hbg, x=0.13, y=0.66, title = '2D signif.: ', thr = 2., tst = 0.045):
    sarea = 0.
    serrtotsq = 0.
    i0,j0 = FindMaxSignifBins2D(hdata, hbg)
    usedbins = MakeEmptyUsedBins2D(hdata)
    usedbins[i0-1][j0-1] = 1
    changed = True
    while changed:
        changed = False
        for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
            for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
                if usedbins[ibin-1][jbin-1] == 0 and BinNeighboursUsedBins2D(ibin-1,jbin-1, usedbins):
                    tot = hdata.GetBinContent(ibin,jbin)
                    bg = hbg.GetBinContent(ibin,jbin)
                    edata = hdata.GetBinError(ibin,jbin)
                    ebg = hbg.GetBinError(ibin,jbin)
                    etotsq = pow(edata, 2) + pow(ebg, 2)
                    if etotsq > 0. and tot-bg > 0:
                        signif = (tot - bg) / sqrt(etotsq)
                        if signif > thr:
                            usedbins[ibin-1][jbin-1] = 1
                            sarea = sarea + tot - bg
                            serrtotsq = serrtotsq + etotsq
                            changed = True

    # full significance
    nused = 0
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
            if usedbins[ibin-1][jbin-1]:
                nused = nused + 1
    fsignif = 0.
    if serrtotsq > 0.:
        fsignif = sarea / sqrt(serrtotsq)

    print('{} used bins for significance: {}. Signif={}'.format(hdata.GetName(), nused, fsignif))
    stex = ROOT.TLatex(x, y, '')
    stex.SetText(x, y, title + '{:2.1f}'.format(fsignif))
    stex.SetTextSize(tst)
    stex.SetNDC()
    return fsignif, stex, nused, usedbins, i0, j0


##########################################
def ComputeExcessSignificance2D(hdata, hbg, x=0.13, y=0.66, title = '2D sig. signif.: ', sumSignif = True):
    sarea = 0.
    errtotsq = 0.
    for ibin in range(1, hdata.GetXaxis().GetNbins() + 1):
        for jbin in range(1, hdata.GetYaxis().GetNbins() + 1):
            tot = hdata.GetBinContent(ibin,jbin)
            bg = hbg.GetBinContent(ibin,jbin)
            edata = hdata.GetBinError(ibin,jbin)
            ebg = hbg.GetBinError(ibin,jbin)
            if sumSignif:
                if ebg > 0.:
                    sarea = sarea + (tot-bg) / ebg
            else:
                sarea = sarea + tot - bg
                errtotsq = errtotsq + pow(edata, 2) + pow(ebg, 2)
    signif = 0.
    if sumSignif:
        signif = sarea
    else:
        if errtotsq > 0.:
            signif = sarea / sqrt(errtotsq)
            
    stex = ROOT.TLatex(x, y, '')
    stex.SetText(x, y, title + '{:2.1f}'.format(signif))
    stex.SetTextSize(0.05) # 0.055
    stex.SetNDC()
    return signif, stex


##########################################
def DrawUsedBinsGrid(hdata,used2d,i0,j0, col = ROOT.kGreen, hcol = ROOT.kGreen + 2, lst = 1):
    xlines = []
    iline = -1
    lm = -1
    for i in range(0,len(used2d)):
        for j in range(0,len(used2d[i])):
            if used2d[i][j] == 1:
                iline = iline + 1
                x0 = hdata.GetXaxis().GetBinLowEdge(i+1)
                wx = hdata.GetXaxis().GetBinWidth(i+1)
                y0 = hdata.GetYaxis().GetBinLowEdge(j+1)
                wy = hdata.GetYaxis().GetBinWidth(j+1)
                myxlines = []
                myxlines.append(ROOT.TLine(x0,y0,x0+wx,y0))
                myxlines.append(ROOT.TLine(x0,y0,x0,y0+wy))
                myxlines.append(ROOT.TLine(x0,y0+wy,x0+wx,y0+wy))
                myxlines.append(ROOT.TLine(x0+wx,y0,x0+wx,y0+wy))
                for xline in myxlines:
                    xline.SetLineColor(hcol)
                    xline.SetLineStyle(lst)
                    if i+1 == i0 and j+1 == j0:
                        xline.SetLineWidth(4)
                        xline.SetLineColor(col)
                        lm = iline
                    else:
                        xline.SetLineWidth(3)
                    xline.Draw()
                xlines.append(myxlines)
    # redraw the box with the maximal significance:
    if lm >= 0 and lm < len(xlines):
        for xline in xlines[lm]:
            xline.Draw()
    return xlines

##########################################
