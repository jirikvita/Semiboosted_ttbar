#!/usr/bin/python3

# jk 25.2.2018, 7.3.2018, 3.4.2018, 11.6.2019, 10.6.2020

from ReadTools import *


from matplotlib import pyplot as plt

PlotByPyPlot = False

# following:
# http://nbviewer.jupyter.org/github/gerbaudo/fbu/blob/v0.0.2/tutorial.ipynb
# modified by jk 30.11.2017, 1.2.2018, 26.2.2018, 28.2.2018

import fbu

myfbu = fbu.PyFBU()
stuff.append(myfbu)

# increase defaults:
###myfbu.nMCMC = 1000000 # default!
myfbu.nMCMC = 5000000
#myfbu.nwalkers = 1000
#myfbu.nBurn = 500000
#myfbu.nThin = 20
#myfbu.verbose = False
myfbu.verbose = True

# The Regularization class allows to impose an additional prior on the unfolded parameters
#from fbu import Regularization
#myfbu.regularization = Regularization('Tikhonov',parameters=[{'refcurv':0.,'alpha':0.01}])
#myfbu.regularization = Regularization('Tikhonov',parameters=[{'refcurv':0.,'alpha':0.1},{'refcurv':0.2,'alpha':0.1}])
# Supply the input distribution to be unfolded as a 1-dimensional list for N bins, with each entry corresponding to the bin content.



####################################################################################################################################################
def Unfold(filename_pseudodata, filename_mc, dirname, respname, outfilename, woption, smear, nrebin, ValSFdown, ValSFup, ftag, ValAddSFLastFirst):
    
    print(woption, smear, nrebin, ValSFdown, ValSFup, ftag, ValAddSFLastFirst)
    
    # additional SF for first and last bin:
    if ValAddSFLastFirst > 0:
        ValSFLastFirst = ValSFdown
    
    rfile_pseudodata = ROOT.TFile(filename_pseudodata, 'read')
    h_response_pseudodata = rfile_pseudodata.Get(dirname + respname)

    rfile_mc = ROOT.TFile(filename_mc, 'read')
    h_response_mc = rfile_mc.Get(dirname + respname)

    #if 'ttbarMass' in respname:
    #    h_response_pseudodata = RemoveFirstFewBins(h_response_pseudodata, 4)
    #    h_response_mc = RemoveFirstFewBins(h_response_mc, 4)
    
    # REBIN!!!
    if nrebin > 1:
        h_response_pseudodata.Rebin2D(nrebin,nrebin)
        h_response_mc.Rebin2D(nrebin,nrebin)

    h_migra = NormalizeResponse(h_response_mc)

    h_reco = MakeProjectionYFromHisto(h_response_pseudodata)
    h_ptcl = MakeProjectionXFromHisto(h_response_pseudodata)

    # data histogram:
    # closure:
    # h_data = h_reco.Clone(h_reco.GetName() + '_clone')
    # smear optionally by Poisson
    h_data = SmearData(h_reco, smear)

    print ('h_reco:')
    print (h_reco)
    print ('h_ptcl:')
    print (h_ptcl)


    # closure:
    myfbu.data = MakeListFromHisto(h_reco)
    print ('data:')
    print (myfbu.data)

    print ('Migra:')
    print (h_migra)

    myfbu.response = MakeListResponse(h_migra)
    print (myfbu.response)

    # Define the boundaries of the hyperbox to be sampled for each bin.
    # optimize for bin contents of the truth distribution?
    #myfbu.lower = [0] * h_reco.GetNbinsX()
    #myfbu.upper = [nrebin*50.e3] * h_reco.GetNbinsX()

    # Jiri:
    # better this way; it turns out it is very important to bound the
    # number of unfolded events from below by sth larger than 0., otherwise
    # instabilities are observed!!!
    myfbu.lower = []
    myfbu.upper = []
    # sf to restrict the range of truth to be sampled around reco
    nbins = h_reco.GetXaxis().GetNbins()
    for i in range(1, nbins+1):

        # DEFAULT!
        val = h_reco.GetBinContent(i)

        # sf to restrict the range of truth to be sampled around truth
        #for i in range(1, h_reco.GetYaxis().GetNbins()+1):

        ### original lines:
        ###if ValSFup < 1.: ValSFup = 1./ValSFup
        ###if ValSFdown < 1.: ValSFdown = 1./ValSFdown

        # note that now valsw dow is the ratio of parton / detector in bins
        # and that valsfup is an additional relative SHIFT in SF!!
        # change names!!!

        # HACK!
        # hm, let's just search around ptcl val ?! :D JK 29.11.2019
        val = h_ptcl.GetBinContent(i)

        myfbu.lower.append(val/(1. + ValSFup))
        # myfbu.lower.append(0.)
        myfbu.upper.append(val*(1. + ValSFup))

        # Nov 2019:
        #val = h_reco.GetBinContent(i)
        # around the known ratio w.r.t. reco bin content?
        #myfbu.lower.append(val/(ValSFdown[i-1] + ValSFup))
        #myfbu.upper.append(val*(ValSFdown[i-1] + ValSFup))

        ### ?!?!?! myfbu.lower[0] = 0. 
        #if i < 2 or i > nbins - 1 and ValAddSFLastFirst > 0:
        #    #myfbu.lower[i-1] = 0. ### HACK!!! jk 18.9.2019 myfbu.lower[-1] / ValAddSFLastFirst
        #    myfbu.lower[i-1] = myfbu.lower[i-1]  / ValAddSFLastFirst
        #    myfbu.upper[i-1] = myfbu.upper[i-1] * ValAddSFLastFirst
        print('Boundaries: {}--{}'.format(myfbu.lower[-1], myfbu.upper[-1]))
            

    # Run the MCMC sampling (this step might take up to several minutes for a large number of bins).
    myfbu.run()

    # Retrieve the N-dimensional posterior distribution in the form of a list of N arrays.
    trace = myfbu.trace
    print (trace)
    print ('Length of the trace: {:}'.format(len(trace)) )


    outfile = ROOT.TFile(outfilename, woption)
    outfile.cd()

    posteriors_diag = MakeTH1Ds(trace)
    print ('Length of the diag: {:}'.format(len(posteriors_diag)) )
    canp = PlotPosterior(posteriors_diag, 'PosteriorsDiag_' + respname + ftag) # TODO: add a tag here!

    h_unfolded = MakeUnfoldedHisto(h_reco, posteriors_diag)

    if PlotByPyPlot:
        for ibin in range(0,h_reco.GetNbinsX()):
            plt.hist(trace[ibin],
                     bins=20,alpha=0.85,
                     normed=True)
            plt.ylabel('probability')
            plt.show()

    #
    #
    # more on plotting:
    # https://matplotlib.org/users/pyplot_tutorial.html
    #
    #

    h_ptcl.Write()
    h_reco.Write()
    h_response_mc.Write()
    h_migra.Write()
    h_data.Write()

    outfile.Write()


    #ROOT.gApplication.Run()

    return
