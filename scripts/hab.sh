#!/bin/bash

# note on additional running for hab

# 29th May 2021


./python/runIndividual_runs_qg.py 
./python/haddToHalves_gq.py 


# plot the gg/qg/qq fractions
./python/PlotPdfInfo.py analyzed_histos_pp_2tj_incl_NLO_nocuts_14TeV_ATLAS_all.root
