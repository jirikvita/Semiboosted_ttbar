#!/bin/bash

# covid-19 times of 2020
# jk 27.7.2020, 29.10.2020, XI 2020, III/2021

##############################################
export MainSample=analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_all.root


##############################################
# getting files
./scripts/getRootFilesFromIga.sh

##############################################
# plotting

# draw jet phase spaces for the ttbar sample, cartoon:
root -l macros/DrawttJetPhaseSpaces.C+

# compare weighted ttbar samples
./python/CompareWeightedSamples.py

# cutflow
./python/CmpCutFlow.py

# topologies fraction
./python/ComputeTopoFractions.py

# tag effs
./python/PlotTagEffs.py 

# correlations between variables 6.8.2020
./python/PlotCompactCorrelations.py ${MainSample}

# unfolding corrections
./python/DrawCorrections.py

# compare replicas, 3.9.2020
./python/compareReplicas.py ${MainSample}

# make sure the files are latest! 
./python/fitJES.py         data/analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_allforJES.root
./python/plotJESClosure.py data/analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_forClosure_all.root

# JK 11.9.2021:
#./python/plotJESClosure.py analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_weighted_forClosure_all.root
#./python/plotJESClosure.py analyzed_histos_pp_2tj_allhad_NLO_ALL_14TeV_ATLAS_allforJES.root

# main stacking and shapes, also 1D and 2D significances!
./scripts/draw.sh


# getting pdfs from IGA:

# after drawing on IGA:
# cd python/
# ./Tar.sh
# cd ../
# cd tex/
# ./tar.sh
# cd ../


# locally:
# cd python/
# getTarred.sh
# cd ../


# paper draft and note pdf update and compilation
# cd ../tex/draft/semiboosted_allhad
# ./getAllPdf.sh
# cd tables/
# ./get.sh
# cd ../
# ./python/FindMaxBHinTables.py
# ./doall.sh




# depricated:
# unfolding on IGA:
# ./python/run_unfold.py 
# plot unfolded results locally:
# ./scripts/getUnfolded.sh
# ./python/draw_unfolded.py


