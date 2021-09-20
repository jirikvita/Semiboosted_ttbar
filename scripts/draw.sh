#!/bin/bash

#./python/CmpPlotHistos_GenDet.py  lists/list_BS_sel.txt -b
#./python/XsectStack.py  lists/list_ZpDM.txt

#####################################################################
# shapes:
#./python/CmpPlotHistos_GenDet.py  lists/list_BS.txt -n -b -s

./python/CmpPlotHistos_GenDet.py  lists/list_shapes.txt -n -b          >& log_shapes.txt
./python/CmpPlotHistos_GenDet.py  lists/list_shapes.txt -n -b -m       >& log_shapes_m.txt
./python/CmpPlotHistos_GenDet.py lists/list_ttOnly.txt -n -b -m      >& log_shapes_tt_m.txt
./python/CmpPlotHistos_GenDet.py lists/list_ttOnly.txt -n -b         >& log_shapes_tt.txt

./python/CmpPlotHistos_GenDet.py lists/list_zpOnly.txt -n -b -m    >& log_shapes_zp.txt
./python/CmpPlotHistos_GenDet.py lists/list_zpOnly.txt -n -b       >& log_shapes_zp_m.txt


##########################################################################################################################################
# JK 26.4.2021
# the lists for the three main samples y0_1000GeV, zp_1000GeV and xdxdtt_100GeV contain
# additional SFs for enhancing the signal and then also individually in topologies, in this order,
# of 2B0S, 1B1S and 0B2S, e.g.
#
#   analyzed_histos_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS_half0.root
#   y_{0}#rightarrow t#bar{t}, m_{y_{0}} = 1 TeV, #Gamma_{y_{0}} = 10 GeV
#   12 0.5 1.5 4
#   m_{y_{0}} = 1 TeV, #Gamma_{y_{0}} = 10 GeV
#   y0_1000GeV
#
# where 12 is the overall SF for the y0 sample, and 0.5, 1.5 and 4 are additional SF for powers of 2 to scale
# the  2B0S, 1B1S and 0B2S topologies by a factor of 2^0.5 etc.
##########################################################################################################################################


#####################################################################
# pseudodata stacked cmp:

./python/XsectStack.py lists/list_zp_1500GeV.txt
#./python/XsectStack.py lists/list_zp_1250GeV.txt
./python/XsectStack.py lists/list_zp_1000GeV.txt
#./python/XsectStack.py lists/list_zp_900GeV.txt 
#./python/XsectStack.py lists/list_zp_800GeV.txt 
#./python/XsectStack.py lists/list_zp_750GeV.txt 
#./python/XsectStack.py lists/list_zp_700GeV.txt 

./python/XsectStack.py lists/list_y0_1000GeV.txt
./python/XsectStack.py lists/list_y0_1000GeV_width100GeV.txt
./python/XsectStack.py lists/list_y0_1000GeV_width300GeV.txt

./python/XsectStack.py lists/list_xdxdtt_10GeV.txt
./python/XsectStack.py lists/list_xdxdtt_100GeV.txt
./python/XsectStack.py lists/list_xdxdtt_300GeV.txt

#####################################################################
# 1D significance comparison over replicas!
# 3.9.2020

./python/signifXsectStackReplicas1D.py lists/list_zp_1500GeV.txt
#./python/signifXsectStackReplicas1D.py lists/list_zp_1250GeV.txt
./python/signifXsectStackReplicas1D.py lists/list_zp_1000GeV.txt
#./python/signifXsectStackReplicas1D.py lists/list_zp_900GeV.txt 
#./python/signifXsectStackReplicas1D.py lists/list_zp_800GeV.txt 
#./python/signifXsectStackReplicas1D.py lists/list_zp_750GeV.txt 
#./python/signifXsectStackReplicas1D.py lists/list_zp_700GeV.txt 

./python/signifXsectStackReplicas1D.py lists/list_y0_1000GeV.txt
./python/signifXsectStackReplicas1D.py lists/list_y0_1000GeV_width100GeV.txt
./python/signifXsectStackReplicas1D.py lists/list_y0_1000GeV_width300GeV.txt

./python/signifXsectStackReplicas1D.py lists/list_xdxdtt_10GeV.txt
./python/signifXsectStackReplicas1D.py lists/list_xdxdtt_100GeV.txt
./python/signifXsectStackReplicas1D.py lists/list_xdxdtt_300GeV.txt

#####################################################################
# 2D Xsect and BumpHunter log(t) and p-val; NOT over replicas!
# 11.9.2020

./python/XsectStack2D.py lists/list_zp_1500GeV.txt
#./python/XsectStack2D.py lists/list_zp_1250GeV.txt
./python/XsectStack2D.py lists/list_zp_1000GeV.txt
#./python/XsectStack2D.py lists/list_zp_900GeV.txt 
#./python/XsectStack2D.py lists/list_zp_800GeV.txt 
#./python/XsectStack2D.py lists/list_zp_750GeV.txt 
#./python/XsectStack2D.py lists/list_zp_700GeV.txt 

./python/XsectStack2D.py lists/list_y0_1000GeV.txt
./python/XsectStack2D.py lists/list_y0_1000GeV_width100GeV.txt
./python/XsectStack2D.py lists/list_y0_1000GeV_width300GeV.txt

./python/XsectStack2D.py lists/list_xdxdtt_10GeV.txt
./python/XsectStack2D.py lists/list_xdxdtt_100GeV.txt
./python/XsectStack2D.py lists/list_xdxdtt_300GeV.txt


#####################################################################
# pseudodata stacked cmp, with both Ratio and Significance sub-pads;)
# 15.9.2020
./python/XsectStackRatioSignif.py lists/list_zp_1500GeV.txt
#./python/XsectStackRatioSignif.py lists/list_zp_1250GeV.txt
./python/XsectStackRatioSignif.py lists/list_zp_1000GeV.txt
#./python/XsectStackRatioSignif.py lists/list_zp_900GeV.txt 
#./python/XsectStackRatioSignif.py lists/list_zp_800GeV.txt 
#./python/XsectStackRatioSignif.py lists/list_zp_750GeV.txt 
#./python/XsectStackRatioSignif.py lists/list_zp_700GeV.txt 

./python/XsectStackRatioSignif.py lists/list_y0_1000GeV.txt
./python/XsectStackRatioSignif.py lists/list_y0_1000GeV_width100GeV.txt
./python/XsectStackRatioSignif.py lists/list_y0_1000GeV_width300GeV.txt

./python/XsectStackRatioSignif.py lists/list_xdxdtt_10GeV.txt
./python/XsectStackRatioSignif.py lists/list_xdxdtt_100GeV.txt
./python/XsectStackRatioSignif.py lists/list_xdxdtt_300GeV.txt


# BumpHunter algorithm 1D and 2D
### ./python/runAllBH2.py
### ./scripts/jobStarter.sh --jobs 8 _run_all_BH2.sh

