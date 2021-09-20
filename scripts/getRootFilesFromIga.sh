#!/bin/bash

# April 2020:

~/bin/myget.py  iga /home/qitek/work/mg5analysis/boosted_gacr/ "*half*.root"

~/bin/myget.py  iga /home/qitek/work/mg5analysis/boosted_gacr/ "*_all.root"

# just pseudodata:
# ~/bin/myget.py  iga /home/qitek/work/mg5analysis/boosted_gacr/ "*pseudo*.root"




####################################
# obsolete:

#rm analyzed_histos_4j_all_14TeV.root
#hadd analyzed_histos_4j_all_14TeV.root analyzed_histos_4j_14TeV.root analyzed_histos_4j_50GeVcut_14TeV.root

#rm analyzed_histos_zp_ttbarj_allhad_allM_14TeV_ATLAS.root
#hadd analyzed_histos_zp_ttbarj_allhad_allM_14TeV_ATLAS.root analyzed_histos_zp_ttbarj_allhad_*GeV_14TeV_ATLAS.root

#rm analyzed_histos_zp_ttbarj_allhad_AllSamples_14TeV_ATLAS.root
#hadd analyzed_histos_zp_ttbarj_allhad_AllSamples_14TeV_ATLAS.root analyzed_histos_zp_ttbarj_allhad_allM_14TeV_ATLAS.root analyzed_histos_4j_50GeVcut_14TeV.root

# 29.10.2020
# rename NEW files
./scripts/renameNEW.sh

