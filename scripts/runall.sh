#!/bin/bash

# 4t incl.:
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 4t_incl_13TeV_ATLAS      ejets Standard >& log_4t_incl_13TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 4t_incl_test_14TeV_ATLAS ejets Standard >& log_4t_incl_test_14TeV_ATLAS.txt &

# 4t allhad
#./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 4t_allhad_14TeV_ATLAS ejets Standard >& log_4t_allhad_14TeV_ATLAS.txt &
# split:
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ "run_0*_4t_allhad_14TeV_ATLAS" ejets Standard >& log_run_0X_4t_allhad_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ "run_1*_4t_allhad_14TeV_ATLAS" ejets Standard >& log_run_1X_4t_allhad_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ "run_2*_4t_allhad_14TeV_ATLAS" ejets Standard >& log_run_2X_4t_allhad_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ "run_3*_4t_allhad_14TeV_ATLAS" ejets Standard >& log_run_3X_4t_allhad_14TeV_ATLAS.txt &

# QCD jets:
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 4j_50GeVcut_14TeV ejets Standard >& log_4j_50GeVcut_14TeV.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 4j_14TeV          ejets Standard >& log_4j_14TeV.txt &

# 2t allhad:
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 2t_allhad_14TeV_LO_ATLAS                ejets Standard >& log_2t_allhad_14TeV_LO_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 2tj_allhad_NLO_ATLAS                    ejets Standard >& log_2tj_allhad_NLO_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS ejets Standard >& log_2tj_allhad_NLO_ptheavy50GeV_14TeV_ATLAS.txt &

# 2t ljets

./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 2tj_ljets_NLO_14TeV_ATLAS              ejets Standard >& log_2tj_ljets_NLO_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ 2tj_ljets_NLO_14TeV_ptheavy50GeV_ATLAS ejets Standard >& log_2tj_ljets_NLO_14TeV_ptheavy50GeV_ATLAS.txt &

# Z':
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ zp_ttbarj_allhad_1500GeV_14TeV_ATLAS   ejets Standard >& log_zp_ttbarj_allhad_1500GeV_14TeV_ATLAS.txt
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ zp_ttbarj_allhad_1000GeV_14TeV_ATLAS   ejets Standard >& log_zp_ttbarj_allhad_1000GeV_14TeV_ATLAS.txt
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ zp_ttbarj_allhad_700GeV_14TeV_ATLAS   ejets Standard >& log_zp_ttbarj_allhad_700GeV_14TeV_ATLAS.txt

echo "DONE! Don't forget to run /.scripts/hadd.sh to merge results when all interactive jobs are done;)!"

