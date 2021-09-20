#!/bin/bash

./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/  pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV_ATLAS ejets Standard >& log_pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/  pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV_ATLAS ejets Standard >& log_pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV_ATLAS.txt &

./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS  ejets Standard >& log_pp2y02tt_y0_1000GeV_NLO_14TeV_ATLAS.txt &

./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS  ejets Standard >& log_pp2Wj_pTj1min50GeV_allhad_14TeV_ATLAS.txt &
./bin/AnalyzeBoosted_x /home/public/data/ttbar/data/qitek/DelphesOut/ pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS ejets Standard >& log_pp2Wbb_pTj1min50GeV_allhad_14TeV_ATLAS.txt &

