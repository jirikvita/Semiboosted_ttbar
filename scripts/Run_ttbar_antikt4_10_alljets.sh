#!/bin/bash
# jk 2017--2019

EXPERIMENT=ATLAS
#EXPERIMENT=CMS

now=`date  '+%Y-%m-%d-%H-%M-%S'`

path=""
exepath=""
cardpath=../cards/delphes/

### WHICH SHOWER? ###
#Shower="parton"
#Shower="Py6"
Shower="Py8"

# Py8:
exe=DelphesHepMC
# DEFAULT:
lhebase=tag_1_pythia8_events
# HACK FOR SECOND ROUNDS!!!
#lhebase=tag_2_pythia8_events
# HACK FOR THIRD ROUNDS!!!
#lhebase=tag_3_pythia8_events
format=hepmc

if [ $Shower == "parton" ] ; then
  # no shower:)
  exe=DelphesLHEF
  lhebase=events
  format=lhe
fi 

if [ $Shower == "Py6" ] ; then
  exe=DelphesSTDHEP
  lhebase=tag_1_pythia_events
  format=hep
fi


#### !!!! ####
# force overwrite output from Delphes:
Force=0
# Your dangerous choice!:)
#Force=1


# generic paths:
exepath=/home/${USER}/workdir/mg5analysis/delphes/Delphes-3.4.0/
cardpath=/home/${USER}/workdir/mg5analysis/cards/delphes/
outdir=/home/public/data/ttbar/data/${USER}/DelphesOut/
mkdir -p ${outdir}
outtag=""

#######################################################
#######################################################
#######################################################

if [ `hostname` == "iga2011" ] ; then

  # Jirkovy cesty:
    if [ ${USER} == "qitek" ] ; then
     exepath=/home/qitek/work/delphes/
     cardpath=/home/qitek/work/mg5analysis/cards/delphes/

     # 4t:
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/4t_incl_test/Events/
     #outtag="_4t_incl_13TeV"

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/4t_incl_test_14TeV/Events/
     #outtag="_4t_incl_test_14TeV"

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/4t_allhad_14TeV/Events/
     #outtag="_4t_allhad_14TeV"

     # 4j:
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_4j_50GeVcut_14TeV/Events/
     #outtag="_4j_50GeVcut_14TeV"

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_4j_14TeV/Events/
     #outtag="_4j_14TeV"


     # 2tj:
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2t_allhad_14TeV_LO/Events/
     #outtag=_2t_allhad_14TeV_LO
     
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO/Events/
     #outtag=_2tj_allhad_NLO

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptheavy50GeVcut/Events/
     #outtag=_2tj_allhad_NLO_ptheavy50GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/2tj_ljets_NLO_14TeV/Events/
     #outtag=_2tj_ljets_NLO_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/2tj_ljets_NLO_14TeV_ptheavy50GeVcut/Events/
     #outtag=_2tj_ljets_NLO_14TeV_ptheavy50GeV
     
     # Z':
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbarj_allhad_1000GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_1000GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbarj_allhad_1000GeV_NEW_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_1000GeV_NEW_14TeV

     # 14.9.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbarj_allhad_1250GeV_NEW_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_1250GeV_NEW_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbarj_allhad_700GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_700GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/FromLaptop/pp_zp_ttbarj_allhad_1500GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_1500GeV_14TeV

     # additional events to Z' allhad 1250 GeV, 6.3.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbarj_allhad_1250GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_1250GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/FromLaptop/pp_zp_ttbarj_allhad_800GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_800GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/FromLaptop/pp_zp_ttbarj_allhad_900GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_900GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/FromLaptop/pp_zp_ttbarj_allhad_750GeV_14TeV/Events/
     #outtag=_zp_ttbarj_allhad_750GeV_14TeV
     

     # re-run 21.5.2021:
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_5_5/pp_zp_ttbarj_ljets_700GeV_13TeV/Events/
     #outtag=_zprime_700_ljets_Py8_ATLAS

     # rerun of mg2.5.5 @ 13 TeV!
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_5_5/pp2ttbarj_ljets_both_NLO/Events/
     # WAS: outtag=_ljets_ttj_nocuts_ADDITIONAL
     # re-run 21.5.2021:
     #outtag=_ljets_ttj_nocuts_ADDITIONAL

     # rerun of mg2.5.5 @ 13 TeV! ttj NLO
     # re-run 21.5.2021:
     #path=/home/public/data/ttbar/data/qitek/MG5_aMC_v2_5_5/pp2ttbarj_ljets_both_LO_v3/Events/
     #outtag=_ljets_ttj_nocuts

     # rerun of mg2.5.5 @ 13 TeV! tt LO
     # re-run 21.5.2021:
     #path=/home/public/data/ttbar/data/qitek/MG5_aMC_v2_5_5/pp2ttbar_ljets_both_LO/Events/
     #outtag=_ljets_tt_nocuts


     # ===============
     
     # 25.3.2020, rerun 18.6.2021
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2y02tt_y0_1000GeV_NLO/Events/
     #outtag=_pp2y02tt_y0_1000GeV_NLO_14TeV

     # 29.3.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2Wj_pTj1min50GeV_allhad_14TeV/Events/
     #outtag=_pp2Wj_pTj1min50GeV_allhad_14TeV

     #29.3.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2Wbb_pTj1min50GeV_allhad_14TeV/Events/
     #outtag=_pp2Wbb_pTj1min50GeV_allhad_14TeV

     # 25.3.2020, 13.6.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad/Events/
     #outtag=_pp2xdxdtt_y0_1000GeV_xd_10GeV_allhad_14TeV
     
     # 1.4.2020, 10.6.2020, 23.6.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad/Events/
     #outtag=_pp2xdxdtt_y0_1000GeV_xd_100GeV_allhad_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad/Events/
     #outtag=_pp2xdxdtt_y0_1000GeV_xd_300GeV_allhad_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2y02tt_y0_1000GeV_width100GeV_NLO/Events/
     #outtag=_pp2y02tt_y0_1000GeV_width100GeV_NLO_14TeV

     # OLD BASELINE
     # used also after 28.5.2021 for q/g fractions after enabling storage of PDF infor from pythia by Radek:)
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1min60_ptj2min60/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1min60_ptj2min60_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2y02tt_y0_1000GeV_width300GeV_NLO/Events/
     #outtag=_pp2y02tt_y0_1000GeV_width300GeV_NLO_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2jjjMatched_pTj1j2min_60GeV/Events/
     #outtag=_pp2jjjMatched_pTj1j2min_60GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2jjjjMatched_pTj1j2min_60GeV_14TeV/Events/
     #outtag=_pp2jjjjMatched_pTj1j2min_60GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2jjjj_pTj1min400GeV_pTj2min60GeV_14TeV/Events/
     #outtag=_pp2jjjj_pTj1min400GeV_pTj2min60GeV_14TeV

     # BAD: 30.4.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2jjjj_mjjmin400GeV_14TeV/Events/
     #outtag=_pp2jjjj_mjjmin400GeV_14TeV
     
     # 25.4.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV/Events/
     #outtag=_pp2wbbjjMatched_allhad_pTj1j2min_60GeV_14TeV

     # 20.6.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1min100_ptj2min100/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1min100_ptj2min100_14TeV

     # 23.6.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1min200_ptj2min200/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1min200_ptj2min200_14TeV
     
     # NEW BASELINE (a) 18.10.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1j2min60_ptj1j2max200_14TeV

     # NEW BASELINE (b) 18.10.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1j2min200/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1j2min200_14TeV

     # NEW BASELINE (c) 23.10.2020
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200/Events/
     #outtag=_pp_2tj_allhad_NLO_ptj1min200_ptj2min60max200_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_zp_ttbar_allhad_1000GeV_14TeV/Events/
     #outtag=_pp_zp_ttbar_allhad_1000GeV_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp2wwbb_allhad_pTj1j2min_60GeV_14TeV/Events/
     #outtag=_pp2wwbb_allhad_pTj1j2min_60GeV_14TeV

     # Radek, oth Z offshell?
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_h_zz_4mu_radek_sm-full/Events/
     #outtag=_pp_h_zz_4mu_radek_sm-full_14TeV
     
     # Higgs 4l, still bad
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_h_zz_4l_sm-full/Events/
     #outtag=_pp_h_zz_4l_sm-full_14TeV

     # Higgs, Loop, OK:
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_H_4l_NLO/Events/
     #outtag=_pp_H_4l_NLO_14TeV

     # 24.5.2021
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/gg2ttbarj_allhad_NLO/Events/
     #outtag=_gg2ttbarj_allhad_NLO

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/qg2ttbarj_allhad_NLO/Events/
     #outtag=_qg2ttbarj_allhad_NLO

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/qq2ttbarj_allhad_NLO/Events/
     #outtag=_qq2ttbarj_allhad_NLO

     # 28.5.2021
     # g/q pdf info enabled by trick by Radek

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_allhad_NLO_nocuts/Events/
     #outtag=_pp_2tj_allhad_NLO_nocuts_14TeV

     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2tj_incl_NLO_nocuts/Events/
     #outtag=_pp_2tj_incl_NLO_nocuts_14TeV

     #path=/home/public/Jirka/MGKarel/DPEpomslicesttbar/
     #outtag=

     # 29.6.2021: QCD samples
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_4j_ptj1min200_ptj2min60max200/Events/
     #outtag=_pp_4j_ptj1min200_ptj2min60max200_14TeV
     
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_4j_ptj1j2min200/Events/
     #outtag=_pp_4j_ptj1j2min200_14TeV
     
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_4j_ptj1j2min60_ptj1j2max200/Events/
     #outtag=_pp_4j_ptj1j2min60_ptj1j2max200_14TeV
     
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV/Events/
     #outtag=_pp_2b2j_LO_matched_ptj1j2min60_ptj1j2max200_14TeV
     
     #path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2b2j_LO_matched_ptj1j2min200_14TeV/Events/
     #outtag=_pp_2b2j_LO_matched_ptj1j2min200_14TeV

     path=/home/public/data/ttbar/dataNTFS/qitek/MG5_aMC_v2_6_4/pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV/Events/
     outtag=_pp_2b2j_LO_matched_ptj1min200_ptj2min60max200_14TeV

fi # JK
  
#######################################################
   if [ ${USER} == "pacalt" ] ; then
     exepath=/home/${USER}/workdir/delphes/Delphes-3.4.0
     cardpath=/home/${USER}/workdir/mg5analysis/cards/delphes
    #Pepovy pokusy
	path = /media/data/pacalt/MGout/ljets_ZP1000/Events/
	#path=/media/data/qitek/MG5_aMC_v2_3_3/pp2ttbarj_ljets_tbarlep_LO/Events/
	outtag="_jp_ljets_ZP1000"
    #Pepovy pokusy
    #path = /media/data/qitek/MG5_aMC_v2_3_3/pp2ttbarj_ljets_tlep_LO/Events/
    #outtag ="_jp_ljets_tlep"
    # 
   fi
 
#######################################################
   if [ ${USER} == "pebaron" ] ; then
     exepath=/home/${USER}/delphes/Delphes-3.4.0/
     cardpath=/home/${USER}/mg5analysis/cards/delphes/
     export LD_LIBRARY_PATH=${exepath}:$LD_LIBRARY_PATH
     export LD_LIBRARY_PATH=/opt/root/lib/:$LD_LIBRARY_PATH
     #Petrovy cesty path ukazuje na to co chci analyzovat
     path=/media/data/pebaron/MG5_aMC_v2_4_3/TT_zp700_semilepton5/Events/
     outtag="_pb_zp700_semilep"
   fi
fi
#######################################################



#######################################################
# remember current directory:
local=`pwd`

tmpdir=/media/data/${USER}/DelphesOut/tmp-${now}
DoCopy=1
if [ -w $path ] ; then
  echo "OK, we don't have to copy as we have write access to the directory with generated events!;)"
  DoCopy=0
  tmpdir=${path}/tmp-${now}
else 
  echo "Sigh, we will have to to copy as we do not have write access to the directory with generated events!:-("
fi

echo "Making temp diretory ${tmpdir}"
mkdir -p ${tmpdir}

# delphes card tag:
for ctag in _AtKt4_and_10 ; do
# for ctag in _PileUp_AtKt10 _PileUp_AtKt4 ; do
# for ctag in _AtKt4_and_10_ljetsUPGRADE ; do

echo "Running $ctag ... "
#for run in `cd  $path ; ls | egrep "run_14|run_15|run_16"` ; do
#for run in `cd  $path ; ls | egrep "run_32|run_33"` ; do
#for run in `cd  $path ; ls | egrep "run_34|run_35|run_36"` ; do
#    for run in `cd  $path ; ls | egrep "run_37|run_38|run_39"` ; do
# special hacks:
#for run in `cd $path ; ls | grep run_13` ; do
#for run in `cd $path ; ls | egrep "run_73"` ; do
#for run in `cd  $path ; ls | egrep -v "run_00|run_21|run_22|run_23|run_24"` ; do
### DEFAULT:
for run in `cd  $path ; ls | grep run_` ; do
    
    if [ -d ${path}/${run} ] ; then

    echo "*** Processing $path $run ***"

    rfile=${outdir}/out_boosted${ctag}_${run}${outtag}_${EXPERIMENT}.root
    echo "Will try to produce ${rfile}"


    if [ $Force == 1 ] ; then
      echo "Forced to remove output file ${rfile}!"
      rm -f $rfile
    fi
    if ! [ -f $rfile ] ; then

      dir=${tmpdir}/unzipped
      echo "Will unzip in ${dir}"
      mkdir -p ${dir}
      if [ ${DoCopy} == 1 ] ; then
	  echo "Copying generated zip file..."
	  cp ${path}/${run}/${lhebase}.${format}.gz ${dir}/
      fi
      cd ${dir}/
      unzipped=${path}/${run}/${lhebase}.${format}
      if [ ! -f $unzipped ] ; then
        echo "Unzipping into ${dir}/${lhebase}.${format}"
        zipped=${path}/${run}/${lhebase}.${format}.gz
        gunzip ${path}/${run}/${lhebase}.${format}.gz -c > ${dir}/${lhebase}.${format}
      else
        ln -s $unzipped ${dir}/
      fi
      cd $local
      cmd="${exepath}/${exe} ${cardpath}/delphes_card_${EXPERIMENT}${ctag}.tcl ${rfile} ${dir}/${lhebase}.${format}"
      echo "Running ${cmd}"
      ${cmd}
      echo "Cleaning..."
      rm ${dir}/*.${format}
      rmdir ${dir}
      rmdir ${tmpdir}

    else
      echo "NOT Running ${cmd} as out file $rfile already exist!"
    fi
    fi # directory with the run exists
  done
done
