#!/bin/bash

# jk 4.6.2021

for tag in ptj1j2min60_ptj1j2max200 ptj1j2min200_14TeV ptj1min200_ptj2min60max200 ; do
  echo "=== $tag ==="
  orig=`ls _run_pp*.sh | grep $tag | grep weigh`
  echo $orig
  script=_to_re_run_${tag}.sh
  rm $script
  for i in `cat weighted_runs_to_rerun_${tag}.sh` ; do
      cat $orig | grep $i >> $script
  done
done
