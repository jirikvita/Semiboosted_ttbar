#!/bin/bash

for run in 07 17 ; do

  for i in `cd runs/ ; ls *run_${run}*.root` ; do
      j=`echo $i | sed "s|${run}|\*|"`
      n=`cd runs/ ; ls ${j} | wc -l`
      echo "${i} : ${n}"
  done
done

