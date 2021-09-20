#!/bin/bash

## based on Peter Hamal's /home/public/navody/debugROOTMacro/setupAthenaGCC.sh
#PATH=/athena/gcc481/gcc-alt-481/x86_64-slc6-gcc48-opt/bin:$PATH
#export PATH
#C_INCLUDE_PATH=/athena/gcc481/gcc-alt-481/x86_64-slc6-gcc48-opt/include
#export C_INCLUDE_PATH
#CPLUS_INCLUDE_PATH=/athena/gcc481/gcc-alt-481/x86_64-slc6-gcc48-opt/include
#export CPLUS_INCLUDE_PATH
#LIBRARY_PATH=/athena/gcc481/gcc-alt-481/x86_64-slc6-gcc48-opt/lib:$LIBRARY_PATH
#export LIBRARY_PATH
#LD_LIBRARY_PATH=/athena/gcc481/gcc-alt-481/x86_64-slc6-gcc48-opt/lib64:$LIBRARY_PATH
#export LD_LIBRARY_PATH

#export delphespath=/home/qitek/install/delphes/

if [ `hostname` == "iga2011" ] ; then 
  #export LD_LIBRARY_PATH=/opt/ROOT/root-6.14.04-CERN/lib/:$LD_LIBRARY_PATH
  #mkdir -p /media/data/${USER}/DelphesOut
  ###delphespath=/home/qitek/work/delphes/
  #delphespath=/home/${USER}/workdir/delphes/Delphes-3.4.0
  delphespath=/home/${USER}/work/delphes/
  if [ ${USER} == "xzy" ] ; then
    delphespath=xyz
  fi

fi



export PYTHONPATH=`pwd`/python:$PYTHONPATH
export LD_LIBRARY_PATH=${delphespath}:$LD_LIBRARY_PATH

