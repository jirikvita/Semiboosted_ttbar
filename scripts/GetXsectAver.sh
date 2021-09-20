#!/bin/bash


for i in `ls | grep pp ` ; do
    if [ -d $i ] ; then
	xs=`cat ${i}/crossx.html | grep run_ | grep result | cut -d \> -f 4 | cut -d \< -f 1 | awk 'BEGIN{n=0; sum=0.;};{sum = sum + $1; n=n+1;};END{printf("%f\n", sum/n);};'`
	evtsPy8=`cat ${i}/crossx.html | grep pythia8 | grep -v log | cut -d \> -f 3 | cut -d \< -f 1 | awk 'BEGIN{};{printf("%d, ", $1);};END{printf("\n")};'`
	evtsPart=`cat ${i}/crossx.html | grep parton\ mad | grep -v log | cut -d \> -f 3 | cut -d \< -f 1 | awk 'BEGIN{};{printf("%d, ", $1);};END{printf("\n")};'`

	echo "----------------------------------------------"
	if [ "${evtsPart}" != "" ] ; then
    		echo "'${i}' : evtsPart: $evtsPart"
	fi
	echo "'$i' : $xs pb,   evtsPy8 : $evtsPy8 ${evtsPart}"
    fi
done



