#!/bin/bash

echo "----------------------------------------------"

for i in `ls | grep pp ` ; do
    if [ -d $i ] ; then
	echo $i
	for j in `ls ${i}/Events/run_*/tag_?_merged_xsecs.txt` ; do xs=`cat $j  | head -n 3 | tail -n 1 | awk '{print $2};'` ; echo $xs ; done | \
	    awk 'BEGIN{sum = 0.; n = 0};{n = n+1; sum = sum + $1};END{printf("%f pb / %i = %f pb\n", sum, n, sum/n)};'
    fi
    echo "----------------------------------------------"
done



