#!/bin/bash

echo "Renaming NEW files..."
for i in `ls *_NEW*.root` ; do
    j=`echo $i | sed "s|_NEW||g"`
    echo $i $j
    mv $i $j
done
