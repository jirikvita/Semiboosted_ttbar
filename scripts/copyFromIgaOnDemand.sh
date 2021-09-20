#!/bin/bash

for item in $* ; do
    ~/bin/myget.py iga  /home/qitek/work/mg5analysis/boosted_gacr/ $item
done

