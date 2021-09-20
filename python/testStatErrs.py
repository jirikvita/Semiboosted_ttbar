#!/usr/bin/python3

from math import *

data = [
    ['Zp\rightarrow t\bar{t}, m_{Zp} = 1 \,\mathrm{TeV} (\times37500*2^{2.0})' , 2300 , 48 , 293300 , 6200 ,  2.1],  
    ['t\bar{t}, p_\mathrm{T}^{j1} \geq 200 \,\mathrm{GeV} p_\mathrm{T}^{j2} \in (60,200) \,\mathrm{GeV}' , 12300 , 110 , 252300 , 2300 ,  0.9],  
    ['t\bar{t}, p_\mathrm{T}^{j1,j2} \geq 200 \,\mathrm{GeV}' , 21100 , 150 , 105400 , 730 ,  0.7],  
    ['t\bar{t}, p_\mathrm{T}^{j1,j2} \in (60, 200) \,\mathrm{GeV}' , 1300 , 36 , 232700 , 6400 ,  2.8],  
    ['WWbb, p_\mathrm{T}^{j1,j2} > 60 \,\mathrm{GeV}' , 3300 , 57 , 411700 , 7200 ,  1.7],  
    ['Wbbjj, p_\mathrm{T}^{j1,j2} > 60 \,\mathrm{GeV}' , 250 , 16 , 25800 , 1600 ,  6.3],  
    ['bbjj, p_\mathrm{T}^{j1} \geq 200 \,\mathrm{GeV} p_\mathrm{T}^{j2} \in (60,200) \,\mathrm{GeV}' , 960 , 31 , 1420600 , 45900 ,  3.2],  
    ['bbjj, p_\mathrm{T}^{j1,j2} \geq 200 \,\mathrm{GeV}' , 5400 , 73 , 1419700 , 19400 ,  1.4],
    ]

sume2 = 0.
sumy = 0.
for d in data:
    print(d)
    y = d[-3]
    e = d[-2]
    sumy = sumy + y
    sume2 = sume2 + e*e

print(sumy, sqrt(sume2), 100*sqrt(sume2) / sumy)

