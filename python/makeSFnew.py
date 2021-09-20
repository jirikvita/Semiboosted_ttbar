#!/usr/in/python

from __future__ import print_function

from math import pow, log

from collections import OrderedDict

# jk 20.11.2020

samples = {'y0' : 256,
           'y0gamma300' : 1024,
           'y0gamma100' : 512,
           'zp1000' : 65536,
           'zp1500' : 131072,
           'zp1250' : 65536,
           'zp900' : 16384,
           'zp800' : 16384,
           'zp750' : 8192,
           'zp700' : 8192,
           'xdxd300' : 65536,
           'xdxd100' : 65536,
           'xdxd10' : 65536,
           
           
}

oldDivPowers = OrderedDict()
oldDivPowers['2B0S'] = 4
oldDivPowers['1B1S'] = 3
oldDivPowers['0B2S'] = 0

reftopo = '2B0S'
for sample in samples:
    globsf = samples[sample]

    refoldpow = oldDivPowers[reftopo]
    refSF = globsf / pow(2, refoldpow)
    print('{} BaseSF: {} New powers: '.format(sample, refSF), end='')    
    for topo in oldDivPowers:
        
        
        oldpow = oldDivPowers[topo]
        oldSF = globsf / pow(2, oldpow)

        newSF = oldSF/refSF

        newpow = log(newSF) / log(2)
        print(' {}'.format(newpow), end='')
    print('')
