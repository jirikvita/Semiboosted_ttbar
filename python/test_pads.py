#!/usr/bin/python

from Tools import *

can = ROOT.TCanvas()
pads,pad_inset =  MakeMultiSubPads(can, [0.7, 0.15,0.15])

#print(update(5))


ROOT.gApplication.Run()
