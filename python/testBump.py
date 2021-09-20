#!/usr/bin/python

from __future__ import print_function
import ROOT
from BumpSignifTools import *

hdata = ROOT.TH2D("a", "a", 10, 0, 10, 7, 0, 7)   

usedbins = MakeEmptyUsedBins(hdata)               

usedbins = MarkAsUsedWindow(usedbins, 3, 2, 2, 2)
print(usedbins)
print(CountUsedBins(usedbins))
usedbins[3][3] = 1           
usedbins[2][3] = 1
print(usedbins)
print(CountUsedBins(usedbins))

wx,wy = GetWidthXY(usedbins)
print(wx,wy)
