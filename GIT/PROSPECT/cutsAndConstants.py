#!/usr/bin/env python
'''
define cuts and constants used for 227Ac setup and processing
20161224
'''
import math
import sys
import numpy
import time,os


class cutsAndConstants():
    def __init__(self):
        self.Po215halflife = 0.001781 # from NNDC
        self.Po215lifetime = self.Po215halflife/math.log(2)
        self.tOffset = 10.*self.Po215lifetime

        self.Qmax = 0.05
        self.psdCut = 0.35
        self.lifeRange = [1,2,3,5]
        self.lowChargeCut = 0.012
        self.promptChargeCut = [0.03,0.04]
        nsig = self.nsigmaCutPo215Peak = 3.0
        mean,sigma = 0.03839, 0.00265 # these values from fit to run74 data
        self.delayChargeCut  = [mean-nsig*sigma,mean+nsig*sigma]

        # derived quantities
        self.totalAlphas = 5.
        self.lowChargeCutEffy =  0.999998600339 # cut range(pC) 12.0 50.0, from detMC
        self.promptChargeCutEffy =  0.947929001184 # 219Rn cut range(pC) 30.0 40.0 from detMC
        
        print 'cutsAndConstants.__init__ Initialized'
        return
if __name__ == '__main__':
    cAC = cutsAndComments()
