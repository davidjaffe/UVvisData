#!/usr/bin/env python
'''
estimate size needed for patch panel
20170316
'''
import math
import sys


#import datetime
import os

#import ROOT
#import graphUtils
#import gfit


class patchpanel():
    def __init__(self):
        self.sep = 1.12  # separation between centers in inches
        #self.sep = 1.0
        
        return
    def go(self):
        sep = self.sep
        border = 2
        for pairs,close in zip([1,7],[10,20]):
        
            rows = [5,6,7,8] # range(2,8)
            cols = range(2,100)
            ncol = 2*pairs
            cellpcol = 11
            pmtpcell = 2
            sphv     = 2
            optical = 0 #6 * pairs
            PMTconnectors = ncol * cellpcol * pmtpcell *sphv
            connectors = optical + PMTconnectors
            print '\n PAIRS',pairs,'Connectors total',connectors,'for PMTs',PMTconnectors,'for optical',optical,'border(in)',border/2.*sep

            for n in rows:
                for m in cols:
                    t = m*n + (m-1)*(n-1)
                    w = (m+border)*sep
                    h = (n+border)*sep
                    area = w*h
                    if abs(t-connectors)<close and t>=connectors:
                        print 'rows,columns,total,ht(in),wid(in),area(in2)',n,m,t,h,w,area
        return

if __name__ == '__main__' :
    pp = patchpanel()
    pp.go()    
