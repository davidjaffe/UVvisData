#!/usr/bin/env python
'''
process HV log files for 227Ac setup
20161229
'''
import math
import sys
#import random
import numpy
#import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
#from scipy.integrate import quad
#import scipy.stats
#import matplotlib
#import matplotlib.pyplot as plt
import datetime,os

import ROOT
import graphUtils
#import gfit
#import readLogFile
#import cutsAndConstants

class procHVlog():
    def __init__(self):
        self.filename = '/Users/djaffe/work/WaveDumpData/HVfiles/CAENGECO2020.log'
        self.figdir = 'Figures/procHVlog/'
        self.gU = graphUtils.graphUtils()
        return
    def process(self):
        '''
        process data from log file
        return timestamp, voltage

        snippet of file
[2016-12-15T12:08:13]: [Geco] bd [0] ch [3] par [IMonH] val [493.05]; 
[2016-12-15T12:08:15]: [Geco] bd [0] ch [3] par [VMon] val [1499.5]; 
[2016-12-15T12:08:15]: [Geco] bd [0] ch [3] par [VMon] val [1499.3]; 
[2016-12-15T12:08:25]: [Geco] bd [0] ch [3] par [VMon] val [1499.5]; 
[2016-12-15T12:08:25]: [Geco] bd [0] ch [3] par [VMon] val [1499.3]; 
[2016-12-15T12:08:39]: [Geco] bd [0] ch [3] par [VMon] val [1499.5]; 
[2016-12-15T12:08:40]: [Geco] bd [0] ch [3] par [VMon] val [1499.3]; 
[2016-12-15T12:08:42]: [Geco] bd [0] ch [3] par [VMon] val [1499.5]; 
        '''
        if not os.path.isfile(self.filename):
            sys.exit( 'procHVlog.process INVALID FILE NAME '+self.filename)
        f = open(self.filename,'r')
        print 'procHVlog.process Reading',self.filename
        voltage = []
        current = []
        ON = True
        for l in f:
            if 'ch [3]' in l:
                line = l.replace(']','').replace('[','').replace(':','').replace(';','')
                s = line.split()
                if 'Status' in line:
                    if int(s[-1])==8197 : ON = False
                    if int(s[-1])==1    : ON = True
                if 'VMon' or 'IMonH' in line:
                    ts = s[0]
                    dt = self.ts2dt(ts)
                    d  = float(s[-1])
                    if s[7]=='VMon':
                        voltage.append( (dt,d, ON) )
                        #if ON and d<1400.: print l[:-1]
                    if s[7]=='IMonH': current.append( (dt,d, ON) )
        f.close()
        hists,graphs = [],[]
        print 'len(voltage),len(current)',len(voltage),len(current)
        print 'voltage[:3]',voltage[:3]
        for Q,W in zip([voltage,current],['Voltage','Current']):
            U,V = [],[]
            for triple in Q:
                u,v,ON = triple

                U.append(u)
                V.append(v)
            name = W[0].lower()+'t'
            title = W+' vs time'
            g = self.gU.makeTGraph(U,V,title,name)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True)
            graphs.append(g)
            m,s = numpy.average(V),numpy.std(V)
            #print 'W,mean,stddev',W,m,s
            h = self.gU.makeTH1D(V,W,W,nx=100) # automatically choose limits
            hists.append(h)
            if W=='Voltage': hh = self.gU.makeTH1D(V,W+'m',W+'m',nx=100,xmi=1490.,xma=1510.)
            if W=='Current': hh = self.gU.makeTH1D(V,W+'m',W+'m',nx=100,xmi=480.,xma=500.)
            hists.append(hh)
        self.gU.drawMultiHists(hists,fname='hv',figdir=self.figdir,setLogy=True,biggerLabels=False)
        outfn = self.figdir + 'hv.root'
        outrf = ROOT.TFile(outfn,'RECREATE')
        for h in hists: outrf.WriteTObject(h)
        for g in graphs:outrf.WriteTObject(g)
        print 'procWaveDump.main Wrote',len(hists),'hists and',len(graphs),'graphs to',outfn
        outrf.Close()
        
        return
    
    def ts2dt(self,ts):
        '''
        return TDatime object for text timestamp ts
        '''
        fmt = '%Y-%m-%dT%H%M%S'
        dt = self.gU.getTDatime(ts,fmt=fmt)
        return dt
if __name__ == '__main__' :
    pHV = procHVlog()
    pHV.process()
