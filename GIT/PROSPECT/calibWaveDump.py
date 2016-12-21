#!/usr/bin/env python
'''
calibrate wavedump histogram files for 227Ac setup
20161220
'''
import math
import sys
import random
import numpy
import scipy
from scipy.stats.mstats import chisquare
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import datetime,os

import ROOT
import graphUtils
import gfit

class calibWaveDump():
    def __init__(self):


        self.gU = graphUtils.graphUtils()
        self.gfit = gfit.gfit()

        self.rootFileDir = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        self.logFileDir  = self.rootFileDir.replace('rootfiles','logfiles')
        self.outFileDir  = 'Output/'
        self.figdir      = 'Figures/calibWaveDump/'
        return
    def subAndFit(self,hs,hb,dname,nsig=3.):
        '''
        return fitted difference of signal hs and background hb histograms with gaussian
        '''
        hd = hs.Clone(dname)
        hd.Add(hb,-1.)
        

        GoodFit,mean,emean, sgm, prob = self.gfit.fit(hd,nsigplus=nsig,nsigminus=nsig,debug=False,quiet=True)
        xmi,xma = mean-nsig*sgm,mean+nsig*sgm
        S,E = self.getIntegral(hd,xmi,xma)
        return hd,mean,emean,sgm,S,E
    def getIntegral(self,hd,xmi,xma):
        '''
        Return integral and uncertainty of 1d hist hd in range (xmi,xma)
        Use fractional contents and error for first,last bin in range
        '''
        i1 = hd.GetXaxis().FindBin(xmi)
        i2 = hd.GetXaxis().FindBin(xma)
        nx = hd.GetNbinsX()
        err = ROOT.Double(-1.)
        S = hd.IntegralAndError(i1+1,i2-1,err)
        err = err*err
        for i,x in zip([i1,i2],[xmi,xma]):
            xc = hd.GetBinCenter(i)
            c  = hd.GetBinContent(i)
            e  = hd.GetBinError(i)
            w  = hd.GetXaxis().GetBinWidth(i)
            if i==i1: f = (xc+w/2.-x)/w
            if i==i2: f = (x-(xc-w/2.))/w
            S += f*c
            err += (f*e)*(f*e)
        err = math.sqrt(err)
        #print 'S,err',S,err
        return S,err
    def loop(self,fn='Output/run00073.root'):
        '''
        loop over relevant hists in input root file and fit them,
        or whatever
        '''
        if not os.path.isfile(fn): return
        hf = ROOT.TFile.Open(fn)
        ROOT.gStyle.SetOptFit(1111)
        ROOT.gStyle.SetStatX(.5)
        hlist = []
        for l in range(1,6):
            clife = str(l)
            life = float(l)
            sname = 'dc'+clife
            bname = 'oc'+clife
            if hf.FindKey(sname):
                hs = hf.Get(sname)
                hb = hf.Get(bname)
                dname = 'd-o'+clife
                hd,mean,emean,sgm,S,E = self.subAndFit(hs,hb,dname)
                hlist.append(hd)
                f = 1./(1.-math.exp(-life))
                Scorr = S*f
                Ecorr = E*f
                print fn,'l,mean,emean,sgm,Scorr,Ecorr {0} {1:.5f} {2:.5f} {3:.5f} {4:.1f} {5:.1f}'.format(l,mean,emean,sgm,Scorr,Ecorr)


        dn = (fn.split('/')[1]).split('.')[0]
        self.gU.drawMultiHists(hlist,fname=dn,figdir=self.figdir,dopt='pe',statOpt=0)
        return
if __name__ == '__main__' :
    cWD = calibWaveDump()
    for i in range(31,90):
        fn='Output/run'+str(i).zfill(5)+'.root'
        cWD.loop(fn=fn)
