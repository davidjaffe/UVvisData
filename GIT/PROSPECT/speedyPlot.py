#!/usr/bin/env python
'''
plot quantities form speedy tree made by calibWaveDump
20170319
'''
import math
import sys
import random
import numpy
import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import datetime,os

import ROOT
import graphUtils
import gfit
import readLogFile
import cutsAndConstants
import Logger
import livetime
import get_filepaths
import twod

import makeTTree

class speedyPlot():
    def __init__(self,Log=True,Speedy=False):


        self.gU = graphUtils.graphUtils()
        self.mTT= makeTTree.makeTTree()
        self.cAC= cutsAndConstants.cutsAndConstants()
        self.figdir = 'Figures/Speedy/'
        return
    def loop(self,fn=None):
        '''
        main processing loop
        '''
        tdict = self.mTT.getTTree(ttName='sWD',fn=fn)

        snList = [] # list of sample numbers
        samples= {}
        for sn,sample in zip(tdict['sn'],tdict['sample']):
            if sn>0 and sn not in snList:
                snList.append( sn )
                samples[sn] = sample
        snList.sort()

        # define useful local variables and calculate expected rates
        Timestamps = tdict['ts']
        Runs       = tdict['run']
        tdict['expect'] = []
        tdict['dexpect']= []
        tdict['mme'] = [] # measured - expected
        tdict['dmme']= [] # unc
        nAlpha = self.cAC.totalAlphas
        for I,sample in enumerate(tdict['sample']):
            ts = tdict['ts'][I]
            Ac227rate = self.cAC.expectAc227Rate(inputDay=int(ts/1000.),sampleName=sample)

            tdict['expect'].append( Ac227rate[0] )
            tdict['dexpect'].append( Ac227rate[1] )
            m = tdict['evtsN'][I]/nAlpha
            dm= tdict['rmsN'][I]/nAlpha
            if 0 : print 'speedyPlot.loop sample,ts,Ac227rate,d,dm',sample,ts,Ac227rate,m,dm
            tdict['mme'].append( m - Ac227rate[0] )
            tdict['dmme'].append( math.sqrt( dm*dm ) )#+ Ac227rate[1]*Ac227rate[1] ) )


        # collect, plot per-sample quantites
        graphs = {}
        hists = []
        TMG = {}

        for cEvts,dEvts in zip(['evts','evtsN','expect','mme',],['rms','rmsN','dexpect','dmme']):
            nA = 1.
            yLimits = None
            if cEvts=='evtsN':
                nA = nAlpha
                yLimits = [399./nAlpha,481./nAlpha]
        
            for j,sn in enumerate(snList):
                Y,dY,X,T = [],[],[],[]
                for i,evt in enumerate(tdict[cEvts]):
                    if tdict['sn'][i]==sn:
                        Y.append( evt/nA )
                        dY.append( tdict[dEvts][i]/nA )
                        X.append( Runs[i] )
                        T.append( Timestamps[i] )
                dX = [0. for x in X]

                sample = samples[sn]

                title = 'Total '+cEvts+' v run '+sample
                name = title.replace(' ','_').replace('#','')
                if len(Y)<1: print 'speedyPlot.loop ERROR title',title,'NO ENTRIES'
                g = self.gU.makeTGraph(X,Y,title,name,ex=dX,ey=dY)
                self.gU.color(g,j,j,setMarkerColor=True)
                self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False,yLimits=yLimits)
                graphs[name] = g

                title = 'Total '+cEvts+' v time '+sample
                name = title.replace(' ','_').replace('#','')
                g = self.gU.makeTGraph(T,Y,title,name,ex=dX,ey=dY)
                self.gU.color(g,j,j,setMarkerColor=True)
                self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True,option="AP",yLimits=yLimits)
                graphs[name] = g

        # histograms of measured/expected & measured - expected
        for sn in snList:
            M,E = [],[]
            for i,samnum in enumerate(tdict['sn']):
                if sn==samnum:
                    M.append( tdict['evtsN'][i]/nAlpha )
                    E.append( tdict['expect'][i] )
            M = numpy.array(M)
            E = numpy.array(E)
            name = 'meas_minus_expect_LiLS'+str(sn)
            title = 'Measured - Expected LiLS#'+str(sn)
            h = self.gU.makeTH1D(M-E,title,name,xmi=-10.,xma=0.)
            hists.append( h ) 
            name = 'meas_by_expect_LiLS'+str(sn)
            title= 'Measured/Expected LiLS#'+str(sn)
            h = self.gU.makeTH1D(M/E,title,name,xmi=0.9,xma=1.0)
            hists.append( h )
        self.gU.drawMultiObjects(hists, fname='hists',figdir=self.figdir,abscissaIsTime=False,Grid=True,
                                 addLegend=False,statOpt=1110,biggerLabels=True)
            
        

        #print 'snList',snList
        for ax in ['time','run']:
            sn1 = 2
            prefix = 'Total_evtsN_v_'+ax+'_LiLS'
            g1 = graphs[prefix+str(sn1)]
            for sn in snList:
                if sn!=sn1:
                    g2 = graphs[prefix+str(sn)]
                    ggname = prefix+str(sn1)+'_and_'+str(sn)
                    self.gU.drawMultiObjects([[g1,g2]],fname=ggname,figdir=self.figdir,abscissaIsTime=(ax=='time'),Grid=True,
                                             addLegend=True,statOpt=0,biggerLabels=False,debug=False)

            expfix = 'Total_expect_v_'+ax+'_LiLS'
            for sn in snList:
                g1 = graphs[prefix+str(sn)]
                g2 = graphs[expfix+str(sn)]
                ggname = prefix+str(sn)+'_and_expected'
                self.gU.drawMultiObjects([[g1,g2]],fname=ggname,figdir=self.figdir,abscissaIsTime=(ax=='time'),Grid=True,
                                             addLegend=True,statOpt=0,biggerLabels=False,debug=False)

        return
        
if __name__ == '__main__' :
    '''
    arguments
    '''
    sP = speedyPlot()
    fn="speedyTrees/tree_20170319_124642_615646.root"
    fn = 'speedyTrees/tree_20170324_181234_373694.root'
    fn = "speedyTrees/tree_20170327_201637_890698.root"
    sP.loop(fn=fn)
