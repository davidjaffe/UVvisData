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

        Timestamps = tdict['ts']
        Runs       = tdict['run']

        graphs = {}
        TMG = {}

        yLimits = None
        for cEvts,dEvts in zip(['evts','evtsN'],['rms','rmsN']):
            if cEvts=='evtsN': yLimits = [399.,471.]
        
            for j,sn in enumerate(snList):
                Y,dY,X,T = [],[],[],[]
                for i,evt in enumerate(tdict[cEvts]):
                    if tdict['sn'][i]==sn and evt>0:
                        Y.append( evt )
                        dY.append( tdict[dEvts][i] )
                        X.append( Runs[i] )
                        T.append( Timestamps[i] )
                dX = [0. for x in X]

                sample = samples[sn]

                title = 'Total '+cEvts+' v run '+sample
                name = title.replace(' ','_').replace('#','')
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



        print 'snList',snList
        for ax in ['time','run']:
            sn1 = 2
            prefix = 'Total_evtsN_v_'+ax+'_LiLS'
            g1 = graphs[prefix+str(sn1)]
            for sn in snList:
                if sn!=sn1:
                    g2 = graphs[prefix+str(sn)]
                    ggname = prefix+str(sn1)+'_and_'+str(sn)
                    self.gU.drawMultiObjects([[g1,g2]],fname=ggname,figdir=self.figdir,abscissaIsTime=True,Grid=True,
                                             addLegend=True,statOpt=0,biggerLabels=False,debug=False)
        
        return
        
if __name__ == '__main__' :
    '''
    arguments
    '''
    sP = speedyPlot()
    fn="speedyTrees/tree_20170319_124642_615646.root"
    fn = 'speedyTrees/tree_20170324_181234_373694.root'
    sP.loop(fn=fn)
