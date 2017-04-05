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

import collections

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

        now = datetime.datetime.now()
        fmt = '%Y%m%d_%H%M%S_%f'
        self.start_time = cnow = now.strftime(fmt)

        self.outRootFile    =  'Speedy/Results/speedy_'+cnow+'.root'
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

        noPopUp = False
        
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
                if sn==2 and cEvts=='evtsN':
                    S2Y = numpy.array(Y)
                    S2DY= numpy.array(dY)
                    S2T = numpy.array(T)

                sample = samples[sn]

                title = 'Total '+cEvts+' v run '+sample
                name = title.replace(' ','_').replace('#','')
                if len(Y)<1: print 'speedyPlot.loop ERROR title',title,'NO ENTRIES'
                g = self.gU.makeTGraph(X,Y,title,name,ex=dX,ey=dY)
                self.gU.color(g,j,j,setMarkerColor=True)
                self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False,yLimits=yLimits,noPopUp=noPopUp)
                graphs[name] = g

                title = 'Total '+cEvts+' v time '+sample
                name = title.replace(' ','_').replace('#','')
                g = self.gU.makeTGraph(T,Y,title,name,ex=dX,ey=dY)
                self.gU.color(g,j,j,setMarkerColor=True)
                self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True,option="AP",yLimits=yLimits,noPopUp=noPopUp)
                graphs[name] = g

        # histograms of measured/expected, measured - expected, (meas-exp)/unc, (meas-meas[S2])
        
        for j,sn in enumerate(snList):
            M,E,D,T = [],[],[],[]
            for i,samnum in enumerate(tdict['sn']):
                if sn==samnum:
                    M.append( tdict['evtsN'][i]/nAlpha )
                    D.append( tdict['rmsN'][i]/nAlpha )
                    E.append( tdict['expect'][i] )
                    T.append( tdict['ts'][i] )
            M = numpy.array(M)
            D = numpy.array(D)
            E = numpy.array(E)
            T = numpy.array(T)
            dX= numpy.array([0. for x in E])
            Mby2,Dby2 = self.interp(T,S2T,S2Y,S2DY)
            MmM = M-Mby2
        
            name = 'meas_minus_S2_LiLS'+str(sn)
            title= 'Measured(LiLS'+str(sn)+') - Measured(LiLS2)'
            h = self.gU.makeTH1D(M-Mby2,title,name)
            hists.append( h )
            name += '_v_time'
            title += ' vs time'
            
            eY = numpy.sqrt(D*D+Dby2*Dby2)
            g = self.gU.makeTGraph(T,MmM,title,name,ex=dX,ey=eY)
            self.gU.color(g,j,j,setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True,option="AP",noPopUp=noPopUp)
            graphs[name] = g
            self.fitGraph(g,Draw=True)
            
            name = 'meas_minus_expect_LiLS'+str(sn)
            title = 'Measured - Expected LiLS#'+str(sn)
            h = self.gU.makeTH1D(M-E,title,name,xmi=-10.,xma=0.)
            hists.append( h )
            name = 'residual_meas_expect_LiLS'+str(sn)
            title = 'Residual(meas - expected) LiLS#'+str(sn)
            h = self.gU.makeTH1D((M-E)/D,title,name)
            hists.append( h )
            name = 'meas_by_expect_LiLS'+str(sn)
            title= 'Measured/Expected LiLS#'+str(sn)
            h = self.gU.makeTH1D(M/E,title,name,xmi=0.9,xma=1.0)
            hists.append( h )
        self.gU.drawMultiObjects(hists, fname='hists',figdir=self.figdir,abscissaIsTime=False,Grid=True,
                                 addLegend=False,statOpt=1110,biggerLabels=True,forceNX=4,noPopUp=noPopUp)
            
        

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
                                             addLegend=True,statOpt=0,biggerLabels=False,debug=False,noPopUp=noPopUp)

            expfix = 'Total_expect_v_'+ax+'_LiLS'
            for sn in snList:
                g1 = graphs[prefix+str(sn)]
                g2 = graphs[expfix+str(sn)]
                ggname = prefix+str(sn)+'_and_expected'
                self.gU.drawMultiObjects([[g1,g2]],fname=ggname,figdir=self.figdir,abscissaIsTime=(ax=='time'),Grid=True,
                                             addLegend=True,statOpt=0,biggerLabels=False,debug=False,noPopUp=noPopUp)

        
        f = ROOT.TFile(self.outRootFile,'RECREATE')
        for g in graphs: f.WriteTObject( graphs[g] )
        for h in hists: f.WriteTObject( h )
        f.Close()
        print 'speedyPlot.loop Wrote',len(graphs),'graphs and',len(hists),'hists to',self.outRootFile
        return
    def fitGraph(self,g,Draw=False):
        '''
        fit graph g with polynomial(s)
        optionally draw result
        '''
        x0 = g.GetXaxis().GetXmin()
        
        sp0 = "[0]"
        POL0= ROOT.TF1("POL0",sp0)
        POL0.SetParName(0,'Const')
        
        sp1 = "[0]+[1]*(x-"+str(x0)+")"
        POL1= ROOT.TF1("POL1",sp1)
        POL1.SetParName(0,'Offset')
        POL1.SetParName(1,'Slope')

        opt = 'NFEM'
        POL0.SetParameter(0,.1)
        g.Fit("POL0",opt)
        c,dc = POL0.GetParameter(0),POL0.GetParError(0)
        chi0,ndf0 = POL0.GetChisquare() , POL0.GetNDF()

        POL1.SetParameter(0,.1),POL1.SetParameter(0,-1.e-4)
        g.Fit("POL1",opt)
        off,doff = POL1.GetParameter(0),POL1.GetParError(0)
        slope,dslope = POL1.GetParameter(1),POL1.GetParError(1)
        chi1,ndf1 = POL1.GetChisquare() , POL1.GetNDF()

        if Draw:
            #ROOT.gROOT.ProcessLine("gROOT->SetBatch(0)")
            c2 = ROOT.TCanvas("c2","c2",1000,800)
            c2.Divide(1,1)
            c2.cd(1)
            xlo,xhi = g.GetXaxis().GetXmin(), g.GetXaxis().GetXmax()
            g.Draw("AIP")
            #raw_input('g.Draw')
            POL0.Draw("same")
            POL1.Draw("same")

            if 1:
                L0 = ROOT.TLine(xlo,c,xhi,c)
                L0.SetLineColor(ROOT.kRed)
                L0.SetLineStyle(2)
                L0.Draw("")
            #raw_input('L0.Draw')

            
                ylo,yhi = off+slope*(xlo-x0),off+slope*(xhi-x0)
                print 'off,slope,x0,xlo,xhi,ylo,yhi',off,slope,x0,xlo,xhi,ylo,yhi
                L1 = ROOT.TLine(xlo,ylo, xhi,yhi)
                L1.SetLineColor(ROOT.kBlue)
                L1.SetLineStyle(3)
                L1.Draw()
            #raw_input('L1.Draw')


            c2.Draw()
            #raw_input('c2.Draw')
            c2.cd()
            c2.Modified()
            c2.Update()
            
            wait = raw_input()
            if wait!='' : sys.exit('DONE')
            c2.IsA().Destructor(c2)
 #           del c2
        return

    def interp(self,t0,T,Y,dY):
        '''
        return estimate of Y and uncertainty at t0 give T,Y,dY 
        assume T is monotonically increasing
        '''
        debug = False

        if debug : print 'speedyPlot.interp t0',t0
        
        AY,AdY = numpy.array(Y), numpy.array(dY)
        AT = numpy.array(T)
        if isinstance(t0,collections.Container):
            y0,dy0 = [],[]
            for i,t in enumerate(t0):
                if i==0 and debug: print 'speedyPlot.interp i,t',i,t
                q =  numpy.interp(t, AT, AY)
                dq=  numpy.interp(t,AT, AdY) 
                y0.append(  q ) 
                dy0.append( dq)
                if debug : print 'speedyPlot.interp t,q,dq',t,q,dq
            y0,dy0 = numpy.array(y0), numpy.array(dy0)
        else:
            y0 = numpy.interp(t0, AT, AY )
            dy0= numpy.interp(t0, AT, AdY)
            if debug : print 'speedyPlot.interp t0,y0,dy0',t0,y0,dy0
        return y0,dy0
        
        
if __name__ == '__main__' :
    '''
    arguments
    '''
    sP = speedyPlot()
    fn="speedyTrees/tree_20170319_124642_615646.root"
    fn = 'speedyTrees/tree_20170324_181234_373694.root'
    fn = "speedyTrees/tree_20170327_201637_890698.root"
    fn = 'speedyTrees/tree_20170328_131720_516281.root'
    fn = 'speedyTrees/tree_20170403_084215_030514.root'
    fn = 'speedyTrees/tree_20170403_123150_108072.root'
    sP.loop(fn=fn)
