#!/usr/bin/env python
'''
estimate alpha rates in LiLS#1, LiLS#9 files
20170829
'''
import math
import sys
import random
import numpy

import datetime,os
import time
import ROOT
import graphUtils

import readLogFile

class nine20170829():
    def __init__(self):
        self.rootprefix = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        self.rfiles = [ ['/Users/djaffe/work/WaveDumpData/rootfiles/run01134_ts1499445083.root','0deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01137_ts1499451494.root','180deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01138_ts1499458474.root','90deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01139_ts1499460963.root','0deg']
                        ]
        self.gU = graphUtils.graphUtils()
        self.figdir = 'Figures/nine20170829/'
        self.rLF = readLogFile.readLogFile()
        return
    def makeFileList(self,Number=9):
        '''
        make file list of LiLS#X root files and a list of starting timestamps where X=Number
        '''
        allowedRuns = []#['run01222','run01228']
        bnList,TSlist = self.rLF.numberNine(Number=Number)
        #print 'nine20170829.makeFileList Number',Number,'bnList',bnList
        fnList,tsList = [],[]
        for bn,ts in zip(bnList,TSlist):
            Reject = len(allowedRuns)>0
            if len(allowedRuns)>0:
                for r in allowedRuns:
                    if r in bn: Reject = False

            if not Reject:
                fn = self.rootprefix + bn + '.root'
                if os.path.isfile(fn):
                    fnList.append(fn)
                    tsList.append(ts)
                    #print 'nine21070829.makeFileList ABORT AFTER ONE FILE FOUND'
                    #break
                
        return fnList,tsList
    def coagulateFileLists(self):
        '''
        fill 2 dicts
        first dict with names of root files to analyze
        second dict with starting timestamps
        '''


        self.LiLSFiles = {}
        self.LiLSstamps= {}
        for Number in [1,9]:
            self.LiLSFiles[Number],self.LiLSstamps[Number] = self.makeFileList(Number=Number)
            print 'nine20170829.coagulateFileLists',len(self.LiLSFiles[Number]),'files for sample',Number
        return
    def main(self):
        self.coagulateFileLists()
        self.loop()
        return
    def goodQandPSD(self,Q,PSD):
        '''
        rejects events in high rate regions. based on run1030 analysis
        '''
        return (Q>0.001 and PSD<0.6)
    def goodTime(self,runnum,t):
        '''
        rejects events in high instantaneous rate regions based on runnum(str) and abs_time
        '''
        if runnum=='run01030':
            return t>1. and abs(t-942.)>11. and abs(t-48129.5577105)>1. and abs(t-57056.1945872)>1.
        if runnum=='run01114':
            return t>1. and abs(t-3.573856e5)>1. and abs(t-619753.)>1. and abs(t-7.11004e5)>1.
        if runnum=='run01154':
            return t>1. and abs(t-1421.8)>1.
        if runnum=='run01068':
            return t>1. and abs(t-3667.88)>1.
        if runnum=='run01081':
            return t>1. and abs(t-15971.1659668)>1.

        return t>1.
    def loop(self):
        '''
        loop over trees in multiple root files, leave files open to make sure hists stay in scope
        '''
        debug = False
        cadence = 1
        if debug: cadence = 100000

        abs_time_cut = 0.004
        psdlo,psdhi,psdCut,psdlimit = 0., 1., 0.3, 0.6
        cpsdCut = ' PSD>{0:.2f} '.format(psdCut)
        Qlo,Qhi,Qcut = 0.0, 0.10, 0.01
        nx,ny = 100,100

        
        

        runStats = {} # {runnum: [livetime, timestamp, Tend]} 
        
        
        h0,h1,h2,h3 = {},{},{},{}
        bookHists=fillHists=True
        if bookHists:
            normedH = []
            #print self.LiLSFiles
            for Number in self.LiLSFiles:
                for fn in self.LiLSFiles[Number]:
                    runnum = fn.split('/')[-1].split('_')[0]
                    Title1 = 'Q '+runnum
                    Name1 = Title1.replace(' ','_')
                    Title1 += cpsdCut
                    H1 = ROOT.TH1D(Name1,Title1,nx,Qlo,Qhi)
                    Title0 = 'NoCuts PSD v Q '+runnum
                    Name0 = Title0.replace(' ','_')
                    H0 = ROOT.TH2D(Name0,Title0,nx,Qlo,Qhi,ny,psdlo,psdhi)
                    Title2 = 'BaseCuts PSD v Q '+runnum
                    Name2 = Title2.replace(' ','_')
                    H2 = ROOT.TH2D(Name2,Title2,nx,Qlo,Qhi,ny,psdlo,psdhi)
                    Title3 = 'PSD '+runnum
                    Name3 = Title3.replace(' ','_')
                    Title3 += ' {0:.2f} < Q < {1:.2f}'.format(Qcut,Qhi)
                    H3 = ROOT.TH1D(Name3,Title3,ny,psdlo,psdlimit)
                    for h in [H0,H1,H2,H3]: h.Sumw2()
                    h1[runnum] = H1
                    h2[runnum] = H2
                    h3[runnum] = H3
                    h0[runnum] = H0

        instRateThres = 100.
        debugLivetime = False

            
        for Number in self.LiLSFiles:
            for fn,ts in zip(self.LiLSFiles[Number],self.LiLSstamps[Number]):
                rf = ROOT.TFile.Open(fn,'r')
                runnum = fn.split('/')[-1].split('_')[0]
                print 'nine20170829.loop Start',runnum
                if fillHists:
                    H0 = h0[runnum]
                    H1 = h1[runnum]
                    H2 = h2[runnum]
                    H3 = h3[runnum]
                tree = rf.Get('tree')
                entries = tree.GetEntriesFast()
                t0,Nevt = None,0
                t1,N1   = None,100
                firstT,lastT = None,None
                validLivetime = [None]
                for jentry in xrange(0,entries+1,cadence):
                    ientry = tree.LoadTree(jentry)
                    if ientry<0:
                        break
                    nb = tree.GetEntry(jentry)
                    Q = tree.QtotalCh0
                    PSD=tree.psdCh0
                    abt=tree.abs_time

                    
                    if firstT is None: firstT = abt
                    lastT = abt
                    if t0 is None:
                        t0 = abt
                        t1 = abt
                    else:
                        dt = abt - t0
                    Nevt += 1
                    if abt!=t0: rate = float(Nevt)/dt
                    rate1 = None
                    if Nevt%N1==0:
                        dt1 = abt-t1
                        rate1 = float(N1)/dt1
                        t1 = abt

                        if abt>abs_time_cut and self.goodTime(runnum,abt):
                            if rate1>instRateThres:
                                txt = ''
                                if validLivetime[0] is not None: txt = 'Rate spike?'
                                if debugLivetime or validLivetime[0] is not None: print 'nine20170829.loop',txt,'Nevt,integralRate,instRate,t',Nevt,rate,rate1,abt
                            else:
                                if validLivetime[0] is None:
                                    validLivetime[0] = abt
                                    validLivetime.append(abt)
                                else:
                                    validLivetime[1] = abt
                    if debug: print 'abt,Q,PSD',abt,Q,PSD
                    if abt>abs_time_cut and self.goodTime(runnum,abt):
                        if fillHists:
                            H0.Fill(Q,PSD)
                            if self.goodQandPSD(Q,PSD):
                                H2.Fill(Q,PSD)
                                if PSD>psdCut and PSD<psdhi: H1.Fill(Q)
                                if Qcut<Q and Q<Qhi: H3.Fill(PSD)

                delT = 0.001
                if lastT-firstT>10000. : delT = 0.01
                if lastT-firstT>100000.: delT = 0.1
                f = 1./delT
                live,dead = 0.,0.
                for t in [x/f for x in range(int(f*firstT),int(f*lastT+delT))] :
                    if t>abs_time_cut and self.goodTime(runnum,t):
                        live += delT
                    else:
                        dead += delT
                totalTime = live+dead
                liveTime  = live
                ltratio = None
                if totalTime>0.: ltratio = liveTime/totalTime
                print 'nine20170829.loop finished {0} liveTime {1:.2f} totalTime {2:.2f} live/total {3:.5f} firstT {4:.4f} lastT {5:.4f}'.format(fn,liveTime,totalTime,ltratio,firstT,lastT)
                runStats[runnum] = [liveTime, ts, lastT]

                    
                                                 
        if fillHists:
            HDICTS = [h0,h1,h2,h3]
            for runnum in h1:
                liveTime,ts,lastT = runStats[runnum]
                for h in HDICTS:
                    hist = h[runnum]
                    words = ' per sec Tstart={0} Tend={1:.2f}'.format(ts,lastT)
                    newhist = self.gU.normHist(hist,liveTime,titleSuffix=words)
                    normedH.append(newhist)
            
            cnow = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            newfn = self.figdir + 'nine_' + cnow + '.root'
            new0 = ROOT.TFile(newfn,'RECREATE')
            #print 'new0',new0,'h1,h2',h1,h2

            for runnum in h1:
                for h in HDICTS: new0.WriteTObject(h[runnum])
            for h in normedH:
                print h
                new0.WriteTObject(h)
            new0.Close()
            print 'wrote',len(h1)+len(h2)+len(normedH),'hists to',newfn

        if debug : print 'h1.values()',h1.values(),  ' [ h1.values() ]',[ h1.values() ]
        
        return
    def overlay(self,Number=9,Hname='PSD_runXXXXXN'):
        '''
        overlay all distributions consistent with histogram name Hname from vial LiLS#N N = Number
        '''
        rfn = 'Figures/nine20170829/nine_20170901_155149.root'
        figdir = 'Figures/nine20170829/'
        rf = ROOT.TFile(rfn,'r')
        print 'nine20170829.overlay Opened',rfn
        
        basenames, timestamps = self.rLF.numberNine(Number=Number)
        runnums = []
        for bn in basenames:
            runnums.append( bn.split('_')[0] )

        hists = []

        for runnum in runnums:
            name = Hname.replace('runXXXXX',runnum)
            h = rf.Get(name)
            if h:
                hists.append(h)
        if len(hists)==0: return # nothing to plot
        legendX = 0.1
        if Number==1: legendX = 0.5
        for setLogy in [True,False]:
            self.gU.drawMultiObjects( [hists] ,figdir=figdir,fname=Hname+str(Number),statOpt=0,setLogy=setLogy,Grid=True,changeColors=True,addLegend=not setLogy,legendColumns=1,legendX=legendX,legendY=0.45,legendDY=0.5,biggerLabels=False)
        return
if __name__ == '__main__' :
    OP = nine20170829()
    if sys.argv>1:
        Number = int(sys.argv[1])
        print 'Number',Number
        OP.overlay(Number=Number)
    else:
        OP.main()                
    
