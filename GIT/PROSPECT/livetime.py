#!/usr/bin/env python
'''
compare last time from root file with reported run time in log file
20170103
'''
import math
import sys
import numpy
import datetime,os

import readLogFile
import get_filepaths
import ROOT,graphUtils

class livetime():
    def __init__(self):
        self.procLog = 'Logfiles/procWaveDump/'
        self.daqLog  = '/Users/djaffe/work/WaveDumpData/logfiles/'
        self.cWD = get_filepaths.get_filepaths()
        self.rLF = readLogFile.readLogFile()
        self.gU  = graphUtils.graphUtils()
        self.procLogList = self.cWD.get_filepaths(self.procLog)
        self.daqLogList  = self.cWD.get_filepaths(self.daqLog)
        self.figdir = 'Figures/livetime/'
        
        print 'livetime Initialized'
        return
    def getProcLog(self,inputrn=None):
        debug = False
        rn = None
        if type(inputrn) is str: rn = inputrn
        if type(inputrn) is int: rn = 'run'+str(inputrn).zfill(5)
        if debug : print 'livetime.getProcLog inputrn',inputrn,'rn',rn
        l = []
        for fn in self.procLogList:
            if rn in fn: l.append(fn)
        if debug : print 'livetime.getProcLog l',l
        l.sort()
        if len(l)==0: return None
        return l[-1]
    def getLastTime(self,fn=None):
        '''
        get last time from log file from latest procWaveDump log file
        input can be either logfilename or run number as string or as integer
        '''
        if not os.path.isfile(fn) : fn = self.getProcLog(fn)
        if fn is None : return None
        f = open(fn,'r')
        lastTime = None
        #print 'livetime.getLastTime fn',fn
        for l in f:
            #print l[:-1]
            if 'Last time' in l:
                lastTime = float(l.split()[-1])
                break
        f.close()
        return lastTime
    def main(self):

        procLogList = self.procLogList
        daqLogList  = self.daqLogList

        # get latest daq log file
        newList = []
        L = len(daqLogList)
        for i,fn in enumerate(daqLogList):
            rn = os.path.basename(fn).split('_')[0]
            dups = [fn]
            for j in range(i+1,L):
                fnj = daqLogList[j]
                if rn in fnj: dups.append(fnj)
            if len(dups)==1:
                newList.append(fn)
            else:
                dups = sorted(dups)
                bestfn = dups[-1]
                if bestfn not in newList: newList.append(bestfn)
        print 'livetime.main len(daqLogList),len(newList)',len(daqLogList),len(newList)
        daqLogList = newList

        run,rtime,ltime,lot = [],[],[],[]
        for dfn in daqLogList:
            rn = os.path.basename(dfn).split('_')[0]
            runNum = int(rn.replace('run',''))
            pfn = None
            for fn in procLogList:
                if rn in fn:
                    pfn = fn
                    break
            if pfn is not None:
                timestamp,sources,sample,runtime = self.rLF.readFile(dfn)
                lastTime = self.getLastTime(pfn)
                #print runNum,runtime,lastTime
                run.append(runNum)
                rtime.append(runtime)
                ltime.append(lastTime)
                lot.append(lastTime/runtime)

        xmi = float(run[0])-0.5
        xma = float(run[-1])+0.5
        nx = int(xma-xmi+.001)
        hists = []
        name = 'lot'
        title = 'lastTime/runTime vs run                                                               .'
        h = ROOT.TH1D(name,title,nx,xmi,xma)
        for x,y in zip(run,lot):
            h.Fill(x,y)
        hists.append(h)
        name = 'rtime'
        title = 'runTime(s) vs run'
        h = ROOT.TH1D(name,title,nx,xmi,xma)
        for x,y in zip(run,rtime):
            h.Fill(x,y)
        hists.append(h)
        name = 'ltime'
        title = 'lastTime(s) vs run'
        h = ROOT.TH1D(name,title,nx,xmi,xma)
        for x,y in zip(run,ltime):
            h.Fill(x,y)
        hists.append(h)


        self.gU.drawMultiHists(hists,'livetime',figdir=self.figdir,statOpt=0,abscissaIsTime=False,biggerLabels=True,dopt='hist',Grid=True)

        fn = self.figdir+'live.root'
        rfn = ROOT.TFile.Open(fn,'RECREATE')
        for h in hists: rfn.WriteTObject(h)
        rfn.Close()
        print 'livetime.main Wrote',len(hists),'hists to',fn
        
        
        return

if __name__ == '__main__' :
    lt = livetime()
    lt.main()
 
