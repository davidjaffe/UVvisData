#!/usr/bin/env python
'''
process wavedump root files for 227Ac setup
20161219
'''
import math
import sys
import random
import numpy
#import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
#import matplotlib
#import matplotlib.pyplot as 
import datetime,os

import ROOT
import graphUtils
import cutsAndConstants

class procWaveDump():
    def __init__(self):
        # input, output files
        self.rootFileDir = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        self.logFileDir  = self.rootFileDir.replace('rootfiles','logfiles')
        self.outFileDir  = 'Output/'
        
        # graphing using root
        self.gU = graphUtils.graphUtils()
        # manipulation of root tree contents
        self.eventDict = {}
        self.tree = None
        self.entries = None
        self.miniTree = None
        self.miniOrder = None
        self.miniOrderIndices = None
        self.usemOI = False

        # initialize cuts, constraints
        self.cAC = cutsAndConstants.cutsAndConstants()
        
        self.Po215halflife = self.cAC.Po215halflife 
        self.Po215lifetime = self.cAC.Po215lifetime 
        self.tOffset = self.cAC.tOffset 

        self.Qmax = self.cAC.Qmax 
        self.psdCut = self.cAC.psdCut 
        self.lifeRange = self.cAC.lifeRange 
        self.lowChargeCut = self.cAC.lowChargeCut 
        self.promptChargeCut = self.cAC.promptChargeCut 
        self.delayChargeCut  = self.cAC.delayChargeCut  

        # defined for each file
        self.runTime = None

        # dict of all histograms
        self.hists = {}
        
        name = 'procWaveDump'
        self.figdir = 'Figures/'+name+'/'
        makeSubDir = False
        if makeSubDir:
            now = datetime.datetime.now()
            fmt = '%Y%m%d_%H%M%S_%f'
            self.start_time = cnow = now.strftime(fmt)
            self.figdir = self.figdir + cnow+'/'
            if os.path.isdir(self.figdir):
                pass
            else:
                try:
                    os.makedirs(self.figdir)
                except IOError,e:
                    print name+'__init__',e
                else:
                    print name+'__init__ created',self.figdir
        return
    def main(self,inputfn='run00066_ts1481922242.root',maxE=99999999999):
        '''
        main routine
        Open input root file, get run time from corresponding log file
        Copy useful contents of root file tree into dict for quicker access
        Book histograms
        Loop over events and fill hists
        Normalize a copy of all hists to running time
        Write out hists
        '''
        fn = os.path.basename(inputfn)
        self.open(fn=fn)
        x = self.getRunTime(fn)

        self.preFill()

        self.book()
        self.eventLoop(maxE=maxE)
        
        self.normHistsByRunTime(fn)

        #self.rf.Close() # close input root file
        outfn = self.outFileDir + fn.split('_')[0]+'.root'
        outrf = ROOT.TFile(outfn,'RECREATE')
        for h in self.hists:
            #print h,self.hists[h].GetName()
            outrf.WriteTObject( self.hists[h] )
        print 'procWaveDump.main Wrote',len(self.hists),'hists to',outfn
        
        return
    def normHistsByRunTime(self,rfn=None):
        '''
        make copies of all hists normalized by run time in seconds
        '''
        runTime = self.getRunTime(rfn)
        print 'procWaveDump.normHistsByRunTime runTime',runTime
        if runTime>0:

            newhists = {}
            for x in self.hists:
                h = self.hists[x]
                cn = h.ClassName()
                if 'TH1' in cn or 'TH2' in cn:
                    title = h.GetTitle()
                    name  = h.GetName()
                    newname = name + 'N'
                    hnew = h.Clone(newname)
                    hnew.SetTitle(title + ' per second')
                    hnew.Scale(1./runTime)
                    newhists[newname] = hnew
                    #print 'newname',newname,'hnew.GetName()',hnew.GetName()
            for h in newhists: self.hists[h] = newhists[h]

        return
    def getRunTime(self,rfn=None):
        '''
        get actual run time in seconds from logfile given root file name rfn
        use global variable if it has already been initialized
        '''
        if self.runTime is not None: return self.runTime
            
        if rfn is None:
            print 'procWaveDump.getRunTime FAILED No input file given'
            return -1.
        lfn = rfn.replace('root','log')
        fulllfn = self.logFileDir + lfn
        f = open(fulllfn,'r')
        request, actual = -1., -1.
        for l in f:
            if 'Request run time' in l: request = float(l.split()[-1])
            if 'Actual run time'  in l: actual  = float(l.split()[-1])
        f.close()
        self.runTime = request 
        if actual>-1. :
            self.runTime = actual
        
        return self.runTime
    def eventLoop(self,maxE=1000,debug=False):
        '''
        loop over events, makes some plots
        '''
        freq = 1000
        Nmax = self.entries
        if maxE is not None: Nmax = min(Nmax,maxE)
        print 'procWaveDump.eventLoop Process',Nmax,'events'

        
        if self.usemOI: # use indices into array
            Iabs_time, IpsdCh0, IQtotalCh0, IgoodCh0 = self.miniOrderIndices
        else:   # use names into dict
            Iabs_time, IpsdCh0, IQtotalCh0, IgoodCh0 = self.miniOrder


        
        maxTimeWindow = float(max(self.lifeRange))*self.Po215lifetime
        
        # main event loop
        for event in range(Nmax):
            if debug : print 'procWaveDump.eventLoop start event',event,
            if event%freq==0 and not debug:
                print '\r',event,
                sys.stdout.flush()

            # get tree for current event. Require 'good' channel data and PSD on prompt.
            # don't look for prompt,delayed candidates too close to end of run
            d =  self.getEvent(event)
            self.hists['goodCh0'].Fill(d[IgoodCh0])
            if d[IgoodCh0]:
                if debug : print 'passed goodCh0',
                name = 'PSD_vs_Charge'
                self.hists[name].Fill(d[IQtotalCh0],d[IpsdCh0])
                if d[IQtotalCh0]>self.lowChargeCut:
                    self.hists['psdc'].Fill(d[IpsdCh0])

                name = 'Charge_no_cut'
                self.hists[name].Fill(d[IQtotalCh0])
                if d[IpsdCh0]>self.psdCut :
                    if debug : print 'passed psd',
                    name = 'Charge_PSD_cut'
                    self.hists[name].Fill(d[IQtotalCh0])

                    tP = d[Iabs_time] # prompt time
                    QP = d[IQtotalCh0]# prompt charge
                    tmax = tP + maxTimeWindow
                    if tmax<self.lastTime:
                        if debug : print 'QP,tP,tmax',QP,tP,tmax,
                        name = 'p'
                        self.hists[name].Fill(QP)

                        # look for delayed candidates.
                        # Apply 'goodCh0' and PSD cuts for delayed cand.
                        # Plot charge and multiplicity of delayed cands for different
                        # multiples of Po215 lifetime.
                        Delayed = {}
                        for l in self.lifeRange: Delayed[l] = []
                        e1,e2 = event-1,event+1 
                        for devent in range(event+1,self.entries):
                            e2 = max(e2,devent)
                            d2 = self.getEvent(devent)
                            tD = d2[Iabs_time]
                            if debug : print 'devt,tD',devent,tD,
                            if tD>tmax:
                                break
                            if d2[IgoodCh0] and d2[IpsdCh0]>self.psdCut:
                                QD = d2[IQtotalCh0]
                                for l in self.lifeRange:
                                    tend = tP + float(l)*self.Po215lifetime
                                    if tD < tend:
                                        Delayed[l].append([tD,QD])
                        if debug : print 'final devt,tD',devent,tD,
                        for l in self.lifeRange:
                            clife = str(l)
                            name = 'dm'+clife
                            self.hists[name].Fill(float(len(Delayed[l])))
                            name2 = 'pvd' + clife
                            name = 'd' + clife
                            dcname = 'dc'+clife
                            pcname = 'pc'+clife
                            for t,q in Delayed[l]:
                                self.hists[name].Fill(q)
                                self.hists[name2].Fill(q,QP)
                                if self.inside(q,self.delayChargeCut):
                                    self.hists[pcname].Fill(QP)
                                if self.inside(QP,self.promptChargeCut):
                                    self.hists[dcname].Fill(q)


                        # look for fake delayed candidates in randomly offset time window
                        needed,tries = True,0
                        while needed:
                            tries += 1
                            e0 = numpy.random.randint(0,self.entries)
                            #print 'tries,e0,e1,e2',tries,e0,e1,e2
                            if e0<e1 or e2<e0:
                                d2 = self.getEvent(e0)
                                tFake = tD = d2[Iabs_time]
                                after = tP+self.tOffset<tD and tD<self.lastTime-maxTimeWindow
                                before= tD<tP-maxTimeWindow
                                #print 'tries,tD,tP+self.tOffset,self.lastTime-maxTimeWindow,tP-maxTimeWindow,before,after',tries,tD,tP+self.tOffset,self.lastTime-maxTimeWindow,tP-maxTimeWindow,before,after
                                if before or after:
                                    name = 'tP_vs_tO'
                                    self.hists[name].Fill(tD,tP)
                                    Delayed = {}
                                    for l in self.lifeRange: Delayed[l] = []
                                    tmax = tFake+maxTimeWindow
                                    needed = False
                                    if before: elast = e1
                                    if after:  elast = self.entries
                                    for devent in range(e0+1,elast):
                                        d2 = self.getEvent(devent)
                                        tD = d2[Iabs_time]
                                        if tD>tmax:
                                            break
                                        if d2[IgoodCh0] and d2[IpsdCh0]>self.psdCut:
                                            QD = d2[IQtotalCh0]
                                            for l in self.lifeRange:
                                                tend = tFake + float(l)*self.Po215lifetime
                                                if tD < tend :
                                                    Delayed[l].append([tD,QD])
                                    for l in self.lifeRange:
                                        clife = str(l)
                                        name = 'om'+clife
                                        #print 'l,name',l,name
                                        self.hists[name].Fill(float(len(Delayed[l])))
                                        name2= 'pvo'+clife
                                        name = 'o' + clife
                                        ocname = 'oc'+clife
                                        for t,q in Delayed[l]:
                                            self.hists[name].Fill(q)
                                            self.hists[name2].Fill(q,QP)
                                            if self.inside(QP,self.promptChargeCut):
                                                self.hists[ocname].Fill(q)
                                                
            if debug : print '...done'
                    
        return
    def book(self):
        '''
        book histograms

        parameters for cuts and histogram limits have been moved to initialization
        
        '''

        ROOT.TH1.SetDefaultSumw2(True) # proper calculation of bin errors

        # calibration data
        nx,xmi,xma = 100,0.,self.Qmax
        for l in self.lifeRange:
            clife = str(l)
            title = 'Prompt with Delayed Charge, PSD, '+clife+'lifetime cuts'
            name = 'pc'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
            title = 'Delayed with Prompt Charge, PSD, '+clife+'lifetime cuts'
            name = 'dc'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
            title = 'Doffset with Prompt Charge, PSD, '+clife+'lifetime cuts'
            name = 'oc'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

        
        nx,xmi,xma = 100, 0., self.Qmax
        ny,ymi,yma = 100, 0., 1.
        title = 'PSD vs Charge'
        self.TH2D(title,nx,xmi,xma,ny,ymi,yma)

        nx,xmi,xma = 100, 0., 1.
        title = 'PSD Charge $>$'+str(self.lowChargeCut)
        name = 'psdc'
        self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

        nx,xmi,xma = 100,0.,self.Qmax
        title = 'Charge no cut'
        self.TH1D(title,nx,xmi,xma)
        title = 'Charge PSD cut'
        name = title.replace(' ','_')
        self.TH1D(title,nx,xmi,xma)

        nx = int(self.runTime)+1
        xmi,xma = 0.,float(nx)
        ny,ymi,yma = nx,xmi,xma
        title = 'tP vs tO'
        self.TH2D(title,nx,xmi,xma,ny,ymi,yma)

        nx,xmi,xma = 100, 0., self.Qmax
        title = 'Prompt PSD cut'
        name = 'p'
        self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        for ilife in self.lifeRange:
            clife = str(ilife)
            nx,xmi,xma = 100, 0., self.Qmax
            title = 'Delay PSD cut ' + clife + ' lifetimes'
            name = 'd'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
            name = 'o'+clife
            title = 'Doffset PSD cut ' + clife + ' lifetimes'
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

            ny,ymi,yma = nx,xmi,xma
            title = 'Prompt v Delay PSD cut ' + clife + ' lifetimes'
            name  = 'pvd'+clife
            self.hists[name] = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
            title = 'Prompt v Doffset PSD cut ' + clife + ' lifetimes'
            name  = 'pvo'+clife
            self.hists[name] = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
            
            nx,xmi,xma = 21,-0.5,20.5
            title = 'Delay PSD cut Mult ' + clife + ' lifetimes'
            name = 'dm'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
            title = 'Doffset PSD cut Mult ' + clife + ' lifetimes'
            name = 'om'+clife
            self.hists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
            #print 'name,clife,title',name,clife,title

            
        
        nx,xmi,xma = 2,-0.5,1.5
        name = title = 'goodCh0'
        self.TH1D(title,nx,xmi,xma)

        print 'procWaveDump.book Booked',len(self.hists),'histograms'
        
        return

    def TH1D(self,title,nx,xmi,xma):
        name = title.replace(' ','_')
        self.hists[name] =  ROOT.TH1D(name,title,nx,xmi,xma)
        return
    def TH2D(self,title,nx,xmi,xma,ny,ymi,yma):
        name = title.replace(' ','_')
        self.hists[name] =  ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
        return
    def inside(self,x,l):
        return l[0]<=x and x<=l[1]

    def getEvent(self,event):
        '''
        load tree variables into event dict OR NUMPY ARRAY given event number
        
        '''
        debug = False
        if debug : print 'procWaveDum.getEvent self.tree',self.tree,'self.entries',self.entries

        if event<0 or event>self.entries:
            print 'procWaveDump.getEvent Invalid event#',event
            return None

        if self.miniTree is not None:
            if self.usemOI:
                d = self.miniTree[event]
            else:
                d = {}
                for i,v in enumerate(self.miniOrder):
                    d[v] = self.miniTree[event][i]
            

        else:

            self.tree.GetEntry(event)

            d = {}
            d['event'] = event
            d['abs_time'] = self.tree.abs_time
            d['dt'] = self.tree.dt
            d['nSamples'] = self.tree.nSamples
            d['averageCh0'] = self.tree.averageCh0
            d['psdCh0'] = self.tree.psdCh0
            d['QtailCh0'] = self.tree.QtailCh0
            d['QtotalCh0'] = self.tree.QtotalCh0
            d['goodCh0'] = self.tree.goodCh0
            d['minCh0']  = self.tree.minCh0
            d['minIdxCh0'] = self.tree.minIdxCh0
            d['halfMinIdxCh0'] = self.tree.halfMinIdxCh0
            d['maxCh0'] = self.tree.maxCh0
            d['maxIdxCh0'] = self.tree.maxIdxCh0
            d['pedCh0'] = self.tree.pedCh0
            d['Ch0_wform'] = self.tree.Ch0_wform

        return d
    def preFill(self):
        '''
        fill global array with useful variables from entire tree
        '''

        print 'procWaveDump.preFill Initializing for',self.entries,'events',   

        self.miniOrder = ['abs_time', 'psdCh0', 'QtotalCh0', 'goodCh0']
        
        Iabs_time = self.miniOrder.index('abs_time')
        IpsdCh0   = self.miniOrder.index('psdCh0')
        IQtotalCh0= self.miniOrder.index('QtotalCh0')
        IgoodCh0  = self.miniOrder.index('goodCh0')
        self.miniOrderIndices = [Iabs_time, IpsdCh0, IQtotalCh0, IgoodCh0]
        self.usemOI = True
        
        a = []
        for event in range(self.entries):
            self.tree.GetEntry(event)
            a.append( [self.tree.abs_time, self.tree.psdCh0, self.tree.QtotalCh0, self.tree.goodCh0] )
        self.miniTree = numpy.array(a)
        print '....Done!'
        self.close()
        return 
    def open(self,fn=None):
        '''
        open a root file
        '''
        debug = False
        
        if fn is None: sys.exit('procWaveDump.open ERROR No input file name')
        fullfn = self.rootFileDir + fn
        self.rf = ROOT.TFile.Open(fullfn) # must keep file in scope to keep tree in scope
#        ROOT.SetOwnership(rf,0)
        self.tree = self.rf.Get('tree')
        self.entries = self.tree.GetEntriesFast()
        nb = self.tree.GetEntry(self.entries-1)
        self.lastTime = self.tree.abs_time
        print 'procWaveDump.open Opened',fullfn,'with',self.entries,'entries. Last time',self.lastTime
                
        if debug: # for debug
            for jentry in xrange( 0,self.entries+1,self.entries/3 ):
                ientry = self.tree.LoadTree(jentry)
                print 'jentry,ientry',jentry,ientry
                if ientry<0:
                    break
                nb = self.tree.GetEntry(jentry)
                dt = self.tree.dt
                print 'nb',nb,'dt',dt

        
        return
    def close(self):
        ''' close the current root file '''
        print 'procWaveDump.close Close',self.rf.GetName()
        self.rf.Close()
        return
if __name__ == '__main__' :
    maxE = 9999999999
    if len(sys.argv)>2: maxE = int(sys.argv[2])
    
    pWD = procWaveDump()
    if len(sys.argv)>1:
        fn = sys.argv[1]
        pWD.main(inputfn=fn,maxE=maxE)
    else:
        pWD.main(maxE=maxE)

