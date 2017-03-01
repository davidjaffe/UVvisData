#!/usr/bin/env python
'''
analysis of scope data for LS QA
adapted from Lindsey's PSDscripts
20161003
'''
import sys
# the ONETON dir contains wfanal.py 
sys.path.append('/Users/djaffe/work/GIT/ONETON')
# the PROSPECT dir contains get_filepaths.py
sys.path.append('/Users/djaffe/work/GIT/PROSPECT')


import wfanal,graphUtils
import get_filepaths
import h5py
import matplotlib.pyplot as plt
import numpy
import ROOT
import os
import math

import gzip,shutil
import datetime,Logger
from scipy import signal as scipy_signal # used by getPeaks
from scipy import interpolate as scipy_interpolate # used by getFOMat6Li
import convertLabviewTime

class lsqa():
    def __init__(self,fn=None,SIMPLE=False):
        self.wfa = wfanal.wfanal()
        self.gU  = graphUtils.graphUtils()
        self.gfp = get_filepaths.get_filepaths()
        self.cLT = convertLabviewTime.convertLabviewTime()

        # file to be updated with summary results for each job
        self.resultsSummaryFile = 'Samples/resultsSummaryFile.txt'
                    
        # job start time
        now = datetime.datetime.now()
        fmt = '%Y%m%d_%H%M%S_%f'
        self.cnow = cnow = now.strftime(fmt)

        if SIMPLE:
            print 'lsqa__init__ SIMPLE processing only'
            return
        
        # process input filename to create new directories for figures and logfile
        dirname,rn = '',''
        if fn is not None:
            bn = os.path.basename(fn)
            rn = bn.split('.')[0]
            sd = fn.split('LSQAData')[1]
            dirname = 'Samples' +sd.replace(bn,rn) + '/'
            #print 'bn,rn,sd,dirname',bn,rn,sd,dirname

        self.sampledir = dirname
        self.figdir = dirname+'Figures/'
        self.logdir = dirname+'Logfiles/'
        for dn in [self.figdir, self.logdir]:
            if os.path.isdir(dn):
                pass
            else:
                try:
                    os.makedirs(dn)
                except IOError,e:
                    print 'lsqa.__init__',e
                else:
                    print 'lsqa.__init__ created',dn
        
        # produce logfile name and route stdout to logfile and terminal
        lfn = self.logdir + rn + '_' + cnow + '.log'
        sys.stdout = Logger.Logger(fn=lfn)
        print 'lsqa__init__ Output directed to terminal and',lfn
        print 'lsqa__init__ Job start time',cnow
                                
        # fitting functions
        self.GG = ROOT.TF1("GG","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])")
        G = self.G = ROOT.TF1("G","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])")
        G.SetParName(0,"N")
        G.SetParName(1,"Mean")
        G.SetParName(2,"Sigma")
        self.CERF = ROOT.TF1("CERF",self.cerf,0.,100.,3)
        self.CERF.SetParName(0,'N')
        self.CERF.SetParName(1,'Mean')
        self.CERF.SetParName(2,'Sigma')

        # definition of constants,etc.
        self.Cs137edge = 478. # keV
        self.LiCaptureEnergy = 540. # keV. From doc-1245-v1 "P50D ambient data stability"
        self.QforAvg = [15.,25.] # values of charge to use when averaging pulses
        self.startoff=100
        self.numavg  = 200
        
        self.PSDcut = 0.1
        self.sampledt = 0.2 # 1 sample/0.2ns
        self.favFast = 100
        self.favTotal= 800
        self.fastValues = numpy.linspace(60,140,5) #20,120,6)
        self.totalValues = numpy.linspace(200,1000,5)
        scan = False
        if scan:
            self.fastValues = numpy.linspace(60,140,5+4)
            self.totalValues= numpy.linspace(200,1000,5+4)
        
        return
    def cerf(self,v,p):
        ''' error function complement for fitting compton edge '''
        x = (v[0]-p[1])/math.sqrt(2)/p[2]
        return p[0]*math.erfc(x)
    def getWFMs(self,fn):
        '''
        get waveforms from hdf5 file
        '''
        f = h5py.File(fn)
        ph,ct = [],[] # array of scope readout, clock time
        for w in  f['Waveforms']:
            a,b = w[0]
            if len(a)>0:
                ph.append(a)
                ct.append(b)
        return ph,ct
    def getQ(self,wfm,fastValues,totalValues):
        '''
        efficiently process input waveforms to return sets of Qfast,Qtot pairs
        given 'fast' and 'total' integration parameters
        '''
        startoff = self.startoff
        numavg   = self.numavg
        minidx = numpy.argmin(wfm)
        BL = self.getBL(wfm,minidxInput=minidx,startoff=startoff,numavg=numavg)
        i1 = int(minidx-startoff)
        results = []
        for fastint in fastValues:
            i2 = int(minidx+fastint)
            fast = sum(wfm[i1:i2])
            fastBL = BL*float(startoff+fastint)
            for totint in totalValues:
                i3 = int(minidx+totint)
                total = fast + sum(wfm[i2:i3])
                totalBL = BL*float(startoff+totint)
                Qfast = fast-fastBL
                Qtot  = total-totalBL
                results.append( (Qfast,Qtot) )
        return results
    def runType(self,F):
        '''
        given filename or open hdf5 file,
        return type of file from length. ~20k = neutron only = 'N', ~5k = gamma and neutron source = 'G'
        if type is ambiguous or unknown, then check in LSQA_run_source.txt
        returns None if cannot determine runType
        '''
        f = None
        closeIt = False
        if type(F) is h5py._hl.files.File: f = F
        if type(F) is str:
            f = h5py.File(F,'r')
            closeIt = True
        
        wfms = f['Waveforms']
        L = len(wfms)
        gamma = abs(L-5000)<10
        neutron = abs(L-20000)<10
        if not (gamma or neutron):
            r = str(f.file).split()[2].replace('.h5','').replace('\"','') # should be 'runxxxxx' where 'xxxxx' is the timestamp
            q = open('LSQA_run_source.txt','r')
            for line in q:
                if r in line:
                    if 'gamma' in line: gamma = True
                    if 'neutron' in line: neutron = True
                    break
            if not (gamma or neutron): print 'lsqa.runType runtype unknown for',f.file

        if closeIt: f.close()
        if gamma: return 'G'
        if neutron: return 'N'
        return None
    def wfanaFile(self,fn=None):
        '''
        Return events containing Qfast,Qtot pairs for different fast, total definitions
        given input file.
        Also return pairs of (leadingEdge, waveform) for waveforms with a total charge
        within self.QforAvg identified as n-like or g-like
        '''

        fastValues = self.fastValues 
        totalValues= self.totalValues 
        if fn is None :
            sys.exit('lsqa.wfanaFile ERROR No input file specified')

        idebug = 0
        IntValues = [[x,y] for x in fastValues for y in totalValues]
        idxFav = IntValues.index( [self.favFast,self.favTotal] )
        print 'lsqa.wfanaFile idxFav',idxFav,'IntValues[idxFav]=',IntValues[idxFav]
        f = h5py.File(fn)
        print 'lsqa.wfanaFile Opened',fn
        Events = []
        gWFToAvg, nWFToAvg= [],[]

        ievt,freq = 0,500
        for w in f['Waveforms']:
            wfm,ct = w[0]
            if len(wfm)>0:
                results = self.getQ(wfm,fastValues=fastValues,totalValues=totalValues)
                Events.append( results )
                
                Qfast,Qtot = results[idxFav]
                #print 'lsqa.wfanaFile Qfast,Qtot,self.QforAvg',Qfast,Qtot,self.QforAvg
                if self.QforAvg[0]<=abs(Qtot) and abs(Qtot)<=self.QforAvg[1]:
                    PSD = (Qtot-Qfast)/Qtot
                    leadingEdge,lE1,lE2 = self.getLeadingEdge(wfm)
                    if PSD<self.PSDcut:
                        gWFToAvg.append( [leadingEdge, wfm] )
                    else:
                        nWFToAvg.append( [leadingEdge, wfm] )
                
                if idebug>0:
                    for iPair,rPair in zip(IntValues,results):
                        fastint,totint = iPair
                        Qfast,Qtot     = rPair
                        print '{0} {1} {2:.2f} {3:.2f}'.format(fastint,totint,Qfast,Qtot)
                    if ievt>3: break
                if ievt%freq==0:
                    print '\r',ievt,
                    sys.stdout.flush()
                ievt+=1

        f.close()
        print '\rlsqa.wfanaFile Processed',fn
        return Events,gWFToAvg,nWFToAvg
    def gammaCalib(self,Events,draw=True):
        '''
        return estimates of compton edge give sets of Qfast,Qtot
        Also return histogram for favorite value of total
        compton edge fitted using complement of error function
        gamma calibration only depends on Qtot
        '''
        debug = False
        fastValues = self.fastValues 
        totalValues= self.totalValues 
        IntValues = [[x,y] for x in fastValues for y in totalValues]
        nx,xmi = 100,5.+2.
        xma = xmi+float(nx)/2
        fopt = "SQ"
        EperQ = {}
        E = numpy.array(Events)
        G = self.CERF
        hists = []
        hfav  = None
        c1 = ROOT.TCanvas() # open, so it can be closed after fitting
        for i,vv in enumerate(IntValues):
            ifast,itot = vv
            key = 't'+str(int(itot))
            if key not in EperQ:
                name = 'Gamma_Calib_'+key
                title = name.replace('_',' ')
                h = ROOT.TH1D(name,title,nx,xmi,xma)
                hists.append(h)
                if hfav is None: hfav = h
                if itot==self.favTotal : hfav = h
                for pair in E[:,i]:
                    Qfast,Qtot = pair
                    h.Fill(abs(Qtot))
                peak,ipeak = h.GetMaximum(),h.GetMaximumBin()
                xpeak = float(ipeak)*(xma-xmi)/float(nx)+xmi
                G.SetParameters(peak,xpeak,2.5)
                if debug : print 'lsqa.gammaCalib peak,xpeak',peak,xpeak
                ptr = h.Fit(G,fopt,"",xmi,4.*xpeak)
                m,s = G.GetParameter(1),G.GetParameter(2)
                if debug : print 'lsqa.gammaCalib m,s',m,s
                ptr = h.Fit(G,fopt,"",xmi,m+2.5*s)
                m,s = G.GetParameter(1),G.GetParameter(2)                
                if debug : print 'lsqa.gammaCalib m,s',m,s
                HWHM = math.sqrt(2.*math.log(2.))*s
                Qedge = m + HWHM
                print 'lsqa.gammaCalib Compton edge at',Qedge,'for',name
                EperQ[key] = -1.
                if Qedge>0: EperQ[key] = self.Cs137edge/Qedge
        if draw : self.gU.drawMultiHists(hists,'gammaCalib',figdir=self.figdir,Grid=True,fitOpt=1111)
        c1.Close()
        writeROOT = False
        if writeROOT:
            fn = self.figdir + 'gammaCalib.root'
            rfn = ROOT.TFile.Open(fn,'RECREATE')
            for h in hists: rfn.WriteTObject(h)
            rfn.Close()
            print 'lsqa.gammaCalib Wrote',len(hists),'to',fn
        return EperQ,hfav
    def timestampFromFilename(self,fn):
        ''' extract time stamp as integer from filename  '''
        bn = os.path.basename(fn)
        ts = int((bn.split('.')[0]).replace('run',''))
        return ts
    def findPreviousOutput(self,fn=None):
        '''
        return filename of root file produced by previous analysis of current file.
        If not, try to find a root file by analysis of run close in time to the current file
        '''
        if fn is None: sys.exit('lsqa.findPreviousOutput ERROR No input filename')
    
        debug = True

        bn = os.path.basename(fn).split('.')[0] # should be runxxxxxx where xxxxxx = run#
        runnum = int(bn.replace('run','')) # integer, current run number
        drun = 10000000000000
        rbn = 'optimize_'+bn+'.root'
        rfn = None
        if bn in self.figdir:  # directory for this run exists, look for root file
            rfn = self.figdir + rbn

        if not os.path.isfile(rfn): # look in other directories
            listoffiles = self.gfp.get_filepaths(self.sampledir)
            sampleName = self.sampledir.split('/')[-2]

            oldrun = None
            for name in listoffiles:
                if sampleName in name:
                    s = name.split('/')
                    i = s.index(sampleName)
                    if 'run' in s[i+1]: 
                        newrun = int(s[i+1].replace('run',''))
                        if newrun!=oldrun:
                            if abs(newrun-runnum)<drun:
                                rbn = 'optimize_'+s[i+1]+'.root'
                                rfn = self.sampledir + s[i+1] + '/Figures/'+rbn
                                if os.path.isfile(rfn):
                                    drun = abs(newrun-runnum)
                                    oldrun = newrun
                            
        if os.path.isfile(rfn): return rfn
        return None
                
    def findGammaRun(self,fn=None):
        '''
        find file with run with gamma source closest in time to input file
        '''
        if fn is None: sys.exit('lsqa.findGammaRun ERROR No input filename')

        debug = 0
        
        h5size = float(408498176) # size of run3566813028.h5
        gzsize = float(30126638)  # size of run3567072982.h5.gz
        sizes = {'.h5':h5size, '.gz':gzsize}
            
        bn = os.path.basename(fn)
        ts0 = self.timestampFromFilename(fn)
        filepath = fn.replace(bn,'')
        if debug>0: print 'lsqa.findGammaRun fn,bn,ts0,filepath',fn,bn,ts0,filepath
        listoffiles = self.gfp.get_filepaths(filepath)
        if debug>1 : print 'lsqa.findGammaRun listoffiles',listoffiles
        gammaFile,dt = None,None
        for n in listoffiles:
            filename,extension = os.path.splitext(n)
            size = float(os.path.getsize(n))
            if debug>0 : print 'lsqa.findGammaRun path,extension,size',n,extension,size
            if extension in sizes:
                if abs(size/sizes[extension]-1.)<0.5:
                    ts = self.timestampFromFilename(n)
                    if debug>0 : print 'lsqa.findGammaRun ts,ts-ts0,dt,size/sizes[extension]',ts,ts-ts0,dt,size/sizes[extension]
                    if dt is None:
                        dt = abs(ts-ts0)
                        gammaFile = n
                    else:
                        if abs(ts-ts0)<dt:
                            dt = abs(ts-ts0)
                            gammaFile = n

        return gammaFile
    def gunzip(self,fn):
        '''
        if necessary, gunzip input file, return output file.
        if output file already exists, just return output file
        '''
        h5f = fn.replace('.gz','')
        if os.path.isfile(h5f): return h5f
        print 'lsqa.gunzip input',fn,'output',h5f,
        sys.stdout.flush()
        with gzip.open(fn,'rb') as fin, open(h5f,'wb') as fout:
            shutil.copyfileobj(fin,fout)
        print 'DONE'
        return h5f
    def main(self,filename=None):
        '''
        process neutron file
        find gamma source file that is nearest in time to neutron file, based on timestamp, to obtain
        charge to keVee conversion

        turn waveforms into sets of Qfast,Qtot pairs for different definitions of fast and tot
        then evaluate FOM as function of Qtot
        Best global values appear to be fast=100,total=600.
        '''
        idebug = 0
        fastValues = self.fastValues 
        totalValues= self.totalValues 
        IntValues = [[x,y] for x in fastValues for y in totalValues]

        PSDcut = self.updatePSDcut(fn=filename)
        print 'lsqa.main Update PSDcut from',self.PSDcut,'to',PSDcut
        self.PSDcut = PSDcut 
        
        allhists = []
        opsobjects,opsLogy = [],[] # One Page Summary objects and logy specifications
        summaryInfo = {}
    
        
        nx = len(fastValues)
        dx = float(fastValues[1]-fastValues[0])/2.
        xmi,xma = fastValues[0]-dx,fastValues[-1]+dx
        ny = len(totalValues)
        dy = float(totalValues[1]-totalValues[0])/2.
        ymi,yma = totalValues[0]-dy,totalValues[-1]+dy
        name = 'Opt_CumFOM'
        title= 'Optimize cumulative FOM'
        hcum = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
        allhists.append(hcum)
        name = 'Opt_6LiFOM'
        title= 'Optimize 6Li FOM'
        h6li = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
        allhists.append(h6li)
        name = 'Opt_6LiZ'
        title= 'Optimize 6Li Z'
        hZ  = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
        allhists.append(hZ)

        # uncompress input file, if needed
        fn = filename
        if filename[-3:]=='.gz': fn = self.gunzip(filename)
        bn = os.path.basename(fn).replace('.h5','')

        ### gamma calibration
        ### identify gamma file appropriate for this neutron only run
        ### create histograms of average pulse shape for gamma and neutron band
        ### then do gamma calibration to get energy/charge for each definition of total integral as well as
        ### histogram of charge distribution for favorite total integral definition
        gammaFile = ZgammaFile = self.findGammaRun(fn)
        print 'lsqa.main fn',fn,'gammaFile',ZgammaFile
        if ZgammaFile[-3:]=='.gz':
            gammaFile = self.gunzip(ZgammaFile)
        print 'lsqa.main gammaFile',gammaFile
        summaryInfo['gammaFile'] = gammaFile
        summaryInfo['gammaRunTime'] = self.fileTime(gammaFile)
        gEvents,gWF,nWF = self.wfanaFile(fn=gammaFile)
        hgg1,hgg2 = self.avgWF(gWF,title='Cs137_gamma')
        hng1,hng2 = self.avgWF(nWF,title='Cs137_neutron')
        avWFhists = [ [hgg1,hng1], [hgg2,hng2] ]
        allhists.extend( [hgg1,hng1,hgg2,hng2])
        EperQ,hfav = self.gammaCalib(gEvents,draw=False)
        opsobjects.append(hfav)
        opsLogy.append(False)
        print 'lsqa.main EperQ',EperQ
        key = 't'+str(self.favTotal)
        
        summaryInfo['fav'] = [self.favFast, self.favTotal]
        summaryInfo[key] = EperQ[key]
        sKeys = ['fav',key,'gammaFile','neutronFile','gammaRunTime','neutronRunTime']
        
        # turn waveforms into sets of Qfast,Qtot pairs
        Events,gWF,nWF = self.wfanaFile(fn=fn)
        summaryInfo['neutronFile'] = fn
        summaryInfo['neutronRunTime'] = self.fileTime(fn)
        hg1,hg2 = self.avgWF(gWF,title='gamma')
        hn1,hn2 = self.avgWF(nWF,title='neutron')
        avWFhists.extend( [ [hg1,hn1], [hg2,hn2] ] )
        allhists.extend( [hg1,hn1, hg2,hn2] )
        opsobjects.append( [hg2,hn2] )
        opsLogy.append(True)
        words = self.runType(fn)  # should be either 'G' (for gamma+neutron source) or 'N' (for neutron source only)

        if idebug>0:
            print 'lsqa.main Events',Events
        Qlo = numpy.linspace(5,65,7)
        Qhi = numpy.array([x+10. for x in Qlo])

        sumFOM = {}
        E = numpy.array(Events)
        hists,qhists,graphs = [],[],[]
        tmg = self.gU.makeTMultiGraph('FOMerific_'+bn)
        tmgZ= self.gU.makeTMultiGraph('Zerific_'+bn)
        nx,xmi,xma = 100,0.,100.
        ny,ymi,yma = 100,0.,0.5
        for i,vv in enumerate(IntValues):
            ifast,itot = vv
            favorite = ifast==self.favFast and itot==self.favTotal
            eKey = 't'+str(int(itot))
            PSDdef = name = 'f'+str(int(ifast))+'_t'+str(int(itot))
            title = 'PSD vs Qtot '+name.replace('_',' ')
            h = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
            qname = 'Etot_'+name+'_wPSDcut'
            qtitle = 'Etot ' + name + 'PSD>'+ '{0:.4f}'.format(self.PSDcut)
            qh= ROOT.TH1D(qname,qtitle,nx,EperQ[eKey]*xmi,EperQ[eKey]*xma)
            if favorite:
                opsobjects.extend(  [h,qh] )
                opsLogy.extend( [False for i in range(len([h,qh]))] )
            
            aQtot,aPSD = [],[]
            for pair in E[:,i]:
                Qfast,Qtot = pair
                PSD = (Qtot-Qfast)/Qtot
                Qtot = abs(Qtot)
                h.Fill(Qtot,PSD)
                if PSD>self.PSDcut: qh.Fill(EperQ[eKey]*Qtot)
                aQtot.append(Qtot)
                aPSD.append(PSD)
            hists.append(h)
            qhists.append(qh)
            aQtot = numpy.array(aQtot)
            aPSD  = numpy.array(aPSD)
            FOM,dFOM,FOMhists,FOMgraph,QAFOM,Zgraph,ZFOM = self.beerFOM(aQtot,aPSD,Qlo,Qhi,EperQ[eKey],PSDdef+words,bn,Draw=True,debug=False)  # QAFOM, ZFOM are 2 element arrays

            if favorite: 
                summaryInfo['FOM'] = QAFOM
                summaryInfo['Z']   = ZFOM
                sKeys.append('FOM')
                sKeys.append('Z')
            
            sFOM,sdFOM = numpy.sum(FOM),numpy.sqrt(numpy.sum([x*x for x in dFOM]))
            sumFOM[PSDdef+words] = (sFOM,sdFOM)+QAFOM
            hcum.Fill(float(ifast),float(itot),sFOM)
            h6li.Fill(float(ifast),float(itot),float(QAFOM[0]))
            hZ.Fill(float(ifast),float(itot),float(ZFOM[0]))
            #print 'lsqa.main {0} sum {1:.2f}({2:.2f})'.format(PSDdef+words,sFOM,sdFOM) 
            if favorite:
                opsobjects.extend( FOMhists )
                opsLogy.extend( [True for i in range(len(FOMhists))] )
                opsobjects.append( FOMgraph )
                opsobjects.append( Zgraph )
                opsLogy.append( False ) # FOMgraph
                opsLogy.append( False ) # Zgraph
            allhists.extend( FOMhists )
            graphs.append( FOMgraph )
            graphs.append( Zgraph )
            self.gU.color(FOMgraph,i,i,setMarkerColor=True)
            self.gU.color(Zgraph,i,i,setMarkerColor=True)
            tmg.Add( FOMgraph )
            tmgZ.Add( Zgraph )
            


        print " PSD_defn      Cum_FOM   6LiFOM     Sorted by Cum_FOM  lsqa.main"
        for g in sorted( sumFOM.items(), key=lambda x: x[1]):
            print ' {0:12} {1:.2f}({2:.2f}) {3:.3f}({4:.3f})'.format(g[0],g[1][0],g[1][1],g[1][2],g[1][3])

        for i,h in enumerate(opsobjects):
            logy = opsLogy[i]
            if type(h) is list:
                for hh in h: print 'lsqa.main opsobjects multi,logy',i,hh.GetName(),logy
            else:
                print 'lsqa.main opsobjects,logy',i,h.GetName(),logy

        self.gU.drawMultiObjects(opsobjects,setLogy=opsLogy,fname='summary_'+bn,figdir=self.figdir,statOpt=0,abscissaIsTime=False,Grid=True,changeColors=True,addLegend=True,biggerLabels=1.5,dopt='histfunc',gopt='AP',forceNX=3)

        self.gU.drawMultiHists(hists,'optimize_'+bn,figdir=self.figdir,forceNX=4,Grid=True)
        self.gU.drawMultiHists(qhists,'keVee_with_PSDcut_'+bn,figdir=self.figdir,forceNX=4,Grid=True)
        self.gU.drawMultiHists([h6li,hcum,hZ],'scan_'+bn,figdir=self.figdir,Grid=False,statOpt=0,dopt='colztextcont3')
        tmg.SetMinimum(0.)
        tmg.SetMaximum(3.0)
        canvas = self.gU.drawMultiGraph(tmg,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='keVee',yAxisLabel='FOM',NLegendColumns=3)
        canvasZ = self.gU.drawMultiGraph(tmgZ,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='keVee',yAxisLabel='Z',NLegendColumns=3)
        graphs.append( tmg )
        graphs.append( tmgZ )


        self.gU.drawMultiHists(avWFhists,fname='avergeWF_'+bn,figdir=self.figdir,statOpt=0,setLogy=True,abscissaIsTime=False,Grid=True,changeColors=True,addLegend=True)



        rfn = self.figdir + 'optimize_'+bn+'.root'
        f = ROOT.TFile.Open(rfn,'RECREATE')
        for h in hists: f.WriteTObject(h)
        for h in allhists: f.WriteTObject(h)
        for g in graphs: f.WriteTObject(g)
        f.Close()
        print 'lsqa.main Wrote',len(hists)+len(allhists),'hists and',len(graphs),'graphs to',rfn


        print 'sKeys',sKeys
        for key in sKeys:
            if key=='fav':
                print 'lsqa.main Summary favorite fast,total integral {0:.1f} {1:.1f}'.format(summaryInfo[key][0],summaryInfo[key][1])
            elif key=='FOM' or key=='Z':
                print 'lsqa.main Summary {0} at 6Li peak {1:.4f}({2:.4f})'.format(key,summaryInfo[key][0],summaryInfo[key][1])
            elif 'File' in key or 'Run' in key:
                print 'lsqa.main Summary',key,summaryInfo[key]
            else:
                print 'lsqa.main Summary keV/charge {0:.3f}'.format(summaryInfo[key])
        if 0: ##### TESTING
            import pickle
            pfn = 'pickleLSAQfile.p'
            pickle.dump(summaryInfo,open(pfn,'wb'))
            print 'lsqa.main pickled summaryInfo to',pfn
    
        self.writeRSF(summaryInfo)
        
        return
    def writeRSF(self,summaryInfo):
        '''
        append summary information to results summary file
        '''
        fast,total = summaryInfo['fav']
        Keys = sorted(summaryInfo.keys())
        Keys.remove('fav')

        run = summaryInfo['neutronFile'].split('/')[-1].split('.')[0]
        a = ' ' + run
        h = '*run#'
        L = max(len(a),len(h))
        aline = a.ljust(L)
        header= h.ljust(L)
        
        key = 't'+str(total)
        Keys.remove(key)
        h = ' ' + 'keV/charge'
        a = ' {0:.2f}'.format(summaryInfo[key])
        L = max(len(a),len(h))
        aline += a.ljust(L)
        header += h.ljust(L)

        for key in ['FOM','Z']:
            Keys.remove(key)
            h = ' ' + key
            a = ' {0:.4f}({1:.4f})'.format(summaryInfo[key][0],summaryInfo[key][1])
            L = max(len(a),len(h))
            aline += a.ljust(L)
            header += h.ljust(L)

        h = ' fast total'
        a = ' {0} {1}'.format(fast,total)
        L = max(len(a),len(h))
        aline += a.ljust(L)
        header += h.ljust(L)

        RKeys,FKeys = [],[]
        for key in Keys:
            if 'Run' in key: RKeys.append(key)
            if 'File' in key:FKeys.append(key)
        RKeys.extend(FKeys)
        for key in RKeys:
            Keys.remove(key)
            a = ' ' + summaryInfo[key].replace('../','')
            h = ' ' + key
            L = max(len(a),len(h))
            aline += a.ljust(L)
            header+= h.ljust(L)
        a = ' ' + self.cnow
        h = ' ' + 'Job time'
        L = max(len(a),len(h))
        aline += a.ljust(L)
        header+= h.ljust(L)
        aline += '\n'
        header+= '\n'
        with open(self.resultsSummaryFile,'a') as rsf:
            rsf.write(header)
            rsf.write(aline)
            rsf.close()
        print 'lsqa.writeRSF Wrote header',header, 'lsqa.writeRSF Wrote line',aline, 'lsqa.writeRSF to file',self.resultsSummaryFile
        return
    def fileTime(self,fn):
        '''
        return time of start of run YYYYMMDD HH:MM:SS using timestamp (aka run number) from filename
        '''
        ts = fn.split('/')[-1].split('.')[0].replace('run','')
        its = int(ts)
        fT = self.cLT.convert(its)
        return fT
    def avgWF(self,lEWF,title='Average waveform'):
        '''
        return histograms of average waveform given pairs of leading edge and waveforms
        1st hist is vs samples, 2d hist is vs time
        '''
        startoff,numavg=self.startoff,self.numavg
        i1,i2 = 0-startoff-numavg,0+self.favTotal
        xmi = float(i1)-0.5
        xma = float(i2)+0.5
        nx = int(xma-xmi)
        nWF = len(lEWF)
        htit = title + ' vs samples'
        name = htit.replace(' ','_')
        htit += ' ' + str(nWF) + ' waveforms averaged'
        h1 = ROOT.TH1D(name,htit,nx,xmi,xma)
        htit = title + ' vs time in ns'
        name = htit.replace(' ','_')
        htit += ' ' + str(nWF) + ' waveforms averaged'
        h2 = ROOT.TH1D(name,htit,nx,xmi*self.sampledt,xma*self.sampledt)

        x = numpy.arange(float(i2-i1))+float(i1)
        for pair in lEWF:
            leadingEdge,WF = pair
            BL = self.getBL(WF,startoff=startoff,numavg=numavg)
            X = numpy.arange(float(len(WF))) - leadingEdge
            f = scipy_interpolate.interp1d(X,WF)
            y = []
            for v in x:
                y.append(abs(f(v)-BL))   # invert and subtract baseline
            y = numpy.array(y)
            y = y/numpy.sum(y)/float(nWF) # normalize
            for u,v in zip(x,y):
                h1.Fill(u,v)
                h2.Fill(u*self.sampledt,v)
        return h1,h2
    def onlPlot(self,fn=None,dname=None):
        '''
        get Qfast,Qtot from online analysis and plot 'em
        Vestigial, used for testing, debugging
        '''
        Cs137edge = self.Cs137edge  # keV
        FOMfitThres = 25. # minimum number of events required to fit for FOM

        CF = -1.e9 # arbitrary factor to change measured charge
        Qmi,Qma = 0.1,20.1 # charge with CF applied
        Qlo = [float(i)+1. for i in range(15)]
        Qhi = [x+1. for x in Qlo]
        if Qma > Qhi[-1]:
            Qlo.append( Qhi[-1] )
            Qhi.append(Qma)
        PSDcut = 0.15

        EperQ = {} # keV/charge

        allBN = []
        
        if type(fn) is not list: fn = [fn]
        hists = []
        graphs = []
        icol = 0
        for filename,w in zip(fn,dname):
            print 'lsaq.onPlot Open',filename
            f = h5py.File(filename)
            basename = os.path.basename(filename)
            bn = basename.replace('.h5','')
            allBN.append(bn)
            words = ''
            if w is not None: words = w
            OD =  f['OnlinePSD']
            print 'lsqa.onlPlot Opened',filename,'with',len(OD),'events and type',words
            plt.clf()
            plt.grid()
            Qfast,Qtot = CF*OD[:,0], CF*OD[:,1]  #NOTE CONVERSION
            Qslow = Qtot - Qfast
            PSD = Qslow/Qtot
            A = Qtot>0.1
            nx,xmi,xma = 100,Qmi,Qma
            ny,ymi,yma = 110,-.1,1.
            title = 'Qslow/Qtot vs Qtot '+words
            name  = 'psd_v_q' + '_' + words.replace(' ','_')
            h = self.gU.makeTH2D(Qtot[A],PSD[A],title,name,nx=nx,xmi=xmi,xma=xma,ny=ny,ymi=ymi,yma=yma)
            hists.append(h)
            name = title = 'Qtot'+ '_' + words.replace(' ','_')
            h = self.gU.makeTH1D(Qtot[A],title,name,nx=nx,xmi=xmi,xma=xma)            
            hists.append(h)
            B = (Qtot>Qmi)*(PSD<PSDcut)
            name = title = 'Qtot_loPSD'+ '_' + words.replace(' ','_')
            h = self.gU.makeTH1D(Qtot[B],title,name,nx=nx,xmi=xmi,xma=xma)            
            hists.append(h)
            G = self.G 
            G.SetParName(0,"N")
            G.SetParName(1,"Mean")
            G.SetParName(2,"Sigma")
            peak = h.GetMaximum()
            G.SetParameters(peak,1.5,.8)
            ptr = h.Fit("G","SLQ","",0.,4.)
            m,s = G.GetParameter(1),G.GetParameter(2)
            HWHM = math.sqrt(2.*math.log(2.))*s
            Qedge = m + HWHM
            print 'lsqa.onlPlot Compton edge at',Qedge,'for',name
            EperQ[words] = Cs137edge/Qedge
            
            PSDdef = 'default'
            FOM,dFOM,FOMhists,FOMgraph,QAFOM = self.beerFOM(Qtot,PSD,Qlo,Qhi,EperQ,PSDdef+words,bn,Draw=True)
            
            f.close()

        graphs.append( FOMgraph )
        hists.extend( FOMhists )

        nameForFile = '_'.join(allBN)

        RFN = self.figdir + nameForFile + '.root'
        rfn = ROOT.TFile.Open(RFN,'RECREATE')
        for h in hists: rfn.WriteTObject(h)
        for g in graphs:rfn.WriteTObject(g)
        rfn.Close()
        print 'lsqa.onlPlot Wrote',len(hists),'hists and',len(graphs),'graphs to',RFN
            
        return
    def beerFOM(self,Qtot,PSD,Qlo,Qhi,EperQ,PSDdef,bn,Draw=True,debug=True):
        '''
        return arrays of figure-of-merit and uncertainty, a list of fitted histograms, a graph of FOM vs Q (or E)
               and the figure-of-merit and uncertainty at the nominal 6Li capture energy, if EperQ is given.
        given paired Qtot,PSD arrays, lower and upper charge limits Qlo,Qhi.
        EperQ, if available, converts charge to energy.
        PSDdef is string defining the PSD and bn is basename.
        PSDdef & bn are used for histogram, graph and filenames if Draw is True
        '''

        FOMfitThres = 25. # minimum number of events required to fit for FOM

        GG = self.GG 
        GG.SetParName(0,"gN")
        GG.SetParName(1,"gMean")
        GG.SetParName(2,"gSigma")
        GG.SetParName(3,"nN")
        GG.SetParName(4,"nMean")
        GG.SetParName(5,"nSigma")
        hists = []
        Qmi,Qma = min(0.,min(Qlo)),max(Qhi)
        nx,xmi,xma = 100,Qmi,Qma
        ny,ymi,yma = 120,-.1,0.5
        ybinsize = (yma-ymi)/float(ny)

        FOMmin,FOMmax = 0.,3.
        Zmin,Zmax     = 0.,0.2

        gMean,nMean,gSigma,nSigma = None,None,0.02,0.02  # initial guesses
        nysig = int(0.02/ybinsize + 0.01)  # 1 sigma in bins
        c1 = ROOT.TCanvas() # open, so it can be closed
        X,dX,Y,dY = [],[],[],[]  # X = charge or energy, Y = FOM
        Z,dZ = [],[]         # Z = difference in means for nuclear and electronic recoils
        for lo,hi in zip(Qlo,Qhi):
            A = (lo<=Qtot)*(Qtot<hi)
            title = name = 'PSD_Q'+str(int(lo)).zfill(2)+'_'+str(int(hi)).zfill(2) + PSDdef
            h = self.gU.makeTH1D(PSD[A],title,name,nx=ny,xmi=ymi,xma=yma)
            hists.append(h)

            Nevt = h.GetEntries()
            FOM = None
            if Nevt>FOMfitThres:

                # use peak finder to estimate means. Use sigmas from previous fit if available
                peak = h.GetMaximum()
                Peaks = self.getPeaks(h,ny,nysig)
                if len(Peaks)==2:
                    gMean,nMean = Peaks
                if len(Peaks)==1:
                    step = .1
                    if gMean is not None: step = abs(gMean-nMean)
                    gMean = Peaks[0]
                    nMean = gMean + step
                if gMean is None:
                    gMean,nMean = 0.05,0.15
                if debug or len(Peaks)<2: print 'lsqa.beerFOM:',name,'Peaks,gMean,nMean',' '.join(['{0:.3f}'.format(p) for p in Peaks]),'{0:.3f} {1:.3f}'.format(gMean,nMean)
                GG.SetParameters(peak, gMean, gSigma, peak, nMean, nSigma)
                ptr = h.Fit(GG,"SLQ")

                #ptr.Print("V") # print everything about fit
                gMean,gSigma = GG.GetParameter(1),GG.GetParameter(2)
                nMean,nSigma = GG.GetParameter(1+3),GG.GetParameter(2+3)
                FWHM = 2.*math.sqrt(2.*math.log(2.))*math.sqrt(gSigma*gSigma + nSigma*nSigma)
                FOM = abs(gMean-nMean)/FWHM

                dgMean,dnMean = GG.GetParError(1),GG.GetParError(1+3)

                par = []
                for ii in range(6): par.append( GG.GetParameter(ii) )
                par = numpy.array(par)
                #print 'par',par
                dFOM  =self.getFOMunc(ptr,par,FOM)

            Elo,Ehi,Units = lo,hi,'charge'
            if EperQ is not None: Elo,Ehi,Units = EperQ*lo,EperQ*hi,'keV'

            if debug: print 'lsqa.beerFOM {0} lo,hi,FOM {1:.1f} {2:.2f} {3}'.format(PSDdef,Elo,Ehi,Units),
            X.append(0.5*(Ehi+Elo))
            dX.append(0.5*(Ehi-Elo))
            if FOM is None:
                Y.append(0.)
                dY.append(10.)
                Z.append(0.)
                dZ.append(10.)
                if debug: print FOM
            else:
                Y.append(FOM)
                dY.append(dFOM)
                Z.append(nMean-gMean)
                dz = math.sqrt(dnMean*dnMean+dgMean*dgMean)
                dZ.append(dz)
                FOMmin = min(FOM-dFOM,FOMmin)
                FOMmax = max(FOM+dFOM,FOMmax)
                #Zmin = min(nMean-gMean-dz,Zmin)
                Zmax = max(nMean-gMean+dz,Zmax)
                if debug: print ' {0:.2f}({1:.2f})'.format(FOM,dFOM)
            
        c1.Close() 
        name = 'FOM_v_Q_'+PSDdef+bn
        nameZ= 'Z_v_Q_'+PSDdef+bn
        QA_FOM,QA_Z = None,None
        if EperQ is not None:
            name,nameZ = 'FOM_v_keVee_'+PSDdef+bn, 'Z_v_keVee_'+PSDdef+bn
            QA_FOM = self.getFOMat6Li(numpy.array(X),numpy.array(Y),numpy.array(dX),numpy.array(dY))
            QA_Z   = self.getFOMat6Li(numpy.array(X),numpy.array(Z),numpy.array(dX),numpy.array(dZ))
        title,titleZ = name.replace('_',' '),nameZ.replace('_',' ')
        g,gZ = self.gU.makeTGraph(X,Y,title,name,ex=dX,ey=dY),self.gU.makeTGraph(X,Z,titleZ,nameZ,ex=dX,ey=dZ)
        self.gU.color(g,0,0,setMarkerColor=True)
        self.gU.color(gZ,0,0,setMarkerColor=True)
        if Draw:
            self.gU.drawGraph(g,figDir=self.figdir,option='AP',yLimits=[FOMmin,FOMmax])
            self.gU.drawGraph(gZ,figDir=self.figdir,option='AP',yLimits=[Zmin,Zmax])


        
        if Draw:
            nameForFile = bn+'_'+PSDdef

            hists1d,hists2d = [],[]
            for h in hists:
                if h.GetDimension()==1: hists1d.append(h)
                if h.GetDimension()==2: hists2d.append(h)

            if debug : print 'lsqa.beerFOM len(hists),len(hists1d),len(hists2d)',len(hists),len(hists1d),len(hists2d)
            if len(hists1d)>0: self.gU.drawMultiHists(hists1d,nameForFile+'_1d',figdir=self.figdir,forceNX=4,Grid=False,biggerLabels=True,fitOpt=111,statOpt=0)
            if len(hists2d)>0: self.gU.drawMultiHists(hists2d,nameForFile+'_2d',figdir=self.figdir,forceNX=2,Grid=True)

        
        return Y,dY,hists,g,QA_FOM,gZ,QA_Z
    def getFOMat6Li(self,E,FOM,dE,dFOM):
        '''
        estimate FOM and uncertainty at 6Li capture energy in keVee given FOM as function of energy
        '''
        debug = False
        Q = {}
        if debug : print 'lsqa.getFOMat6Li E,FOM,dE,dFOM',E,FOM,dE,dFOM
        for i in [-1.,0.,1.]:
            f = scipy_interpolate.interp1d(E,FOM+i*dFOM)
            Q[i] = float(f(self.LiCaptureEnergy))
        FOM6Li = Q[0.]
        dFOM6Li= 0.5*( abs(Q[1.]-FOM6Li) + abs(Q[-1.]-FOM6Li) )
        if debug : print 'lsqa.getFOMat6Li',Q,FOM6Li,dFOM6Li,'FOM6Li {0:.3f}({1:.3f})'.format(FOM6Li,dFOM6Li)
        return (FOM6Li,dFOM6Li)
    def getFOMunc(self,ptr,pars,FOM):
        '''
        return uncertainty in FOM given pointer to fit result and best fit parameters
        '''

        N = len(pars)
        # get covariance matrix as 2-d numpy array by brute force
        cov = []
        for ii in range(N):
            for jj in range(N):
                cov.append( ptr.CovMatrix(jj,ii) )
        #print 'len(cov)',len(cov)
        cov = numpy.array(cov)
        cov = numpy.reshape(cov, (N,N))
        cov = numpy.asmatrix( cov )
        #print 'cov.shape',cov.shape

        # compute array of derivatives evaluated at best fit value
        C = 2.*math.sqrt(2.*math.log(2.))
        Q = FOM
        s2 = pars[2]*pars[2] + pars[5]*pars[5]
        s = math.sqrt(s2)
        A = numpy.array( [ 0.,-1./C/s, -pars[2]*Q/s2, 0., 1./C/s, -pars[5]*Q/s2] )
        mA = numpy.mat(A)
        mAT= mA.T
        #print 'A.shape,A.T.shape,mA.shape,mAT.shape',A.shape,A.T.shape,mA.shape,mAT.shape
        v = numpy.dot(numpy.dot(mA,cov),mAT)
        u = -1.
        if v>=0.: u = math.sqrt(v)
        #print 'FOMunc variance',v,'unc',u
        

        return u
    def getPeaks(self,h,nx,sig,maxN=2):
        '''
        return list of peaks with width sig in input 1d hist h with nx bins
        maximum number of peaks that can be found is maxN
        uses scipy.signal.find_peaks_cwt
        '''
        debug = False
        c = []
        for i in range(nx):
            c.append( h.GetBinContent(i+1) )
        c = numpy.array(c)
        n = max(2,sig)
        ipeaks,oldpeaks = None,None
        while ipeaks is None or len(ipeaks)>maxN:
            oldpeaks = ipeaks
            ipeaks = scipy_signal.find_peaks_cwt(c,numpy.arange(1,n))
            if debug : print 'lsqa.getPeaks: ipeaks,n',ipeaks,n
            n += 1
        if len(ipeaks)==1 and oldpeaks is not None and len(oldpeaks)==3:
            ipeaks = [oldpeaks[0],oldpeaks[2]]
            if debug : print 'lsqa.getPeaks: ipeaks',ipeaks,'fixed'
        peaks = []
        for i in ipeaks: peaks.append( h.GetXaxis().GetBinCenter(i) )
        return peaks
    def getLeadingEdge(self,wf,frac=0.5):
        '''
        returning leading edge in floating bin number for input waveform
        also return 1st,last bin number used to compute leading edge
        '''
        startoff = self.startoff
        numavg   = self.numavg
        minidx = numpy.argmin(wf)
        H = wf[minidx]
        BL = self.getBL(wf,minidxInput=minidx,startoff=startoff,numavg=numavg)
        hLE = frac*(H-BL)+BL
        leadBins,leadWF = [],[]
        i = minidx
        while wf[i]<BL and i>0:
            leadBins.append(i)
            leadWF.append(wf[i])
            i-=1
        leadBins = numpy.array(leadBins,'float')
        leadWF   = numpy.array(leadWF,'float')
        f = scipy_interpolate.interp1d(leadWF,leadBins)
        leadingEdge = f(hLE)
        #print 'minidx,leadingEdge,leadBins[-1],leadBins[0]',minidx,leadingEdge,leadBins[-1],leadBins[0]
        return leadingEdge,leadBins[-1],leadBins[0]
    def getBL(self,wf,minidxInput=None,startoff=None,numavg=None):
        '''
        return baseline for input waveform
        '''
        minidx = minidxInput
        if minidxInput is None : minidx = numpy.argmin(wf)
        BL = sum(wf[int(minidx-startoff-numavg):int(minidx-startoff)])/float(numavg)
        return BL
    def simple(self,fn = '/Users/djaffe/work/LSQAData/LiLS01/run3566934507.h5'):
        '''
        used for testing, debugging, visualizing
        '''
        bn = os.path.basename(fn).replace('.h5','')
        usePulseAnal = False
        startoff, fastint, totint, numavg = self.startoff, 120, 600, self.numavg
        nsd = 3.0
        ph,ct = self.getWFMs(fn)
        print 'len(ph),len(ct)',len(ph),len(ct)
        for ip in range(0,len(ph),233):
            wf = ph[ip]
            if usePulseAnal:
                ped,pedsd,iPulse,subPperP,pArea,pTime = self.wfa.pulseAnal(wf,'LSQA',debug=1,nsd=nsd)
                print 'ped,pedsd,iPulse,subPperP,pArea,pTime',ped,pedsd,iPulse,subPperP,pArea,pTime
                words = 'pulseAnal'
            else:
                iPulse = []
                minidx = numpy.argmin(wf)
                windows = [ [minidx-startoff-numavg,minidx-startoff], [minidx-50,minidx+totint], [minidx-50,minidx+fastint] ]
                lE,b0,b1 = self.getLeadingEdge(wf)
                
                windows.append( [int(b0),int(b1)])
                windows.append( [int(lE),int(lE+1.)] )
                title = 'Evt '+str(ip) + ' ' + bn
                self.drawPulse(wf,windows=windows,title=title)
            if len(iPulse)>0:
                for j,p in enumerate(iPulse):
                    title = words+ ' Evt '+str(ip)+ ' pulse#'+str(j)+' of '+str(len(iPulse)-1)+' t '+str(pTime[j])
                    self.drawPulse(wf,window=iPulse[j],title=title)
        return       
    def drawPulse(self,wf,windows=None,title='waveform'):
        '''
        draw a waveform. overlay windows used for integration
        windows = [ [minBL,maxBL], [minT,maxT], [minF,maxF] ]
        where [minX,maxX] is the min,max sample number corresponding
        to X=BL,T,F = baseline,total,fast
        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        colorpoint =  ['b.', 'go', 'r.', 'co', 'mD', 'kD']
        i1,i2 = windows[0]
        for w in windows:
            i1 = min(i1,w[0])
            i2 = max(i2,w[1])
        i1 = max(0,i1-100)
        i2 = min(len(wf),i2+100)
        W = [ [i1,i2] ]
        W.extend(windows)
        for i,w in enumerate(W):
            i1,i2 = w
            x = numpy.array(range(i1,i2))
            y = numpy.array(wf[i1:i2])
            plt.plot(x,y,colorpoint[i])
        plt.show()
        #plt.savefig(figpdf)
        #print 'lsqa.drawPulse Wrote',figpdf
        return
    def updatePSDcut(self,fn = None):
        '''
         #'Samples/P50-1/run3568649727/Figures/optimize_run3568649727.root'):

        return optimal PSD cut based on past (or near in time) processing of same sample
        open fitted hist, extract function, determine position of valley between
        two peaks
        '''
        debug = False
        PSDcut = self.PSDcut 
        rfn = self.findPreviousOutput(fn=fn)
        if rfn is None: return PSDcut # no file found
        print 'lsqa.updatePSDcut For update, use',rfn
        f = ROOT.TFile(rfn,'r')
        favFast = self.favFast
        favTotal= self.favTotal
        Qlo,Qhi = self.QforAvg
        hname = 'PSD_Q'+ str(int(Qlo)).zfill(2) + '_' + str(int(Qhi)).zfill(2)
        hname+= 'f' + str(int(favFast)) + '_t' + str(int(favTotal)) + 'N'
        h = f.Get(hname)
        ff = None
        for a in h.GetListOfFunctions():
            if 'TF1' in a.ClassName(): ff = a
        if ff is None:
            print 'lsqa.updatePSDcut WARNING PSD cut NOT updated. Could not find function in hist',hname
            return PSDcut
        if debug: print 'lsqa.updatePSDcut hist',h.GetName(),'function',ff.GetName(),'Npar',ff.GetNpar()

        if ff.GetName()=='GG' and ff.GetNpar()==6:
            xmi,xma = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
            xmi,xma = ff.GetParameter(1),ff.GetParameter(1+3) # fitted means of gaussians
            nx = 100
            dx = (xma-xmi)/float(nx)
            xlo,FX = xmi,ff.Eval(xmi)
            for i in range(nx):
                x = xmi + (float(i)+0.5)*dx
                fx = ff.Eval(x)
                if i%10==0 and debug : print 'ls.updatePSDcut scan x',x,'func(x)',fx
                if fx<FX : FX,xlo = fx,x
            PSDcut = xlo
            if debug: print 'lsqa.updatePSDcut minimum of',FX,'at updated PSDcut',xlo
                
        return PSDcut

        
if __name__ == '__main__':
    if sys.argv[1]=='SIMPLE':
        L = lsqa(SIMPLE=True)
        L.simple()
        sys.exit('')

    if sys.argv[1]=='PICKLE':
        L = lsqa(SIMPLE=True)
        pfn = 'pickleLSAQfile.p'
        import pickle
        summaryInfo =     pickle.load(open(pfn,'rb'))
        print 'read pickled summaryInfo from',pfn
        L.writeRSF(summaryInfo)
        sys.exit('END OF lsqa.writeRSF test')


    
    if len(sys.argv)<=1:
        sys.exit('lsqa.py ERROR! Must specify input filename with full path')
    fn = sys.argv[1]
    L = lsqa(fn=fn)
    L.main(filename=fn)


    if 0:
        L = lsqa()
        simple = True
        if simple:
            fn='/Users/djaffe/work/LSQAData/P20/run3566844679.h5' #'/Users/djaffe/work/LSQAData/P20/run3566553226.h5'
            fn='/Users/djaffe/work/LSQAData/LiLS01/run3567106894.h5' 
            L.main(fn=fn)
            #L.simple(fn=fn)
        else:

            nfn = '/Users/djaffe/work/LSQAData/P20/run3566553226.h5'
            gfn = '/Users/djaffe/work/LSQAData/P20/run3566498398.h5'
            L.onlPlot(fn=[gfn,nfn],dname=['G','N'])#dname=['AmBe','AmBe_137Cs'])

            nfn = '/Users/djaffe/work/LSQAData/P20/run3566576160.h5'
            L.onlPlot(fn=[gfn,nfn],dname=['G','N'])#dname=['AmBe','AmBe_137Cs'])

    
