#!/usr/bin/env python
'''
analysis of scope data for LS QA
adapted from Lindsey's PSDscripts
20161003
'''
import sys
# the ONETON dir contains wfanal.py 
sys.path.append('/Users/djaffe/work/GIT/ONETON')


import wfanal,graphUtils
import h5py
import matplotlib.pyplot as plt
import numpy
import ROOT
import os
import math

class lsqa():
    def __init__(self):
        self.wfa = wfanal.wfanal()
        self.gU  = graphUtils.graphUtils()
        self.figdir = 'Figures/'

        self.GG = ROOT.TF1("GG","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])")
        self.G = ROOT.TF1("G","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])")
        return
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
    def getQ(self,wfm,fastValues,totalValues,startoff=100,numavg=200):
        minidx = numpy.argmin(wfm)
        BL = sum(wfm[int(minidx-startoff-numavg):int(minidx-startoff)])/float(numavg)
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
        determine type of file from length. ~20k = neutron only, ~5k = gamma and neutron source
        if type is ambiguous or unknown, then check in LSQA_run_source.txt
        '''
        f = None
        if type(F) is h5py._hl.files.File: f = F
        if type(F) is str: f = h5py.File(F,'r')
        
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
            if not (gamma or neutron): 
                print 'lsqa.runType runtype unknown for',f.file
                return None
        if gamma: return 'G'
        if neutron: return 'N'
        return None
    def wfanaFile(self,fn=None,fastValues=None,totalValues=None):
        '''
        Return events containing Qfast,Qtot pairs for different fast, total definitions
        given input file and ranges of values
        '''
        if fn is None or fastValues is None or totalValues is None:
            sys.exit('lsqa.wfanaFile ERROR Invalid input')

        idebug = 0
        IntValues = [[x,y] for x in fastValues for y in totalValues]
        f = h5py.File(fn)
        print 'lsqa.wfanaFile Opened',fn
        Events = []

        ievt,freq = 0,500
        for w in f['Waveforms']:
            wfm,ct = w[0]
            if len(wfm)>0:
                results = self.getQ(wfm,fastValues=fastValues,totalValues=totalValues)
                Events.append( results )
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
        return Events
    def main(self,fn=None):
        '''
        processing file
        determine type of file from length. ~20k = neutron only, ~5k = gamma and neutron source
        turn waveforms into sets of Qfast,Qtot pairs for different definitions of fast and tot
        then evaluate FOM as function of Qtot
        Best global values appear to be fast=100,total=600.
        '''
        idebug = 0
        fastValues = numpy.linspace(100,140,3) #80,140,4) #60,140,5) #20,120,6)
        totalValues= numpy.linspace(400,1000,4) #200,1000,5)
        IntValues = [[x,y] for x in fastValues for y in totalValues]
        bn = os.path.basename(fn).replace('.h5','')
        
        # waveforms to sets of Qfast,Qtot pairs
        Events = self.wfanaFile(fn=fn,fastValues=fastValues,totalValues=totalValues)
        words = self.runType(fn)  # should be either 'G' (for gamma+neutron source) or 'N' (for neutron source only)

        if idebug>0:
            print 'lsqa.main Events',Events
        Qlo = numpy.linspace(5,65,7)
        Qhi = numpy.array([x+10. for x in Qlo])
        EperQ = {}
        sumFOM = {}
        E = numpy.array(Events)
        hists,allhists,graphs = [],[],[]
        tmg = self.gU.makeTMultiGraph('FOMerific')
        nx,xmi,xma = 100,0.,100.
        ny,ymi,yma = 100,0.,1.
        for i,vv in enumerate(IntValues):
            ifast,itot = vv
            PSDdef = name = 'f'+str(int(ifast))+'_t'+str(int(itot))
            title = 'PSD vs Qtot '+name.replace('_',' ')
            h = ROOT.TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
            aQtot,aPSD = [],[]
            for pair in E[:,i]:
                Qfast,Qtot = pair
                PSD = (Qtot-Qfast)/Qtot
                h.Fill(abs(Qtot),PSD)
                aQtot.append(abs(Qtot))
                aPSD.append(PSD)
            hists.append(h)
            aQtot = numpy.array(aQtot)
            aPSD  = numpy.array(aPSD)
            FOM,dFOM,FOMhists,FOMgraph = self.beerFOM(aQtot,aPSD,Qlo,Qhi,EperQ,PSDdef+words,bn,Draw=True)
            sFOM,sdFOM = numpy.sum(FOM),numpy.sqrt(numpy.sum([x*x for x in dFOM]))
            sumFOM[PSDdef+words] = (sFOM,sdFOM)
            print 'lsqa.main {0} sum {1:.2f}({2:.2f})'.format(PSDdef+words,sFOM,sdFOM) 
            allhists.extend( FOMhists )
            graphs.append( FOMgraph )
            self.gU.color(FOMgraph,i,i,setMarkerColor=True)
            tmg.Add( FOMgraph )

        print 'lsqa.main Sorted FOMs'
        for g in sorted( sumFOM.items(), key=lambda x: x[1]):
            print 'lsqa.main {0} {1:.2f}({2:.2f})'.format(g[0],g[1][0],g[1][1])
        self.gU.drawMultiHists(hists,'optimize',figdir=self.figdir,forceNX=5,Grid=True)
        tmg.SetMinimum(0.)
        tmg.SetMaximum(3.0)
        canvas = self.gU.drawMultiGraph(tmg,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Charge',yAxisLabel='FOM',NLegendColumns=3)
        graphs.append( tmg )
        
        rfn = self.figdir + 'optimize.root'
        f = ROOT.TFile.Open(rfn,'RECREATE')
        for h in hists: f.WriteTObject(h)
        for h in allhists: f.WriteTObject(h)
        for g in graphs: f.WriteTObject(g)
        f.Close()
        print 'lsqa.main Wrote',len(hists)+len(allhists),'hists and',len(graphs),'graphs to',rfn
        return

        
    def onlPlot(self,fn=None,dname=None):
        '''
        get Qfast,Qtot from online analysis and plot 'em
        '''
        Cs137edge = 478. # keV
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
            FOM,dFOM,FOMhists,FOMgraph = self.beerFOM(Qtot,PSD,Qlo,Qhi,EperQ,PSDdef+words,bn,Draw=True)
            
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
        return arrays of figure-of-merit and uncertainty, a list of fitted histograms and a graph of FOM vs Q (or E)
        given paired Qtot,PSD arrays, lower and upper charge limits Qlo,Qhi.
        EperQ, if available, converts charge to energy.
        PSDdef is string defining the PSD and bn is basename.
        PSDdef & bn are used for histogram, graph and filenames if Draw is True
        '''
        #FOM,dFOM,FOMhists,FOMgraph = self.beerFOM(Qtot,PSD,Qlo,Qhi,EperQ,words,Draw=True)

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
        ny,ymi,yma = 110,-.1,1.

        FOMmin,FOMmax = 0.,3.
        
        X,dX,Y,dY = [],[],[],[]            
        for lo,hi in zip(Qlo,Qhi):
            A = (lo<=Qtot)*(Qtot<hi)
            title = name = 'PSD_Q'+str(int(lo)).zfill(2)+'_'+str(int(hi)).zfill(2) + PSDdef
            h = self.gU.makeTH1D(PSD[A],title,name,nx=ny,xmi=ymi,xma=yma)
            hists.append(h)

            Nevt = h.GetEntries()
            FOM = None
            if Nevt>FOMfitThres:

                peak = h.GetMaximum()
                mean = h.GetMean()
                gMean = mean - 0.05
                nMean = mean + 0.05
                GG.SetParameters(peak, gMean, .02, peak, nMean, .02)
                ptr = h.Fit(GG,"SLQ")

                #ptr.Print("V") # print everything about fit
                gMean,gSigma = GG.GetParameter(1),GG.GetParameter(2)
                nMean,nSigma = GG.GetParameter(1+3),GG.GetParameter(2+3)
                FWHM = 2.*math.sqrt(2.*math.log(2.))*math.sqrt(gSigma*gSigma + nSigma*nSigma)
                FOM = abs(gMean-nMean)/FWHM

                par = []
                for ii in range(6): par.append( GG.GetParameter(ii) )
                par = numpy.array(par)
                #print 'par',par
                dFOM  =self.getFOMunc(ptr,par,FOM)

            Elo,Ehi,Units = lo,hi,'charge'
            if 'G' in EperQ: Elo,Ehi,Units = EperQ['G']*lo,EperQ['G']*hi,'keV'

            if debug: print 'lsqa.beerFOM {0} lo,hi,FOM {1:.1f} {2:.2f} {3}'.format(PSDdef,Elo,Ehi,Units),
            X.append(0.5*(Ehi+Elo))
            dX.append(0.5*(Ehi-Elo))
            if FOM is None:
                Y.append(0.)
                dY.append(10.)
                if debug: print FOM
            else:
                Y.append(FOM)
                dY.append(dFOM)
                FOMmin = min(FOM-dFOM,FOMmin)
                FOMmax = max(FOM+dFOM,FOMmax)
                if debug: print ' {0:.2f}({1:.2f})'.format(FOM,dFOM)
            
        name = 'FOM_v_QE_'+PSDdef+bn
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,Y,title,name,ex=dX,ey=dY)
        self.gU.color(g,0,0,setMarkerColor=True)
        if Draw: self.gU.drawGraph(g,figDir=self.figdir,option='AP',yLimits=[FOMmin,FOMmax])


        
        if Draw:
            nameForFile = bn+'_'+PSDdef

            hists1d,hists2d = [],[]
            for h in hists:
                if h.GetDimension()==1: hists1d.append(h)
                if h.GetDimension()==2: hists2d.append(h)

            print 'len(hists),len(hists1d),len(hists2d)',len(hists),len(hists1d),len(hists2d)
            self.gU.drawMultiHists(hists1d,nameForFile+'_1d',figdir=self.figdir,forceNX=4,Grid=False,biggerLabels=True,fitOpt=111,statOpt=0)
            self.gU.drawMultiHists(hists2d,nameForFile+'_2d',figdir=self.figdir,forceNX=2,Grid=True)

        
        return Y,dY,hists,g
    def getFOMunc(self,ptr,pars,FOM):
        '''
        return uncertainty in FOM given pointer to fit result and best fit parameters '''

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
    def simple(self,fn = '/Users/djaffe/work/LSQAData/Test1/run3554994172.h5'):
        usePulseAnal = False
        startoff, fastint, totint, numavg = 100, 200, 500, 200
        nsd = 3.0
        ph,ct = self.getWFMs(fn)
        print 'len(ph),len(ct)',len(ph),len(ct)
        for ip in range(0,len(ph),281):
            wf = ph[ip]
            if usePulseAnal:
                ped,pedsd,iPulse,subPperP,pArea,pTime = self.wfa.pulseAnal(wf,'LSQA',debug=1,nsd=nsd)
                print 'ped,pedsd,iPulse,subPperP,pArea,pTime',ped,pedsd,iPulse,subPperP,pArea,pTime
                words = 'pulseAnal'
            else:
                minidx = numpy.argmin(wf)
                iPulse = [[minidx-startoff, minidx+ totint]]
                pTime = [float(minidx)]
                words = 'simple'
            if len(iPulse)>0:
                for j,p in enumerate(iPulse):
                    title = words+ ' Evt '+str(ip)+ ' pulse#'+str(j)+' of '+str(len(iPulse)-1)+' t '+str(pTime[j])
                    self.drawPulse(wf,window=iPulse[j],title=title)
        return       
    def drawPulse(self,wf,window=None,title='waveform'):
        '''
        draw a waveform
        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        colorpoint = {0: 'b.', 1:'ro'}
        i1,i2 = 0,len(wf)
        loop = 1
        if window is not None: loop = 2
        for l in range(loop):
            if window is not None:
                i1,i2 = window
                if l==0:
                    d1 = max(200,i2-i1)
                    d2 = max(400,i2-i1)
                    i1 = max(0,i1-d1)
                    i2 = min(len(wf),i2+d2)
            #print l,i1,i2
            x = numpy.array(range(i1,i2))
            y = numpy.array(wf[i1:i2])
            plt.plot(x,y,colorpoint[l])
        plt.show()
        return
        
        
if __name__ == '__main__':
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
    
    
