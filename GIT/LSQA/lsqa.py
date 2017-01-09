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
        FOMmin,FOMmax = 0.,3.
        
        if type(fn) is not list: fn = [fn]
        hists = []
        graphs = []
        icol = 0
        for filename,w in zip(fn,dname):
        
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
            G = ROOT.TF1("G","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])")
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

            X,dX,Y,dY = [],[],[],[]
            
            for lo,hi in zip(Qlo,Qhi):
                A = (lo<=Qtot)*(Qtot<hi)
                title = name = 'PSD_Q'+str(int(lo)).zfill(2)+'_'+str(int(hi)).zfill(2) + words
                h = self.gU.makeTH1D(PSD[A],title,name,nx=ny,xmi=ymi,xma=yma)
                hists.append(h)

                Nevt = h.GetEntries()
                FOM = None
                if Nevt>FOMfitThres:
                
                    GG = ROOT.TF1("GG","[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])")
                    GG.SetParName(0,"gN")
                    GG.SetParName(1,"gMean")
                    GG.SetParName(2,"gSigma")
                    GG.SetParName(3,"nN")
                    GG.SetParName(4,"nMean")
                    GG.SetParName(5,"nSigma")
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
                
                print 'lsqa.onlPlot {0} lo,hi,FOM {1:.1f} {2:.2f} {3}'.format(words,Elo,Ehi,Units),
                if FOM is None:
                    print FOM
                else:
                    X.append(0.5*(Ehi+Elo))
                    dX.append(0.5*(Ehi-Elo))
                    Y.append(FOM)
                    dY.append(dFOM)
                    FOMmin = min(FOM-dFOM,FOMmin)
                    FOMmax = max(FOM+dFOM,FOMmax)
                    print ' {0:.2f}({1:.2f})'.format(FOM,dFOM)
            
            f.close()

        name = 'FOM_vs_keVee_'+bn
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,Y,title,name,ex=dX,ey=dY)
        self.gU.color(g,0,0,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir,option='AP',yLimits=[FOMmin,FOMmax])
        graphs.append(g)

        nameForFile = '_'.join(allBN)

        hists1d,hists2d = [],[]
        for h in hists:
            if h.GetDimension()==1: hists1d.append(h)
            if h.GetDimension()==2: hists2d.append(h)
        
        
        self.gU.drawMultiHists(hists1d,nameForFile+'_1d',figdir=self.figdir,forceNX=4,Grid=False,biggerLabels=True,fitOpt=111,statOpt=0)
        self.gU.drawMultiHists(hists2d,nameForFile+'_2d',figdir=self.figdir,forceNX=2,Grid=True)

        RFN = self.figdir + nameForFile + '.root'
        rfn = ROOT.TFile.Open(RFN,'RECREATE')
        for h in hists: rfn.WriteTObject(h)
        for g in graphs:rfn.WriteTObject(g)
        rfn.Close()
        print 'lsqa.onlPlot Wrote',len(hists),'hists and',len(graphs),'graphs to',RFN
            
        return
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
    def simple(self):
        fn = '/Users/djaffe/work/LSQAData/Test1/run3554994172.h5'
        nsd = 3.0
        ph,ct = self.getWFMs(fn)
        print 'len(ph),len(ct)',len(ph),len(ct)
        for ip in range(0,len(ph),281):
            wf = ph[ip]
            ped,pedsd,iPulse,subPperP,pArea,pTime = self.wfa.pulseAnal(wf,'LSQA',debug=1,nsd=nsd)
            print 'ped,pedsd,iPulse,subPperP,pArea,pTime',ped,pedsd,iPulse,subPperP,pArea,pTime
            if len(iPulse)>0:
                for j,p in enumerate(iPulse):
                    title = 'Evt '+str(ip)+ ' pulse#'+str(j)+' of '+str(len(iPulse)-1)+' t '+str(pTime[j])
                    self.drawPulse(wf,window=iPulse[j],title=title)
        return       
    def drawPulse(self,wf,window=None,title='waveform'):
        '''
        draw a waveform
        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        colorpoint = {0: 'bo', 1:'ro'}
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
    #    L.simple()
    nfn = '/Users/djaffe/work/LSQAData/P20/run3566553226.h5'
    gfn = '/Users/djaffe/work/LSQAData/P20/run3566498398.h5'
    L.onlPlot(fn=[gfn,nfn],dname=['G','N'])#dname=['AmBe','AmBe_137Cs'])

    nfn = '/Users/djaffe/work/LSQAData/P20/run3566576160.h5'
    L.onlPlot(fn=[gfn,nfn],dname=['G','N'])#dname=['AmBe','AmBe_137Cs'])
    
    
