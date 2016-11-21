#!/usr/bin/env python
'''
simulation to study how reference sample can be used to estimate adsorption in material sample
20161026
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

class adsorpMC():
    def __init__(self,decayCon=100000.):
        self.dt = 10. # minutes per run
        self.dtsec = self.dt*60. # seconds per run
        self.rate = 100. # Hz

        self.refRate0 = 100.
        self.samRate0 = 96.8/100.*self.refRate0
        self.samDecay = 1./decayCon # minutes decay const

        adConInvDays = self.samDecay * 24.*60.
        print 'adsorpMC.__init__ sample adsorption constant (1/days) ',adConInvDays
        self.adConInvDays = adConInvDays

        self.RefX,self.RefY = None,None

        self.Interpolate = True
        
        self.samOscAmp = 0.01

        self.plotFits = False
        self.theFits = {}
        self.Duration = None

        return
    def refFcn(self,t):
        f = self.refRate0
        evts = numpy.random.poisson(f * self.dtsec)
        return float(evts)
    def samFcn(self,t,nfcn=0):
        '''
        time-dependence of sample
        default is exponential decrease for nfcn<=0
        '''
        f = self.samRate(nfcn)*math.exp(-t*self.samDecay * float(nfcn) )

        if 0: print nfcn,t,f    
        evts = numpy.random.poisson( f  * self.dtsec )
        return float(evts)
    def samRate(self,nfcn=0):
        '''
        set overall rate of sample
        '''
        samR = self.samRate0 * (1.0 - float(nfcn)/100.)
        return samR
    def collect(self,duration):
        t0 = 0.
        t =  t0
        et = self.dt/2.
        i = 0
        ref,sam = [],[]
        while (t-t0<duration):
            if i%2==0:
                r = self.refFcn(t)
                ref.append( (t+et,r) )
            else:
                s = self.samFcn(t)
                sam.append( (t+et,s) )
            t += self.dt
            i += 1
        return ref,sam
    def collectMany(self,duration,Nsamples=1,t0=0.,tstart=0.,Measurements={}):
        '''
        Measurements = {} # key,value=sample#, (time,counts). sample#0 == reference
        '''
        #t = t0 = 0.
        t = tstart
        et = self.dt/2.
        if len(Measurements)==0:
            for i in range(Nsamples+1): Measurements[i] = []
        i,j = 0,0
        while (t-tstart<duration):
            if i%2==0:
                r = self.refFcn(t-t0)
                Measurements[0].append( (t+et-t0,r) )
            else:
                j += 1
                if j>Nsamples: j = 1
                s = self.samFcn(t-t0,nfcn=j-1)
                Measurements[j].append( (t+et-t0,s) )
            i += 1
            t += self.dt
        if 0:
            for key in sorted(Measurements.keys()):
                print 'adsorpMC.collectMany key',key,'Msmt',Measurements[key]
        return Measurements
    def fitRef(self,ref,order=1):
        '''
        fit reference data to polynomial of order=order
        '''
        x,y = numpy.array([a for a,b in ref]),numpy.array([b for a,b in ref])
        par = numpy.polyfit(x,y,order)
        if self.plotFits:
            xf = numpy.linspace(0.,max(self.RefX),1000) #self.Duration)
            yf = numpy.polyval(par,xf)
            self.theFits[0] = [ xf,yf ]
        return par
    def func1(self,t,A0):
        if self.Interpolate:
            y = numpy.array([])
            for x in t:
                y = numpy.append(y, self.interpRef(x) * A0)
            f = y
            
        else:
            f = numpy.polyval(self.baseF,t)*A0
        #print 'adsorpMC.func1 f',f
        return f
    def interpRef(self,t):
        '''
        return interpolated value from reference data for points neighboring input data point t.
        If t has only a single neighbor, then use ordinate value of that neighbor
        '''
        X,Y = self.RefX,self.RefY
        
        j = None
        for i,x in enumerate(X):
            if t<x and i==0:
                j = i
                break
            if t>x and i+1==len(X):
                j = i
                break
            if x<=t and t<=X[i+1]:
                j = i + 1
                break
        if j is None:
            sys.exit('adsorpMC.interpRef ERROR no range found. t '+str(t)+' i '+str(i))
        if i==j:
            y = Y[i]
        else:
            y = ((Y[j]-Y[i])*t + X[j]*Y[i]-X[i]*Y[j])/(X[j]-X[i])
        return y
    def chi2Sam(self,sam,baseF,debug = False,nfcn=-1):
        '''
        determine chi2 for sample given base function baseF
        Fit function specified by baseF to sample data to extract constant, then determine chisquare of fit
        '''
        
        x,y = numpy.array([a for a,b in sam]),numpy.array([b for a,b in sam])
        if debug : print 'adsorpMC.chi2Sam sam',sam,'x',x,'y',y
        self.baseF = baseF

        # do the fit
        [a0,covA] = curve_fit(self.func1,x,y)
        if debug : print 'adsorpMC.chis2Sam a0',a0

        if self.plotFits:
            if nfcn<=0 : sys.exit('adsorpMC.chi2Sam ERROR nfcn='+str(nfcn)+' should be >0')
            xf = numpy.linspace(0.,max(self.RefX),1000) #self.Duration,100)
            yf = self.func1(xf,a0)
            self.theFits[nfcn] = [ xf,yf ]

        # evaluate chisquare
        observed = y
        expected = self.func1(x,a0) 
        if debug : print 'adsorpMC.chis2Sam observed',observed
        if debug : print 'adsorpMC.chis2Sam expected',expected
        ndf = len(x)-1
        chisq,pvalue = chisquare(observed, expected, -1 )
        if debug : print 'adsorpMC.chis2Sam chisq,pvalue,ndf',chisq,pvalue,ndf
        return chisq,ndf,pvalue
    def expt(self,duration=120.,faker=False):
        '''
        perform a toy experiment duration minutes long
        '''
        ref,sam = self.collect(duration)
        if faker :   
            ref2,sam2 = self.collect(duration)
            sam = ref2
        self.RefX,self.RefY = numpy.array([a for a,b in ref]),numpy.array([b for a,b in ref])
        par = aMC.fitRef(ref)
        baseF = par/par[1]
        chi2,ndf,pvalue = aMC.chi2Sam(sam,baseF)
        return chi2,pvalue
    def exptMany(self,Nsamples=1,duration=120.,days=5,plot=False):
        '''
        perform a toy experiment duration minutes long with Nsamples samples each day for days
        '''
        Msmts = {}
        
        for day in range(days):
            tstart = float(day)*24.*60.
            Msmts = self.collectMany(duration,Nsamples=Nsamples,Measurements=Msmts,tstart=tstart)
        Results = {}
        ref = Msmts[0]
        self.RefX,self.RefY = numpy.array([a for a,b in ref]),numpy.array([b for a,b in ref])
        #print 'adsorpMC.exptMany min(self.RefX),max(self.RefX)',min(self.RefX),max(self.RefX)
        par = self.fitRef(ref)
        baseF = par/par[1]
        for key in sorted(Msmts.keys()):
            if key>0:
                same = Msmts[key]
                chi2,ndf,pvalue = self.chi2Sam(same,baseF,nfcn=key,debug=False)
                #print 'adsorpMC.exptMany key,chi2,ndf,pvalue',key,chi2,ndf,pvalue
                Results[key] = [chi2,pvalue]
        if plot: self.plotMany(Msmts,Results=Results)
        return Results
    def toy(self,Nexpt=100,duration=120.,faker=False):
        '''
        perform Nexpt toy experiments, each experiment duration minutes long
        '''
        chi,p = [],[]
        for iexpt in range(Nexpt):
            chi2,pvalue = self.expt(duration=duration,faker=faker)
            chi.append(chi2)
            p.append(pvalue)
        chi2 = sum(chi)/float(len(chi))
        pvalue = sum(p)/float(len(p))
        print 'adsorpMC.toy',Nexpt,'expts. Mean chi2',chi2,'pvalue',pvalue,'faker',faker
        return
    def plot(self,ref,sam):
        '''
        show data in a figure
        '''
        x,y = numpy.array([a for a,b in ref]),numpy.array([b for a,b, in ref])
        u,v = numpy.array([a for a,b in sam]),numpy.array([b for a,b, in sam])
        plt.plot(x,y,'bo', label='Reference')
        plt.plot(u,v,'r+',label='Sample')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Counts')
        plt.grid(True)
        
        plt.savefig('/Users/djaffe/Desktop/example.pdf')
        return
    def plotMany(self,Msmts,Results=None):
        '''
        show data from many samples in figure
        '''
        colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        points = {0:'o', 1:'s',2:'D'} # filled circle, square, diamond
        lines  = {0:'--', 1:'-.',2:':'}
        ref = Msmts[0]
        x,y = numpy.array([a for a,b in ref]),numpy.array([b for a,b, in ref])
        xmi,xma,ymi,yma = min(x),max(x),min(y),max(y)
        plt.clf() # clear figure needed?
        lab = 'Reference'
        plt.plot(x,y,'ko',label=lab)
        for key in sorted(Msmts.keys()):
            if key>0:
                sam = Msmts[key]
                u,v = numpy.array([a for a,b in sam]),numpy.array([b for a,b, in sam])
                lab = 'Sample'+str(key)
                if Results is not None:
                    chi2,pvalue = Results[key]
                    #print 'key,chi2,pvalue',key,chi2,pvalue
                    if math.isnan(pvalue) : pvalue = 0.
                    lab = 'S'+str(key)+ ' ' + str( int(0.5+pvalue*10000.)/100. )
                c = colors[key%len(colors)]
                p = points[key/len(colors)]
                #print 'c+p',c+p
                xmi,xma,ymi,yma = min(xmi,min(u)),max(xma,max(u)),min(ymi,min(v)),max(yma,max(v))
                plt.plot(u,v,c+p,label=lab)
        nrow = 2
        if len(Msmts)>10: nrow = 3
        fontsize = 3*10/nrow
        ncol = (max(Msmts.keys())+1)/nrow+1
        #print 'len(Msmts)',len(Msmts),'max(Msmts.keys())',max(Msmts.keys()),'ncol',ncol
        plt.legend(bbox_to_anchor=(0.03, 1.02-.02, 1., .102), handletextpad=0., loc=3, ncol=ncol, mode="expand", borderaxespad=0., numpoints=1, fontsize=fontsize)

        if self.plotFits:
            for key in sorted(self.theFits.keys()):
                UV = self.theFits[key]
                #print 'key',key,'UV',UV
                u,v = UV
                c = colors[key%len(colors)]
                l = lines[key/len(colors)]
                plt.plot(u,v,c+l)
        
        plt.xlabel('Time (minutes)')
        plt.ylabel('Counts')
        words = 'Input adsorption constants (SampleNumber/(' + str(1./self.adConInvDays) + ' days)'
        plt.text(0.01*(xma-xmi)+xmi,1.05*(yma-ymi)+ymi,words)
        plt.grid(True)
        
        plt.savefig('/Users/djaffe/Desktop/exampleMany.pdf')
        return
    def toyMany(self,Nexpt=1, Nsamples=11, Times=4, days=1):
        '''
        main routine for run numerous toy experiments with many samples/experiment
        '''
        print '\n',Nexpt,'expts, ',Nsamples,'samples with time/msmt',self.dt,'and',Times,'msmts/sample/day for',days,'days'
        self.Duration = duration = 2*(Nsamples)*Times*self.dt
        self.plotFits = True
        allResults = {}
        for iexpt in range(Nexpt):
            plot = iexpt==0
            allResults[iexpt] = self.exptMany(Nsamples=Nsamples,duration=self.Duration,days=days,plot=plot)


        samChi2 = {}
        sampval = {}
        for isam in range(1,1+Nsamples):
            achi,ap = numpy.array([]),numpy.array([])
            for iexpt in range(Nexpt):
                results = allResults[iexpt]
                chi2,pvalue = results[isam]
                achi = numpy.append(achi,chi2)
                ap   = numpy.append(ap,pvalue)
            samChi2[isam] = achi
            sampval[isam] = ap
        font = {'size'   : 8}
        matplotlib.rc('font', **font)

        nrows,ncols = 4,3
        for words,samDict in zip(['chi2','pvalue'],[samChi2,sampval]):
            fig,ax = plt.subplots(nrows=nrows,ncols=ncols)
            nbin = 100
            for isam in range(1,1+Nsamples):
                jsam = Nsamples+1-isam
                icol = 2-((jsam)%ncols)
                irow = jsam/ncols
                print 'adsorpMC.toyMany isam,jsam,icol,irow',isam,jsam,icol,irow

                ax[irow,icol].hist(samDict[isam],nbin,histtype='step')
                ax[irow,icol].set_title(words+' sample'+str(isam),y=0.8)
            pdf = '/Users/djaffe/Desktop/toyMany_'+words+'.pdf'
            plt.savefig(pdf)
            print 'adsorpMC.toyMany Wrote',pdf
                
        
            
        return
   

            
if __name__ == '__main__' :
    if len(sys.argv)>1:
        aMC = adsorpMC( float(sys.argv[1]) )
    else:
        aMC = adsorpMC()

    Nexpt = 1
    if len(sys.argv)>2: Nexpt = int(sys.argv[2])
    Nsamples = 11
    Times = 4
    days = 1
    aMC.toyMany(Nexpt=Nexpt,Nsamples=Nsamples,Times=Times,days=days)

    
    if 0:
        n = 1000
        aMC.toy(Nexpt=n)
        aMC.toy(Nexpt=n,faker=True)
    if 0:
        print '\nsingle expt'
        duration = 120.
        ref,sam = aMC.collect(duration)
        print 'ref',ref
        print 'sam',sam
        par = aMC.fitRef(ref)
        print 'par',par
        baseF = par/par[1]
        print 'baseF',baseF
        chi2,ndf,pvalue = aMC.chi2Sam(sam,baseF,debug=True)
        print 'chi2,ndf,pvalue,decayConstant',chi2,ndf,pvalue, aMC.samDecay
        aMC.plot(ref,sam)
    if 0:
        print '\nsingle expt, many samples'
        Nsamples = 11
        Times = 4
        aMC.Duration = duration = 2*(Nsamples)*Times*aMC.dt
        aMC.plotFits = True
        if 0:
            Msmts = aMC.collectMany(duration,Nsamples=Nsamples)
            print 'duration(min)',duration,' ',str(int(duration/60.))+'h'+str(int(duration-float(60.*int(duration/60.))))+'m'
            for i in sorted(Msmts.keys()):
                print 'msmt',i,'#',len(Msmts[i])
            aMC.plotMany(Msmts)
        else:
            results = aMC.exptMany(Nsamples=Nsamples,duration=aMC.Duration)
    
        