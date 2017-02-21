#!/usr/bin/env python
'''
calibrate wavedump histogram files for 227Ac setup
20161220
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

import makeTTree

class calibWaveDump():
    def __init__(self,Log=True):


        self.gU = graphUtils.graphUtils()
        self.gfit = gfit.gfit()
        self.rLF = readLogFile.readLogFile()
        self.cAC = cutsAndConstants.cutsAndConstants()
        self.lt = livetime.livetime()
        self.gfp = get_filepaths.get_filepaths()
        self.mTT = makeTTree.makeTTree()

        self.psdCut = self.cAC.psdCut
        self.lowChargeCut = self.cAC.lowChargeCut

        
        name = 'calibWaveDump'
        self.rootFileDir = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        self.logFileDir  = self.rootFileDir.replace('rootfiles','logfiles')
        self.outFileDir  = 'Output/'
        self.figdir      = 'Figures/'+name+'/'
        self.perrunfigdir= self.figdir + 'perRun/'

        self.logdir = 'Logfiles/'+name+'/'
        now = datetime.datetime.now()
        fmt = '%Y%m%d_%H%M%S_%f'
        self.start_time = cnow = now.strftime(fmt)

        self.outTreeDir = 'CalibTrees/'
        self.outTree    = self.outTreeDir + 'tree_'+cnow+'.root'
        
        if Log:
            lfn = self.logdir + cnow + '.log'
            sys.stdout = Logger.Logger(fn=lfn)
            print 'procWaveDump__init__ Output directed to terminal and',lfn
            print 'procWaveDump__init__ Job start time',self.start_time

        
        return
    def get_filepaths(self,directory,exclude=None):
        return self.gfp.get_filepaths(directory,exclude=exclude)
    def subAndFit(self,hs,hb,dname,nsig=3.):
        '''
        return fitted difference of signal hs and background hb histograms with gaussian
        '''
        hd = hs.Clone(dname)
        hd.Add(hb,-1.)
        

        GoodFit,mean,emean, sgm, prob = self.gfit.fit(hd,nsigplus=nsig,nsigminus=nsig,debug=False,quiet=True)
        xmi,xma = mean-nsig*sgm,mean+nsig*sgm
        S,E = self.getIntegral(hd,xmi,xma)
        return hd,mean,emean,sgm,S,E
    def getIntegral(self,hd,xmi,xma):
        '''
        Return integral and uncertainty of 1d hist hd in range (xmi,xma)
        Use fractional contents and error for first,last bin in range
        '''
        i1 = hd.GetXaxis().FindBin(xmi)
        i2 = hd.GetXaxis().FindBin(xma)
        nx = hd.GetNbinsX()
        err = ROOT.Double(-1.)
        S = hd.IntegralAndError(i1+1,i2-1,err)
        err = err*err
        for i,x in zip([i1,i2],[xmi,xma]):
            xc = hd.GetBinCenter(i)
            c  = hd.GetBinContent(i)
            e  = hd.GetBinError(i)
            w  = hd.GetXaxis().GetBinWidth(i)
            if i==i1: f = (xc+w/2.-x)/w
            if i==i2: f = (x-(xc-w/2.))/w
            S += f*c
            err += (f*e)*(f*e)
        err = math.sqrt(err)
        #print 'S,err',S,err
        return S,err
    def getStats(self,h,xmi=None,xma=None):
        '''
        return mean,std. dev., median and number of entries for hist h in range [xmi,xma]
        If xmi or xma is None, then use first or last bin, resp.
        '''
        i1 = 1
        if xmi is not None: i1 = h.GetXaxis().FindBin(xmi)
        i2 = h.GetNbinsX()
        if xma is not None: i2 = h.GetXaxis().FindBin(xma)
        x = numpy.array([])
        N = 0
        for i in range(i1,i2+1):
            xc = h.GetBinCenter(i)
            c = h.GetBinContent(i)
            x = numpy.append(x, numpy.full(c,xc))
        mean = numpy.mean(x)
        median = numpy.median(x)
        stddev = numpy.std(x)
        nument = len(x)
        return mean,stddev,median,nument
    def getPstats(self,fn='Output/run00073.root'):
        '''
        return mean,std dev, median & number of entries in prompt charge histogram
        in input file
        '''
        if not os.path.isfile(fn): return
        hname = 'Charge_PSD_cut'
        hf = ROOT.TFile.Open(fn,'r')
        h = hf.Get(hname)
        mean,sd,med,n = self.getStats(h,xmi=self.lowChargeCut)
        hf.Close()
        return mean,sd,med,n
    def getGstats(self,fn='Output/run00037.root'):
        '''
        return estimate of effy & unc. of goodCh0 requirement
        '''
        if not os.path.isfile(fn): return
        hname = 'goodCh0'
        hf = ROOT.TFile.Open(fn,'r')
        h = hf.Get(hname)
        nogood = h.GetBinContent(1)
        good   = h.GetBinContent(2)
        effy,err = -1.,-1.
        if good+nogood>0:
            effy = float(good)/float(good+nogood)
            err  = 1.-effy
        return effy,err
    def getLife(self,fn='Output/run00073',withPromptCut=False,fixLifetime=False):
        '''
        return normalization (# of events) and fitted lifetime (and uncertainties)
        from delayed energy vs tDelay-dPrompt distribution using off-time window
        background subtraction
        At present delayed energy vs deltaT runs from deltaT = (0,5*Po215lifetimes)
        from max(lifeRange in cutsAndConstants.py).
        20170105 return normalization when lifetime is fixed to nominal lifetime,
        but fitted lifetime when both normalization and lifetime are free fit
        '''
        if not os.path.isfile(fn): return None

        hf = ROOT.TFile.Open(fn,'r')
        bn = os.path.basename(fn).replace('.root','')
        dname,oname = 'dvdt','ovdt'
        if withPromptCut: dname,oname = 'dvdtc','ovdtc'
        hd,ho = hf.Get(dname),hf.Get(oname)
        hs = hd.Clone("hs")
        hs.SetTitle("Qdelay_vs_dt_bkgd_sub")
        hs.Add(ho,-1.)
        hist2d = [hd,ho,hs]

        nx = hs.GetXaxis().GetNbins()
        hs_px = hs.ProjectionX("hs_px",1,nx,"e")
        hs_px.SetTitle("Po215 events vs dt")
        w = hs.GetXaxis().GetBinWidth(2)
        hist1d = [hs_px]

        # starting value of N for fit to N/tau * exp(-t/tau)
        xmi = hs.GetXaxis().GetBinLowEdge(1)
        xma = hs.GetXaxis().GetBinLowEdge(nx+1) #
        sum,dsum = self.getIntegral(hs_px,xmi,xma)
        #print 'nx,xmi,xma,sum,dsum',nx,xmi,xma,sum,dsum

        E = ROOT.TF1("E","[0]*exp(-x/[1])/[1]*"+str(w))
        E.SetParName(0,"N")
        E.SetParName(1,"Lifetime")
        E.SetParameter(0,sum)
        E.SetParameter(1,self.cAC.Po215lifetime)
        
        # first fit with lifetime fixed to nominal. Q=quiet, 0 = do not plot, S = tfitresultptr should be returned
        E.FixParameter(1,self.cAC.Po215lifetime)
        fr0 = hs_px.Fit(E,"QS0")
        N0,dN0 = E.GetParameter(0),E.GetParError(0)

        # second fit with lifetime free (if requested)
        E.ReleaseParameter(1)
        E.SetParameter(1,self.cAC.Po215lifetime)
        if fixLifetime : E.FixParameter(1,self.cAC.Po215lifetime)

        ROOT.gStyle.SetOptFit(1111)
        ROOT.gStyle.SetOptStat(111110)
        # Q=quiet, L=likelihood, I=integral of function in bin, instead of value at bin center, S=should return TFitResultPtr , 0 = do not plot result of fit
        fr = hs_px.Fit(E,"QS")  # fit projection
        N,dN = E.GetParameter(0),E.GetParError(0)
        tau,dtau = E.GetParameter(1),E.GetParError(1)
        print 'calibWaveDump.getLife',bn,hs_px.GetTitle(),'N0',N0,'(',dN0,')','N',N,'(',dN,') lifetime',tau*1.e3,'(',dtau*1.e3,') ms'

        self.gU.drawMultiHists(hist1d,fname=bn+'_getLife1d',figdir=self.perrunfigdir,abscissaIsTime=False,dopt="hist e0 func",biggerLabels=False,statOpt=1001111,fitOpt=1111)
        twoD = False
        if twoD: self.gU.drawMultiHists(hist2d,fname=bn+'_getLife2d',figdir=self.perrunfigdir,abscissaIsTime=False,dopt="COLZ",statOpt=0,biggerLabels=False)


            

        hf.Close()
        return N0,dN0,tau,dtau
        
            
    def loop(self,fn='Output/run00073.root',withPromptCut=True):
        '''
        loop over relevant delayed hists in input root file and
        fit them to extract mean, error on mean, sigma of gaussian, 
        and calculate area of peak (Scorr) and uncertainty (Ecorr)
        return dict with key =lifetime cut, value = [mean,emean,sgm, Scorr,Ecorr]
        '''
        results = {}
        if not os.path.isfile(fn): return None
        hf = ROOT.TFile.Open(fn,'r')
        ROOT.gStyle.SetOptFit(1111)
        statX = ROOT.gStyle.GetStatX() # default
        ROOT.gStyle.SetStatX(.5)
        hlist = []
        for l in range(1,6):
            clife = str(l)
            life = float(l)
            if withPromptCut:
                sname = 'dc'+clife
                bname = 'oc'+clife
            else:
                sname = 'd'+clife
                bname = 'o'+clife
            if hf.FindKey(sname):
                hs = hf.Get(sname)
                hb = hf.Get(bname)
                dname = 'd-o'+clife
                if withPromptCut: dname = 'd-oc'+clife
                
                hd,mean,emean,sgm,S,E = self.subAndFit(hs,hb,dname)
                hd,hs,hb = self.commonOrdinate(hd,hs,hb)
                hlist.append([hd,hs,hb])
                f = 1./(1.-math.exp(-life))
                Scorr = S*f
                Ecorr = E*f
                #print fn,'l,mean,emean,sgm,Scorr,Ecorr {0} {1:.5f} {2:.5f} {3:.5f} {4:.1f} {5:.1f}'.format(l,mean,emean,sgm,Scorr,Ecorr)
                results[l] = [mean,emean,sgm,Scorr,Ecorr]


        dn = (fn.split('/')[1]).split('.')[0]
        if not withPromptCut : dn += 'withNoPromptCut'
        self.gU.drawMultiHists(hlist,fname=dn,figdir=self.perrunfigdir,statOpt=0,biggerLabels=False,dopt='hist e0 func')
        ROOT.gStyle.SetStatX(statX) # restore to default
        hf.Close()
        return results
    def commonOrdinate(self,h1,h2,h3):
        up,lo = h1.GetMaximum(),h1.GetMinimum()
        for h in [h1,h2,h3]:
            up = max(up,h.GetMaximum())
            lo = min(lo,h.GetMinimum())
        up = 1.1*up
        if lo<0.: lo = min(1.1*lo,-up/10.)
        for h in [h1,h2,h3]:
            h.SetMaximum(up)
            h.SetMinimum(lo)
        h1.SetLineColor(ROOT.kBlue)
        h2.SetLineColor(ROOT.kBlack)
        h3.SetLineColor(ROOT.kRed)
        return h1,h2,h3
    def fitPSD(self,fn='Output/run00073.root'):
        '''
        fit PSD dist with two gaussians and estimate efficiency of psd cut with two methods.
        First method: Fractional area above PSD cut based on fit to n-recoil peak
        Second method: Fraction histogram area above PSD cut by subtracting fitted e-recoil peak

        Final result is average of two methods with uncertainty taken as sum in quad. of two methods plus the difference between the two methods

        Notes:
        gfit.twoG returns fit parameters par,ga,ga
        where par = c1,mean1,emean1,sgm1, c2, mean2,emean2,sgm2, prob, gsumstat
                  = constant, mean, unc. on mean, sigma for each gaussian,
                    prob of chi2 of fit, status (='CONVERGED ' or not)
            ga,ga = TF1 objects with best fit parameters
        '''
        debug = False
        
        hname = 'psdc'
        psdCut = self.psdCut
        if not os.path.isfile(fn): return None
        hf = ROOT.TFile.Open(fn,'r')
        h = hf.Get(hname)
        par,ga,gb = self.gfit.twoG(h,x1min=0.,x1max=psdCut,x2max=1.,sgm=0.02,Unconstrained=True)

        # estimate fractional area of n-recoil peak above cut from gaussian fit
        G1 = lambda x: math.exp(-0.5*((x-par[1])*(x-par[1])/par[3]/par[3]))        
        G2 = lambda x: math.exp(-0.5*((x-par[1+4])*(x-par[1+4])/par[3+4]/par[3+4]))
        num = quad(G2,psdCut,1.)
        den = quad(G2,0.,1.)
        effy,erry = -1.,-1.
        if den[0]>0:
            effy = num[0]/den[0]
            if num[1]>0: erry = num[1]*num[1]/num[0]/num[0]
            erry = effy*math.sqrt(erry + den[1]*den[1]/den[0]/den[0])

        # estimate fractional area of n-recoil peak by subtracting gaussian fit of e-recoil peak
        hd = h.Clone('cloned_'+hname)
        hd.Add(ga,-1.)
        num = self.getIntegral(hd,psdCut,1.)
        den = self.getIntegral(hd,0.,1.)
        effy2,erry2 = -1.,-1.
        if den[0]>0:
            effy2 = num[0]/den[0]
            if num[1]>0: erry2 = num[1]*num[1]/num[0]/num[0]
            erry2 = effy2*math.sqrt(erry2 + den[1]*den[1]/den[0]/den[0])

        # overall estimate of effy is average, uncertainty is total in quadrature of two methods plus the difference between the two
        ef3,er3,nf = 0.,0.,0.
        for x,dx in zip( [effy,effy2], [erry,erry2]):
            if x>-1.:
                nf += 1.
                ef3 += x
            if dx>-1.: er3 += dx*dx
        if nf>1.:
            er3 += (effy-effy2)*(effy-effy2)
            er3 = math.sqrt(er3)
        if nf>0.: ef3 = ef3/float(nf)
            

        if debug: print 'calibWaveDump.fitPSD effy,erry',effy,erry,'effy2,erry2',effy2,erry2,'ef3,er3',ef3,er3
        hf.Close()      
        return effy,erry, effy2,erry2, ef3,er3
    def main(self,withPromptCut=True,useLastTime=True):
        '''
        main routine
        generate list of input files, find corresponding logfiles to get run info,
        get stats for a run,
        plot stats vs run # and vs time
        '''
        listOfHistFiles = self.get_filepaths(self.outFileDir,exclude='summary')
        listOfLogFiles  = self.get_filepaths(self.logFileDir)

        

        Threshold = 10. # minimum number of coincidences for a GOOD run
        Ac227only = True # only process runs with Ac227 (no Cs137, for example)
        if Ac227only:
            goodSources = ['Ac-227']
            badSources  = ['Cs-137','Ba-133']
        else:
            goodSources = ['Ac-227','Cs-137']
            badSources  = []
        print 'calibWaveDump.main Threshold',Threshold,'goodSources',goodSources,'badSources',badSources,'withPromptCut',withPromptCut,'useLastTime',useLastTime
            
        
        X,Y,dY,PoPeak,dPoPeak,PoS = [],{},{},{},{},{}
        PSD = [],[],[],[],[],[]
        Prompt = [],[],[],[]
        dPrompt= [],[],[],[]
        PromptName = ['Prompt_Mean','Prompt_StdDev','Prompt_Median','Prompt_Rate_Hz']
        LIFE  = [],[],[],[] # N,dN,tau,dtau
        LIFEName = ['Po215Hz','dPo215Hz','Fitted_Tau','dFitted_Tau']
        ecLIFE = [],[]
        ecLIFEName = ['effcPo215Hz','deffcPo215Hz']
        T = []  # time of start of run
        ER,dER = [],[] # expected rate and uncertainty of sample at start of run
        GoodCut = [],[]
        Run = []
        sampleName = []
            
            

        listOfRunsAS = []     # AS=AsStrings
        listOfGoodRunsAS = []
        Quick = False #True
        nQuick = 10
        
        for fn in listOfHistFiles:
            rn = os.path.basename(fn).split('.')[0] # run00xxx
            runNum = int(rn.replace('run','')) # integer run number
            listOfRunsAS.append(rn)
            lf = None
            for lfn in listOfLogFiles:
                if rn in lfn:
                    lf = lfn
                    break
            if lf is not None:
                timestamp,sources,sample,runtime = self.rLF.readFile(lf)
                if useLastTime:
                    lastTime = self.lt.getLastTime(rn)
                    #print rn,'lastTime',lastTime,'runtime',runtime,
                    if lastTime is not None: runtime = lastTime
                    #print 'fixed runtime',runtime
                GOOD = False
                if sources is not None:
                    for src in goodSources:
                        if src in sources: GOOD = True
                    for src in badSources:
                        if src in sources: GOOD = False
                if not GOOD: print 'calibWaveDump.main',rn,'no good sources'
                if runNum in self.cAC.badRuns :
                    GOOD = False
                    print 'calibWaveDump.main',rn,'bad run'
                if GOOD:
                    # loop returns results of fit to bkgd-subtracted delayed hist
                    # results[key] are mean,dmean,sigma, area,darea for
                    #   key=integer cut on coincidence window in number of lifetimes
                    results = self.loop(fn=fn,withPromptCut=withPromptCut)
                    for l in results:
                        if results[l][3]<Threshold : GOOD = False
                    if not GOOD: print 'calibWaveDump.main',rn,'too few coincidences'
                if GOOD:
                    listOfGoodRunsAS.append(rn)

                    nQuick -= 1
                    Run.append( runNum )
                    sampleName.append( sample )
                    
                    # estimate effy of GoodCh0 cut
                    g = self.getGstats(fn=fn)
                    for i,p in enumerate(GoodCut):
                        p.append(g[i])

                    # get statistics on prompt distribution, then normalize # of prompt
                    # events by run time
                    pStats  = self.getPstats(fn=fn) #mean,stddev,median,N
                    norm = (1.,1.,1.,runtime)
                    for i,p in enumerate(Prompt):
                        p.append( float(pStats[i])/norm[i] )
                        if i==0:
                            dPrompt[i].append( pStats[1]/math.sqrt(float(pStats[3]) ) )
                        if i==1 or i==2 : dPrompt[i].append( 0. )
                        if i==3:
                            dPrompt[i].append( math.sqrt(float(pStats[3]))/norm[i] )
                            
                    # fit PSD distribution to estimate PSD cut effy
                    # fitPSD returns  eff_est1, deff_est1, eff_est2, deff_est2, eff_combine,deff_combine
                    # eff_est1 = fractional fitted gaussian area above PSD cut
                    # eff_est2 = fractional bkgd-subtracted hist area above PSD cut
                    # eff_combine = average eff_est1 & eff_est2
                    psdEffy = self.fitPSD(fn=fn) 
                    for i,psd in enumerate(PSD):
                        psd.append( psdEffy[i] )
                        
                    # fit putative Po215 lifetime distribution to extract
                    # fitted lifetime, number and rate of Po215 events
                    # getLife returns N,dN, tau,dtau
                    # tau,dtau = fitted lifetime, error for 2 parameter (N,tau) fit
                    # N,dN = fitted # events, error for 1 parameter fit, tau fixed at known Po215 lifetime
                    lf = self.getLife(fn=fn,withPromptCut=withPromptCut)
                    norm = (runtime,runtime,1.,1.)
                    for i,life in enumerate(LIFE):
                        life.append( lf[i]/norm[i] )

                    # get per run results (start timestamp, sources used, sample name, runtim
                    # add to results dict with key=rn=string run number =runXXXXX
                    results[rn] = [timestamp,sources,sample,runtime]
                    # normalize Po215 peak area and uncertainty for each lifetime cut
                    norm = self.normResults(results)

                    # expected Ac227 rate at start of run given timestamp(converted from ms) and sample name
                    expectedRate,der = self.cAC.expectAc227Rate(inputDay=timestamp/1000,sampleName=sample)
                    ER.append(expectedRate)
                    dER.append(der)
                    words = rn + ' rates(Hz) '
                    stuff = 'livetime={0:.1f} Evts(unc)'.format(runtime)
                    X.append(float(rn.replace('run','')))
                    T.append(timestamp/1000) # convert from ms
                    for x in sorted(norm.keys()):
                        if x is not rn:
                            #print rn,'x,norm[x]',x,norm[x]
                            words += ' {0:.1f}({1:.1f})'.format(norm[x][3],norm[x][4])
                            stuff += ' {0:.1f}({1:.1f})'.format(results[x][3],results[x][4])
                            if x not in Y:
                                Y[x],dY[x],PoPeak[x],dPoPeak[x],PoS[x] = [],[],[],[],[]
                            Y[x].append(norm[x][3])
                            dY[x].append(norm[x][4])
                            PoPeak[x].append(results[x][0]) #mean
                            dPoPeak[x].append(results[x][1])#emean
                            PoS[x].append(results[x][2])    #sigma
                    print words,' ',stuff
                    if Quick and nQuick<=0:
                        print '.................QUICK job'
                        break


        '''
        End of processing
        '''

        #print 'X',X
        #print 'Y',Y
        #print 'dY',dY
        X = numpy.array(X)
        dX = numpy.array([0. for x in X])
        T = numpy.array(T)

        tdict = {}
        # fill dict for ttree
        for i,V in enumerate(Prompt):
            n = PromptName[i]
            dn = 'd'+n
            tdict[n] = V
            tdict[dn]= dPrompt[i]
        for n,V in zip(LIFEName,LIFE):
            tdict[n] = V
        tdict['ER'] = ER
        tdict['dER']= dER
        tdict['ts'] = T
        tdict['run']= Run
        tdict['sample'] = sampleName


        tmg = self.gU.makeTMultiGraph('Po215_Rate_Hz_vs_run')
        tmgT= self.gU.makeTMultiGraph('Po215_Rate_Hz_vs_time')
        tmgEC = self.gU.makeTMultiGraph('EffCorr_Po215_Rate_Hz_vs_run')
        tmgTEC= self.gU.makeTMultiGraph('EffCorr_Po215_Rate_Hz_vs_time')
        tmgP= self.gU.makeTMultiGraph('Po215_peak_vs_run')
        tmgS= self.gU.makeTMultiGraph('Po215_sigma_vs_run')
        tmgPSD = self.gU.makeTMultiGraph('PSD_effy_vs_run')
        tmgN= self.gU.makeTMultiGraph('HzPo215_from_lifetimeFit_vs_run')
        tmgL= self.gU.makeTMultiGraph('Fitted_lifetime_vs_run')
        tmgAll = self.gU.makeTMultiGraph('Ac227_rate_Hz_vs_run')
        tmgAllT = self.gU.makeTMultiGraph('Ac227_rate_Hz_vs_time')
        tmgEff = self.gU.makeTMultiGraph('Efficiency_vs_run')
        graphs = [tmg,tmgT,tmgEC,tmgTEC,tmgP,tmgS,tmgPSD,tmgN,tmgL,tmgAll,tmgAllT, tmgEff]
        hists = []

        title = 'Expected sample rate vs run'
        name  = title.replace(' ','_')
        y,dy = numpy.array(ER),numpy.array(dER)
        g = self.gU.makeTGraph(X,y,title,name,ey=dy,ex=dX)
        self.gU.color(g,3,3,setMarkerColor=True)
        tmg.Add(g)
        tmgEC.Add(g)
        tmgN.Add(g)
        tmgAll.Add(g)
        graphs.append(g)

        title = 'Expected sample rate vs time'
        name  = title.replace(' ','_')
        g = self.gU.makeTGraph(T,y,title,name,ey=dy,ex=dX)
        self.gU.color(g,3,3,setMarkerColor=True)
        tmgT.Add(g)
        tmgTEC.Add(g)
        tmgAllT.Add(g)
        graphs.append(g)
        

        for l in sorted(Y.keys()):
            y,dy = numpy.array(Y[l]),numpy.array(dY[l])
            title = 'Cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            delta = (float(l)-2.5)/5.
            g = self.gU.makeTGraph(X+delta,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmg.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir)
            graphs.append(g)

            name = name+'_vs_timestamp'
            g = self.gU.makeTGraph(T,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgT.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True)
            graphs.append(g)

            y,dy = numpy.array(PoPeak[l]),numpy.array(dPoPeak[l])
            title='Po215 peak cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            g = self.gU.makeTGraph(X+delta,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgP.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            graphs.append(g)
            
            y    = numpy.array(PoS[l])
            title='Po215 sigma cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            g = self.gU.makeTGraph(X+delta,y,title,name)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgS.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            graphs.append(g)

        for i in range(3):
            y,dy = numpy.array( PSD[2*i] ), numpy.array( PSD[2*i+1] )
            title = 'PSD effy'+str(i)+' vs run'
            name  = title.replace(' ','_')
            g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,i,i,setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            tmgPSD.Add(g)
            graphs.append(g)

        for i in range(2):
            y,dy = numpy.array( LIFE[2*i] ), numpy.array( LIFE[2*i+1] )
            title= LIFEName[2*i] + ' vs run'
            name = title.replace(' ','_')
            g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,i,i,setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            graphs.append(g)
            if i==0: tmgN.Add(g)
            if i==1:
                tmgL.Add(g)
                title = 'Fitted lifetime (s)'
                name = 'Fitted_lifetime_s'
                xma,xmi,dq = max(y),min(y),(max(y)-min(y))/10.
                hists.append( self.gU.makeTH1D(y,title,name,xmi=xmi-dq,xma=xma+dq,nx=50) )

                title = 'Fitted lifetime residual'
                name = title.replace(' ','_')
                q = (y-self.cAC.Po215lifetime)/dy
                xma,xmi,dq = max(q),min(q),(max(q)-min(q))/10.
                hists.append( self.gU.makeTH1D(q,title,name,xmi=xmi-dq,xma=xma+dq,nx=50) )
                    
            

        for i,p in enumerate(Prompt):
            y,dy = numpy.array(p),numpy.array(dPrompt[i])
            name = PromptName[i]+'_vs_run'
            title = name.replace('_',' ')
            g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,i,i,setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            graphs.append(g)

        y,dy = numpy.array(GoodCut[0]),numpy.array(GoodCut[1])
        name = 'GoodCh0Effy_vs_run'
        title= name.replace('_',' ')
        g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
        self.gU.color(g,1,1,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        graphs.append(g)
        tmgEff.Add(g)

        # apply effy corrections to Po215 rate and prompt spectrum (entire alpha spectrum)
        Ycorr,dYcorr = [],[] # will become corrected Po215 rate
        Pcorr,dPcorr = [],[] # will become prompt corrected rate
        # effPeak is efficiency of selecting Po215 in bkgd-subtracted delayed spectrum
        effPeak = quad(scipy.stats.norm.pdf,-self.cAC.nsigmaCutPo215Peak,self.cAC.nsigmaCutPo215Peak)[0]
        # effPrompt is putative effy of cut on prompt energy distribution to isolate Po215 in delayed dist.
        effPrompt = 1.
        if withPromptCut : effPrompt = self.cAC.promptChargeCutEffy
        effLo     = self.cAC.lowChargeCutEffy
        nAlpha    = self.cAC.totalAlphas
        effPSD = numpy.array(PSD[4])
        errPSD = numpy.array(PSD[5])
        effGood= numpy.array(GoodCut[0])
        errGood= numpy.array(GoodCut[1])
        #print 'effPSD.shape,errPSD.shape,effGood.shape,errGood.shape',effPSD.shape,errPSD.shape,effGood.shape,errGood.shape
        tdict['effGood'] = effGood
        tdict['errGood'] = errGood
        tdict['effPSD']  = effPSD
        tdict['errPSD']  = errPSD
        tdict['effPrompt']=numpy.array([effPrompt for x in effPSD]) # same value correct number of times
        tdict['effLo']    =numpy.array([effLo     for x in effPSD]) # same value correct number of times

        name = 'PSDEffy_vs_run'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,effPSD,title,name,ex=dX,ey=errPSD)
        self.gU.color(g,2,2,setMarkerColor=True)
        graphs.append(g)
        tmgEff.Add(g)

        name = 'LowEnCutEffy_vs_run'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,tdict['effLo'],title,name,ex=dX,ey=dX) # no uncertainty on this effy
        self.gU.color(g,3,3,setMarkerColor=True)
        graphs.append(g)
        tmgEff.Add(g)
        

        # apply efficiency correction to Po215 rate from lifetime fit
        # effy corr is square of product of effy of goodCh0, PSD and lowChargeCut
        #  because same cuts are applied to both prompt and delayed candidates.
        #  Assume no uncertainty on lowChargeCut
        for i,Nlf in enumerate(LIFE[0]):
            dNlf = LIFE[1][i] # fitted error on N of events from lifetime fit
            fG,dG = effGood[i],errGood[i]
            fP,dP = effPSD[i],errPSD[i]
            ec = math.pow(fG*fP*effLo,2)
            R  = Nlf/ec
            dR= math.sqrt( dG*dG/fG/fG + dP*dP/fP/fP ) * 2. * R
            dR= math.sqrt( dR*dR/R/R + dNlf*dNlf/Nlf/Nlf ) * R
            ecLIFE[0].append( R )
            ecLIFE[1].append( dR)
        for word,value in zip(ecLIFEName,ecLIFE):
            tdict[word] = value
            
        
        
        
        for l in sorted(Y.keys()):
            y = numpy.array(Y[l])
            #print 'l',l,'y,effPSD,effGood,effPeak,effPrompt',y,effPSD,effGood,effPeak,effPrompt
            ycorr = y/effPSD/effGood/effPeak/effPrompt
            dycorr= ycorr*numpy.sqrt( errPSD*errPSD/effPSD/effPSD + errGood*errGood/effGood/effGood)
            tdict['corrPo215Hz'+str(l)]  =  ycorr
            tdict['dcorrPo215Hz'+str(l)] = dycorr

            
            name = 'EffCorr_Po215Rate_Hz_vs_run_with_'+str(l)+'lifetime_cut'
            title = name.replace('_',' ')
            g = self.gU.makeTGraph(X,ycorr,title,name,ex=dX,ey=dycorr)
            self.gU.color(g,int(l),int(l),setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir)
            tmgEC.Add(g)

            graphs.append(g)

            name = 'EffCorr_Po215Rate_Hz_vs_time_with_'+str(l)+'lifetime_cut'
            title = name.replace('_',' ')
            g = self.gU.makeTGraph(T,ycorr,title,name,ex=dX,ey=dycorr)
            self.gU.color(g,int(l),int(l),setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir)
            tmgTEC.Add(g)

            graphs.append(g)
            

        # effy-corrected and Nalpha-corrected prompt rate
        p = numpy.array(Prompt[3])
        dp= numpy.array(dPrompt[3])
        pcorr = p/effPSD/effLo/effGood/nAlpha
        #print 'pcorr[0],p[0],effPSD[0],effLo,effGood[0],nAlpha',pcorr[0],p[0],effPSD[0],effLo,effGood[0],nAlpha
        dpcorr=  pcorr*numpy.sqrt( errPSD*errPSD/effPSD/effPSD + errGood*errGood/effGood/effGood)
        tdict['corrPromptHz'] =  pcorr
        tdict['dcorrPromptHz']= dpcorr
        name = 'EffCorr_PromptRate_Hz_vs_run'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,pcorr,title,name,ex=dX,ey=dpcorr)
        self.gU.color(g,4,4,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        tmgAll.Add(g)
        graphs.append(g)

        name = 'EffCorr_PromptRate_Hz_vs_time'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(T,pcorr,title,name,ex=dX,ey=dpcorr)
        self.gU.color(g,4,4,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        tmgAllT.Add(g)
        graphs.append(g)
        
        name = 'EffCorr_PromptRate_Hz'
        title = name.replace('_',' ')
        xma,xmi,dq = max(pcorr),min(pcorr),(max(pcorr)-min(pcorr))/10.
        h = self.gU.makeTH1D(pcorr,title,name,xmi=xmi-dq,xma=xma+dq,nx=50)
        hists.append(h)

        er,der = numpy.array(ER),numpy.array(dER)
        name = 'EffCorr_PromptRate_Residual'
        title = name.replace('_',' ')
        q = (pcorr-er)/numpy.sqrt(der*der + dpcorr*dpcorr)
        xma,xmi,dq = max(q),min(q),(max(q)-min(q))/10.
        h = self.gU.makeTH1D(q,title,name,xmi=xmi-dq,xma=xma+dq,nx=50)
        hists.append(h)
        
        

        # effy-corrected Po215 rate from lifetime fit
        y,dy = numpy.array(ecLIFE[0]),numpy.array(ecLIFE[1])
        name = 'EffCorr_Po215_rate_from_lifetime_fit_vs_run'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        tmgAll.Add(g)
        graphs.append(g)

        name = 'EffCorr_Po215_rate_from_lifetime_fit_vs_time'
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(T,y,title,name,ex=dX,ey=dy)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        tmgAllT.Add(g)
        graphs.append(g)

        # difference between effy-correcgted Po215 rate from lifetime fit and expected rate
        name = 'Po215_rate_from_lifetime_fit_minus_expected_rate_vs_run'
        title = name.replace('_',' ')
        v1,v2 = numpy.array(ecLIFE[1]),numpy.array(dER)
        Q,dQ = numpy.subtract(ecLIFE[0],ER),numpy.sqrt( v1*v1 + v2*v2 )
        g = self.gU.makeTGraph(X,Q,title,name,ex=dX,ey=dQ)
        self.gU.color(g,1,1,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        graphs.append(g)

        name = name.replace('run','time')
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(T,Q,title,name,ex=dX,ey=dQ)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        graphs.append(g)

        hQ = self.gU.makeTH1D(Q/dQ,'(Po215rate(lifetime fit) - expected rate)/uncertainty','Po215_er_by_unc',nx=50,xmi=-5.,xma=5.)
        hists.append( hQ )
            
        # difference between effy-correcgted Ac227 rate from prompt rate  and expected rate
        name = 'Ac227_rate_from_prompt_rate_minus_expected_rate_vs_run'
        title = name.replace('_',' ')
        Q,dQ = numpy.subtract(pcorr,ER),numpy.sqrt( dpcorr*dpcorr + v2*v2 )
        g = self.gU.makeTGraph(X,Q,title,name,ex=dX,ey=dQ)
        self.gU.color(g,3,3,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        graphs.append(g)

        name = name.replace('run','time')
        title = name.replace('_',' ')
        g = self.gU.makeTGraph(T,Q,title,name,ex=dX,ey=dQ)
        self.gU.color(g,4,4,setMarkerColor=True)
        self.gU.drawGraph(g,figDir=self.figdir)
        graphs.append(g)

        hQ = self.gU.makeTH1D(Q/dQ,'(Ac227rate(prompt events) - expected rate)/uncertainty','Ac227_er_by_unc',nx=50,xmi=-5.,xma=5.)
        hists.append( hQ )
            

        self.gU.drawMultiGraph(tmg,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgT,figdir=self.figdir,abscissaIsTime=True,xAxisLabel='Time',yAxisLabel='Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgEC,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Effy corrected Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgTEC,figdir=self.figdir,abscissaIsTime=True,xAxisLabel='Time',yAxisLabel='Effy corrected Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgP,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Fitted mean of 215Po peak (nC)')
        self.gU.drawMultiGraph(tmgS,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Fitted sigma of 215 Po peak (nC)')
        self.gU.drawMultiGraph(tmgPSD,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='PSD effy estimate')
        self.gU.drawMultiGraph(tmgN,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Po215 rate (Hz) from lifetime fit')
        self.gU.drawMultiGraph(tmgL,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Fitted lifetime (s)')

        self.gU.drawMultiGraph(tmgAll,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Ac227 rate (Hz)')
        self.gU.drawMultiGraph(tmgAllT,figdir=self.figdir,abscissaIsTime=True,xAxisLabel='Run number',yAxisLabel='Ac227 rate (Hz)')

        self.gU.drawMultiGraph(tmgEff,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Efficiency')
        
        rl = sorted(listOfGoodRunsAS)
        print 'calibWaveDump.main Opened',len(listOfRunsAS),'root files with',len(listOfGoodRunsAS),'good runs. First,last good runs',rl[0],rl[-1]
        
        fn = self.outFileDir + 'summary_' + rl[0] + '_' + rl[-1] + '.root'
        rfn = ROOT.TFile.Open(fn,'RECREATE')
        for g in graphs: rfn.WriteTObject(g)
        for h in hists: rfn.WriteTObject(h)
        rfn.Close()
        print 'calibWaveDump.main Wrote',len(graphs),'graphs,',len(hists),'hists to',fn



        self.mTT.makeTTree(tdict,fn=self.outTree,treename='cWD',debug=False)
          
        return
    def normResults(self,results):
        '''
        return results as rates where applicable
        '''
        norm = {}
        rn = None
        for x in results:
            if type(x) is str: rn = x
        timestamp,sources,sample,runtime = results[rn]
        norm[rn] = results[rn]
        for x in results:
            if x is not rn:
                #print rn,'x',x
                mean,emean,sgm, S,E = results[x]
                norm[x] = [mean,emean,sgm, S/runtime,E/runtime]
        return norm
if __name__ == '__main__' :
    '''
    arguments
    1 = if True or true or 1, then apply prompt cut
    2 = if false or no, then no logging
    '''
    Log = True
    withPromptCut = False
    if len(sys.argv)>1 : withPromptCut = sys.argv[1].lower()=='true' or sys.argv[1]=='1'
    if len(sys.argv)>2 : Log = sys.argv[2].lower()=='false' or sys.argv[1].lower()=='no'
    
    cWD = calibWaveDump(Log=Log) 
    cWD.main(withPromptCut=withPromptCut)


    
    if 0:
        for i in range(31,90):
            fn='Output/run'+str(i).zfill(5)+'.root'
            cWD.loop(fn=fn)
