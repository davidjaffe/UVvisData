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
import matplotlib
import matplotlib.pyplot as plt
import datetime,os

import ROOT
import graphUtils
import gfit
import readLogFile
import cutsAndConstants

class calibWaveDump():
    def __init__(self):


        self.gU = graphUtils.graphUtils()
        self.gfit = gfit.gfit()
        self.rLF = readLogFile.readLogFile()
        self.cAC = cutsAndConstants.cutsAndConstants()

        self.psdCut = self.cAC.psdCut
        self.lowChargeCut = self.cAC.lowChargeCut

        

        self.rootFileDir = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        self.logFileDir  = self.rootFileDir.replace('rootfiles','logfiles')
        self.outFileDir  = 'Output/'
        self.figdir      = 'Figures/calibWaveDump/'
        return
    def get_filepaths(self,directory):
        '''
        20160906 taken from http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
        This function will generate the file names in a directory 
        tree by walking the tree either top-down or bottom-up. For each 
        directory in the tree rooted at directory top (including top itself), 
        it yields a 3-tuple (dirpath, dirnames, filenames).
        '''
        file_paths = []  # List which will store all of the full filepaths.

        # Walk the tree.
        for root, directories, files in os.walk(directory):
            for filename in files:
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)  # Add it to the list.

        return file_paths  # Self-explanatory.
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
    def loop(self,fn='Output/run00073.root'):
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
        ROOT.gStyle.SetStatX(.5)
        hlist = []
        for l in range(1,6):
            clife = str(l)
            life = float(l)
            sname = 'dc'+clife
            bname = 'oc'+clife
            if hf.FindKey(sname):
                hs = hf.Get(sname)
                hb = hf.Get(bname)
                dname = 'd-o'+clife
                hd,mean,emean,sgm,S,E = self.subAndFit(hs,hb,dname)
                hlist.append(hd)
                f = 1./(1.-math.exp(-life))
                Scorr = S*f
                Ecorr = E*f
                #print fn,'l,mean,emean,sgm,Scorr,Ecorr {0} {1:.5f} {2:.5f} {3:.5f} {4:.1f} {5:.1f}'.format(l,mean,emean,sgm,Scorr,Ecorr)
                results[l] = [mean,emean,sgm,Scorr,Ecorr]


        dn = (fn.split('/')[1]).split('.')[0]
        self.gU.drawMultiHists(hlist,fname=dn,figdir=self.figdir,dopt='pe',statOpt=0,biggerLabels=False)
        hf.Close()
        return results
    def fitPSD(self,fn='Output/run00073.root'):
        '''
        fit PSD dist with two gaussians and estimate efficiency of psd cut with two methods.
        First method: Fractional area above PSD cut based on fit to n-recoil peak
        Second method: Fraction histogram area above PSD cut by subtracting fitted e-recoil peak

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

        if debug: print 'calibWaveDump.fitPSD effy,erry',effy,erry,'effy2,erry2',effy2,erry2
        hf.Close()      
        return effy,erry, effy2,erry2
    def main(self):
        '''
        main routine
        generate list of input files, find corresponding logfiles to get run info,
        get stats for a run,
        plot stats vs run # and vs time
        '''
        listOfHistFiles = self.get_filepaths(self.outFileDir)
        listOfLogFiles  = self.get_filepaths(self.logFileDir)

        Threshold = 10. # minimum number of coincidences for a GOOD run
        Ac227only = True # only process runs with Ac227 (no Cs137, for example)
        if Ac227only:
            goodSources = ['Ac-227']
            badSources  = ['Cs-137']
        else:
            goodSources = ['Ac-227','Cs-137']
            badSources  = []
        print 'calibWaveDump.main Threshold',Threshold,'goodSources',goodSources,'badSources',badSources
            
        
        X,Y,dY,PoPeak,dPoPeak,PoS = [],{},{},{},{},{}
        PSD = PSD1,dPSD1,PSD2,ePSD2     = [],[],[],[]
        Prompt = [],[],[],[]
        dPrompt= [],[],[],[]
        PromptName = ['Prompt_Mean','Prompt_StdDev','Prompt_Median','Prompt_Rate_Hz']
        T = []
        
        for fn in listOfHistFiles:
            rn = os.path.basename(fn).split('.')[0] # run00xxx
            lf = None
            for lfn in listOfLogFiles:
                if rn in lfn:
                    lf = lfn
                    break
            if lf is not None:
                timestamp,sources,sample,runtime = self.rLF.readFile(lf)
                GOOD = False
                for src in goodSources:
                    if src in sources: GOOD = True
                for src in badSources:
                    if src in sources: GOOD = False
                if not GOOD: print 'calibWaveDump.main',rn,'no good sources'
                if GOOD:
                    results = self.loop(fn=fn)
                    for l in results:
                        if results[l][3]<Threshold : GOOD = False
                    if not GOOD: print 'calibWaveDump.main',rn,'too few coincidences'
                if GOOD:
                    pStats  = self.getPstats(fn=fn) #mean,stddev,median,N
                    norm = (1.,1.,1.,runtime)
                    for i,p in enumerate(Prompt):
                        p.append( float(pStats[i])/norm[i] )
                        if i==0:
                            dPrompt[i].append( pStats[1]/math.sqrt(float(pStats[3]) ) )
                        if i==1 or i==2 : dPrompt[i].append( 0. )
                        if i==3:
                            dPrompt[i].append( math.sqrt(float(pStats[3]))/norm[i] )
                        
                    psdEffy = self.fitPSD(fn=fn)
                    for i,psd in enumerate(PSD):
                        psd.append( psdEffy[i] )
                        
                        
                    results[rn] = [timestamp,sources,sample,runtime]
                    norm = self.normResults(results)
                    #print rn,'results',results
                    #print rn,'norm',norm
                    words = rn + ' '
                    stuff = '{0:.1f}'.format(runtime)
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
        #print 'X',X
        #print 'Y',Y
        #print 'dY',dY
        X = numpy.array(X)
        dX = numpy.array([0. for x in X])
        T = numpy.array(T)
        tmg = self.gU.makeTMultiGraph('Rate_vs_run')
        tmgT= self.gU.makeTMultiGraph('Rate_vs_time')
        tmgP= self.gU.makeTMultiGraph('Po215_peak_vs_run')
        tmgS= self.gU.makeTMultiGraph('Po215_sigma_vs_run')
        tmgPSD = self.gU.makeTMultiGraph('PSD_effy_vs_run')

        for l in sorted(Y.keys()):
            y,dy = numpy.array(Y[l]),numpy.array(dY[l])
            title = 'Cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            delta = (float(l)-2.5)/5.
            g = self.gU.makeTGraph(X+delta,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmg.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir)

            name = name+'_vs_timestamp'
            g = self.gU.makeTGraph(T,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgT.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=True)

            y,dy = numpy.array(PoPeak[l]),numpy.array(dPoPeak[l])
            title='Po215 peak cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            g = self.gU.makeTGraph(X+delta,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgP.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            
            y    = numpy.array(PoS[l])
            title='Po215 sigma cut at '+str(l)+' lifetimes'
            name = title.replace(' ','_')
            g = self.gU.makeTGraph(X+delta,y,title,name)
            self.gU.color(g,l,l,setMarkerColor=True)
            tmgS.Add(g)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)

        for i in range(2):
            y,dy = numpy.array( PSD[2*i] ), numpy.array( PSD[2*i+1] )
            title = 'PSD effy'+str(i)+' vs run'
            name  = title.replace(' ','_')
            g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,i,i,setMarkerColor=True)
            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)
            tmgPSD.Add(g)

        print 'meanPrompt',
        for v,dv in zip(Prompt[0],dPrompt[0]): print '{0:.5f}({1:.5f})'.format(v,dv),
        print ''
        print 'PromptRateHz',
        for v,dv in zip(Prompt[3],dPrompt[3]): print '{0:.1f}({1:.1f})'.format(v,dv),
        print ''
        for i,p in enumerate(Prompt):
            y,dy = numpy.array(p),numpy.array(dPrompt[i])
            name = PromptName[i]+'_vs_run'
            title = name.replace('_',' ')
            g = self.gU.makeTGraph(X,y,title,name,ex=dX,ey=dy)
            self.gU.color(g,i,i,setMarkerColor=True)

            self.gU.drawGraph(g,figDir=self.figdir,abscissaIsTime=False)

        self.gU.drawMultiGraph(tmg,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgT,figdir=self.figdir,abscissaIsTime=True,xAxisLabel='Time',yAxisLabel='Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgP,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Fitted mean of 215Po peak (nC)')
        self.gU.drawMultiGraph(tmgS,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Fitted sigma of 215 Po peak (nC)')
        self.gU.drawMultiGraph(tmgPSD,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='PSD effy estimate')

        
            
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
    cWD = calibWaveDump()
    #cWD.fitPSD()
    #sys.exit('no more')
    cWD.main()
    if 0:
        for i in range(31,90):
            fn='Output/run'+str(i).zfill(5)+'.root'
            cWD.loop(fn=fn)
