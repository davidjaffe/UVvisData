#!/usr/bin/env python
'''
investigate chi2 behavior for pmt calibration algorithm
20191026
'''
import math
import sys,os

import datetime
import numpy
from scipy.stats import poisson,norm,relfreq
#import copy

import re
import glob # used in __init__

import matplotlib.pyplot as plt

class pmtcal():
    def __init__(self,makeFigures=False):

        self.debug = 1

        self.currentConfig = None
        
        self.figdir = 'FIGURES/'
        self.makeFigures = makeFigures
        if self.makeFigures:
            print 'pmtcal.__init__ will create new sub-directories in',self.figdir
            for dirpath, dirnames, filenames in os.walk(self.figdir):
                if dirpath!=self.figdir:
                    try:
                        os.rmdir(dirpath)
                        print 'pmtcal.__init__ Deleted',dirpath
                    except OSError as ex:
                        print ex
            

        # enables use of latex 
        from matplotlib import rc
        rc('text', usetex=True)
            
        print 'pmtcal.__init__ Done'
        return
    def takeData(self,mu,nData,smear=0.5,tailFrac=0.01,tailMu=100.,debug=0):
        '''
        generate nData events drawn from two poisson distributions with 
        means mu and tailMu with frequency 1-tailFrac, tailFrac
        with gaussian smearing with sigma=smear
        non negative results required
        '''
        if tailFrac<0. or tailFrac>1. :
            sys.exit('pmtcal.takeData ERROR tailFrac='+str(tailFrac)+' must be in [0,1]')
        n1 = int(float(nData)*(1.-tailFrac))
        n2 = nData - n1
        if debug>0: print 'pmtcal.takeData tailFrac,n1,n2',tailFrac,n1,n2
        pile1 = poisson.rvs(mu,size=n1) + smear*norm.rvs(size=n1)
        pile2 = poisson.rvs(tailMu,size=n2) + smear*norm.rvs(size=n2)
        pile = numpy.append(pile1,pile2)
        data = numpy.maximum( numpy.zeros(len(pile)), pile )
        return data
    def generateMC(self,mu,nMC,smear=0.5):
        '''
        generate nMC MC events drawn from poisson distribution, mean mu
        '''
        pile = poisson.rvs(mu,size=nMC) + smear*norm.rvs(size=nMC)
        MC = numpy.maximum( numpy.zeros(len(pile)), pile )
        return MC
    def makeConfigs(self):
        '''
        return dict that defines various configurations
        dict[key] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        nData = # of data events
        muData= mean of poisson for data
        tailFrac = fraction of data that has poisson with mean=tailMu
        analogous definitions for MC
        text = text defining the configuration
        '''
        ic = 0
        config = {}
        nData = 10000
        nMC   = 10*nData
        muData = 8.3
        muMC   = 8.3
        tailMu = 40.
        tailFrac = 0.
        text = '\n{6}: nD={0}, nM={1}, $\mu_D=${2:0.2f}, $\mu_M=${3:0.2f}, tailF={4:0.2f}, $\mu_t={5:0.2f}'.format(nData,nMC,muData,muMC,tailFrac,tailMu,ic)
        config[ic] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        print 'pmtcal.makeConfigs STOP WITH',ic+1,'CONFIGURATIONS'
        return config

    
        for tailFrac in [0., 0.01, 0.05]:
            for muMC in [6.0, 10.0]:
                ic += 1
                text = '\n{6}: nD={0}, nM={1}, $\mu_D=${2:0.2f}, $\mu_M=${3:0.2f}, tailF={4:0.2f}, $\mu_t$={5:0.2f}'.format(nData,nMC,muData,muMC,tailFrac,tailMu,ic)
                config[ic] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        return config
            
    def main(self):
        '''
        main routine
        generate data and mc
        bin data and mc
        plot data and mc
        '''

        debug = 0
        
        nData = 10000
        nMC = 10*nData
        muData = 8.3
        tailMu = 40.
        muMC   = 8.3

        configs = self.makeConfigs()

        for ic in sorted(configs):
            self.currentConfig = ic
            nData,nMC, muData,muMC, tailMu,tailFrac, text = configs[ic]
            print 'pmtcal.main Configuration#',ic,text

            words = text
            data = self.takeData(muData,nData,tailFrac=tailFrac,tailMu=tailMu)
            MC   = self.generateMC(muMC,nMC)

            origData,origMC,nmax,limits = self.binData(data,MC,nThres=-99)
            self.drawHist(nmax,limits,origData,'NPE','Events','Original Data'+words)
            self.drawHist(nmax,limits,origMC,'NPE','Events','Original MC'+words)


            
            # bin the data, require >nThres events per bin. Overflows in bin failing overflow criterion
            nThres = 10
            histData, histMC, nbins, limits = self.binData(data,MC,nThres=nThres)

            self.drawHist(nbins,limits,histData,'NPE','Events','Data'+words)
            self.drawHist(nbins,limits,histMC  ,'NPE','Events','MC'+words)
            sMC = sum(histData)/sum(histMC)*histMC
            self.draw2Hist(nbins,limits,histData,sMC,'NPE','Events','Data and scaled MC'+words,label1='Data',label2='Scaled MC')

            
            bChi,Chi2,MChists = {},{},{}
            nf = 11
            fmi,fma = 0.5,1.5
            df = (fma-fmi)/float(nf-1)
            while df>0.001:
                fcalValues = numpy.linspace(fmi,fma,nf)
                minChi,fatmin = 1.e20, fmi
                for fcal in fcalValues:
                    bChi[fcal],MChists[fcal] = self.binnedChi(fcal, histData, MC, nbins, limits)
                    Chi2[fcal] = sum(bChi[fcal]*bChi[fcal])
                    if debug>0: print 'pmtcal.main fcal,Chi2',fcal,Chi2[fcal]
                    if Chi2[fcal]<minChi:
                        minChi = Chi2[fcal]
                        fatmin = fcal
                fmi,fma = fatmin-df,fatmin+df
                nf = (nf-1)*2 + 1
                df = (fma-fmi)/float(nf-1)           

            fcal = fatmin
            title = 'Best fit at calibration factor of {0:.3f}'.format(fcal)+words
            self.drawComp(nbins,limits,histData,MChists[fcal],bChi[fcal],'NPE','Events',title,label1='Data',label2='Scaled MC')
                
            x = sorted(Chi2)
            for fcal in [ x[3],x[-3] ]:
                title = 'Fit result for calibration factor of {0:.3f}'.format(fcal)+words
                self.drawComp(nbins,limits,histData,MChists[fcal],bChi[fcal],'NPE','Events',title,label1='Data',label2='Scaled MC')
            
            y = []
            for fcal in x: y.append( Chi2[fcal] )
            self.drawIt(x,y,'calibration factor','Chi2','Chi2 vs calib factor'+words)
            ix = y.index(min(y))
            xmi = x[ix]-0.1
            xma = x[ix]+0.1
            ymi = min(y)-1.
            yma = ymi + 11.
            self.drawIt(x,y,'calibration factor','Chi2','Chi2 vs calib factor Zoom on chi2min'+words,ylims=[ymi,yma],xlims=[xmi,xma])

        
        return
    def binnedChi(self,fcal, histData, MC, nbins, limits):
        '''
        return array of (data - MC)/sigma for each bin and MC hist
        given the calibration scale factor fcal, histData = histogram of data, MC = MC events, nbins, limits = #bins, limits of histogram 
        MC is scaled by the calibration factor fcal
        '''
        nData = sum(histData)
        nMC   = float(len(MC))
        rdmc  = nData/nMC
        sMC = fcal*MC
        histMC = self.histMaker(sMC, nbins, limits)
        sig = numpy.sqrt(histData + rdmc*rdmc*histMC)
        histMC *= rdmc
        histChi = (histData - histMC)/sig
        return histChi, histMC
    def binData(self,data,MC,nThres=10,debug=0):
        '''
        return histData,histMC,nbins,(lower limit, upper limit)
        given data,MC and threshold counts for data histogram

        First make a histogram of data with unit bins that contains all data.
        Then find the first bin with content<nThres = thresBin
        Then make a new histogram of data with contents of bins>thresBin in thresBin
        Also make a hist with same limits as data for MC

        if threshold is negative, just return hists with binning that contains all data
        '''
        maxdata = max(data)
        
        nmax = int(maxdata)+1
        limits = (0., float(nmax))
        if debug>1: print 'pmtcal.binData mindata,maxdata,nmax,limits,nThres',min(data),maxdata,nmax,limits,nThres
        
        if nThres<0 :   # JUST RETURN HISTS BINNED WITHOUT DATA OVERFLOWS
            histData = self.histMaker(data,nmax,limits,debug=debug,caller='pmtcal.binData original data')
            histMC   = self.histMaker(MC,nmax,limits,debug=debug,caller='pmtcal.binData original MC')
            return histData,histMC,nmax,limits
    
        
        hist, lowerlimit, binsize, extrapoints = relfreq(data,nmax,limits)

        if debug>1: print 'pmtcal.binData histogram #bins,limits,lowerlimit, binsize, extrapoints=',nmax-1,limits,lowerlimit,binsize,extrapoints
        if extrapoints>0:
            print 'pmtcal.binData ERROR extrapoints',extrapoints,'. It should be zero!'
        hist *= len(data)

        joverflow = numpy.where(hist-nThres<0.)[0][0]
        if debug>1: print 'pmtcal.binData joverflow',joverflow
        
        ioverflow = joverflow+1
        if debug>1:
            print 'pmtcal.binData threshold,overflow bin,bin contents',nThres,ioverflow,hist[ioverflow]
            print 'pmtcal.binData bin#,content',
            for i,x in enumerate(hist): print '{0:},{1:.0f} '.format(i,x),
            print ''

        upperlimit = lowerlimit + float(ioverflow)*binsize
        limits = (lowerlimit,upperlimit)
        nmax = ioverflow
        histData = self.histMaker(data,nmax,limits,debug=debug,caller='pmtcal.binData data')
        histMC   = self.histMaker(MC  ,nmax,limits,debug=debug,caller='pmtcal.binData MC')

        if debug>0 :
            for hist,words in zip([histData,histMC],['histData','histMC']):
                print 'pmtcal.binData',words,'bin#,content',
                for i,x in enumerate(hist): print '{0:},{1:.0f} '.format(i,x),
                print ''
        
        return histData, histMC, nmax, limits
    def histMaker(self,data, nmax, limits, debug=0, caller=None):
        '''
        return a histogram of input data with nmax bins between limits.
        Overflows are put in the uppermost bin
        '''
        upperlimit = limits[1]
        trunc = numpy.minimum( data, numpy.ones(len(data))*upperlimit)
        hist, lowerlimit, binsize, extrapoints = relfreq(trunc,nmax,limits)
        hist *= len(data)
        words = 'pmtcal.histMaker'
        if caller is not None: words = caller
        if debug>0:
            print words,'rebinned result #bins,limits,lowerlimit, binsize, extrapoints',nmax,limits,lowerlimit,binsize,extrapoints
        if extrapoints>0 :
            print words,'ERROR rebinning, extrapoints',extrapoints,'. It should be zero!'
        return hist
    def drawIt(self,x,y,xtitle,ytitle,title,figDir=None,ylog=False,xlims=None,ylims=None):
        '''
        draw graph defined by x,y

        '''
        plt.clf()
        plt.grid()
        plt.title(r''+title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        X = numpy.array(x)
        Y = numpy.array(y)
        plt.plot(X,Y,'o-')
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : plt.yscale('log')
        if xlims is not None: plt.xlim(xlims)
        if ylims is not None: plt.ylim(ylims)

        self.showOrPrint(plt,title)
        return    
    def drawHist(self,nbins,limits,histData,xtitle,ytitle,title,debug=0):
        '''
        draw histogram defined by nbins, limits, histData
        '''
        if debug>0: print 'pmtcal.drawHist nbins',nbins,'limits',limits,'histData',histData
        
        plt.clf()
        plt.grid()
        plt.title(r''+title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        dx = (limits[1]-limits[0])/float(nbins)
        x  = [limits[0]+i*dx for i in range(nbins)]
        plt.bar(x,histData,width=1,align='edge',color='white')

        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        #plt.show()
        self.showOrPrint(plt,title)
        return
    def draw2Hist(self,nbins,limits,hist1,hist2,xtitle,ytitle,title,label1='hist1',label2='hist2',debug=0):
        '''
        draw two histograms 
        error bars for first hist will be taken as sqrt of content
        no error bars draw for second hist
        '''
        if debug>0: print 'pmtcal.draw2Hist nbins',nbins,'limits',limits,'hist1',hist1,'hist2',hist2
        yma1,yma2 = max(hist1),max(hist2)
        yma = 1.025*max(yma1,yma2)
        if debug>0: print 'pmtcal.draw2Hist yma1,yma2,yma',yma1,yma2,yma
            
        plt.clf()
        plt.grid() # doesn't work?
        plt.title(r''+title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        dx = (limits[1]-limits[0])/float(nbins)
        x1  = [limits[0]+i*dx for i in range(nbins)]
        x2  = [limits[0]+dx/2.+i*dx for i in range(nbins)]
        yerr1 = numpy.sqrt(hist1)
        
        plt.bar(x1,hist1,width=1,align='edge',color='yellow',label=label1,yerr=yerr1,ecolor='black')
        plt.bar(x2,hist2,width=0.6,align='center',color='blue',label=label2,alpha=0.9)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
#        plt.ylim( [0., yma] )
        plt.legend(loc='best')
        #plt.show()
        self.showOrPrint(plt,title)
        return
    def drawComp(self,nbins,limits,hist1,hist2,y,xtitle,ytitle,title,label1='hist1',label2='hist2',
                     ylim=[-5.,5],debug=1):
        '''
        two unequal size panels
        top panel: plot two hists
        bottom panel : y = (data-MC)/sigma 
        '''
        plt.clf()
        plt.grid() # doesn't work?
        
        n = 6
        m = n-1
        ax1 = plt.subplot2grid((n, 1), (0, 0), rowspan=m)
        ax2 = plt.subplot2grid((n, 1), (m, 0))

        ax1.set_title(r''+title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        dx = (limits[1]-limits[0])/float(nbins)
        x1  = [limits[0]+i*dx for i in range(nbins)]
        x2  = [limits[0]+dx/2.+i*dx for i in range(nbins)]
        yerr1 = numpy.sqrt(hist1)
        
        ax1.bar(x1,hist1,width=1,align='edge',color='yellow',label=label1,yerr=yerr1,ecolor='black')
        ax1.bar(x2,hist2,width=0.6,align='center',color='blue',label=label2,alpha=0.9)
        ax2.set_xlabel(xtitle)
        ax2.set_ylabel('(Data-MC)/$\sigma$')
        ax1.set_ylabel(ytitle)
        ax1.set_xticklabels([]) # no tick labels on top plot

        ax1.legend(loc='best')

        ax2.grid()
        #ax2.plot(x2,y,'o')
        ax2.bar(x2,y,width=0.5,align='center',color='red')
        ax2.set_ylim(ylim)
        plt.subplots_adjust(hspace=0.) # zero height space between subplots
        
        #plt.show()
        self.showOrPrint(plt,title)
        
        return
    def renderTitle(self,title):
        i = title.find('\n')
        if i>-1: title = title[:i]
        figpdf = title.replace(' ','_')
        figpdf = 'FIG_' + figpdf + '.pdf'
        return figpdf
    def showOrPrint(self,plot,title):
        if self.makeFigures:
            if self.currentConfig is None:
                sys.exit('pmtcal.showOrPrint ERROR currentConfig is None')
            C = '{:>02d}'.format(self.currentConfig)
            figpdf = self.renderTitle(title)
            
            directory = self.figdir + C + '/'
            if not os.path.exists(directory):
                os.makedirs(directory)
                print 'pmtcal.showOrPrint Created',directory
                
            fn = directory + figpdf
            plot.savefig(fn)
        else:
            plot.show()
        return
if __name__ == '__main__' :

    makeFigures = False
    if len(sys.argv)>1:
        if sys.argv[1].lower()=='print': makeFigures=True
    pc = pmtcal(makeFigures=makeFigures)
    pc.main()
