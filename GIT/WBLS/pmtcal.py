#!/usr/bin/env python
'''
investigate chi2 behavior for pmt calibration algorithm
20191026
'''
import math
import sys,os,shutil

import datetime
import numpy
from scipy.stats import poisson,norm,relfreq
#import copy

import random

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
                    shutil.rmtree(dirpath)
                    print 'pmtcal.__init__ Deleted',dirpath

            

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
        onlyOne = False #True
        shortLists = False # True
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
        if onlyOne:
            print 'pmtcal.makeConfigs STOP WITH',ic+1,'CONFIGURATIONS'
            return config

        for muData in [8.3, 16.6]:
            tailFracs = [0., 0.01, 0.05]
            muMCs     = [float(i) for i in range(int(muData)-2,int(muData)+3) ]
            if shortLists:
                tailFracs = [tailFracs[0]]
                muMCs     = [muMCs[0]]
                print 'pmtcal.makeConfigs Only use single-entry lists for loops over tailFrac and muMC'


            for tailFrac in tailFracs:
                for muMC in muMCs:
                    ic += 1
                    text = '\n{6}: nD={0}, nM={1}, $\mu_D=${2:0.2f}, $\mu_M=${3:0.2f}, tailF={4:0.2f}, $\mu_t$={5:0.2f}'.format(nData,nMC,muData,muMC,tailFrac,tailMu,ic)
                    config[ic] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        return config
            
    def main(self):
        '''
        main routine
        for each experiment configuration:
           generate data and mc
           bin data and mc
           plot data and mc
           scan chi2 and get best fit value
        summarize results
        '''

        debug = 0


        configs = self.makeConfigs()
        results = {}

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

            # iterative scan to find minimum
            bChi,Chi2,MChists = {},{},{}
            nf = 11
            fmi,fma = 0.1,1.9
            df = (fma-fmi)/float(nf-1)
            dfmin = 0.0001 # 0.0001
            while df>dfmin:
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

            # show data, MC comparison for best fit and 2 random fits
            fcal = fatmin
            results[ic] = [fcal, Chi2[fcal], nbins]
            title = 'Best fit at calibration factor of {0:.3f}'.format(fcal)+words
            self.drawComp(nbins,limits,histData,MChists[fcal],bChi[fcal],'NPE','Events',title,label1='Data',label2='Scaled MC')
                
            x = sorted(Chi2)
            rx = random.sample(x,2) # select and draw a couple results at random
            for fcal in rx:
                title = 'Fit result for calibration factor of {0:.3f}'.format(fcal)+words
                self.drawComp(nbins,limits,histData,MChists[fcal],bChi[fcal],'NPE','Events',title,label1='Data',label2='Scaled MC')
                
            # show chi2 over full range and around minimum
            y = []
            for fcal in x: y.append( Chi2[fcal] )
            self.drawIt(x,y,'calibration factor','Chi2','Chi2 vs calib factor'+words)
            ix = y.index(min(y))
            ymi = min(y)-1.
            yma = ymi + 11.
            xmi,xma = None,None
            for j,u in enumerate( y[ix:]):
                if u>yma :
                    xma = x[j+ix]
                    xmi = x[max(ix-j,0)]
                    break
            for j,u in enumerate(y[:ix]):
                if u<yma :
                    xmi = x[j]
                    break
            if xmi is None: xmi = x[ix]-0.05
            if xma is None: xma = x[ix]+0.05
            self.drawIt(x,y,'calibration factor','Chi2','Chi2 vs calib factor Zoom on chi2min'+words,ylims=[ymi,yma],xlims=[xmi,xma])

        # show best fit vs expectation 
        self.drawBias(configs, results)

        # summary table and summary figures for latex
        self.makeSummary(configs, results)
        self.makeFigureFile(configs)
            
        return
    def drawBias(self,configs, results):
        '''
        for configuration ic
        configs[ic] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        results[ic] = [fcal, Chi2[fcal], nbins]

        compare best fit calibration factor with expectation as a function of tailFrac

        '''
        markers = ['o',   '<',    '^',  '>',    'v']
        colors  = ['blue','green','red','black','orange']

        
        x,y = {},{}
        xmi,xma,ymi,yma = 0.,2.,0.,2.
        for ic in sorted(configs):
            tailFrac=configs[ic][5]
            muData = configs[ic][2]
            muMC   = configs[ic][3]
            fexp   = muData/muMC
            fbest  = results[ic][0]
            xma = yma = max(xma,fexp,fbest)
            
            if tailFrac not in x:
                x[tailFrac],y[tailFrac] = [],[]
            x[tailFrac].append(fexp)
            y[tailFrac].append(fbest)

        plt.clf()
        plt.grid()
        title = 'Compare best fit vs expectation as a function of tail fraction'
        plt.title(title)
        plt.xlabel('Expected calibration factor')
        plt.ylabel('Best fit calibration factor')
        plt.xlim( [xmi,xma] )
        plt.ylim( [ymi,yma] )
        plt.plot( [xmi,xma], [ymi,yma],'--') # diagonal dashed line
        for i,tailFrac in enumerate(sorted(x)):
            X = numpy.array( x[tailFrac] )
            Y = numpy.array( y[tailFrac] )
            lg = 'Tail fraction = {0:.3f}'.format(tailFrac)
            plt.plot(X,Y,'s',marker=markers[i],color=colors[i],label=lg)

        plt.legend(loc='best')
        self.showOrPrint(plt,title)
        return
    def makeFigureFile(self,configs):
        '''
        use figurePacker to make an latex input file for figures for each configuration
        then create a single latex input for all the files
        '''
        confs = sorted(configs)
        c1 = '{0:>02d}'.format(confs[0])
        c2 = '{0:>02d}'.format(confs[-1])
        fnall = 'TEX/all_figint_'+c1+'_'+c2+'.tex'
        f = open(fnall,'w')
        for c in confs:
            figint = self.figurePacker(c)
            f.write('\\input{' + figint + '} \n')
        f.close()
        print 'pmtcal.makeFigureFile Wrote',len(confs),'input files into',fnall
        return
    def makeSummary(self,configs, results):
        '''
        for configuration ic
        configs[ic] = [nData,nMC, muData,muMC, tailMu,tailFrac, text]
        results[ic] = [fcal, Chi2[fcal], nbins]
        '''
        latexFile = self.figdir + 'results' + '_table.tex'
        latex = open(latexFile,'w')
        
        label = 'tab:' + 'results'
        caption = 'Different configurations and results. '
        caption += '$\mu_d = $ mean PE in data, $\mu_M =$ mean PE in MC, $\mu_t = $ mean PE in the tail, '
        caption += 'tailF = tail fraction, $f_{exp}=$ expected calibration factor, $f_{best} =$ best fit calibration factor, '
        caption += '$\chi^2_{min} = $ value of $\chi^2$ at minimum and nBin = number of bins in histogram.'

        ds = ' \\\\ \n' 
        
        latex.write('\\begin{table}[htp] \n')
        #latex.write( '\\setlength{\\tabcolsep}{2pt} \n')
        latex.write('\\begin{center} \n')
        latex.write('\\begin{tabular}{|r|rr| rr| rr| rr| rr|} \n' )
        fmt = '{0:>6} {1:>6} {2:>6} {3:>6} {4:>6} {5:>6} {6:>6} {7:>6} {8:>6} {9:>6} {10:>6}'
        fmt = fmt.replace('} {','}&{')
        latex.write( '\\hline \n' )
        latex.write(fmt.format('config','nData','nMC','$\mu_d$','$\mu_M$','$\mu_t$','tailF', '$f_{exp}$', '$f_{best}$', '$\chi^2_{min}$','nBin')+ds)
        latex.write( '\\hline \n' )
        fmt = '{9:>6} {0:>6} {1:>6} {2:>6.2f} {3:>6.2f} {4:>6.2f} {5:>6.2f} {6:>6.2f} {7:>6.2f} {8:>6.2f} {10:>6}'
        fmt = fmt.replace('} {','}&{')
        oldTF = None
        for ic in sorted(configs):
            nData,nMC, muData,muMC, tailMu,tailFrac, text = configs[ic]
            if oldTF is not None and oldTF!=tailFrac: latex.write( '\\hline \n' )
            fexp = muData/muMC
            fcal, chi2min, nbins = results[ic]
            latex.write(fmt.format(nData,nMC, muData,muMC, tailMu,tailFrac, fexp,fcal, chi2min, ic, nbins) + ds )
            oldTF = tailFrac

        latex.write( '\\hline \n' )
        latex.write('\\end{tabular} \n')
        latex.write('\\label{'+label+'} \n')
        latex.write('\\caption{'+caption+'} \n')
        latex.write('\\end{center} \n')
        latex.write('\\end{table} \n') 

            
        latex.close()
        print 'pmtcal.makeSummary wrote',latexFile
        return
        
    def binnedChi(self,fcal, histData, MC, nbins, limits):
        '''
        return array of (data - MC)/sigma for each bin and MC hist
        given the calibration scale factor fcal, histData = histogram of data, MC = MC events, nbins, limits = #bins, limits of histogram 
        MC is scaled by the calibration factor fcal

        handle the case where both hists have zero entries
        '''
        nData = sum(histData)
        nMC   = float(len(MC))
        rdmc  = nData/nMC
        sMC = fcal*MC
        histMC = self.histMaker(sMC, nbins, limits)
        sig = numpy.sqrt(histData + rdmc*rdmc*histMC)
        histMC *= rdmc
        if not numpy.all(sig):
            while not numpy.all(sig):
                i = numpy.argmin(sig)
                if histData[i]==0. and histMC[i]==0. :
                    sig[i] = 1.
                else:
                    print 'pmtcal.binnedChi i',i,'sig[i]',sig[i],\
                      'histData[i]',histData[i],'histMC[i]',histMC[i],\
                      ' THIS SHOULD NOT HAPPEN. ALL VALUES BESIDES INDEX SHOULD BE ZERO'
            for i,s in enumerate(sig):
                if s==0. and histData[i]==0. and histMC[i]==0.: sig[i] = 1.
            
        histChi = (histData - histMC)/sig
        return histChi, histMC
    def binData(self,data,MC,nThres=10,debug=0):
        '''
        return histData,histMC,nbins,(lower limit, upper limit)
        given data,MC and threshold counts for data histogram

        First make a histogram of data with unit bins that contains all data.
        Then find the first bin past the bin containing the maximum content with content<nThres = thresBin
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
    
        # make a hist of all data
        hist, lowerlimit, binsize, extrapoints = relfreq(data,nmax,limits)

        if debug>1: print 'pmtcal.binData histogram #bins,limits,lowerlimit, binsize, extrapoints=',\
            nmax-1,limits,lowerlimit,binsize,extrapoints
        if extrapoints>0:
            print 'pmtcal.binData ERROR extrapoints',extrapoints,'. It should be zero!'
        hist *= len(data)

        # add a 'bump' to hist to make sure that first value found below threshold is past the peak of the hist
        max_index = numpy.argmax(hist)
        max_value = hist[max_index]
        bump = numpy.array( [max_value*float(i<max_index) for i in range(len(hist)) ] )

        joverflow = numpy.where(hist+bump-nThres<0.)[0][0]
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
        figpdf = title.replace(' ','_').replace('.','_')
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
    def figurePacker(self,ic):
        '''
        create generic packages that to include all figures produced in a latex file

        return the filename of a file that contains something like the following:

        \begin{figure}[htbp]
        \begin{center}
        \includegraphics[width=0.45\textwidth]{../FIGURES/00/FIG_Original_Data.pdf}
        \includegraphics[width=0.45\textwidth]{../FIGURES/00/FIG_Original_MC.pdf}
        \includegraphics[width=0.45\textwidth]{../FIGURES/00/FIG_Data.pdf}
        \includegraphics[width=0.45\textwidth]{../FIGURES/00/FIG_MC.pdf}
        \caption{NPE histograms for data and MC for configuration 00. }
        \label{fig:dmc_00}
        \end{center}
        \end{figure}

        '''
        debug = 0
        fnprefix = '../'
        conf = '{0:>02d}'.format(ic)
        fmatch = 'FIGURES/'+conf+'/*.pdf'
        filelist = glob.glob(fmatch)
        if debug>0 :
            print 'pmtcal.figurePacker files matched to',fmatch,'by glob'
            for f in filelist: print f
        captions = ['Data compared to nominal MC, MC scaled by the best fit calibration factor,'+
                        ' scans of $\chi^2$ over a large range and about the minimum for configuration '+conf+
                        '. Data compared to MC scaled by two randomly chosen calibration factors.',
                        'NPE histograms for data and MC for configuration '+conf+'. Top are original hists.'+
                        ' Bottom are hists after truncation at bin containing less than 10 entries with that bin containing overflows.']
        labels = ['tab:best_'+conf,'tab:npe_'+conf]
        keys1 = ['_Data_and_scaled_MC','Best_fit','Chi2', 'Fit_result']
        keys2 = ['Original', '_Data.', '_MC.']
        list1,list2 = [],[]
        keys = [ keys1, keys2 ]
        lists= [ list1, list2 ]
        for k,KL in enumerate(zip(keys,lists)):
            key,l = KL
            for akey in key:
                if debug>1: print 'pmtcal.figurePacker akey',akey
                deletionList = []
                for i,f in enumerate(filelist):
                    if debug>1 : print 'pmtcal.figurePacker f',f
                    if akey in f:
                        l.append(fnprefix+f)
                        deletionList.append(i)
                deletionList.sort(reverse=True)
                if debug>0 : print 'pmtcal.figurePacker deletionList',deletionList
                for i in deletionList:
                    if debug>0 : print 'pmtcal.figurePacker --------> Delete',filelist[i]
                    del filelist[i]
        if debug>0:
            print 'pmtcal.figurePacker Here are the lists'
            for i,l in enumerate(lists):
                print 'pmtcal.figurePacker List#',i
                for x in l: print x
        if debug>0 or len(filelist)>0: 
            print 'pmtcal.figurePacker There are',len(filelist),'files left in filelist:',
            L2 = []
            for f in filelist:
                print f
                L2.append(fnprefix+f)
            lists.append(L2)
            labels.append('tab:extra_'+conf)
            captions.append('Comparison of best fit with expectation as a function of tail fractions')
                
                

        
        fname = 'TEX/figint_'+conf+'.tex'
        f = open(fname,'w')
        hdr = '\n \\begin{figure}[htbp] \\begin{center} \n'
        for i,figlist in enumerate(lists):
            fig = self.packFigures(figlist)
            trl = self.getTrailer(captions[i],labels[i])

            f.write(hdr)
            f.write(fig)
            f.write(trl)
        f.write('\\clearpage\n ')
        f.close()
        print 'pmtcal.figurePacker Wrote',fname

        fname = fname.replace('TEX/','')
        return fname
    def getTrailer(self,caption,label):
        trl = '\\caption{'+caption+'} \n'
        trl += '\\label{'+label+'} \n'
        trl += '\\end{center} \\end{figure} \n'
        return trl
    def packFigures(self,figlist):
        figpre = '\includegraphics[width=0.45\\textwidth]{'
        if len(figlist)==1 : figpre = figpre.replace('0.45','1.00')
        figpost= '} \n'
        f = ''
        for fig in figlist:
            f += figpre + fig + figpost
        return f
if __name__ == '__main__' :

    makeFigures = False
    if len(sys.argv)>1:
        if sys.argv[1].lower()=='print': makeFigures=True
    pc = pmtcal(makeFigures=makeFigures)

    #pc.figurePacker(0)
    pc.main()
