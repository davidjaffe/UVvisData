#!/usr/bin/env python
'''
eventually combine E949/E787 with NA62 for Br(K+ => pi+,nu,nubar)
20191023
20191205 clean up. copy original to old_combiner.py and remove unused code
'''
import math
import sys,os

import datetime
import numpy
#import copy

import re
#import glob # used in __init__

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,  AutoMinorLocator)

class combiner():
    def __init__(self,debug=0,drawEach=True,drawToFile=False):

        self.debug = debug
        self.drawEach = drawEach
        self.drawToFile = drawToFile
        print 'combiner.__init__ debug',self.debug,'drawEach',self.drawEach,'drawToFile',self.drawToFile

        self.AssumedBr = self.AssumedBR = 1.73e-10 # PRD79, 092004 (2009) Table VIII

        print 'combiner.__init__ self.AssumedBr',self.AssumedBr,' *****************'

        rmi,rma,steps = 0.,10.,10001
        if self.debug>1 : steps = 11
        dr = (rma-rmi)/float(steps-1)
        self.ratioRange = [rmi+float(i)*dr for i in range(steps)]

        self.figDir = 'FIGURES/'
        
        print 'combiner.__init__ Did something'
        return
    def main(self):
        '''
        cleverly named main routine for loading E787/E949 data and computing -2*loglikelihood
        '''
        debug = self.debug
        drawEach = self.drawEach

        # load data, report it, correct it to have same assumed branching fraction, report that,
        # then calculate m2ll = -2*loglike for each dataset
        cands = self.loadData()
        self.reportData(cands,mode='raw')
        cands = self.setAssBr(cands)
        self.reportData(cands,mode='same_assumed_Br')
        cands = self.fillM2LL(cands)

        # plot and combine likelihoods for all datasets
        globalM2LL = None
        x = numpy.array(self.ratioRange)
        M2LL = {}
        xtitle = 'Br(K+ => pi+,nu,nubar)/'+str(self.AssumedBr)
        ytitle = '-2*loglikelihood'
        for dataset in sorted(cands):
            cand = cands[dataset]
            m2ll = numpy.array(cand['m2ll'])
            M2LL[dataset] = m2ll-min(m2ll)
            ratmin = ', min at '+str(x[numpy.argmin(m2ll)])
            print 'combiner.main dataset',dataset,ratmin
            if debug>1 : print 'combiner.main dataset,len(x),len(m2ll)',dataset,len(x),len(m2ll)
            if debug>2 : print 'combiner.main dataset',dataset,[str(a)+'/{0:.2f}'.format(b) for a,b in zip(x,m2ll)]
            if drawEach : self.drawIt(x,m2ll,xtitle,ytitle,dataset+ratmin,mark='-')
            if globalM2LL is None:
                globalM2LL = numpy.array(m2ll)
            else:
                globalM2LL += numpy.array(m2ll)
                
        M2LL['all'] = globalM2LL-min(globalM2LL)
        ratmin = ', min at '+str(x[numpy.argmin(globalM2LL)])
        title = 'All E787/E949 data'+ratmin
        print 'combiner.main',title
        if drawEach: self.drawIt(x,globalM2LL,xtitle,ytitle,title,mark='-')
            
        title = '-2*loglikelihood'
        loc = 'upper right'
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title,loc=loc)
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title+' restricted y range',ylims=[0.,10.],loc=loc)
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title+' restricted x and y ranges',ylims=[0.,10.],xlims=[0.,2.],loc=loc)
        return
    def m2loglike(self,cand,ratio):
        '''
        calculate -2 * log likelihood from NK,Atot,[s/b], given ratio = BR/self.AssumedBr
        '''
        NK = cand['NK']
        Atot=cand['Atot']
        soverb = cand['soverb']
        like = ratio*self.AssumedBr*NK*Atot
        for x in soverb:
            like -= math.log(1. + ratio*x)
        like *= 2.
        return like
    def fillM2LL(self,cands):
        '''
        loop over datasets and add array of -2*loglike(ratio) for ratio in self.ratioRange to dict cands
        Note that input dict cands is modified by this module.
        '''
        debug = self.debug
        for dataset in sorted(cands):
            cand = cands[dataset]
            if debug>0: print 'combiner.fillM2LL dataset,soverb',dataset,cand['soverb']
            if 'm2ll' in cand:
                sys.exit('combiner.fillM2LL ERROR key `m2ll` already exists for dataset '+dataset+', perhaps due to multiple calls to this routine?')
            m2ll = []
            for ratio in self.ratioRange:
                x = self.m2loglike(cand,ratio)
                m2ll.append(x)
            cands[dataset]['m2ll'] = m2ll
        return cands
    def loadData(self):
        '''
        return dict loaded with all E787/E949 data
        For each dataset we have 
        dataset name
        journal reference
        NK = stopped kaons
        Atot = total acceptance
        Ncand = number of candidates
        AssumedBr = assumed B(K+ => pi+,nu,nubar) used for s_i/b_i ratio
        s_i/b_i = signal to background ratio in ith cell (cells containing candidates)
        '''
        cands = {}
        
        dataset = 'pnn1_E787_95-7'
        journal = 'PRL88_041803'
        NK = 3.2e12
        Atot = 2.1e-3
        Ncand = 1
        soverb = [35.]
        AssumedBr = 7.5e-11
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn1_E787_98'
        journal = 'PRL88_041803'
        NK = 2.7e12
        Atot = 1.96e-3
        Ncand = 1
        soverb = [3.6]
        AssumedBr = 7.5e-11
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn1_E949'
        journal = 'PRD77_052003'
        NK = 1.77e12
        Atot = 1.694e-3
        AssumedBr = self.AssumedBr
        Ncand = 1
        b = 5.7e-5
        s = 3.628e5*AssumedBr
        soverb = [s/b]
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn2_E787_96'
        journal = 'PLB537_2002_211'
        NK = 1.12e12
        Atot = 0.765e-3
        Ncand = 1
        b = 0.734 # +- 0.117
        AssumedBr = self.AssumedBr
        s = NK*Atot*AssumedBr*float(Ncand)
        soverb = [s/b]
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn2_E787_97'
        journal = 'PRD70_037102'
        NK = 0.61e12
        Atot = 0.97e-3
        Ncand = 0
        AssumedBr = self.AssumedBr
        soverb = []
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn2_E949'
        journal = 'PRD79_092004'
        NK = 1.71e12
        Atot = 1.37e-3 #(+-0.14e-3)
        Ncand = 3
        AssumedBr = 1.73e-10
        soverb = [0.47, 0.42, 0.20]
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}
        
        return cands
    def setAssBr(self,cands):
        '''
        return cands with s/b using a common assumed branching fractions self.AssumedBr
        '''
        for dataset in cands:
            cand = cands[dataset]
            factor = self.AssumedBr/cand['AssumedBr']
            soverb = cand['soverb']
            newsb = [factor*x for x in soverb]
            cand['soverb'] = newsb
            cand['AssumedBr'] = self.AssumedBr
        print 'combiner.setAssBr set assumed Br to',self.AssumedBr,'for s/b'
        return cands
    def reportData(self,cands,mode=None):
        '''
        report contents of all data in cands
        mode = 'raw' ==> no change to data
        mode = 'sameBr' ==> give all s/b at self.AssumedBr
        '''
        units = {'NK':1e12, 'Atot':1e-3, 'AssumedBr':1e-10, 'SES':1.e-10}

        
        print 'combiner.reportData mode is',mode,'. `AssBr` means Assumed Branching fraction of K+ => pi,nu,nubar for sig/bkg column'
        print '{0:<15} {1:>5}({2:5.0e}) {3:>5}({4:5.0e}) {5:>5}({6:5.0e}) {7:>5}({8:5.0e}) {9:>5} {10:>15} {11:>15}'.format('dataset','NK',units['NK'],'Atot',units['Atot'],'SES',units['SES'],'AssBr',units['AssumedBr'],'Ncand','sig/bkg','journal')
        #           0      1    2           3       4            5      6            7      8                 9        10        11
        for dataset in sorted(cands):
            cand = cands[dataset]
            NK = cand['NK']/units['NK']
            Atot = cand['Atot']/units['Atot']
            SES = 1./cand['NK']/cand['Atot']/units['SES']
            AssBr = cand['AssumedBr']/units['AssumedBr']
            Ncand = cand['Ncand']
            soverb= cand['soverb']
            journal = cand['journal']
            wsb = ''
            for x in soverb:
                wsb += '{0:>5.2f}'.format(x)
            print '{0:<15} {1:>12.2f} {2:>12.2f} {3:>12.2f} {4:>12.2f} {5:>5} {6:>15} {7:>15}'.format(dataset,NK,Atot,SES,AssBr,Ncand,wsb,journal)
            
        return
    def titleAsFilename(self,title):
        '''
        return ascii suitable for use as a filename
        list of characters to be replaced is taken from https://stackoverflow.com/questions/4814040/allowed-characters-in-filename
        '''
        r = {'_': [' ', ',',  '\\', '/', ':', '"', '<', '>', '|'], 'x': ['*']}
        filename = title
        filename = ' '.join(filename.split()) # substitutes single whitespace for multiple whitespace
        for new in r:
            for old in r[new]:
                if old in filename : filename = filename.replace(old,new)
        return filename
    def drawIt(self,x,y,xtitle,ytitle,title,ylog=False,xlims=None,ylims=None,mark='o-'):
        '''
        draw graph defined by x,y

        '''
        plt.clf()
        plt.grid()
        plt.title(title)

        X = numpy.array(x)
        Y = numpy.array(y)
        plt.plot(X,Y,mark)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : plt.yscale('log')
        if xlims is not None: plt.xlim(xlims)
        if ylims is not None: plt.ylim(ylims)

        if self.drawToFile:
            fn = self.titleAsFilename(title)
            figpdf = 'FIG_'+fn + '.pdf'
            figpdf = self.figDir + figpdf
            plt.savefig(figpdf)
            print 'combiner.drawIt wrote',figpdf
        else:
            plt.show()
        return    
    def drawMany(self,x,y,xtitle,ytitle,title,ylog=False,xlims=None,ylims=None,loc='best'):
        '''
        draw many graphs with same abscissa and different ordinate values on same plot defined by x,y
        y = dict
        ytitle = keys of dict

        '''
        fig,ax = plt.subplots()
        plt.grid()
        plt.title(title)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))


        ls = ['-','--','-.',':','-','--',':']
        c  = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
        
        X = numpy.array(x)
        for i,key in enumerate(ytitle):
            Y = numpy.array(y[key])
            ax.plot(X,Y,linestyle=ls[i],color=c[i],label=key)

        plt.xlabel(xtitle)
#        plt.ylabel(ytitle)
        if ylog : ax.yscale('log')
        if xlims is not None: plt.xlim(xlims)
        if ylims is not None: plt.ylim(ylims)

            
        ax.legend(loc=loc)
        if self.drawToFile : 
            fn = self.titleAsFilename(title)
            figpdf = 'FIG_'+fn + '.pdf'
            figpdf = self.figDir + figpdf
            plt.savefig(figpdf)
            print 'combiner.drawMany wrote',figpdf
        else:
            plt.show()
        return    
if __name__ == '__main__' :

    drawEach = False # draw loglikelihood for each dataset?
    debug    = 0 # >0 gives output
    drawToFile = False # plots go to file instead of to terminal (use savefig() instead of show())
    if len(sys.argv)>1:
        if sys.argv[1].lower()=='draweach': drawEach=True
        if sys.argv[1].lower()=='help' : sys.exit( 'usage: python combiner drawEach debug drawToFile' )
    if len(sys.argv)>2:
        debug = int(sys.argv[2])
    if len(sys.argv)>3:
        if sys.argv[3].lower()=='drawtofile' : drawToFile = True
        
    cb = combiner(debug=debug,drawEach=drawEach,drawToFile=drawToFile)
    cb.main()
