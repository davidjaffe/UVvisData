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
import random
#import copy

import re
#import glob # used in __init__

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,  AutoMinorLocator)

class combiner():
    def __init__(self,debug=0,drawEach=True,drawToFile=False,turnOnSyst=False,studyVar=False):

        self.debug      = debug
        self.drawEach   = drawEach
        self.drawToFile = drawToFile
        self.turnOnSyst = turnOnSyst
        self.studyVar   = studyVar
        print 'combiner.__init__ debug',self.debug,'drawEach',self.drawEach,'drawToFile',self.drawToFile,'turnOnSyst',self.turnOnSyst,'studyVar',self.studyVar

        self.systOn = False
        self.systAcc = 0.10
        self.systN   = 10000
        if self.systOn : print 'combiner.__init__ SYSTEMATIC VARIATIONS APPLIED. systAcc,systN',self.systAcc,self.systN
        
        self.AssumedBr = self.AssumedBR = 1.73e-10 # PRD79, 092004 (2009) Table VIII

        print 'combiner.__init__ self.AssumedBr',self.AssumedBr,' *****************'

        self.Groups = None
        self.dataSets = None

        # define the binning to scan for the ratio wrt the assumed BR
        r = numpy.arange(0.,0.1,0.05)
        r = numpy.append( r, numpy.arange(0.1,0.8,0.05)  )
        r = numpy.append( r, numpy.arange(0.8,1.2,0.005) )
        r = numpy.append( r, numpy.arange(1.2,2.5,0.1)   )
        r = numpy.append( r, numpy.arange(2.5,10.,0.5)   )
        # smaller steps near minima of fits to subsets of data
        dx,step = 0.2,0.001
        for x0 in [1.50, 2.70, 4.95, 7.80, 1.70]:
            r = numpy.append( r, numpy.arange(x0-dx,x0+dx,step) )

        r = numpy.unique( numpy.sort(r) )
        self.ratioRange = r

        self.figDir = 'FIGURES/'
        
        print 'combiner.__init__ Did something'
        return
    def reportSyst(self):
        '''
        return string with report on systematics parameters
        '''
        s = 'NO SYSTEMATIC VARIATIONS. Flag is False'
        if self.systOn: s = 'SYSTEMATIC VARIATION systAcc {0:.2f} systN {1}'.format(self.systAcc,self.systN)
        return s
    def studyVariations(self,cands):
        '''
        make a deepcopy of the input cands and monkey around
        to understand dependency of the BR that minimizes -2*loglike on various parameters
        '''
        import copy

        x = numpy.array(self.ratioRange)
        x = numpy.append(x, numpy.arange(0.9,1.1,.0001) )
        x = numpy.sort(x)
        x = numpy.unique(x)
        ytitle = 'Br(K+ => pi+,nu,nubar)/'+str(self.AssumedBr)


    
        local_cands = self.collate(cands)


        vary = {'NK'    :numpy.arange( 0.90, 1.10, 0.01),
                'soverb':numpy.arange( 0.80, 1.20, 0.02) }

        for key in vary:
            br = []


            for v in vary[key]:
                title = 'Variation is ' + key + ' times ' + str(v)
                allc = copy.deepcopy(local_cands)
                CAND = allc['all']
                if type(CAND[key]) is list:
                    l = [v*z for z in CAND[key]]
                    CAND[key] = l
                else:
                    CAND[key] = v*CAND[key]
                allc = self.fillM2LL(allc,x)
                if self.debug>1: print 'combiner.studyVariations allc.keys()',allc.keys()
                if self.debug>1: print 'combiner.studyVariations allc[all].keys()',allc['all'].keys()
                m2ll = numpy.array(allc['all']['m2ll'])
                m2ll = m2ll-min(m2ll)
                if self.debug>1: print 'combiner.studyVariations',title,'len(x),len(m2ll)',len(x),len(m2ll)
                xatmin = x[numpy.argmin(m2ll)]
                if self.debug>0: print 'combiner.studyVariations',title,'minimized at',xatmin
                br.append( xatmin )
                #self.drawIt(x,m2ll,xtitle,ytitle,title,mark='-')
                del allc
            i = numpy.argmin(abs(numpy.array(br)-1.0))
            i1= max(i-1,0)
            i2= min(i+1,len(br)-1)
            slope = (br[i2]-br[i1])/(vary[key][i2]-vary[key][i1])
            y1,y2= br[i1],br[i2]
            x1,x2= vary[key][i1],vary[key][i2]
            best = ((x2-x1) - (x2*y1-y2*x1))/(y2-y1)
            print 'combiner.studyVariations {0}, slope of scale factor {1:.3f}, best scale factor {2:.3f}'.format(key,1./slope,best)
            self.drawIt(vary[key],br,key+' scale factor',ytitle,'Variation of '+key,mark='o-')
            
        return
    def main(self):
        '''
        cleverly named main routine for loading E787/E949 data and computing -2*loglikelihood
        '''
        debug = self.debug
        drawEach = self.drawEach
        x = numpy.array(self.ratioRange)
        xtitle = 'Br(K+ => pi+,nu,nubar)/'+str(self.AssumedBr)
        ytitle = '-2*loglikelihood'

        # load data, report it, correct it to have same assumed branching fraction, report that,
        # then, if requested, study fitted Br for variations in input parameters
        # then calculate m2ll = -2*loglike for each dataset
        cands = self.loadData()
        self.reportData(cands,mode='raw')
        cands = self.setAssBr(cands)
        self.reportData(cands,mode='same_assumed_Br')
        if self.studyVar : self.studyVariations(cands)
            
        # group candidates, calculate minimum of chi2, report results by group
        groupCands = {}
        for group in sorted(self.Groups.keys()):
            groupCands[group] = self.collate(cands,keyw=group)[group]
            if debug>0: print 'combiner.main group',group,'groupCands[group].keys()',groupCands[group].keys(),'groupCands[group]',groupCands[group]
        if debug>0: print 'combine.main groupCands.keys()',groupCands.keys()
        groupCands = self.fillM2LL(groupCands)
        Results = {}
        gLL = {}
        for group in sorted(self.Groups.keys()):
            m2ll = numpy.array(groupCands[group]['m2ll'])
            gLL[group] = m2ll = m2ll-min(m2ll)
            xatmin = x[numpy.argmin(m2ll)]
            Results[group] = xatmin*self.AssumedBR
            if debug>0: print 'combine.main {0} minimized at BF {1:.2e}'.format(group,xatmin*self.AssumedBR)

        self.reportGroups(Results)
        title = 'Groups'
        loc = 'best'
        self.drawMany(x,gLL,xtitle,gLL.keys(),title,loc=loc)
        self.drawMany(x,gLL,xtitle,gLL.keys(),title+' restricted x and y ranges',ylims=[0.,4.],xlims=[0.,2.],loc=loc)
            
        cands = self.fillM2LL(cands)


        
        M2LL = {}
        # combined candidates and likelihood
        if self.turnOnSyst:
            self.systOn = True
            if self.systOn : print 'combiner.main Systematics calculation enabled for combined likelihood.',self.reportSyst()
            allcands = self.collate(cands)
            allcands = self.fillM2LL(allcands)
            m2ll = numpy.array(allcands['all']['m2ll'])
            m2ll = m2ll-min(m2ll)
            if debug>0: print 'combiner.main allcands minimized at',x[numpy.argmin(m2ll)]
            title = '-2*loglikelihood with systematics'
            self.drawIt(x,m2ll,xtitle,ytitle,title,mark='-')
            self.drawIt(x,m2ll,xtitle,ytitle,title+' restrict ranges',mark='-',xlims=[0.5,1.5],ylims=[0.,0.1])
            M2LL['all_with_syst'] = m2ll

            self.systOn = False

        if 0: # this stuff has been mostly replaced
            
            # plot and combine likelihoods for all datasets
            globalM2LL = None
            for dataset in sorted(cands):
                cand = cands[dataset]
                m2ll = numpy.array(cand['m2ll'])
                M2LL[dataset] = m2ll-min(m2ll)
                xatmin = x[numpy.argmin(m2ll)]
                BFatmin = ' BF={0:.2e}'.format(xatmin*self.AssumedBR)
                ratmin = ', min at '+str(xatmin)+BFatmin
                if debug>0: print 'combiner.main dataset',dataset,ratmin
                if debug>1 : print 'combiner.main dataset,len(x),len(m2ll)',dataset,len(x),len(m2ll)
                if drawEach : self.drawIt(x,m2ll,xtitle,ytitle,dataset+ratmin,mark='-')
                if globalM2LL is None:
                    globalM2LL = numpy.array(m2ll)
                else:
                    globalM2LL += numpy.array(m2ll)

            M2LL['all'] = globalM2LL-min(globalM2LL)
            ratmin = ', min at '+str(x[numpy.argmin(globalM2LL)])
            title = 'All E787/E949 data'+ratmin
            if debug>0: print 'combiner.main',title
            if drawEach: self.drawIt(x,globalM2LL,xtitle,ytitle,title,mark='-')

            title = '-2*loglikelihood'
            loc = 'upper right'
            self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title,loc=loc)
            self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title+' restricted y range',ylims=[0.,10.],loc=loc)
            self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title+' restricted x and y ranges',ylims=[0.,4.],xlims=[0.,2.],loc=loc)
            self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title+' fanatical x and y ranges',ylims=[0.,0.2],xlims=[0.8,1.2],loc=loc)
        return
    def m2loglike(self,cand,RATIO):
        '''
        calculate -2 * log likelihood from NK,Atot,[s/b], given ratio = BR/self.AssumedBr
        optionally include averaging over systematic variation of global acceptance 
        '''
        if type(cand['NK']) is list:
            NKlist = cand['NK']
            Atotlist = cand['Atot']
        else:
            NKlist = [cand['NK']]
            Atotlist = [cand['Atot']]
        soverb = cand['soverb']

        v = [1.]
        if self.systOn:
            v = numpy.random.normal(1.,self.systAcc,self.systN)

        like = 0.
        for f in v:
            ratio = f*RATIO
            for NK,Atot in zip(NKlist,Atotlist):
                like += ratio*self.AssumedBr*NK*Atot
            for x in soverb:
                like -= math.log(1. + ratio*x)
        like = like/float(len(v))
        like *= 2.
        return like
    def fillM2LL(self,cands,ratRange=None):
        '''
        loop over datasets and add array of -2*loglike(ratio) for ratio in ratRange to dict cands
        Note that input dict cands is modified by this module.
        ratRange defaults to self.ratioRange if no input is provided
        '''
        debug = self.debug
        ratioRange = ratRange
        if ratRange is None : ratioRange = self.ratioRange
        if debug>0 : print 'combiner.fillM2LL cands',cands
        for dataset in sorted(cands):
            cand = cands[dataset]
            if debug>0: print 'combiner.fillM2LL dataset,cand',dataset,cand
            if 'm2ll' in cand:
                sys.exit('combiner.fillM2LL ERROR key `m2ll` already exists for dataset '+dataset+', perhaps due to multiple calls to this routine?')
            m2ll = []
            for ratio in ratioRange:
                x = self.m2loglike(cand,ratio)
                m2ll.append(x)
            cands[dataset]['m2ll'] = m2ll
        return cands
    def collate(self,cands,keyw='all'):
        '''
        create a dict with candidates from datasets specified by keyw
        AssumedBr must be the same for all combined datasets
        Cannot perform collation if a dataset in cands contains -2*loglike array.
        combination defined by keyw must already be defined

        '''
        allcands = {}
        allcands[keyw] = {'NK':[], 'Atot':[], 'soverb':[], 'AssumedBr':None}
        groups = self.Groups
        AssBr = None

        if keyw=='all':
            setList = sorted(cands.keys())
        elif keyw in groups:
            setList = groups[keyw]
        else:
            sys.exit('combiner.collate ERROR Invalid keyw '+keyw)

        
        for dataset in setList: 
            if 'm2ll' in cands[dataset]:
                sys.exit('combiner.collate ERROR key `m2ll` in input dict cands for dataset '+dataset)
            cand = cands[dataset]
            NK = cand['NK']
            Atot = cand['Atot']
            soverb = cand['soverb']
            if AssBr is None: AssBr = cand['AssumedBr']
            if AssBr!=cand['AssumedBr']:
                print 'combiner.collate ERROR dataset,AssumedBr',dataset,cand['AssumedBr'],'is not equal to',AssBr,'found for first dataset'
                sys.exit('combiner.collate ERROR inconsistent assumed Br')
            allcands[keyw]['NK'].append( NK )
            allcands[keyw]['Atot'].append( Atot )
            allcands[keyw]['soverb'].extend( soverb )
            allcands[keyw]['AssumedBr'] = AssBr
        return allcands
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

        self.dataSets = []
        
        dataset = 'pnn1_E787_95-7'
        self.dataSets.append( dataset )
        journal = 'PRL88_041803'
        NK = 3.2e12
        Atot = 2.1e-3
        Ncand = 1
        #print '\n combiner.loadData TEMPORARY CHANGE S/B FOR PNN1 EVENT 95A ********************'
        soverb = [35.] ## Changing to 25.7 has no significant effect on combined BFs
        #print 'combiner.loadData TEMPORARY CHANGE S/B FOR PNN1 EVENT 95A\n'
        AssumedBr = 7.5e-11
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn1_E787_98'
        self.dataSets.append( dataset )
        journal = 'PRL88_041803'
        NK = 2.7e12
        Atot = 1.96e-3
        Ncand = 1
        soverb = [3.6]
        AssumedBr = 7.5e-11
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn1_E949'
        self.dataSets.append( dataset )
        journal = 'PRD77_052003'
        NK = 1.77e12
        Atot = 2.22e-3
        AssumedBr = self.AssumedBr
        Ncand = 1
        b = 5.7e-5
        s = 3.628e5*AssumedBr
        soverb = [s/b]
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn2_E787_96'
        self.dataSets.append( dataset )
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
        self.dataSets.append( dataset )
        journal = 'PRD70_037102'
        NK = 0.61e12
        Atot = 0.97e-3
        Ncand = 0
        AssumedBr = self.AssumedBr
        soverb = []
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}

        dataset = 'pnn2_E949'
        self.dataSets.append( dataset )
        journal = 'PRD79_092004'
        NK = 1.71e12
        Atot = 1.37e-3 #(+-0.14e-3)
        Ncand = 3
        AssumedBr = 1.73e-10
        soverb = [0.47, 0.42, 0.20]
        cands[dataset] = {'NK':NK, 'Atot':Atot, 'Ncand':Ncand, 'soverb':soverb, 'AssumedBr':AssumedBr, 'journal':journal}
        groups = {'pnn1_pub': ['pnn1_E787_95-7','pnn1_E787_98','pnn1_E949'],
                    'All E787': ['pnn1_E787_95-7','pnn1_E787_98','pnn2_E787_96','pnn2_E787_97'],
                    'All E949': ['pnn1_E949', 'pnn2_E949'],
                    'All pnn1': ['pnn1_E787_95-7','pnn1_E787_98','pnn1_E949'],
                    'All pnn2': ['pnn2_E787_96','pnn2_E787_97','pnn2_E949'],
                    'all'     : self.dataSets,
                    'E949 pnn1': ['pnn1_E949'],
                    'E949 pnn2': ['pnn2_E949']
                      }
        self.Groups = groups
        # E949 Technote K074.v1 Table 89 for fitted BR in 1.e-10 units
        # labelled 'Joss' in this file for reasons that cannot be revealed
        self.Joss  = {'pnn1_pub': 1.47, 
                    'All E787'  : 1.49, 
                    'All E949'  : 2.80, 
                    'All pnn1'  : 1.46,
                    'All pnn2'  : 5.11, 
                    'all'       : 1.73, 
                    'E949 pnn1' : 0.96,
                    'E949 pnn2' : 7.89
                      }

        return cands
    def reportGroups(self,Results):
        '''
        report content of self.Groups with fitted BR and Joss's fitted BR
        '''
        if self.Groups is None: sys.exit('combiner.reportGroups ERROR self.Groups not initialized')
        print '{0:^12} {1:^12} {2:^6} {3:<15} | {4}'.format('BF(1e-10)','Joss(1e-10)','BF/J','Group name','data sets')
        for group in sorted(self.Groups.keys()):
            mine = Results[group]/1.e-10
            Joss = self.Joss[group]
            r    = mine/Joss
            print '{0:>12.2f} {1:<12.2f} {2:^6.2f} {3:<15} |'.format(mine,Joss,r,group),'%s' % ' '.join(map(str,self.Groups[group]))
        return
            
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
        major = 1.
        if xlims is not None:
            if xlims[1]-xlims[0]<2: major = (xlims[1]-xlims[0])/10
        ax.xaxis.set_major_locator(MultipleLocator(major))
        minor = major/5.
        ax.xaxis.set_minor_locator(MultipleLocator(minor))
        if self.debug>1: print 'combiner.drawMany major,minor',major,minor,'xlims',xlims


        ls = ['-','--','-.',':','-','--',':']
        ls.extend(ls[::-1])
        c  = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        c.extend(c)
        
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
    turnOnSyst = False # include systematics in BR determination
    studyVar   = False # variation of inputs for alternate BR determinations
    if len(sys.argv)>1:
        if sys.argv[1].lower()=='draweach': drawEach=True
        if sys.argv[1].lower()=='help' : sys.exit( 'usage: python combiner.py drawEach debug drawToFile turnOnSyst studyVar' )
    if len(sys.argv)>2:
        debug = int(sys.argv[2])
    if len(sys.argv)>3:
        if sys.argv[3].lower()=='drawtofile' : drawToFile = True
    if len(sys.argv)>4:
        if sys.argv[4].lower()=='turnonsyst' : turnOnSyst = True
    if len(sys.argv)>5:
        if sys.argv[5].lower()=='studyvar'   : studyVar   = True
        
    cb = combiner(debug=debug,drawEach=drawEach,drawToFile=drawToFile,turnOnSyst=turnOnSyst,studyVar=studyVar)
    cb.main()
