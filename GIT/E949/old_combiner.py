#!/usr/bin/env python
'''
eventually combine E949/E787 with NA62 for Br(K+ => pi+,nu,nubar)
20191023
20191205 Deprecated
'''
import math
import sys,os

import datetime
import numpy
#import copy

import re
import glob # used in __init__

import matplotlib.pyplot as plt

class old_combiner():
    def __init__(self):

        self.debug = 1

        # directory with Joe's stuff
        self.JoeDir = '/Users/djaffe/work/GIT/E949/Joe/home/sher/test/'
        self.bogus_eclnames = ['ecl9598inside_x.dat','ecl.probbr.dat']
        

        self.AssumedBr = self.AssumedBR = 1.73e-10 # PRD79, 092004 (2009) Table VIII

        print 'combiner.__init__ self.AssumedBr',self.AssumedBr,' *****************'

        rmi,rma,steps = 0.,10.,10001
        #steps = 11
        dr = (rma-rmi)/float(steps-1)
        self.ratioRange = [rmi+float(i)*dr for i in range(steps)]

        
        self.E949data = 'DATA/E949_data_PRD79_092004.txt'
        # acceptance factors for tight cuts from Table IX of PRD79, 092004
        self.E949AccFac = {'KINT':0.812, 'TDT':0.812, 'DCT':0.911, 'PVT':0.522}
        Temp = {}
        for cut in self.E949AccFac: Temp[cut] = self.E949AccFac[cut]
        for cutT in Temp:
            AccT = self.E949AccFac[cutT]
            cutR = cutT[:-1] + 'R'
            self.E949AccFac[cutR] = 1. - AccT
        if 0:
            print '{0:>6} {1:>6} Relative acceptance of each cut combiner.__init__'.format('Cut','Acc')
            for cut in sorted(self.E949AccFac):
                print '{0:>6} {1:>6.3f}'.format(cut,self.E949AccFac[cut])
        
        print 'combiner.__init__ Did something'
        return
    def reproE949(self):
        '''
        try to reproduce the E949/E787 likelihood vs Br
        '''
        E949cells = self.readE949(debug=1)
        if 0: # acceptance numbers not used
            Acc = {}
            for cell in E949cells:
                defn = E949cells[cell][0]
                cutlist = defn.split('.')
                acc = 1.
                for cutn in cutlist:
                    if cutn in self.E949AccFac: acc *= self.E949AccFac[cutn] 
                Acc[cell] = acc
        
        
        print '{0:>6} {1:>18} {2:>6} {3:>6} {4:>6}'.format('cell','Definition','b','s/b','d')
        for cell in sorted(E949cells):
            if cell!='SES':
                print '{0:>6} {1:>18} {2:>6.3f} {3:>6.2f} {4:>6.0f}'.format(cell,E949cells[cell][0],E949cells[cell][1],E949cells[cell][2],E949cells[cell][3])
        if 'SES' in E949cells: print 'SES = {0:.2e} +- {1:.2e}'.format(E949cells['SES'][0],E949cells['SES'][1])

        useE949 = True
        if useE949:
            title = 'E949 combined result'
        else:
            title = 'E949 pnn1 only'
            
        x,y = [],[]
        rmi,rma,steps = 0.,10.,1001
        dr = (rma-rmi)/float(steps-1)
        ratio = rmi
        while ratio<rma: 

            if useE949:
                like = self.likeE949(ratio,E949cells)
            else:
                like = self.likepnn1(ratio)
            x.append(ratio)
            y.append(like)
            ratio += dr 
        ymi = min(y)
        ratmin = x[y.index(ymi)]
        y = [z-ymi for z in y]
        
        if useE949:
            print 'combiner.reproE949 SES={0:.2e} All. AssumedBR/SES={1:.2e}, AssumedBR={2:.2e} BRatmin={3:.2e}'.format(self.SES,self.AssumedBR/self.SES,self.AssumedBR,ratmin*self.AssumedBR)
            for key in sorted(self.SEStable):
                print 'combiner.reproE949 SES={0:.2e} {1}'.format(self.SEStable[key],key)

        xtitle = 'Branching fraction/'+str(self.AssumedBR)
            
        self.drawIt(x,y,xtitle,'-2loglike',title)
            
        return
    def likepnn1(self,ratio):
        '''
        calculate -2 * log likelihood for the single pnn1 event 
        from Section IV.B of PRD77 052003
        
        '''
        NK   = 1.77e12
        Atot = 1.694e-3
        Ai   = 1.21e-4 #  acceptance of cell with signal relative to Atot

        b    = 5.75e-5
        B    = ratio*self.AssumedBR
        s    = NK*Atot*Ai*B
        likeA = 2.*NK*Atot*B
        likeB = -2.*math.log(1. + s/b)
        like = likeA+likeB
        print 'combiner.likepnn1 ratio,BF,likeA,likeB,like',ratio,B,likeA,likeB,like
        return like
    def likeE949(self,ratio,cells):
        '''
        calculate -2 * log likelihood from s/b in cells for ratio = BR/self.AssumedBR
        '''
        debug = 0
        pnn2_E949_only = True
        
        SES = {}
        SES['pnn2_E949'] = 4.28e-10 # PRD79, 092004    # Omit this and likeE949 minimized at ratio=1.!?!?!?!
        if not pnn2_E949_only:
            SES['pnn1_E787_9597'] = 1./3.2e12/2.1e-3 # PRL88,041803 Tables I,II
            SES['pnn1_E787_98']   = 1./2.7e12/1.96e-3# PRL88,041803 Table I,II
            SES['pnn1_E949'] = 2.55e-10 # PRD77, 052003
            SES['pnn2_E787_96'] = 11.7e-10 # PLB537(2002)211
            SES['pnn2_E787_97'] = 16.9e-10 # PRD70(2004)037102

        self.SEStable = SES
        P = 0.
        for x in SES: P += 1./SES[x]
        P = 1./P
        SES = self.SES = P
        likeA = ratio*self.AssumedBR/SES
        likeB = 0.
        for cell in sorted(cells):
            if cell!='SES':
                soverb = cells[cell][2]
                d = cells[cell][3] # events observed
                likeB -= d*math.log(1. + ratio*soverb)
                if debug>1 and ratio==0. : print 'combiner.likeE949 cell,soverb',cell,soverb
        likeA = 2.*likeA
        likeB = 2.*likeB
        like = likeA+likeB
        if debug>0 : print 'combiner.likeE949 ratio,likeA,likeB,like',ratio,likeA,likeB,like
        return like
    def main(self):
        '''
        cleverly named main routine for loading E787/E949 data and computing -2*loglikelihood
        '''
        debug = 0
        drawEach = True
        cands = self.loadData()
        self.reportData(cands,mode='raw')
        cands = self.setAssBr(cands)
        self.reportData(cands,mode='same_assumed_Br')
        cands = self.fillM2LL(cands)
        
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
        if drawEach: self.drawIt(x,globalM2LL,xtitle,ytitle,title)
        title = '-2*loglikelihood'
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title)
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title,ylims=[0.,10.])
        self.drawMany(x,M2LL,xtitle,sorted(M2LL.keys()),title,ylims=[0.,10.],xlims=[0.,2.])
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
        debug = 0
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
    def readE949(self,debug=0):
        '''
        read cell, b, s/b, d for E949/E787 data
        '''
        f = open(self.E949data,'r')
        E949cells = {}
        for line in f:
            if debug>0: print 'combiner.readE949 line',line[:-1]
            if '#' not in line:
                r = line[:-1].split('\t')
                s = []
                for x in r:
                    q = x.strip()
                    if q!='': s.append(q)
                if debug>0: print 'combiner.readE949 s',s
                if s[0]=='SES':
                    ses = float(s[1]) # single event sensitivity
                    dses= float(s[2]) # uncertainty on ses
                    E949cells['SES'] = [ses,dses]
                else:
                    cell = int(s[0])
                    defn = s[1]
                    b    = float(s[2])
                    soverb=float(s[3])
                    d     =float(s[4])
                    if cell in E949cells:
                        sys.exit('combiner.readE949 ERROR Duplicate cell#'+str(cell))
                    E949cells[cell] = [defn,b,soverb,d]
        f.close()
        print 'combiner.readE949 Read',len(E949cells),'cells from',self.E949data
        return E949cells
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
    def readJoeDat(self,name):
        '''
        return dict of contents of ecl<name>.dat in Joe's directory
        '''
        debug = 0

        eclname = name if (name[:3]=='ecl') else 'ecl' + name + '.dat'
        if eclname in self.bogus_eclnames: return None

        fn = self.JoeDir + eclname
        f = open(fn,'r')
        print 'combiner.readJoeDat Opened',fn
        contents = {}
        for line in f:
            s = line[:-1].split(':')
            key = s[0].strip()
            value = float(s[1])
            if key not in contents : contents[key] = []
            contents[key].append(value)
        f.close()
        L = {}
        sameL = None
        wrongL = False
        for key in sorted(contents):
            L[key] = len(contents[key])
            if sameL is None: sameL = L[key]
            if sameL != L[key] : wrongL = True
        if wrongL or debug>0:
            if wrongL: print 'combiner.readJoeDat name',name,'ERROR? found keys with different number of entries'
            for key in sorted(L):
                print 'combiner.readJoeDat name',name,'key',key,'has',L[key],'entries'
        return contents
    def plotJoeDat(self,name,contents,xname,yname):
        '''
        given dict of contents of ecl<name>.dat plot values of yname vs xname
        '''
        
        eclname = name if (name[:3]=='ecl') else 'ecl' + name + '.dat'
        
        if xname not in contents:
            sys.exit('combiner.plotJoeDat name ' + eclname + ' ERROR invalid xname ' + xname)
        if yname not in contents:
            sys.exit('combiner.plotJoeDat name ' + eclname + ' ERROR invalid yname ' + yname)
        X = contents[xname]
        Y = contents[yname]
        yma = max(Y)
        i = Y.index(yma)
        xma = X[i]
        words = name + ':' + yname + ' max at ' + str(yma) + ' at ' + xname + ' of ' + str(xma)
        self.drawIt(X,Y,xname,yname,words)
        return
    def plotJoeDatLoop(self,name,contents,xname):
        '''
        given dict of contents of ecl<name>.dat plot all value of all keys vs values of key=xname
        '''
        for yname in sorted(contents):
            if yname!=xname:
                self.plotJoeDat(name,contents,xname,yname)
        return
    def showAllJoeDat(self,xname='signal scaling',yname='Xobs'):
        '''
        plot Xobs vs signal scaling for all ecl*.dat files
        '''
        filelist = glob.glob(self.JoeDir + 'ecl*.dat')
        for fn in filelist:
            name = os.path.basename(fn)
            contents = self.readJoeDat(name)
            if contents is not None : self.plotJoeDat(name,contents,xname,yname)
        return
    def drawIt(self,x,y,xtitle,ytitle,title,figDir=None,ylog=False,xlims=None,ylims=None,mark='o-'):
        '''
        draw graph defined by x,y

        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        X = numpy.array(x)
        Y = numpy.array(y)
        plt.plot(X,Y,mark)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : plt.yscale('log')
        if xlims is not None: plt.xlim(xlims)
        if ylims is not None: plt.ylim(ylims)

        if figDir is not None:
            figpdf = figDir + figpdf
            plt.savefig(figpdf)
            print 'abslen_model.drawIt wrote',figpdf
        else:
            plt.show()
        return    
    def drawMany(self,x,y,xtitle,ytitle,title,figDir=None,ylog=False,xlims=None,ylims=None):
        '''
        draw many graphs with same abscissa and different ordinate values on same plot defined by x,y
        y = dict
        ytitle = keys of dict

        '''
        from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,  AutoMinorLocator)


        fig,ax = plt.subplots()
        plt.grid()
        plt.title(title)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

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

            
        ax.legend(loc='best')
        if figDir is not None:
            figpdf = figDir + figpdf
            plt.savefig(figpdf)
            print 'abslen_model.drawIt wrote',figpdf
        else:
            plt.show()
        return    
if __name__ == '__main__' :
   
    cb = old_combiner()
    cb.main()
    if 0: 
        cb.reproE949()
    if 0: 
        cb.showAllJoeDat()
    if 0:
        name = '98inside'
        contents = cb.readJoeDat(name)
        xname = 'signal scaling'
        cb.plotJoeDatLoop(name,contents,xname)
