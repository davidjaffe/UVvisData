#!/usr/bin/env python
'''
eventually combine E949/E787 with NA62 for Br(K+ => pi+,nu,nubar)
20191023
'''
import math
import sys,os

import datetime
import numpy
#import copy

import re
import glob # used in __init__

import matplotlib.pyplot as plt

class combiner():
    def __init__(self):

        self.debug = 1

        self.AssumedBR = 1.73e-10 # PRD79, 092004 (2009) Table VIII
        self.E949data = 'DATA/E949_data_PRD79_092004.txt'
        # acceptance factors for tight cuts from Table IX of PRD79, 092004
        self.E949AccFac = {'KINT':0.812, 'TDT':0.812, 'DCT':0.911, 'PVT':0.522}
        Temp = {}
        for cut in self.E949AccFac: Temp[cut] = self.E949AccFac[cut]
        for cutT in Temp:
            AccT = self.E949AccFac[cutT]
            cutR = cutT[:-1] + 'R'
            self.E949AccFac[cutR] = 1. - AccT
        print '{0:>6} {1:>6} Relative acceptance of each cut combiner.__init__'.format('Cut','Acc')
        for cut in sorted(self.E949AccFac):
            print '{0:>6} {1:>6.3f}'.format(cut,self.E949AccFac[cut])
        
        print 'combiner.__init__ Did something'
        return
    def reproE949(self):
        '''
        try to reproduce the E949/E787 likelihood vs Br
        '''
        E949cells = self.readE949()
        Acc = {}
        for cell in E949cells:
            defn = E949cells[cell][0]
            cutlist = defn.split('.')
            acc = 1.
            for cutn in cutlist:
                if cutn in self.E949AccFac: acc *= self.E949AccFac[cutn] 
            Acc[cell] = acc
        
        
        print '{0:>6} {1:>18} {2:>6} {3:>6} {4:>10}'.format('cell','Defintion','b','s/b','Sacc')
        for cell in sorted(E949cells):
            print '{0:>6} {1:>18} {2:>6.3f} {3:>6.2f} {4:>10.4f}'.format(cell,E949cells[cell][0],E949cells[cell][1],E949cells[cell][2],Acc[cell])

        x,y = [],[]
        for iratio in range(0,210,10):
            ratio = float(iratio)/100.
            like = self.likeE949(ratio,E949cells)
            x.append(ratio)
            y.append(like)
        ymi = min(y)
        y = [z-ymi for z in y]

        print 'combiner.reproE949 SES={0:.2e} AssumedBR/SES={1:.2e}'.format(self.SES,self.AssumedBR/self.SES)
        
        self.drawIt(x,y,'Br scale factor','-2loglike','whats the frequency kenneth?')
            
        return
    def likeE949(self,ratio,cells):
        '''
        calculate -2 * log likelihood from s/b in cells for ratio = BR/self.AssumedBR
        '''

        SESpnn2b = 4.28e-10
        SESpnn1  = 2.55e-10
        SESpnn2a = 6.87e-10
        SES = self.SES = 1./( 1./SESpnn2b + 1./SESpnn1 + 1./SESpnn2a )
        likeA = ratio*self.AssumedBR/SES
        likeB = 0.
        for cell in cells:
            soverb = cells[cell][2]
            if cell==1: likeB -= math.log(1. + ratio*soverb)
        likeA = 2.*likeA
        likeB = 2.*likeB
        print ratio,likeA,likeB
        like = likeA+likeB
        return like
    def readE949(self,debug=0):
        '''
        read cell, b, s/b for E949/E787 data
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
                cell = int(s[0])
                defn = s[1]
                b    = float(s[2])
                soverb=float(s[3])
                if cell in E949cells:
                    sys.exit('combiner.readE949 ERROR Duplicate cell#'+str(cell))
                E949cells[cell] = [defn,b,soverb]
        f.close()
        print 'combiner.readE949 Read',len(E949cells),'cells from',self.E949data
        return E949cells
    def drawIt(self,x,y,xtitle,ytitle,title,figDir=None,ylog=False,xlims=None,ylims=None):
        '''
        draw graph defined by x,y

        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        X = numpy.array(x)
        Y = numpy.array(y)
        plt.plot(X,Y,'o-')
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
if __name__ == '__main__' :
   
    cb = combiner()
    cb.reproE949()
