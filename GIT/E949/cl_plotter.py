#!/usr/bin/env python
'''
plot cls curves from arXiv: arXiv:0903.0030 from https://www.phy.bnl.gov/e949/E949Archive/E949_results.tar.gz
20191231
'''
import math
import sys,os

import datetime
import numpy
import random
#import copy
from scipy.interpolate import interp1d
import re
import glob # used in __init__

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,  AutoMinorLocator)

class cl_plotter():
    def __init__(self,debug=0,drawToFile=False):

        self.debug      = debug
        self.drawToFile = drawToFile

        self.supplDir = 'SUPPLEMENTAL/'
        self.filelist = glob.glob(self.supplDir + '*SM*.dat')

        self.figDir = 'FIGURES/'
        
        print 'cl_plotter.__init__ Done'
        return
    def main(self):
        '''
        evocative name!
        '''
        X,Y = {},{}
        for filename in self.filelist:
            basename = os.path.basename(filename)
            name = basename.replace('.dat','')
            br,CLs = self.read_cls(filename)
            X[name] = br
            Y[name] = CLs
            BF2 = br[numpy.argmin(abs(CLs-0.16))]
            BF1 = br[numpy.argmin(abs(CLs-0.84))]
            title = name + ', 68%CL interval ({0:.2f},{1:.2f})e-10'.format(BF1,BF2)
            self.drawIt(br,CLs,'br','CLs',title,mark='.',xlims=[0.,1.])
            self.drawIt(br,CLs,'br','CLs',title,mark='.',xlims=[1.1,1.2])

        self.drawMany(X,Y,sorted(X.keys()),'Br(1e-10)','CLs','',xlims=[0.,6.])
        #self.combine(X,Y)
        return
    def combine(self, X,Y):
        '''
        compare combination of SM_I and SM_II with SM
        '''
        f = {}
        for key in X:
            f[key] = interp1d(X[key],Y[key])
        Ynew = []
        Xnew = []
        keySM = 'cls_SM'
        k1    = 'cls_SM_I'
        k2    = 'cls_SM_II'
        low   = max(min(X[k1]),min(X[k2]))
        high  = min(max(X[k1]),max(X[k2]))
        for br in X[keySM]:
            if low<=br and br<=high:
                c1 = f[k1](br)
                c2 = f[k2](br)
                c = c1*c2*(1. - math.log(c1*c2))
                Ynew.append(c)
                Xnew.append(br)
        X['new'] = numpy.array(Xnew)
        Y['new'] = numpy.array(Ynew)

        
        self.drawMany(X,Y,sorted(X.keys()),'Br(1e-10)','CLs','',xlims=[0.,6.])
        return
    def read_cls(self,filename=None):
        '''
        return arrays of br,CLs from input filename
        column 1: br = K+ -> pi+,nu,nubar branching ratio times 10e10
        column 2: CLs = signal conf. level as defined by Junk's method. All uncertainties included
        Note: file cls_SM.dat has text that must be skipped
        '''
        if filename is None: sys.exit('cl_plotter.read_cls ERROR No input filename given')
        use = os.path.basename(filename)!='cls_SM.dat'
        
        f = open(filename,'r')
        if self.debug>0: print 'cl_plotter.read_cls Opened',filename
        br,CLs = [],[]
        for line in f:
            if use:
                s = line.split()
                br.append( float(s[0] ) )
                CLs.append(float(s[1] ) )
            if '====' in line: use = True
        f.close()
        if self.debug>-1: print 'cl_plotter.read_cls Read',len(br),'lines from',filename
        br = numpy.array(br)
        CLs= numpy.array(CLs)
        return br,CLs
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
            print 'cl_plotter.drawIt wrote',figpdf
        else:
            plt.show()
        return    
    def drawMany(self,x,y,leg,xtitle,ytitle,title,ylog=False,xlims=None,ylims=None,loc='best'):
        '''
        draw many graphs values on same plot defined by dicts x,y with keywords in leg

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
        if self.debug>1: print 'cl_plotter.drawMany major,minor',major,minor,'xlims',xlims


        ls = ['-','--','-.',':','-','--',':']
        ls.extend(ls)
        c  = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        c.extend(c)

        
        for i,key in enumerate(leg):
            X = numpy.array(x[key])
            Y = numpy.array(y[key])
            ax.plot(X,Y,linestyle=ls[i],color=c[i],label=key)

        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : ax.yscale('log')
        if xlims is not None: plt.xlim(xlims)
        if ylims is not None: plt.ylim(ylims)

            
        ax.legend(loc=loc)
        if self.drawToFile : 
            fn = self.titleAsFilename(title)
            figpdf = 'FIG_'+fn + '.pdf'
            figpdf = self.figDir + figpdf
            plt.savefig(figpdf)
            print 'cl_plotter.drawMany wrote',figpdf
        else:
            plt.show()
        return    
if __name__ == '__main__' :

    debug    = 0 # >0 gives output
    drawToFile = False # plots go to file instead of to terminal (use savefig() instead of show())
    if len(sys.argv)>1:
        if sys.argv[1].lower()=='help' : sys.exit( 'usage: python combiner.py drawEach debug drawToFile turnOnSyst studyVar' )
        debug = int(sys.argv[1])
    if len(sys.argv)>2:
        if sys.argv[2].lower()=='drawtofile' : drawToFile = True
        
    clp = cl_plotter(debug=debug,drawToFile=drawToFile)
    clp.main()
