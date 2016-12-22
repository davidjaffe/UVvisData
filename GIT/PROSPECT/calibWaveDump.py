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
from scipy.stats.mstats import chisquare
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import datetime,os

import ROOT
import graphUtils
import gfit
import readLogFile

class calibWaveDump():
    def __init__(self):


        self.gU = graphUtils.graphUtils()
        self.gfit = gfit.gfit()
        self.rLF = readLogFile.readLogFile()

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
    def loop(self,fn='Output/run00073.root'):
        '''
        loop over relevant hists in input root file and
        fit them to extract mean, error on mean, sigma of gaussian, 
        and calculate area of peak (Scorr) and uncertainty (Ecorr)
        return dict with key =lifetime cut, value = [mean,emean,sgm, Scorr,Ecorr]
        '''
        results = {}
        if not os.path.isfile(fn): return None
        hf = ROOT.TFile.Open(fn)
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
        self.gU.drawMultiHists(hlist,fname=dn,figdir=self.figdir,dopt='pe',statOpt=0)
        return results
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
        
        X,Y,dY = [],{},{}
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
                results = cWD.loop(fn=fn)
                GOOD = False
                for l in results:
                    if results[l][3]>Threshold : GOOD = True
                if GOOD:
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
                                Y[x],dY[x] = [],[]
                            Y[x].append(norm[x][3])
                            dY[x].append(norm[x][4])
                    print words,'\n',stuff
        #print 'X',X
        #print 'Y',Y
        #print 'dY',dY
        X = numpy.array(X)
        dX = numpy.array([0. for x in X])
        T = numpy.array(T)
        tmg = self.gU.makeTMultiGraph('Rate_vs_run')
        tmgT= self.gU.makeTMultiGraph('Rate_vs_time')

        for l in Y:
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
            

        self.gU.drawMultiGraph(tmg,figdir=self.figdir,abscissaIsTime=False,xAxisLabel='Run number',yAxisLabel='Coincidence rate(Hz)')
        self.gU.drawMultiGraph(tmgT,figdir=self.figdir,abscissaIsTime=True,xAxisLabel='Time',yAxisLabel='Coincidence rate(Hz)')
        
            
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
    cWD.main()
    if 0:
        for i in range(31,90):
            fn='Output/run'+str(i).zfill(5)+'.root'
            cWD.loop(fn=fn)
