#!/usr/bin/env python
'''
plots of charge for certain runs to study optical effects with viton (LiLS#8, elog91)
20170721
'''
import math
import sys
import random
import numpy

import datetime,os
import time
import ROOT
import graphUtils

class optical20170721():
    def __init__(self):
        self.rfiles = [ ['/Users/djaffe/work/WaveDumpData/rootfiles/run01134_ts1499445083.root','0deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01137_ts1499451494.root','180deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01138_ts1499458474.root','90deg'],
                        ['/Users/djaffe/work/WaveDumpData/rootfiles/run01139_ts1499460963.root','0deg']
                        ]
        self.gU = graphUtils.graphUtils()
        self.figdir = 'Figures/optical20170721/'
        return
    def loop(self):
        '''
        loop over trees in multiple root files, leave files open to make sure hists stay in scope
        '''
        debug = False
        cadence = 1
        if debug: cadence = 100000

        abs_time_cut = 0.004
        psdlo,psdhi,psdCut = 0., 1., 0.3
        Qlo,Qhi = 0.0, 0.05
        nx,ny = 100,100

        
        
        
        h1,h2 = {},{}
        for blob in self.rfiles:
            fn,words = blob[0],blob[1]
            runnum = fn.split('/')[-1].split('_')[0]
            Title1 = 'Q '+runnum
            Name1 = Title1.replace(' ','_')
            Title1 += ' ' + words
            H1 = ROOT.TH1D(Name1,Title1,nx,Qlo,Qhi)
            Title2 = 'PSD v Q '+runnum
            Name2 = Title2.replace(' ','_')
            Title2 += ' ' + words
            H2 = ROOT.TH2D(Name2,Title2,nx,Qlo,Qhi,ny,psdlo,psdhi)
            h1[runnum] = H1
            h2[runnum] = H2
        for blob in self.rfiles:
            fn,words = blob[0],blob[1]
            rf = ROOT.TFile.Open(fn,'r')
            runnum = fn.split('/')[-1].split('_')[0]
            print runnum
            H1 = h1[runnum]
            H2 = h2[runnum]
            tree = rf.Get('tree')
            entries = tree.GetEntriesFast()
            for jentry in xrange(0,entries+1,cadence):
                ientry = tree.LoadTree(jentry)
                if ientry<0:
                    break
                nb = tree.GetEntry(jentry)
                Q = tree.QtotalCh0
                PSD=tree.psdCh0
                abt=tree.abs_time
                if debug: print 'abt,Q,PSD',abt,Q,PSD
                if abt>abs_time_cut:
                    H2.Fill(Q,PSD)
                    if PSD>psdCut and PSD<psdhi: H1.Fill(Q)
        newfn = 'new0.root'
        new0 = ROOT.TFile(newfn,'RECREATE')
        print 'new0',new0,'h1,h2',h1,h2

        for h in h1:
            print h1[h]
            new0.WriteTObject(h1[h])
        for h in h2: new0.WriteTObject(h2[h])
        new0.Close()
        print 'wrote',len(h1)+len(h2),'to',newfn

        if debug : print 'h1.values()',h1.values(),  ' [ h1.values() ]',[ h1.values() ]

        self.gU.drawMultiObjects( h1.values() ,'Q_lils8',figdir=self.figdir,biggerLabels=False)
        self.gU.drawMultiObjects([ h1.values() ],'Q_lils8_overlay',figdir=self.figdir,biggerLabels=False,statOpt=0,addLegend=True,changeColors=True)
        self.gU.drawMultiObjects( h2.values() , 'PSD_v_Q_lils8',figdir=self.figdir,biggerLabels=False,dopt='COLZ')
        
        return
if __name__ == '__main__' :
    OP = optical20170721()
    OP.loop()                
    
