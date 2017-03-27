#!/usr/bin/env python
'''
make a TTree from a dict
20170215
'''
import sys
import numpy
import os

import ROOT
from array import array

class makeTTree():
    def __init__(self):
        return
    def makeTTree(self,C,fn='TTREE.ROOT',treename='tree',debug=False):
        '''
        given dict C create tree treename in file fn

        strings in Branch from http://stackoverflow.com/questions/19691593/python-array-of-strings-in-a-root-ttree#19718029
        floats, ints in Branch from http://wlav.web.cern.ch/wlav/pyroot/tpytree.html
        
        '''
        f = ROOT.TFile( fn, 'recreate' )
        t = ROOT.TTree( treename, treename)


        TD = {}
        orderedkeys = sorted( C.keys() )
        KIND = {}
        LENGTH = {}
        L = 0
        for idx,key in enumerate(orderedkeys):
            D = C[key]
            kind,decl = self.getType(D)
            KIND[key] = kind
            L = max(L,len(D))
            if debug : print 'makeTTree.makeTTree key,kind,decl,D',key,kind,decl,D,[x for x in D]
            if kind=='C':
                # stuff for string
                l = max([len(q) for q in C[key]])
                LENGTH[key] = l
                if debug: print 'self.makeTTree l',l
                TD[key] =   bytearray(l+1) 
                A = key + '['+str(l+1)+']/'+kind
            else:
                TD[key] = array(decl,[D[0]])
                A = key+'/'+kind
            if debug: print 'makeTTree.makeTTree key,idx,TD[key],A',key,idx,TD[key],A
            t.Branch(key,TD[key],A)
        if debug: print 'makeTTree after Branches TD[key]',TD[key]
        if debug: print ''
        if debug: t.Scan()
    
        for i in range(L):
            for idx,key in enumerate(orderedkeys):
                kind = KIND[key]
                D = C[key]
                if kind=='C':
                    for j in range(LENGTH[key]):
                        TD[key][j] = ' '
                    if debug: print 'prefill key,TD[key]',key,TD[key]
                    for j,c in enumerate(D[i]):
                        TD[key][j] = c
                else:
                    if debug: print 'makeTTree.makeTTree key,kind,i,D[i]',key,kind,i,D[i]
                    TD[key][0] =  D[i]
            if debug : print 'makeTTree.makeTTree i,TD',i,TD
            t.Fill()


        if debug : t.Print()
        f.Write()
        f.Close()
        print 'makeTTree.makeTTree Wrote TTree',treename,'to',fn
        return
    def getType(self,a):
        '''
        return type of variable in container a in two formats
        kind is used by TTree
        decl is used by array.array
        '''
        
        kind = None
        decl = None
        vt = type(a[0])
        if vt is int or vt is numpy.int64:
            kind = 'I'
            decl = 'i'
        elif vt is float or vt is numpy.float64:
            kind = 'F'
            decl = 'f'
        elif vt is str:
            kind = 'C'
            decl = 'c'
        return kind,decl
    def a2v(self,a,debug=False):
        '''
        convert input container a into ROOT.vector
        return type of variable in container
        and array type declaration
        '''
        kind,decl = None, None
        vt = type(a[0])
        if debug : print 'makeTTree.a2v type',vt
        if vt is int or vt is numpy.int64:
            V = ROOT.vector('int')()
            kind = 'I'
            decl = 'i'
        elif vt is float or vt is numpy.float64:
            V = ROOT.vector('float')()
            kind = 'F'
            decl = 'f'
        elif vt is str:
            V = ROOT.vector('string')()
            kind = 'C'
            decl = 'c'

        if debug : print 'makeTTree.a2v V',V
        for x in a:
            if kind=='F': y = float(x)
            if kind=='I': y = int(x)
            if kind=='C': y = str(x)
            V.push_back(y)
        return V,kind,decl
    def getTTree(self,fn=None,ttName=None):
        '''
        file dict with contents of TTree ttName in file with name fn
        only works for simple TTrees
        '''
        # bogosity check
        if fn is None: sys.exit('makeTTree.getTTree ERROR No file name specified')
        if ttName is None: sys.exit('makeTTree.getTTree ERROR No tree name given for file '+fn)
        #open the file
        myfile = ROOT.TFile( fn,'r' )

# retrieve the ntuple of interest
        tt = ROOT.gDirectory.Get( ttName )
        LB = tt.GetListOfBranches()
        bn = []
        for b in LB:
            bn.append( b.GetName() )

        tdict = {}
        for b in bn: tdict[b] = []
        
        entries = tt.GetEntriesFast()
        tt.LoadTree(0)

        # special treatment needed for string/character variables.
        isAChar = {}
        for b in bn:
            l = tt.GetLeaf(b)
            if l.GetTypeName()=='Char_t':
                isAChar[b] = bytearray(40)
                tt.SetBranchAddress(b,isAChar[b])


        
        for jentry in xrange( entries ):
 # get the next tree in the chain and verify
            ientry = tt.LoadTree( jentry )
            if ientry < 0:
                break

        # copy next entry into memory and verify
        # take care with character/string and conversion of integers
            nb = tt.GetEntry( jentry )
            if nb <= 0:
                continue
            for b in bn:
                l = tt.GetLeaf(b)
                if b in isAChar:
                    v = str(isAChar[b]).strip().replace('\x00','')
                else:
                    v = l.GetValue()
                    if l.GetTypeName()=='Int_t': v = int(v)  
                if  0 and b=='sample':
                    if jentry%100==0:
                        print 'jentry,b,v,l.GetTypeName()',jentry,b,v,l.GetTypeName()
                tdict[b].append(v)
            #print 'j,evts,evtsN,sample,sn',jentry,sWD.evts,sWD.evtsN,sWD.sample,sWD.sn

        myfile.Close()
        print 'makeTTree.getTTree filled',len(tdict[bn[0]]),'entries for',len(bn),'variables from TTree',ttName,'from file',fn
        return tdict

        
        
 
