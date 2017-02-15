#!/usr/bin/env python
import ROOT
from array import array
import sys


if 0:
    from ROOT import TFile, TTree
    from array import array


    f = TFile( 'test.root', 'recreate' )
    t = TTree( 't1', 'tree with histos' )

    maxn = 10
    TD = {}
    TD['n'] = array( 'i', [ 0 ] )
    TD['d'] = array( 'f', maxn*[ 0. ] )
    t.Branch( 'n', TD['n'], 'n/I' )
    t.Branch( 'd', TD['d'], 'd[n]/F' )

    for i in range(25):
       TD['n'][0] = min(i,maxn)
       for j in range(TD['n'][0]):
          TD['d'][j] = i*0.1+j
       t.Fill()

    f.Write()
    f.Close()
    sys.exit(' >>>>>>>>>>>>>>> AINT NO MORE >>>>>>>>>>')



import makeTTree
mtt = makeTTree.makeTTree()
if 0:
    name = ROOT.vector('string')()
    run  = ROOT.vector('int')()
    x    = ROOT.vector('double')()
    name.push_back('x1')
    name.push_back('x2')
    name.push_back('x3')
    run.push_back(1)
    run.push_back(2)
    run.push_back(3)
    x.push_back(1.1)
    x.push_back(2.2)
    x.push_back(3.3)
if 1:
    name,run,x = [],[],[]
    for i in range(5):
        name.append('ABC'+ ''.join(['x' for q in range(i)]))
        run.append(i)
        x.append(float(i)+float(i)/9.)
        
    print 'name',name,'\nrun',run,'\nx',x,'\n'

C = {'run': run, 'x': x, 'name': name}
mtt.makeTTree(C,fn='TTREE.ROOT',treename='tree')
sys.exit('that is all')

