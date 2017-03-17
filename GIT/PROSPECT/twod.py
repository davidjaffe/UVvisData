#!/usr/bin/env python
'''
fit various slices and projections of a twod histogram
20170316
'''
import math
import sys
from scipy.integrate import quad
import numpy

import datetime,os

import ROOT
import graphUtils
import gfit


class twod():
    def __init__(self):
        self.h = None
        self.f = None # needed to keep histograms in scope
        self.gU = graphUtils.graphUtils()

        self.sr2pi = math.sqrt(2.*math.pi)
        
        
        return
    def defineFF(self,w,xlo=None,xhi=None):
        '''
        define fit functions given bin width w within limits xlo,xhi
        '''
        # fitting functions: 2.5066282746310002 = sqrt(2pi)
        sr2pi = str(math.sqrt(2.*math.pi))
        sw = str(w)
        sGG = "([0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]/"+sr2pi+"+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])/[5]/"+sr2pi+")*"+sw
        if xlo is None or xhi is None:
            GG = self.GG = ROOT.TF1("GG",sGG)
        else:
            GG = self.GG = ROOT.TF1("GG",sGG,xlo,xhi)
        GG.SetParName(0,"gN")
        GG.SetParName(1,"gMean")
        GG.SetParName(2,"gSigma")
        GG.SetParName(3,"nN")
        GG.SetParName(4,"nMean")
        GG.SetParName(5,"nSigma")
        self.sG = sG = "([0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]/"+sr2pi+")*"+sw
        G = self.G = ROOT.TF1("G",sG)
        G.SetParName(0,"N")
        G.SetParName(1,"Mean")
        G.SetParName(2,"Sigma")
        
    def getHist(self,fn,hn):
        '''
        return hist with name hn from file named fn
        '''
        self.f = f = ROOT.TFile.Open(fn,'r')
        self.h = h = f.Get(hn)
        print 'twod.getHist: fn',fn,'hn',hn
        return h
    def anal(self,h):
        '''
        analysis of 2d hist
        work with a copy that gets deleted
        '''
        #print 'twod.anal h',h

        hname = h.GetName()
        newname = 'new_'+hname
        hnew = h.Clone(newname)
        nx,ny = hnew.GetNbinsX(),hnew.GetNbinsY()
        wx,wy = hnew.GetXaxis().GetBinWidth(1),hnew.GetYaxis().GetBinWidth(1)

        upper = [ny,140,120,100]
        lower = [10,12,14,16,18,20]
        
        c1 = ROOT.TCanvas("c1")
        c1.SetGrid()
        c1.SetTicks()

        h.Draw("colz")

        ymi,yma = h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(ny)
        

        line = {}
        for x1 in numpy.concatenate((upper,lower)):
            X1 = h.GetXaxis().GetBinLowEdge(x1)
            #print 'x1,ymi,yma,X1',x1,ymi,yma,X1
            line[x1] = ROOT.TLine(X1,ymi,X1,yma)
            line[x1].SetLineColor(ROOT.kRed)
            line[x1].Draw("same")
        
        c1.Modified()
        c1.Update()
        wait = raw_input()
        c1.Clear()
        if wait!='' : return
        
        ncv = 2
        fo = ["QMR","QMLR"]
        xcut = [[0.2,1.],[0.2,1.]]
        c1.Divide(1,ncv)
        c1.SetGrid()
        c1.SetTicks()
        ROOT.gStyle.SetOptFit(1111)
        ROOT.gStyle.SetOptStat(1)
        
        gN,gMean,gSigma, nN,nMean,nSigma = 10000.,.25,.02, 10000.,.4,.02

        nv = len(upper)*len(lower)*ncv
        iv = 0
        v,dv,x,dx = [],[],[],[]

        WAIT = True#False
        
        x0,x1 = 1,ny # first,last bin to include in projection
        for x1 in upper:
            sx1 = '_'+str(x1)
            for x0 in lower:
                sx0 = '_'+str(x0)

                for i in range(ncv):
                    fitopt = fo[i]
                    hnew_py = hnew.ProjectionY(newname+sx0+sx1+fitopt,x0,x1)

                    c1.cd(i+1)                
                    c1.SetLogy()
                    self.defineFF(wy,xcut[i][0],xcut[i][1])
                    GG = self.GG
                    GG.SetParameters(gN, gMean, gSigma, nN, nMean, nSigma)

                    
                    hnew_py.Fit(GG,fitopt) #####
                    P = []
                    for i in range(6):
                        P.append( GG.GetParameter(i) )
                        
                    PSDcut = self.setPSDcut(hnew_py)
                    n,dn = self.getNevt(hnew_py,P,PSDcut)
                    
                    print sx0+sx1+fitopt,'PSDcut',PSDcut,'nN',n,'+-',dn

                    iv += 1
                    v.append( n )
                    dv.append( dn )
                    x.append( float(iv) )
                    dx.append( 0. )
                    c1.SetLogy()
                    hnew_py.Draw("same")
                    c1.SetLogy()
                    gN,gMean,gSigma, nN,nMean,nSigma = P
                #print 'x0,x1',x0,x1,'gN,gMean,gSigma, nN,nMean,nSigma',gN,gMean,gSigma, nN,nMean,nSigma
                c1.Modified()
                c1.Update()
                if WAIT : wait = raw_input()

        c1.Clear()
        c1.Divide(1,2)
        c1.cd(1)
        ROOT.gStyle.SetOptStat(11111111)
        dv = numpy.array(dv)
        v = numpy.array(v)
        ev = max(max(5.*dv),max(v/100.))
        vlo = min(v-ev)
        vhi = max(v+ev)
        hh = self.gU.makeTH1D(v,'Yield','Yield',nx=40,xmi=vlo,xma=vhi)
        hh.Draw()

        c1.cd(2)
        g = self.gU.makeTGraph(x,v,'yield','yield',ex=dx,ey=dv)
        self.gU.color(g,2,2,setMarkerColor=True)
        g.Draw("AP")
        
        w = sum(1./dv/dv)
        m = sum(v/dv/dv)/w
        s = math.sqrt(1./w)
        rms = math.sqrt( sum( (v-m)*(v-m) )/(float(len(v))-1.) )
        ae = sum(dv)/float(len(dv))
        print 'Weighted mean {0:.1f}({1:.1f}) sqrt(rms) {2:.1f} average unc. {3:.1f}'.format(m,s,rms,ae)
        #print 'len(v),len(dv),len(x),len(dx)',len(v),len(dv),len(x),len(dx)
        #print 'v,dv',[[q,dq] for q,dq in zip(v,dv)]

        L = {}
        for j in [-1,0,1]:
            L[j] = ROOT.TLine(min(x),m+float(j)*s,max(x),m+float(j)*s)
            L[j].SetLineColor(ROOT.kRed)
            if j!=0 : L[j].SetLineStyle(2)
            L[j].Draw("same")
        LL = {}
        for j in [-1,1]:
            LL[j] = ROOT.TLine(min(x),m+float(j)*rms,max(x),m+float(j)*rms)
            LL[j].SetLineColor(ROOT.kBlue)
            LL[j].SetLineStyle(2)
            LL[j].Draw("same")
        
        c1.Modified()
        c1.Update()
        wait = raw_input()

        
                
        del hnew
        return
    def makeG(self,h,P,suffix='_gA'):
        '''
        return TF1 object of gaussian with parameters P nominally a fit to histogram h
        '''
        name = h.GetName() + suffix
        xmi = h.GetXaxis().GetXmin()
        xma = h.GetXaxis().GetXmax()
        gA = ROOT.TF1(name,self.sG,xmi,xma)
        for i in range(3):
            gA.SetParameter(i,P[i])
        return gA
    def getNevt(self,h,P,PSDcut):
        '''
        return estimate number of events above PSDcut
        given fitted histogram h and double gaussian parameters P
        
        '''

        gN,gMean,gSigma, nN,nMean,nSigma = P


        
        # estimate area and fractional area of n-recoil peak above cut from gaussian fit
        G1 = lambda x: gN*math.exp(-0.5*((x-gMean)*(x-gMean)/gSigma/gSigma))/gSigma/self.sr2pi
        G2 = lambda x: nN*math.exp(-0.5*((x-nMean)*(x-nMean)/nSigma/nSigma))/nSigma/self.sr2pi
        num = quad(G2,PSDcut,1.)  #num[0]=area, num[1]=estimate of unc. in area
        den = quad(G2,0.,1.)
        effy,erry = -1.,-1.
        if den[0]>0:
            effy = num[0]/den[0]
            if num[1]>0: erry = num[1]*num[1]/num[0]/num[0]
            erry = effy*math.sqrt(erry + den[1]*den[1]/den[0]/den[0])

        A1 = num[0] # effy corrected area


        # estimate area of n-recoil peak by subtracting gaussian fit of e-recoil peak
        ga = self.makeG(h,P[:3])
        hd = h.Clone('cloned_'+h.GetName())
        hd.Add(ga,-1.)
        A2,dA2 = self.getIntegral(hd,PSDcut,1.)

        # average of two areas is estimate
        A = 0.5*(A1+A2)
        dA= 0.5*abs(A1-A2)
        dA= A/effy*math.sqrt( dA*dA/A/A + erry*erry/effy/effy + dA2*dA2/A2/A2)
        A = A/effy
        #print 'A1,A2,dA2,effy,erry',A1,A2,dA2,effy,erry,'A,dA',A,dA
        return A,dA


        
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
    def setPSDcut(self,h,debug=False):
        '''
        return PSDcut given histogram h fitted with double gaussian GG
        '''

        PSDcut = None
        ff = None
        for a in h.GetListOfFunctions():
            if 'TF1' in a.ClassName(): ff = a
        if ff is None:
            print 'twod.setPSDcut WARNING PSD cut NOT updated. Could not find function in hist',hname
            return PSDcut
        if debug: print 'twod.setPSDcut hist',h.GetName(),'function',ff.GetName(),'Npar',ff.GetNpar()

        if ff.GetName()=='GG' and ff.GetNpar()==6:
            xmi,xma = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
            xmi,xma = ff.GetParameter(1),ff.GetParameter(1+3) # fitted means of gaussians
            nx = 100
            dx = (xma-xmi)/float(nx)
            xlo,FX = xmi,ff.Eval(xmi)
            for i in range(nx):
                x = xmi + (float(i)+0.5)*dx
                fx = ff.Eval(x)
                if i%10==0 and debug : print 'twod.setPSDcut scan x',x,'func(x)',fx
                if fx<FX : FX,xlo = fx,x
            PSDcut = xlo
            if debug: print 'twod.setPSDcut minimum of',FX,'at updated PSDcut',xlo
                
        return PSDcut
    def main(self,fn):
        hn = "PSD_vs_Charge"
        h = self.getHist(fn,hn)
        #print 'twod.main h',h,'self.h',self.h
        self.anal(h)

        self.f.Close()
        return

if __name__ == '__main__' :
    T = twod()
    fn = 'Speedy/Output/run00491.root' # lils#2
    T.main(fn)
    fn = 'Speedy/Output/run00660.root' # lils#8
    T.main(fn)
