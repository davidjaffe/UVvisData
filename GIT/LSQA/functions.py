#!/usr/bin/env python
'''
functions for fitting, etc. 
20170127
'''
import sys
# the ONETON dir contains wfanal.py , graphUtils.py
sys.path.append('/Users/djaffe/work/GIT/ONETON')
# the PROSPECT dir contains get_filepaths.py

import numpy 
import matplotlib.pyplot as plt
import math
import ROOT
import os

import graphUtils

class functions():
    def __init__(self):
        self.gU = graphUtils.graphUtils()
        self.cGwE = ROOT.TF1("cGwE",self.convGwE,-10.,150.,3) # name,func,xmi,xma,Npar
        self.cGwE.SetParName(0,'Amp')
        self.cGwE.SetParName(1,'Lifetime')
        self.cGwE.SetParName(2,'Sigma')

        self.cGwEN= ROOT.TF1("cGwEN",self.convGwEN,-10.,250.,7)
        self.cGwEN.SetParName(0,'Amp')
        self.cGwEN.SetParName(1,'Lifetime0')
        self.cGwEN.SetParName(2,'Sigma')
        self.cGwEN.SetParName(3,'frac1')
        self.cGwEN.SetParLimits(3,0.,1.)
        self.cGwEN.SetParName(4,'Lifetime1')
        self.cGwEN.SetParName(5,'frac2')
        self.cGwEN.SetParLimits(5,0.,1.)
        self.cGwEN.SetParName(6,'Lifetime2')

        

        return
    def convGwEN(self,v,par):
        '''
        convolve N exponentials with gaussian
        '''
        A0   = par[0]
        life0 = par[1]
        sig  = par[2]
        f1   = par[3]
        life1 = par[4]
        f2   = par[5]
        life2 = par[6]
        t = v[0]
        tq= -t
        f0 = (1.-f1-f2)
        p = [A0*f0,life0,sig]
        F0 = self.convGwE(v,p)
        p = [A0*f1,life1,sig]
        F1 = self.convGwE(v,p)
        p = [A0*f2,life2,sig]
        F2 = self.convGwE(v,p)
        return F0+F1+F2
    
    def convGwE(self,v,par):
        '''
        convolve exponential with gaussian over domain (0,+infinity)
        http://math.stackexchange.com/questions/685041/convolute-exponential-with-a-gaussian/872561#872561
        '''
        #print 'functions.convGwE v,par',v,par
        A   = par[0]
        lam = 1./par[1]
        sig = par[2]
        t = v[0]
        #print 'functions.convGwE A,lam,sig,t',A,lam,sig,t
        #tq= -(t-sig*sig/lam) # from ref
        tq= -t  # put at origin
        try:
            F = lam/2. * math.exp(sig*sig*lam*lam/2.)*math.exp(-lam*t)
        except OverflowError,e:
            print 'lsqa.convGwE',e,'sig,lam,t',sig,lam,t
            F = 1.e100
        F = A*F*(1.-math.erf(tq/sig/math.sqrt(2)))
        return F
    def test1(self):
        '''
        draw exponential convolved with gaussian
        '''
        x = (numpy.arange(2000.)-200.)/10.
        lifet = 40.
        lam = 1./40.
        sig = 5.
        A = 1.
        par = [A, lifet,sig]
        y = []
        for v in x:
            y.append( self.convGwE([v],par) )
        y = numpy.array(y)
        plt.plot(x,y,'r.')
        plt.grid()
        plt.show()
        return
    def test4(self):
        '''
        draw 3 exponentials convolved with gaussian
        '''
        x = (numpy.arange(2000.)-200.)/10.
        life0 = 9.3
        life1 = 40.
        life2 = 400.
        A = .16
        f1 = 0.5
        f2 = 0.1
        sig = 1.1

        par = [A, life0,sig,f1,life1,f2,life2]
        y = []
        for v in x:
            y.append( self.convGwEN([v],par) )
        y = numpy.array(y)
        plt.plot(x,y,'r.')
        plt.grid()
        plt.show()
        return
    def test2(self):
        '''
        try to fit decay curve
        '''
        fn = 'Samples/LiLS01/run3566934507/Figures/optimize_run3566934507.root'
        canvas = ROOT.TCanvas("c1")
        f = ROOT.TFile(fn,'r')
        h = f.Get('neutron_vs_time_in_ns')
        FUNC = self.cGwE
        FUNC.SetParameters(.1,40.,.2)
        h.Fit(FUNC,"L","",-3.,20.)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(1111)
        h.Draw()
        canvas.SetLogy(1)
        canvas.SetGrid()
        canvas.Draw()
        canvas.Update()
        raw_input("test2 Waiting....Hit <CR> to continue")
        return
    def test3(self):
        '''
        try to fit decay curve
        iterative approach
        '''
        fn = 'Samples/LiLS01/run3566934507/Figures/optimize_run3566934507.root'
        fn = 'Samples/P50-1/run3568649727/Figures/optimize_run3568649727.root'
        canvas = ROOT.TCanvas("c1")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(1111)
        f = ROOT.TFile(fn,'r')
        horiginal = f.Get('neutron_vs_time_in_ns')
        h = horiginal.Clone('Clone')
        h = self.setHistContent(h,136.,146.) # change errors near blatant reflections
        h.SetMaximum(horiginal.GetMaximum()*1.05)
        xma = h.GetXaxis().GetXmax()
        xmi = -10.
        h.GetXaxis().SetRangeUser(xmi,xma)
        FUNC = self.cGwEN
        FUNC.SetParameters(.16,9.,1.,.5,40,.1,400.)

        fopt = "LB" #weighted bins, likelihood, parameters bounded

        xlo,xhi = -3.,20.
        FUNC.FixParameter(3,0.)
        FUNC.FixParameter(4,40.)
        FUNC.FixParameter(5,0.)
        FUNC.FixParameter(6,400.)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        Lifetime0 = FUNC.GetParameter(1)
        sig = FUNC.GetParameter(2)
        print 'test3 Lifetime0,sig',Lifetime0,sig
        self.drawAndWait(canvas)
        
        FUNC.FixParameter(2,sig)
        xhi += 5.*Lifetime0
        FUNC.SetParameter(3,.5)
        FUNC.SetParLimits(3,0.,1.)
        FUNC.SetParLimits(4,Lifetime0,10*Lifetime0)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        
        Lifetime0, Lifetime1 = FUNC.GetParameter(1),FUNC.GetParameter(4)
        print 'test3 Lifetime0,Lifetime1',Lifetime0,Lifetime1
        self.drawAndWait(canvas)
        
        FUNC.SetParameter(5,.1)
        FUNC.SetParLimits(5,0.,1.)
        FUNC.SetParLimits(6,Lifetime1,20*Lifetime1)
        xhi = min(xhi+10*Lifetime1,xma)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        self.drawAndWait(canvas)

        FUNC.SetParLimits(2,sig/2.,5*sig)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        self.drawAndWait(canvas)

        fopt += "M"
        FUNC.SetParLimits(2,sig/2.,5*sig)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        self.drawAndWait(canvas)

        fopt = ""
        FUNC.SetParLimits(2,sig/2.,5*sig)
        h.Fit(FUNC,fopt,"",xlo,xhi)
        h.Draw("hist")
        FUNC.Draw("same")
        self.drawAndWait(canvas)

        self.drawAndWait(canvas,SetLogy=False)
        a,b = ROOT.Double(0),ROOT.Double(0)
        print 'par# name low value hi'
        for i in range(FUNC.GetNpar()):
            FUNC.GetParLimits(i,a,b)
            print i,FUNC.GetParName(i),a,FUNC.GetParameter(i),b
        
        
        

        return
    def getHistContent(self,h,getErrors=False):
        '''
        return abscissa, ordinate and error on ordinate for input hist
        '''
        x,y,dy = [],[],[]
        for i in range(h.GetNbinsX()):
           y.append( h.GetBinContent(i+1) )
           x.append( h.GetBinCenter(i+1) )
           ey = 0.
           if getErrors: ey = h.GetBinError(i+1)
           dy.append(ey)
        return x,y,dy
    def blowupErrors(self,x,dy,xlo,xhi,factor=1.e20):
        '''
        increase errors for bins in range (xlo,xhi) by factor
        '''
        newdy = []
        for u,dv in zip(x,dy):
            if xlo<=u and u<=xhi:
                newdy.append( dv*factor)
            else:
                newdy.append( dv )
        return newdy
    def setHistContent(self,h,xlo,xhi):
        x,y,dy = self.getHistContent(h,getErrors=True)
        newdy  = self.blowupErrors(x,dy,xlo,xhi)
        for i,du in enumerate(newdy):
            h.SetBinError(i+1,du)
        return h
        
    def drawAndWait(self,canvas,SetLogy=True):
        ''' secret knowledge of how to draw hist and then wait for carriage return '''
        canvas.SetLogy(SetLogy)
        canvas.SetGrid()
        canvas.Draw()
        canvas.Update()
        raw_input("Waiting....Hit <CR> to continue")
        return
    def test5(self,fn = 'Samples/P50-1/run3568649727/Figures/optimize_run3568649727.root'):
        '''
        test drawMultiObjects
        '''
        name_hn = 'neutron_vs_time_in_ns'
        name_hg = 'gamma_vs_time_in_ns'
        f = ROOT.TFile.Open(fn,'r')
        #f.ls()
        bn = os.path.basename(fn).replace('.root','')
        hn = f.Get(name_hn)
        hg = f.Get(name_hg)
        objlist,Logy = [ [hn,hg] ], [True]

        #objlist.extend( [hn, hg] )     # temp for testing # of panels
        #Logy.extend( [False, False] ) # temp for testing # of panel

        objlist.append(f.Get('f100_t800'))
        Logy.append(False)
        objlist.append(f.Get('FOMerific_run3568649727'))
        Logy.append(False)
        objlist.append(f.Get('FOM_v_Q_f100_t800Nrun3568649727'))
        Logy.append(False)
        for i in range(5,65,10):
            name = 'PSD_Q'+str(i).zfill(2)+'_'+str(i+10)+'f100_t800N'
            objlist.append(f.Get(name))
            Logy.append( True)
        self.gU.drawMultiObjects(objlist,fname='test5_'+bn,statOpt=0,setLogy=Logy,abscissaIsTime=False,Grid=True,changeColors=True,addLegend=True,dopt='histfunc',gopt='AP',biggerLabels=False)
        
        return
    def test6(self,fn = 'Samples/P50-1/run3568649727/Figures/optimize_run3568649727.root'):
        '''
        test finding optimal PSD cut
        open fitted hist, extract function, determine position of valley between
        two peaks
        '''
        debug = True
        f = ROOT.TFile.Open(fn,'r')
        favFast = 100
        favTotal= 800
        Qlo,Qhi = 15.,25.
        hname = 'PSD_Q'+ str(int(Qlo)).zfill(2) + '_' + str(int(Qhi)).zfill(2)
        hname+= 'f' + str(int(favFast)) + '_t' + str(int(favTotal)) + 'N'
        h = f.Get(hname)
        ff = None
        for a in h.GetListOfFunctions():
            if 'TF1' in a.ClassName(): ff = a
        print 'hist',h.GetName(),'function',ff.GetName(),'Npar',ff.GetNpar()
        if ff.GetName()=='GG' and ff.GetNpar()==6:
            xmi,xma = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
            xmi,xma = ff.GetParameter(1),ff.GetParameter(1+3) # fitted means of gaussians
            nx = 100
            dx = (xma-xmi)/float(nx)
            xlo,FX = xmi,ff.Eval(xmi)
            for i in range(nx):
                x = xmi + (float(i)+0.5)*dx
                fx = ff.Eval(x)
                if i%10==0: print 'x',x,'func(x)',fx
                if fx<FX : FX,xlo = fx,x
            print 'minimum of',FX,'at',xlo
            text = ROOT.TText(xlo,FX,'MINIMUM')
            c = ROOT.TCanvas()
            h.Draw()
            h.Draw("same func")
            text.SetTextAngle(90.)
            text.Draw()
            self.drawAndWait(c,SetLogy=False)
                
        return
    def test7(self,Quiet=False):
        '''
        20170919
        fit compton edge ratio vs time(minutes) for vials 1-7 given canvases
        fit function is A1*exp(-t/tau1) + A2*exp(-t/tau2)
        '''
        
        for K in range(7):
            fn = '../../PROCESSED_LSQA_data/ComptonEdgeRatio_vs_Time_P'+str(K+1)+'.root'
            bn = os.path.basename(fn)
            print '\n Start fit of',bn

            f = ROOT.TFile(fn,'r')
            c1 = f.Get("c1");
            c1.Draw()
            if not Quiet: raw_input("original " + str(K+1))
                
            g =  c1.GetPrimitive("gComptonEdgeRatio")
            F = None
            for f in g.GetListOfFunctions():
                #print f,f.GetName()
                if f.GetName()!='stats': F = f

            Initial = {}
            for J in range(4):
                #print J,F.GetParameter(J)
                Initial[J] = [F.GetParameter(J), F.GetParError(J)]
            g.Draw()
            c2E = ROOT.TF1("c2E","[0]*exp(-x/[1])+[2]*exp(-x/([3]))")
            c2E.SetParameter(0,.1)
            c2E.SetParameter(1,500.)
            c2E.SetParameter(2,1.)
            c2E.SetParameter(3,1000.*60.*24.)
            ParName = {}
            for i,name in zip([0,1,2,3],['A1new','tau1new','A2new','tau2new']):
                c2E.SetParName(i,name)
                ParName[i] = name
            basic = ""
            if Quiet: basic = "Q"
            for opt in [basic,basic]:
                g.Fit(c2E,opt)
                if 0: raw_input("fit option is "+opt)
                
            c1.Draw()
            c1.Update()
            self.gU.drawMultiObjects([g],fname=str(K+1),figdir='TwoExp/',statOpt=0,fitOpt=1111,gopt='PA',biggerLabels=False)
            Final = {}
            for J in range(4):
                #print J,c2E.GetParameter(J)
                Final[J] = [c2E.GetParameter(J), c2E.GetParError(J)]
            print 'Compare initial and final fit parameters'
            print 'Parameter {0:>20} {1:>30} {2:>30}'.format('Initial','Final','Final-Initial')
            fmt = '{0} {1:>20.3e}({2:.3e}) {3:>20.3e}({4:.3e})  {5:>20.3e}({6:.3e})'
            for J in range(4):
                A,dA = Initial[J]
                B,dB = Final[J]
                print fmt.format(ParName[J],A,dA,B,dB,B-A,math.sqrt(dA*dA+dB*dB))
            A = Initial[0][0]+Initial[2][0]
            dA= math.sqrt(math.pow(Initial[0][1],2)+math.pow(Initial[2][1],2))
            B = Final[0][0]+Final[2][0]
            dB= math.sqrt(math.pow(Final[0][1],2)+math.pow(Final[2][1],2))
            Both = ParName[0]+'+'+ParName[2]
            print fmt.format(Both,A,dA,B,dB,B-A,math.sqrt(dA*dA+dB*dB))

            
            if not Quiet: raw_input("done for "+bn)
        return
        
if __name__ == '__main__':
    F = functions()
    #F.test1()
    #F.test4()
    #F.test2()
    #F.test3()
    #F.test5()
    #F.test6()
    F.test7(len(sys.argv)>1)
