#!/usr/bin/env python
'''
gaussian fitter. Copied from NSRL_run13A 20161221
comment out langaus usage
'''

import ROOT
from array import array
import math
import os
import graphUtils,cutsAndConstants
import numpy

class gfit():
    def __init__(self):
        #ROOT.gROOT.ProcessLine('.L langaus.C++')
        #self.lgfun = ROOT.langaufun
        self.gU = graphUtils.graphUtils()
        self.cAC= cutsAndConstants.cutsAndConstants()
        print 'gfit initialized'
        return

    def getProp(self,hname):
        name = hname.GetName()
        ave = hname.GetMean()
        rms = hname.GetRMS()
        tot = hname.GetEntries()
        under = hname.GetBinContent(0)
        bin1 = hname.GetBinContent(1)
        over  = hname.GetBinContent(hname.GetXaxis().GetNbins()+1)
        ent = tot - under - over 
        xmi = hname.GetXaxis().GetXmin()
        xma = hname.GetXaxis().GetXmax()
        return name,ave,rms,tot,under,bin1,over,ent,xmi,xma

    def getIntegral(self,hname,xmi,xma):
        nbins = hname.GetXaxis().GetNbins()
        i1 = max( 1,hname.GetBinLowEdge(xmi) )
        i2 = min(nbins, hname.GetBinLowEdge(xma)+1 )
        return hname.Integral(i1,i2)
        

    def empAds(self,t,p): # no self?
        '''
        example of external function for TF1 from /Users/djaffe/work/GIT/WATER/scint_time_dist.py

        empirical function for 227Ac adsorption
        t in seconds, t0=p[3]
        f1 = A*(exp(- ((t-t0)/tauA)^1.5 ) -1.) + B
        f = f1*exp(-(t-t0)/life)
        tauA = 22.5days
        life = 22.7/ln(2) years 
        '''
        #print 't',t,'p',p
        oneday = 60.*60.*24.
        t0 = self.T0
        T = float(t[0])-t0
        #print 'T',T,'t0',t0
        T = max(T,0.)/oneday
        A = float(p[0])
        B = float(p[1])
        tau = float(p[2])
        Q = math.pow(T/(22.5),1.5)
        f1 = A*(math.exp(-Q) - 1.) + B
        f  = f1*math.exp(-T/tau)
        return f
    def expF(self,g,pspre=None):
        
        '''
        bunch of fits to time series.
        At least one is an exponential fit of form y = p0 * exp((x-x0)/p1)
        where x0 is the first point
        '''
        x,y = self.gU.getPoints(g)
        print 'min(y),max(y)',min(y),max(y)
        oneday = 60.*60.*24.
        x0 = min(x)
        sx0=str(x0)
        soneday = str(oneday)
        sEXP = "[0]*exp(-(x-"+sx0+")/[1]/"+soneday+")"
        sEXPA= "[0]*(exp(-(x-"+sx0+")/[1]/"+soneday+"))+[2]"

        sP1  = "[1]*(x-"+sx0+")+[0]"
        
        npar = 2

        c1 = self.gU.makeCanvas()

        fEXP = ROOT.TF1('fEXP',sEXP)
        fEXP.SetParName(0,'R0')
        fEXP.SetParName(1,'TauDays')
        fEXP.SetParameter(0,max(y))
        tau = self.cAC.Ac227lifetime * 365.
        fEXP.SetParameter(1,tau)
        print '227Ac lifetime',tau,'days'

        fEXPA = ROOT.TF1('fEXPA',sEXPA)
        fEXPA.SetParName(0,'R0')
        fEXPA.SetParName(1,'TauDays')
        fEXPA.SetParName(2,'Asymptote')
        fEXPA.SetParameter(0,max(y)-min(y))
        fEXPA.SetParameter(1,100.)
        fEXPA.SetParameter(2,min(y))

        self.T0 = x0
        fEMA = ROOT.TF1("empAds",self.empAds,min(x),max(x),3)
        fEMA.SetParName(0,'R0')
        fEMA.SetParName(1,'Asymptote')
        fEMA.SetParName(2,'TauDays')
        fEMA.SetParameter(0,max(y)-min(y))
        fEMA.SetParameter(1,min(y))
        fEMA.SetParameter(2,self.cAC.Ac227lifetime*365.)


        
        
        fP1 = ROOT.TF1('fP1',sP1)
        fP1.SetParName(0,'R0')
        fP1.SetParName(1,'Slope')
        fP1.SetParameter(0,max(y))
        fP1.SetParameter(1,-1./1000.)

        fEXP.FixParameter(1,tau)
        g.Draw("AP")
        g.Fit(fEXP)
        c1.Update()
        if pspre is not None:
            PS = pspre+'Exp_FixTau.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:    
            raw_input('CR to continue')

        fEXP.ReleaseParameter(1)
        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit(fEXP)
        c1.Update()
        if pspre is not None:
            PS = pspre+'Exp.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')

        fEMA.FixParameter(2, self.cAC.Ac227lifetime*365.)
        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit(fEMA)
        g.Fit(fEMA,"M")
        g.Fit(fEMA,"ME")
        c1.Update()
        if pspre is not None:
            PS = pspre+'EmpiricalAdsorption_FixTau.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')

        fEMA.ReleaseParameter(2)
        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit(fEMA)
        g.Fit(fEMA,"M")
        g.Fit(fEMA,"ME")
        c1.Update()
        if pspre is not None:
            PS = pspre+'EmpiricalAdsorption.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')

        

        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit(fEXPA)
        g.Fit(fEXPA,"M")
        g.Fit(fEXPA,"ME")
        c1.Update()
        if pspre is not None:
            PS = pspre+'ExpA.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')

        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit("pol0","F")
        c1.Update()
        if pspre is not None:
            PS = pspre+'Constant.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')
        
        c1 = self.gU.makeCanvas()
        g.Draw("AP")
        g.Fit(fP1,"F")
        #g.Fit(fP1,"F")
        c1.Update()
        if pspre is not None:
            PS = pspre+'Linear.ps'
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        else:
            raw_input('CR to continue')
        
        return
        
    def fourG(self,h,WAIT=False,quietly=True,debug=False,PS=None,DrawFitBox=True):
        '''
        Return area, unc of nuclear recoil peak in PSD distribution in input hist h
        Four gaussian fit to PSD distribution.
        Initial parameters based on Danielle's work.
        '''
        qopt = ''
        if quietly: qopt = 'Q'  # quiet, no output to terminal
        if not debug : qopt += '0' # don't draw
        noPopUp = PS is not None
        if noPopUp : ROOT.gROOT.ProcessLine("gROOT->SetBatch()")

        qopt2 = ''
        if quietly and PS is not None:
            qopt2 = 'Q'

        if debug: print 'gfit.fourG h',h

        w = h.GetXaxis().GetBinWidth(1) # assume equal size bins
        sw = '*'+str(w)
        sr2pi = str(math.sqrt(2.*math.pi))
        sGG = "([0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]/"+sr2pi+  \
            "+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])/[5]/"+sr2pi+  \
            "+[6]*([9]*exp(-0.5*(x-[7])*(x-[7])/[8]/[8])/[8]/"+sr2pi+  \
            "+(1.-[9])*exp(-0.5*(x-[10])*(x-[10])/[11]/[11])/[11]/"+sr2pi+"))"+sw
        if debug: print 'gfit.fourG sGG',sGG

        sG = "([0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]/"+sr2pi+")"+sw

        npar = 12
        xlo,xhi = 0.15,0.55

        c1 = ROOT.TCanvas('c1')
        G4 = ROOT.TF1('G4',sGG)
        ROOT.gStyle.SetOptStat(0)
        if DrawFitBox:
            ROOT.gStyle.SetOptFit(1111)
            ROOT.gStyle.SetTitleX(0.8)
            ROOT.gStyle.SetStatW(0.5) # size of stats box (and text?)
            ROOT.gStyle.SetStatH(0.5)
        else:
            ROOT.gStyle.SetOptFit(0)
        c1.SetGrid(1)
        c1.SetTicks(1)

        # fix and set means, sigmas to best guesses 
        k = 0
        means = [.2629,.2523,.4261,.4165]
        meanLimits = [ [0.2,0.3], [0.2,0.3], [0.3,0.5], [0.3,0.5] ]
        sigmas= [.0222,.0470,.0237,.0310]
        
        for i in [1,4,7,10]:
            j = i+1
            sk = str(k)
            G4.SetParName(i,'Mean'+sk)
            G4.SetParameter(i,means[k])
            G4.FixParameter(i,means[k])
            G4.SetParName(j,'Sigma'+sk)
            G4.SetParameter(j,sigmas[k])
            G4.FixParameter(j,sigmas[k])
            k+=1
        # set amplitudes, including fixed fraction
        G4.SetParName(0,'A0')
        G4.SetParameter(0,100./w)
        G4.SetParName(3,'A1')
        G4.SetParameter(3,100./w)
        iN = 6
        G4.SetParName(iN,'N')
        G4.SetParameter(iN,400./w)
        G4.SetParName(9,'f')
        G4.SetParameter(9,0.4)
        G4.FixParameter(9,0.4)

        fitopt = ''+qopt

        # 3par fit: amplitudes only
        h.Fit(G4,fitopt,"",xlo,xhi)
        c1.Update()
        c1.SetLogy() ###########
        ps = h.GetListOfFunctions().FindObject("stats")
        if ps:
            ps.SetX1NDC(0.11)
            ps.SetX2NDC(0.51)
            ps.SetY1NDC(0.40-0.03)
            ps.SetY2NDC(1.00-0.03)
            ps.SetTextSize(0.01*2*1.5)
        c1.Update()
        if WAIT: wait = raw_input()

        # 4par fit: release fraction
        G4.ReleaseParameter(9)
        G4.SetParameter(9,0.4)
        G4.SetParLimits(9,0.,1.)
        self.setPositive(G4,npar)
        h.Fit(G4,fitopt,"",xlo,xhi)
        c1.Update()
        if WAIT: wait = raw_input()

        # if fraction is at limit, reset it
        f = G4.GetParameter(9)
        if abs(f-0.5)<0.01: f = 0.4
        G4.SetParameter(9,f)
    
            

        # 12par fits: release means, sigmas
        k = 0
        for i in [1,4,7,10]:
            j = i+1
            sk = str(k)
            G4.ReleaseParameter(i)
            G4.SetParameter(i,means[k])
            G4.SetParLimits(i,meanLimits[k][0],meanLimits[k][1])
            G4.ReleaseParameter(j)
            G4.SetParameter(j,sigmas[k])
            k += 1

        N,dN = None,None
        for fitopt in [''+qopt,'I'+qopt,'IM'+qopt,'IME'+qopt2]:

            self.setPositive(G4,npar)
            h.Fit(G4,fitopt,"",xlo,xhi)
            c1.Update()
            ps = h.GetListOfFunctions().FindObject("stats")
            if ps:
                ps.SetX1NDC(0.11)
                ps.SetX2NDC(0.51)
                ps.SetY1NDC(0.40-0.03)
                ps.SetY2NDC(1.00-0.03)
                ps.SetTextSize(0.01*2*1.5)
            c1.Update()
            N = G4.GetParameter(iN)
            dN= G4.GetParError(iN)
            if debug: print 'N {0:.2f}({1:.2f})   dN/N {2:.4f}  1/sqrt(N) {3:.4f}'.format(N,dN,dN/N,1./math.sqrt(N))

            if WAIT: wait = raw_input()

        # draw the 4 gaussians
        par = self.getPar(G4,npar)
        F = {}
        for i in range(4):
            if i<2:
                p = par[i*3:i*3+3]
            elif i==2:
                p = [ par[iN]*par[9], par[i*3+1], par[i*3+2] ]
            else:
                p = [ par[iN]*(1.-par[9]), par[i*3+1], par[i*3+2] ]

            sf = 'f'+str(i)
            F[sf] = ROOT.TF1(sf,sG)
            for j,x in enumerate(p): F[sf].SetParameter(j,x)
            F[sf].SetLineColor(i+2)
            F[sf].Draw("SAME")
            
        c1.Update()
        if PS is not None:
            pdf = PS.replace('.ps','.pdf')
            self.gU.finishDraw(c1,PS,pdf,ctitle=pdf)
        if debug : wait = raw_input('all done. cr to continue')
        
        return N,dN
    def getPar(self,F,npar):
        p = []
        for i in range(npar): p.append( F.GetParameter(i) )
        return p
    def setPositive(self,F,npar):
        for i in range(npar):
            p = F.GetParameter(i)
            F.SetParameter(i,abs(p))
        return
    

         
            
    def twoG(self,hname,debug=False,x1min=-100.,x1max=100.,x2max=300.,sgm=None,Unconstrained=False): # two gaussian fit
        name,ave,rms,tot,under,bin1,over,ent,xmi,xma = self.getProp(hname)

        # check if fit should be performed
        par = {}
        gsumstat = 'FAILED'
        par['gsum'] = [gsumstat]
        if tot<=0 or float(ent)/float(tot)<0.50:
            return par['gsum']

        # fit options: Restricted range, Likelihood, Quiet, do not plot result
        fitopt = "RLQ"
        if debug : fitopt = "RL"


        # define ranges, initial means, sigmas of 2 gaussians
        # this is setup for fitting SPE in NSRL13A run
        x1mean, x1sig =  0.5*(x1min+x1max), 30.
        if sgm is not None: x1sig = sgm
        x2min = x1max
        x2mean, x2sig = 0.5*(x2min+x2max), 50.
        if sgm is not None: x2sig = sgm
        g1 = ROOT.TF1("g1","gaus",x1min,x1max)
        g1.SetParameter(1,x1mean)
        g1.SetParameter(2,x1sig)
        g2 = ROOT.TF1("g2","gaus",x2min,x2max)
        g2.SetParameter(1,x2mean)
        g2.SetParameter(2,x2sig)

        # fit first gaussian. retry once if failed. save parameters
        if debug: print 'gfit.twoG fit g1'
        hname.Fit("g1",fitopt)
        g1stat = ROOT.gMinuit.fCstatu
        if g1stat!='CONVERGED ': # note trailing blank
            if debug : print 'gfit.twoG fit g1, 2d try'
            g1.SetParameter(1,x1mean)
            g1.SetParameter(2,x1sig)
            hname.Fit("g1",fitopt)
            g1stat = ROOT.gMinuit.fCstatu
        c,mean,emean,sgm,prob = self.getg1Par(g1,offset=0)
        par['g1'] = [c,mean,emean,sgm,prob,g1stat]
        if debug: print 'gfit.twoG g1',par['g1']

        # fit second gaussian. retry once if failed. save parameters
        if debug: print 'gfit.twoG fit g2'
        hname.Fit("g2",fitopt)
        g2stat = ROOT.gMinuit.fCstatu
        if g2stat!='CONVERGED ': # note trailing blank
            if debug: print 'gfit.twoG fit g2, 2d try'
            g2.SetParameter(1,x2mean)
            g2.SetParameter(2,x2sig)
            hname.Fit("g2",fitopt)
            g2stat = ROOT.gMinuit.fCstatu
        c,mean,emean,sgm,prob = self.getg1Par(g2,offset=0)
        par['g2'] = [c,mean,emean,sgm,prob,g2stat]
        if debug: print 'gfit.twoG g2',par['g2']

        # fit with both gaussians
        gsum = ROOT.TF1("gsum","gaus(0)+gaus(3)",x1min,x2max)
        gsum.SetParameter(0,par['g1'][0])
        gsum.SetParameter(1,par['g1'][1])
        gsum.SetParameter(2,par['g1'][3])
        gsum.SetParameter(0+3,par['g2'][0])
        gsum.SetParameter(1+3,par['g2'][1])
        gsum.SetParameter(2+3,par['g2'][3])

        
        if debug:
            dbgpar = gsum.GetParameters()
            print 'gfit.twoG fit gsum. initial pars',dbgpar
        hname.Fit("gsum",fitopt)
        gsumRstat = ROOT.gMinuit.fCstatu
        if gsumRstat!='CONVERGED ': # note trailing blank
            if 'CONVERGED' in par['g1'][-1]:
                gsum.SetParameter(0,par['g1'][0])
                gsum.SetParameter(1,par['g1'][1])
                gsum.SetParameter(2,par['g1'][3])
            else:
                gsum.SetParameter(1,x1mean)
                gsum.SetParameter(2,x1sig)
                
            if 'CONVERGED' in par['g2'][-1]:
                gsum.SetParameter(0+3,par['g2'][0])
                gsum.SetParameter(1+3,par['g2'][1])
                gsum.SetParameter(2+3,par['g2'][3])
            else:
                gsum.SetParameter(1+3,x2mean)
                gsum.SetParameter(2+3,x2sig)
            if debug: print 'gfit.twoG gsum, 2d try'
            hname.Fit("gsum",fitopt)
            gsumRstat = ROOT.gMinuit.fCstatu
        c1,mean1,emean1,sgm1,prob = self.getg1Par(gsum,offset=0)
        c2,mean2,emean2,sgm2,prob = self.getg1Par(gsum,offset=3)
        par['gsumR'] = [c1,mean1,emean1,sgm1, c2, mean2,emean2,sgm2, prob, gsumRstat]
        if debug: print 'gfit.twoG gsum fitopt,par',fitopt,par['gsumR']

        #Unconstrained = False
        if Unconstrained : 
            fitopt = fitopt.replace('R','') # no constraints
            if debug: print 'gfit.twoG do unconstrained gsum fit'
            hname.Fit("gsum",fitopt)
            gsumstat = ROOT.gMinuit.fCstatu
            if gsumstat!='CONVERGED ': # note trailing blank
                hname.Fit("gsum",fitopt)
                gsumstat = ROOT.gMinuit.fCstatu
            c1,mean1,emean1,sgm1,prob = self.getg1Par(gsum,offset=0)
            c2,mean2,emean2,sgm2,prob = self.getg1Par(gsum,offset=3)
            par['gsum'] = [c1,mean1,emean1,sgm1, c2, mean2,emean2,sgm2, prob, gsumstat]
        else:
            par['gsum'] = par['gsumR']
            if debug:
                print 'gfit.twoG No unconstrained double-gaussian fit done'
                print 'gfit.twoG gsum',fitopt,par['gsum']

        # make two gaussians fcns TF1 corresponding to the 2 gaussians of the best fit
        ga,gb = self.makeBothG(hname,'gsum')

        #g1.Delete()
        #g2.Delete()
        #gsum.Delete()
        return par['gsum'],ga,gb

    def makeBothG(self,hname,fname):
        # define the two gaussians that were fitted to hname
        f = hname.GetFunction(fname)
        name = hname.GetName()+"_gA"
        xmi = hname.GetXaxis().GetXmin()
        xma = hname.GetXaxis().GetXmax()
        gA = ROOT.TF1(name,"gaus",xmi,xma)
        for i in range(3):
            gA.SetParameter(i,f.GetParameter(i))

        name = hname.GetName()+"_gB"
        gB = ROOT.TF1(name,"gaus",xmi,xma)
        for i in range(3):
            gB.SetParameter(i,f.GetParameter(i+3))

        return gA,gB
        
    def fit(self,hname,nsigplus=1.5,nsigminus=2.0,quiet=False,debug=False,drawFit=False):
        name,ave,rms,tot,under,bin1,over,ent,xmi,xma = self.getProp(hname)

        GoodFit = False

        # reject hists with no entries or with too few entries
        # in the meaningful bins (>1 and <Nbin)
        if tot<=0 or float(ent)/float(tot)<0.50 :
            return GoodFit, ave,rms, -1., -1.

        if debug :
            print 'name,tot,under,bin1,over,ent',name,tot,under,bin1,over,ent

        fit_options = "RL0" # Restricted range,Likelihood,do not draw fit
        if quiet: fit_options += "Q" # quiet
        if drawFit: fit_options.replace("0","")
        nsig = 3.
        mean = ave
        sgm  = rms
        for nsig in [3., 2., -1]:
            xlo = max(xmi,mean - nsig*sgm)
            xhi = min(xma,mean + nsig*sgm)
            if nsig<0:
                xlo = max(xmi,mean - nsigminus*sgm)
                xhi = min(xma,mean + nsigplus*sgm)
            g1 = ROOT.TF1("g1","gaus",xlo,xhi)
            hname.Fit("g1",fit_options)
            c,mean,emean,sgm,prob = self.getg1Par(g1,offset=0)
            if debug : print 'name,nsig,mean,emean,sgm,prob',name,nsig,mean,emean,sgm,prob
            
        emean = g1.GetParError(1)
        GoodFit = True
        return GoodFit,mean,emean, sgm, prob
        
    def getg1Par(self,g1,offset=0):
        # print 'gfit.getg1Par g1,offset',g1,offset
        c = g1.GetParameter(0+offset)
        mean = g1.GetParameter(1+offset)
        sgm  = g1.GetParameter(2+offset)
        emean= g1.GetParError(1+offset)
        prob = g1.GetProb()
        return c,mean,emean,sgm,prob
    
    def lgfit(self,hname,xlo=-1.0):
        # fit hname with landau dist convolved with gaussian
        results = {'GoodFit':False}

        
        # require a minimum number of entries
        name,ave,rms,tot,under,bin1,over,ent,xmi,xma = self.getProp(hname)
        #print 'name,ave,rms,tot,ent,xmi,xma=',name,ave,rms,tot,ent,xmi,xma
        if tot<=10. or ent<=10.: return False,results

        # avoid fitting peak at zero
        if ave<xlo: return False, results

        # guess at parameters of function
        npar = 4
        s = math.sqrt(ave)
        w = s/5.
        startguess = array('d',[w , ave, ent, s])
        #print 'startguess',startguess

        lf = ROOT.TF1("lf",self.lgfun,xmi,xma,npar)
        lf.SetParameters(startguess)
        lf.SetParNames("Width","MostProb","Area","Sigma")

        hname.Fit(lf,"R") # restrict range, chi2 method
        hname.Fit(lf,"RL")# restrict range, likelihood method

        results = {}
        for i in range(npar):
            pname = lf.GetParName(i)
            results[pname] = [lf.GetParameter(i), lf.GetParError(i)]
        results['chi2'] = lf.GetChisquare()
        results['Prob'] = lf.GetProb()
        results['ndf']  = lf.GetNDF()
        results['status'] = ROOT.gMinuit.fCstatu
        results['GoodFit'] = True

        return True, results
    def analyze2D(self):

        fn = 'Speedy/Output/run00164.root'
        f = ROOT.TFile(fn)
        rn = os.path.basename(fn).replace('.root','')
        hists = ["PSD_vs_Charge","PSD_vs_ChargeN"]
        for hn in hists:
            h = f.Get(hn)
            newname = 'new_'+hn
            hnew = h.Clone(newname)
           
            hnew_py = hnew.ProjectionY(newname+'py',1,hnew.GetNbinsY())

            N,dN = self.fourG(hnew_py,WAIT=False,quietly=True)
            if N is not None:
                print '{0} {5} Nalpha {1:.2f}({2:.2f}) dN/Nalpha {3:.4f}  sqrt(1./N) {4:.4f}'.format(rn,N,dN,dN/N,math.sqrt(1./N),hn)

        return
    def analyzePSD(self):
        '''
        main routine for analyzing PSD distributions using 4gaussian fit
        '''


        fn = 'Figures/nine20170829/nine_20170901_111848.root'
        fn = 'Figures/nine20170829/nine_20170901_113548.root'
        fn = 'Figures/nine20170829/nine_20170901_144900.root'
        fn = 'Figures/nine20170829/nine_20170901_155149.root'
        files = [fn]
        figdir = 'Figures/nine20170829/'
        bhn = 'PSD_runXXXXXN'
        r1,r2 = 1000,2000

        Y,dY,T,dT = [],[],[],[]
        
        for fn in files:
            f = ROOT.TFile(fn)
            print 'gfit.analyzePSD Opened',fn
            for r in range(r1,r2):
                hn = bhn.replace('XXXXX','{0:05d}'.format(r))
                h = f.Get(hn)
                if h:
                    t = h.GetTitle()
                    s = t.split(' ')
                    Tstart = float((s[-2]).replace('Tstart=',''))/1000. # convert to seconds from milliseconds
                    Tend   = float((s[-1]).replace('Tend=',''))+Tstart # absolute time
                    ps = figdir + hn + '.ps'
                    N,dN = self.fourG(h,PS=ps,DrawFitBox=False) #None,debug=True,WAIT=True)
                    if N is not None:
                        print '{0} Nalpha {1:.2f}({2:.2f}) dN/Nalpha {3:.4f}  sqrt(1./N) {4:.4f}'.format(hn,N,dN,dN/N,math.sqrt(1./N))
                        Y.append(N)
                        dY.append(dN)
                        T.append((Tstart+Tend)/2.)
                        dT.append((Tend-Tstart)/2.)
            f.Close()

        graphs = []
        Y,dY,T,dT = numpy.array(Y),numpy.array(dY),numpy.array(T),numpy.array(dT)
        g = self.gU.makeTGraph(T,Y,'AlphaRate9','AlphaRate9',ex=dT,ey=dY)
        graphs.append(g)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figdir,abscissaIsTime=True,option='AP',yLimits=[16.0,19.0])
        
        g = self.gU.makeTGraph(T,Y,'AlphaRate9z','AlphaRate9z',ex=dT,ey=dY)
        graphs.append(g)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figdir,abscissaIsTime=True,option='AP',yLimits=[16.5,17.5])

        g = self.gU.makeTGraph(T,Y,'AlphaRate1','AlphaRate1',ex=dT,ey=dY)
        graphs.append(g)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figdir,abscissaIsTime=True,option='AP',yLimits=[-0.1,0.2])

        rateCut = 3. # Hz
        aY,adY,aT,adT = [],[],[],[]
        bY,bdY,bT,bdT = [],[],[],[]
        for i,y in enumerate(Y):
            if y>rateCut:
                aY.append(y)
                adY.append(dY[i])
                aT.append(T[i])
                adT.append(dT[i])
            else:
                bY.append(y)
                bdY.append(dY[i])
                bT.append(T[i])
                bdT.append(dT[i])
        aY,adY,aT,adT = numpy.array(aY),numpy.array(adY),numpy.array(aT),numpy.array(adT)
        g = self.gU.makeTGraph(aT,aY,'AlphaRateHigh','AlphaRateHigh',ex=adT,ey=adY)
        graphs.append(g)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figdir,abscissaIsTime=True,option='AP')
        
        bY,bdY,bT,bdT = numpy.array(bY),numpy.array(bdY),numpy.array(bT),numpy.array(bdT)
        g = self.gU.makeTGraph(bT,bY,'AlphaRateLow','AlphaRateLow',ex=bdT,ey=bdY)
        graphs.append(g)
        self.gU.color(g,2,2,setMarkerColor=True)
        self.gU.drawGraph(g,figdir,abscissaIsTime=True,option='AP')
        
                
        

        gfn = fn.replace('.','_graphs.')
        f = ROOT.TFile(gfn,'RECREATE')
        for g in graphs: f.WriteTObject(g)
        f.Close()
        print 'gfit.analyzePSD Wrote',len(graphs),'graphs to',gfn
        
                        
        return
    def OLDanalyzePSD(self):
        '''
        main routine for analyzing PSD distributions using 4gaussian fit
        '''
        fn = 'Danielle/Histograms_RootTrees/run00701_ts1490121118_Hists.root'


        files = ['Danielle/Histograms_RootTrees/run00697_ts1490118343_Hists.root',
            'Danielle/Histograms_RootTrees/run00698_ts1490118955_Hists.root',
            'Danielle/Histograms_RootTrees/run00699_ts1490119669_Hists.root',
            'Danielle/Histograms_RootTrees/run00700_ts1490120403_Hists.root',
            'Danielle/Histograms_RootTrees/run00701_ts1490121118_Hists.root',
            'Danielle/Histograms_RootTrees/run00702_ts1490121840_Hists.root',
            'Danielle/Histograms_RootTrees/run00703_ts1490122550_Hists.root',
            'Danielle/Histograms_RootTrees/run00704_ts1490123341_Hists.root']

        figdir = 'Danielle/Figures/'

        for fn in files:

            bn = os.path.basename(fn)
            rn = bn.split('_')[0]
            ps = figdir + rn + '.ps'
            f = ROOT.TFile(fn)
            h = f.Get('hPSD')
            N,dN = self.fourG(h,PS=ps)
            if N is not None:
                print '{0} Nalpha {1:.2f}({2:.2f}) dN/Nalpha {3:.4f}  sqrt(1./N) {4:.4f}'.format(rn,N,dN,dN/N,math.sqrt(1./N))
        return
    def fitRates(self):
        '''
        fits of rate vs time
        '''
        fn = 'Figures/nine20170829/nine_20170901_155149_graphs.root'
        pspre = 'Figures/nine20170829/nine_20170901_155149_graphs_fit_'
        f = ROOT.TFile(fn)
        for name in ['AlphaRateHigh', 'AlphaRateLow']:
            g = f.Get(name)
            if name=='AlphaRateLow':
                g.GetYaxis().SetLimits(-0.005,0.085)
                g.GetYaxis().SetRangeUser(-0.005,0.085)
            self.expF(g,pspre=pspre + name + '_')

        
        return
if __name__ == '__main__' :
    T = gfit()
    analyzePSD = False
    fitRates   = True
    if fitRates:   T.fitRates()
    if analyzePSD: T.analyzePSD()
    #T.analyze2D() # for Danielle's plots
    
