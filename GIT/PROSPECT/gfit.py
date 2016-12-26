#!/usr/bin/env python
'''
gaussian fitter. Copied from NSRL_run13A 20161221
comment out langaus usage
'''

import ROOT
from array import array
import math

class gfit():
    def __init__(self):
        #ROOT.gROOT.ProcessLine('.L langaus.C++')
        #self.lgfun = ROOT.langaufun
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
            
    def twoG(self,hname,debug=False,x1min=-100.,x1max=100.,x2max=300.,sgm=None,Unconstrained=False): # two gaussian fit
        name,ave,rms,tot,under,bin1,over,ent,xmi,xma = self.getProp(hname)

        # check if fit should be performed
        par = {}
        gsumstat = 'FAILED'
        par['gsum'] = [gsumstat]
        if tot<=0 or float(ent)/float(tot)<0.50:
            return par['gsum']

        # fit options: Restricted range, Likelihood, Quiet
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
        
    def fit(self,hname,nsigplus=1.5,nsigminus=2.0,quiet=False,debug=False):
        name,ave,rms,tot,under,bin1,over,ent,xmi,xma = self.getProp(hname)

        GoodFit = False

        # reject hists with no entries or with too few entries
        # in the meaningful bins (>1 and <Nbin)
        if tot<=0 or float(ent)/float(tot)<0.50 :
            return GoodFit, ave,rms, -1., -1.

        if debug :
            print 'name,tot,under,bin1,over,ent',name,tot,under,bin1,over,ent

        fit_options = "RL" # Restricted range,Likelihood,
        if quiet: fit_options += "Q" # quiet
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
