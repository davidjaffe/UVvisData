#!/usr/bin/env python
'''
use solution to Bateman equation to understand how breaking of secular
equilibrium will affect adsorption studies.
20161114

add toy MC

'''
import math
import sys
import random
import numpy
import scipy
from scipy.stats.mstats import chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pickle

import cutsAndConstants

class bateman():
    def __init__(self):

        self.cAC = cutsAndConstants.cutsAndConstants()

        
        self.h = 60.*60.
        self.d = 24.*self.h
        self.y = 365.25*self.d
        self.isotopes = ['227Ac','227Th','223Ra','219Rn']
        self.colors   = {'227Ac':'ro','227Th':'bo','223Ra':'gD','219Rn':'co'}
        self.halflives= [21.772*self.y, 18.68*self.d, 11.43*self.d, 3.96]
        self.lifetimes= [math.log(2.)*t for t in self.halflives]
        self.initialAmounts = [100.*self.lifetimes[0], 0., 0., 0.]

        self.draw = False
        
        self.figdir = 'Figures/bateman/'
        return
    def bateman(self,t,initialAmounts,lifetimes):
        '''
        calculate amount of isotopes given time t and initial amounts of isotopes
        use notation of Bateman
        https://archive.org/details/cbarchive_122715_solutionofasystemofdifferentia1843
        '''

        P0,Q0,R0,S0 = initialAmounts
        l1,l2,l3,l4 = [1./x for x in lifetimes]

        P = P0*math.exp(-l1*t)

        Q = l1/(l2-l1) * P0 * math.exp(-l1*t)
        Q+= (l1*P0/(l1-l2) + Q0) * math.exp(-l2*t)

        R = l1*l2*P0/(l2-l1)/(l3-l1)*math.exp(-l1*t)
        R+= (l1*l2*P0/(l1-l2)/(l3-l2) + l2*Q0/(l3-l2))*math.exp(-l2*t)
        R+= (l1*l2*P0/(l1-l3)/(l2-l3) + l2*Q0/(l2-l3) + R0)*math.exp(-l3*t)

        S = l1*l2*l3*P0/(l2-l1)/(l3-l1)/(l4-l1)*math.exp(-l1*t)
        S+= (l1*l2*l3*P0/(l1-l2)/(l3-l2)/(l4-l2) + l2*l3*Q0/(l3-l2)/(l4-l3))*math.exp(-l2*t)
        S+= (l1*l2*l3*P0/(l1-l3)/(l2-l3)/(l4-l3) + l2*l3*Q0/(l2-l3)/(l4-l3) + l3*R0/(l4-l3))*math.exp(-l3*t)
        S+= (l1*l2*l3*P0/(l1-l4)/(l2-l4)/(l3-l4) + l2*l3*Q0/(l2-l4)/(l3-l4) + l3*R0/(l3-l4) + S0)*math.exp(-l4*t)
        return P,Q,R,S
    def toyMC(self,startTime=None):
        '''
        simulate a bunch of decays given conditions at equilibrium
        '''
        if startTime is None: startTime = 2.*self.y
        Amounts = self.bateman(startTime,self.initialAmounts,self.lifetimes)
        Am = {}
        Tot= 0.
        i = 0
        for iso,q in zip(self.isotopes,Amounts):
            s = q/self.lifetimes[i]
            Am[iso] =s
            Tot += s
            i += 1
        Cum = {}
        c = 0.
        for iso in self.isotopes:
            Am[iso] = Am[iso]/Tot
            c += Am[iso]
            Cum[iso] = c
        print 'bateman.toyMC startTime(y)',startTime/self.y,'Total',Tot,'Fractions/isotope$cumulative',
        for iso in self.isotopes:
            print '{0} {1:.4f}${2:.4f}'.format(iso,Am[iso],Cum[iso]),
        print ''

        Ac227Rate = 100. # Hz
        Results = {}
        for Ac227Rate in [0.1, 0.5, 1., 2., 3., 4., 10., 20.]: # Hz
            Results[Ac227Rate] = []
            Ralpha = 4. * Ac227Rate

            N = 25000*4
            t = 0.
            Events = []
            Icount = {}
            for iso in self.isotopes: Icount[iso] = 0
            iso = '215Po'
            Icount[iso] = 0
            while N>0:
                r = random.random()
                for j,iso in enumerate(self.isotopes):
                    if r<Cum[iso]:
                        break
                tdecay = random.expovariate(Ralpha)
                t += tdecay
                N -= 1
                Events.append( (iso,t) )
                Icount[iso] += 1
                if iso=='219Rn': # generate 215Po decay
                    iso = '215Po'
                    Icount[iso] += 1
                    tdecay = random.expovariate(1./self.cAC.Po215lifetime)
                    Events.append( (iso,t+tdecay) )
                    N -= 1
            print '\nbateman.toyMC Produced',len(Events),'events. By isotope',Icount,'Sorting...'
            Events = sorted(Events,key=lambda x: x[1]) # sort by time
            print 'bateman.toyMC Events sorted. Look for disorder'

            oldt = -1.
            oldpair = None
            for pair in Events:
                iso,t = pair
                disorder = ''
                if t<oldt: disorder = 'EVENT OUT OF ORDER'
                if disorder!='': print iso,t,disorder,'oldpair',oldpair
                oldt = t
                oldpair = pair

            tLast = Events[-1][1]
            alphaRate = float(len(Events))/tLast
            print 'bateman.toyMC N(events)',len(Events),'Ac227Rate(Hz)',Ac227Rate,'alphaRate(Hz)',alphaRate,'tLast(s)',tLast,\
              'Look for coincidences'

            for life in [1,2,3,4,5]:
                nLife = float(life)

                signal,bkgd = 0,0


                tWindow = nLife*self.cAC.Po215lifetime 
                for ievt,pair in enumerate(Events):
                    iP,tP = pair
                    if tP+tWindow<=tLast:  # avoid end effect
                        delayed = []
                        if ievt+1<len(Events):
                            for devt,dpair in enumerate(Events[ievt+1:]):
                                if dpair[1]-tP<tWindow:
                                    if devt==0 and dpair[0]=='215Po' and iP=='219Rn':
                                        signal += 1
                                    else:
                                        bkgd += 1
                                        delayed.append(dpair)
                                else:
                                    break
                        if 0: print 'ievt',ievt,'iP,tP',iP,tP,'delayed',delayed


                U = signal + 2.*bkgd
                print 'bateman.toyMC nLife',nLife,'Nsignal',signal,'Nbkgd',bkgd,'U',U,'sqrt(U)/signal',math.sqrt(U)/signal
                Results[Ac227Rate].append( {nLife: [signal,bkgd]} )

        fn = 'pickled_bateman.pickle'
        f = open(fn,'w')
        pickle(Results,f)
        f.close()
        print 'bateman.toyMC wrote',fn
                
        return 
        
    def main(self,amounts=None,name='bateman',Duration=None,specIsotope='All'):
        '''
        main routine
        calculates and plots time evolution of isotopes given initial amounts
        plot to name.pdf
        returns final amounts of isotopes (when t=Duration)
        if specIsotope == 'All' then plot all isotopes
        if specIsotope == isotope name, then plot only named isotope
        '''
        #print 'amounts',amounts
        initialAmounts = self.initialAmounts
        if amounts is not None: initialAmounts = amounts
        #print 'bateman.main',name,'initialAmounts',initialAmounts
        if Duration is None: Duration = self.y/2.
        if type(Duration) is list:
            T = numpy.linspace(Duration[0],Duration[1])
        else:
            T = numpy.linspace(0.,Duration)
        Results = {}
        finalResults = None
        for A in self.isotopes:
            Results[A] = []
        for t in T:
            results  = self.bateman(t,initialAmounts,self.lifetimes)
            if t==T[-1]: finalResults = results
            for A,r in zip(self.isotopes,results):
                Results[A].append(r)
                

        t = T/self.d
        plt.clf() # clear figure needed
        special = None
        for A,color in zip(self.isotopes,['ro','bo','gD','co']):
            b = self.lifetimes[self.isotopes.index(A)]/self.d
            y = numpy.array([a/b for a in Results[A] ])
            if specIsotope=='All' or specIsotope==A:
                #print A,t/self.d,y
                plt.plot(t,y,color,label=A)
                if specIsotope==A: special = [t,y]
        plt.xlabel('Time (days)')
        plt.ylabel('Amount (decays/day)')
        plt.title(name)
        plt.grid()
        plt.legend()#bbox_to_anchor=(0.63, 1.02-.02, 1.-.6, .102), loc=3, ncol=2, mode="expand", borderaxespad=0., numpoints=1)
        if self.draw:
            pdf = self.figdir +name.replace(' ','_')+'.pdf'
            plt.savefig(pdf)
            print 'bateman.main Wrote',pdf
        else:
            plt.show()
        print 'bateman.main Complete, special isotope',specIsotope
        return finalResults, special
    def standard(self,daysDuration=7.,fracAdsorbed=1.0):
        '''
        standard treatment, calculations and plots
        show results over a duration in days of daysDuration
        
        '''
        cpctAdsorb = ' {0:.1f}% '.format(100.*fracAdsorbed)
        print 'bateman.standard daysDuration',daysDuration,cpctAdsorb,'adsorption'
        duration = daysDuration*self.d
        A,specialTY = self.main(specIsotope='All')
        name = 'At equilibrium with no adsorption'
        newA,NoAdsorptionTY = self.main(amounts=A,name=name,Duration=duration,specIsotope='219Rn')
        #print 'newA',newA
        #print 'NoAdsorptionTY',NoAdsorptionTY 
        for specIso in ['All', '219Rn']:
            sI = '_'+specIso
            specialTY = {}
            results,specialTY['StartAtEquil'] = self.main(amounts=A,name='start at equilibrium'+sI,Duration=duration,specIsotope=specIso)

            for i,iso in enumerate(self.isotopes):
                if iso != '219Rn':
                    Z = numpy.array([min(1.,float(j!=i)-fracAdsorbed+1.) for j in range(len(self.isotopes))])
                    name = 'At equilibrium then ' + iso + cpctAdsorb + ' adsorbed' + sI
                    print 'iso,Z',iso,Z,name
                    results,specialTY[iso] = self.main(amounts=A*Z, name=name, Duration=duration, specIsotope=specIso)
                elif 0:

                    fby2 = fracAdsorbed/2.
                    Z = numpy.array( [ 1.-fby2, 1.-fby2, 1., 1.] )
                    c = ' {0:.1f}% '.format(100*fby2)
                    name = 'At equilibrium then '+self.isotopes[0]+ c + 'adsorbed and '+self.isotopes[1]+c+'adsorbed'+sI
                    print 'iso,Z',iso,Z,name
                    R,S = self.main(amounts=A*Z, name=name, Duration=duration, specIsotope=specIso) # makes a plot, but returned data unused
            #print specialTY

            if any( specialTY.values() ): # true if specIso!='All'
                tryFit = True # try to model rate relative to no adsorption as exponent?
                
                print 'bateman.standard specIso is',specIso
                ymax = -1.
                for mode in ['absolute','normed','ratioToNo']:
                    nymi = 1.
                    plt.clf()
                    tNo,yNo = NoAdsorptionTY
                    for iso in specialTY:
                        if iso in self.isotopes:
                            t,y = specialTY[iso]
                            if mode=='absolute':
                                ymax=max(max(y),ymax)
                                #print 'iso,ymax',iso,ymax
                            if mode=='normed'  :
                                y = y/ymax
                            if mode=='ratioToNo':
                                y = y/yNo
                            nymi = min(min(y),nymi)
                            plt.plot(t,y,self.colors[iso],label=iso)
                            if mode=='ratioToNo' and iso=='227Ac':
                                A = max(y)-min(y)
                                B = max(y)
                                tau = 22.5
                                for pw,sym in zip([1.5],['g']): #zip([1.4,1.5,1.6],['k','g','b']):
                                    tp = numpy.power(t/tau,pw)
                                    f = A*(numpy.exp(-tp)-1.)+B
                                    plt.plot(t,f,sym+'-',label='exp(-(t/'+str(tau)+'days)^'+str(pw)+')')

                            
                    plt.xlabel('Time (days)')
                    if mode=='absolute':
                        plt.ylabel('219Rn Rate (decays/day)')
                    elif mode=='normed':
                        plt.ylabel('219Rn relative rate')
                        plt.ylim( (max(0.,nymi-.01),1.01) )
                    elif mode=='ratioToNo':
                        plt.ylabel('219Rn relative to no adsorption')
                        plt.ylim( (max(0.,nymi-.001),1.001) )
                    plt.title('Single isotope'+cpctAdsorb+'adsorption',loc='left')
                    plt.grid()
                    plt.legend()#bbox_to_anchor=(0.6,1.,1.-.6,.1), ncol=2, numpoints=1)
                    if self.draw:
                        pdf = self.figdir +mode +  '_standard' + cpctAdsorb.replace(' ','_').replace('%','percent') + '.pdf'
                        plt.savefig(pdf)
                        print 'batemean.standard Wrote',pdf
                    else:
                        plt.show()

                
                        
        print 'bateman.standard Complete'
        return

    
  
if __name__ == '__main__' :
    B = bateman()
    if 0: # failed attempt to optimize spike rate
        B.toyMC()
    if 0:
        B.main(Duration=[00.*B.d,350.*B.d])
    if 1:
        daysDuration = 100.
        obsF = 1.-72./93.
        for fracAdsorbed in [obsF,0.,0.001,0.01,0.15,0.3,1.]:
            B = bateman()
            B.draw = True
        
            B.standard(daysDuration=daysDuration,fracAdsorbed=fracAdsorbed)
