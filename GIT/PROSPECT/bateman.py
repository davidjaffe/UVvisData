#!/usr/bin/env python
'''
use solution to Bateman equation to understand how breaking of secular
equilibrium will affect adsorption studies.
20161114
'''
import math
import sys
import random
import numpy
import scipy
from scipy.stats.mstats import chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class bateman():
    def __init__(self):

        self.h = 60.*60.
        self.d = 24.*self.h
        self.y = 365.25*self.d
        self.isotopes = ['227Ac','227Th','223Ra','219Rn']
        self.colors   = {'227Ac':'ro','227Th':'bo','223Ra':'gD','219Rn':'co'}
        self.halflives= [21.772*self.y, 18.68*self.d, 11.43*self.d, 3.96]
        self.lifetimes= [math.log(2.)*t for t in self.halflives]
        self.initialAmounts = [100.*self.lifetimes[0], 0., 0., 0.]

        self.figdir = 'Figures/'
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
        plt.legend(bbox_to_anchor=(0.63, 1.02-.02, 1.-.6, .102), loc=3, ncol=2, mode="expand", borderaxespad=0., numpoints=1)

        pdf = self.figdir +name.replace(' ','_')+'.pdf'
        plt.savefig(pdf)
        print 'bateman.main Wrote',pdf
        return finalResults, special
    def standard(self):
        '''
        standard treatment, calculations and plots
        '''
        A,specialTY = self.main(specIsotope='All')
        for specIso in ['All', '219Rn']:
            sI = '_'+specIso
            specialTY = {}
            duration = 7.*self.d
            results,specialTY['StartAtEquil'] = self.main(amounts=A,name='start at equilibrium'+sI,Duration=duration,specIsotope=specIso)

            for i,iso in enumerate(self.isotopes):
                if iso != '219Rn':
                    Z = numpy.array([float(j!=i) for j in range(len(self.isotopes))])
                    name = 'At equilibrium then ' + iso + ' adsorbed' + sI
                    results,specialTY[iso] = self.main(amounts=A*Z, name=name, Duration=duration, specIsotope=specIso)
            #print specialTY

        plt.clf()
        for iso in specialTY:
            if iso in self.isotopes:
                t,y = specialTY[iso]
                plt.plot(t,y,self.colors[iso],label=iso)
        plt.xlabel('Time')
        plt.ylabel('219Rn Rate (decays/day)')
        plt.title('Single isotope disappearance',loc='left')
        plt.grid()
        plt.legend(bbox_to_anchor=(0.6,1.,1.-.6,.1), ncol=2, numpoints=1)
        pdf = self.figdir + 'standard' + '.pdf'
        plt.savefig(pdf)
        print 'batemean.standard Wrote',pdf
        return

    
  
if __name__ == '__main__' :
    B = bateman()
    B.standard()
