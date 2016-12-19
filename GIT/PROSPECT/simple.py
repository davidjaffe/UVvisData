#!/usr/bin/env python
'''
simple attempt to fit adsorption/evolution with exponential
20161201
'''
import math
import sys
#import random
import numpy
#import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class simple():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/simple/'
        return
    def fac10(self,t):
        ''' form for 10% Ac adsorption starting from equilibrium '''
        f = -9.46273*math.pow(10,-9)* math.exp(-0.175037*t) + 72.983 *math.exp(-1.22046*math.pow(10,-6)*t) - 368.71* math.exp(-7.01884*math.pow(10,-7)*t) + 395.174* math.exp(-4.29472*math.pow(10,-7)*t)
        return f
    def main(self):
        '''
        main routine

        '''
        x,y = [],[]
        for ix in range(100):
            t = 10000.*float(ix)
            f = self.fac10(t)
            x.append(t)
            y.append(f)
        X,Y = numpy.array(x),numpy.array(y)    
        plt.clf()
        plt.grid()
        name = 'simple'
        plt.title(name)
        ic,ip = 0,0
        plt.plot(X,Y,self.colors[ic]+self.points[ip])
        
        
        pdf = self.figdir + name + '.pdf'
        plt.savefig(pdf)
        print 'simple.main Wrote',pdf
 
        return

    
  
if __name__ == '__main__' :
    A = simple()
    A.main()
