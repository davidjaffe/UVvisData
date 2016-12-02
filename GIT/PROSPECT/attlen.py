#!/usr/bin/env python
'''
try to estimate change in attenuation length using results shown in doc1487
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

class attlen():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/attlen/'
        return
    def rate(self,x,lam,L):
        f = (math.exp(-x/lam) + math.exp(-(L-x)/lam))/lam/(1.-math.exp(-L/lam))
        return f
    def ave(self,lam,L):
        f = lam*(1.-math.exp(-L/lam))
        return f
    def relrate(self,lam1,lam2,L=1200.):
        '''
        estimate relative rate for two different measurements of response
        over length of the cell of length L with attenuation lengths of lam1, lam2
        Assumes that average response is normalized to same value
        '''

        lx,lr = [],[]
        for ix in range(0,int(L),10):
            x = float(ix)
            lx.append(x)
            f1 = self.rate(x,lam1,L)
            f2 = self.rate(x,lam2,L)
            r = f1/f2
            lr.append(r)

        X,Y = numpy.array(lx),numpy.array(lr)
           
        return X,Y
    
    def main(self):
        '''
        main routine

        '''
        L = 1200. # units are mm
        lam0 = 870.
        for il2 in [0]:
            lam2 = lam0 + 100.*float(il2)
            plt.clf()
            plt.grid()
            name = 'al_lam2_'+str(int(lam2))
            plt.title(name)
            ave2 = self.ave(lam2,L)
            for il1 in range(1,10):
                lam1 = lam2*(1.-float(il1)*0.01)
                X,Y = self.relrate(lam1,lam2,L=L)
                ic = il1%len(self.colors)
                ip = int(il1/len(self.points))
                #print 'ic,ip',ic,ip
                ave1 = self.ave(lam1,L)
                label = 'lam1={0} R={1:.2f}'.format(int(lam1),ave1/ave2)
                plt.plot(X,Y,self.colors[ic]+self.points[ip],label=label)
                print 'lam2',lam2,'lam1',lam1,'ave2',self.ave(lam2,L),'ave1',self.ave(lam1,L)
            pdf = self.figdir + name + '.pdf'
            plt.legend(loc="upper center")
            #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102/2), loc=3, ncol=2, mode="expand", borderaxespad=0.)
            plt.savefig(pdf)
            print 'attlen.main Wrote',pdf
 
        return

    
  
if __name__ == '__main__' :
    A = attlen()
    A.main()
