#!/usr/bin/env python
'''
projection of daya bay sensitivity
20190523

starts with estimate from Chao, Kam-Biu for US portfolio review in 2018

PRL.121.241805 (1958 days) t13 0.0856(0.0029), dm2ee 2.522(+0.068,-0.070)e-3 eV^2

'''
import math
import sys
import numpy
import matplotlib.pyplot as plt
import datetime

class systs():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/'

        self.estimate = {'t13': [ [201507, 0.0027, 0.0019, 0.0033], # yearmonth, stat, syst, sum
                                  [201709, 0.0022, 0.0017, 0.0028],
                                  [202012, 0.0016, 0.0014, 0.0021]],
                         'dm2ee' : [ [201507, 0.058, 0.062, 0.085],
                                     [201709, 0.047, 0.055, 0.072],
                                     [202012, 0.035, 0.046, 0.058]]
                                     }
        self.syst2017 = {}
        for par in self.estimate:
            A = self.estimate[par]
            for B in A:
                dt,st,sy,tot = B
                if dt==201709: self.syst2017[par] = sy
        self.current = {'t13': 0.0841, 'dm2ee':2.50 } # dm2ee in eV^2

        self.PRL1958 = {'t13': [0.0856, 0.0029], 'dm2ee':[2.522, 0.070] }
        self.PRLdate = '201709' # 26Jan2017 + 217 days

        
        return
    def plot(self):
        '''
        
        '''
        #plt.rcdefaults() # why?


        dlim = [datetime.datetime.strptime(d,'%Y%m').date() for d in ['201412','202206']]
        
        pars = ['t13', 'dm2ee']

        for card in ['absolute','relative']:
            wpar = ['$sin^22{\Theta}_{13}$','$\Delta m_{ee}^2$ eV${}^2$']
            if card=='relative': wpar = ['$sin^22{\Theta}_{13}$','$\Delta m_{ee}^2$']

        
            plt.clf()
            fig,ax = plt.subplots(2,1)

            fig.suptitle('Parameter ' + card + ' uncertainties')
        
            for i,par in enumerate(pars):
                cv = self.current[par] # central value of 'current' estimate
                dates, stat, syst, tot = [],[],[],[]
                denom = 1.
                if card=='relative': denom = self.current[par]
                for a in self.estimate[par]:
                    dates.append( a[0] )
                    stat.append(  a[1]/denom )
                    syst.append(  a[2]/denom )
                    tot.append(   a[3]/denom )
                ax[i].set_xlabel('Date of last data')
                ax[i].set_xlim(dlim)
                ax[i].grid()
                X = numpy.array(dates)
                DT = X = [datetime.datetime.strptime(str(d),'%Y%m').date() for d in dates]
                ST = Y = numpy.array(stat)
                ax[i].plot(X,Y,'o-',label='Stat')
                SY = Y = numpy.array(syst)
                ax[i].plot(X,Y,'s-',label='Syst')
                Y = numpy.sqrt( ST*ST + SY*SY)
    #            print 'systs.plot check total for',par,Y
                ax[i].plot(X,Y,'D-',label='Total')
                ax[i].set_ylabel(wpar[i])
                
                # compare to the PRL result (1958 days)
                X = [datetime.datetime.strptime(self.PRLdate,'%Y%m').date()]
                Y = self.PRL1958[par][1]
                if card=='relative': Y = Y/self.PRL1958[par][0]
                ax[i].plot(X,Y,marker='*',color='black',label='PRL 1958days')

                # project using 2017 systematics
                X = DT[-1]
                sy = self.syst2017[par]
                Y = numpy.sqrt( ST[-1]*ST[-1] + sy*sy/denom/denom)
                ax[i].plot(X,Y,marker="v",color='maroon',label='with 2017 syst')
                
                if i==0:
                    b,t = ax[0].get_ylim()
                    ax[0].set_ylim(b,1.2*t)
                    ax[0].legend(loc='upper right',fontsize='small',numpoints=1)

        
            pdf = self.figdir + 'sensitivity_'+card+'.pdf'
            plt.savefig(pdf)
            print 'systs.plot Wrote',pdf
        return
      
if __name__ == '__main__' :
    S = systs()
    S.plot()

