#!/usr/bin/env python
'''
calculations for AD1 spike
20170815
'''
import math
import sys
#import random
import numpy
#import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import cutsAndConstants

class spike():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/spike/'
        self.cAC = cutsAndConstants.cutsAndConstants()
        return
    def main(self,makePDF=False,goalUnc=0.005,rateUnc=0.05):
        '''
        goalUnc = 0.005 # relative mass msmt uncertainty in AD1 spike
        if goalUnc<0, then minimize total relative uncertainty in spike

        rateUnc = expected uncertainty in measured rate of activity remaining in vial in Bq
        Default is 0.05Bq estimated for P50X spike LiLS#9

        
        C = activity concentration of spiked LiLS stock (=B2)
        m1 = mass of LiLS in vial V1 before spike
        m2 = mass of spike from B2
        m3 = mass dispensed from vial after spike
        AR = activity removed from V1 to be added to AD1
        MR = mass remaining in V1 after removal = m1 + m2 - m3
        sM = uncertainty in a single mass measurement
        Desire AR = 1.8 Bq, MR = 10 g
          '''
        AR = 1.8 # Bq
        sM = 0.010 # g uncertainty in single measurement

        coincRate  = 3.39 # Bq doc-1704 Danielle's P50X spike remainder measurement
        dcoincRate = rateUnc # Bq
        coincUnc = dcoincRate/coincRate # estimated relative uncertainty
        precrate = '{0:.2f} Bq'.format(dcoincRate)

        precmg = '{0:.1f} mg'.format(sM*1000.)
        permille = u"\u2030"
        
        Cjune2017 = 9.356 # activity concentration of B2 in Bq/g in June 2017
        Ac227tau = self.cAC.Ac227lifetime # in years
        dt = 0.5 # activity concentration in Dec 2017
        C = Cjune2017*math.exp(-0.5/Ac227tau)
        print 'Activity concentration in stock spiked LiLS',C,'Bq/g'
        mgmi = 500
        dmg  =  50/2
        mgma =2500+dmg+1 
        xmin = float(mgmi)/1000.
        xmax = float(mgma)/1000.
        ymin = 0.
        ymax = 20.+0.
        
        plt.clf()
        plt.grid()
        plt.xlabel('m2 mass of spike from stock added to vial (g)')
        plt.ylabel('mass (g) or activity (Bq) or '+permille+' uncertainty')
        name = 'AD1 {0:.2f} Bq spike Assume {1} mass msmt and {2} rate msmt unc'.format(AR,precmg,precrate)
        plt.title(name)
        axes = plt.gca()
        markersize = 4
        axes.set_xlim([xmin,xmax])
        axes.set_xticks(numpy.arange(xmin,xmax,0.1))
        plt.yticks(numpy.arange(ymin,ymax+1.0,1.0))
        axes.set_yticks(numpy.arange(ymin,ymax+1.0,0.25),minor=True)
        axes.set_ylim([ymin,ymax])
        
        GOAL = None
        for MR in [10.]: # grams
            line = '-'
            if MR>10.: line = ':'

            '''
            AL = activity (Bq) remaining in vial after AD1 spike
            Utot = total relative uncertainty in AD1 spike
                adding relative uncertainties of mass measurement
                and P50X coincidence rate measurement in quadrature
            Uest = estimated relative uncertainty in AD1 spike
                adding relative uncertainty of mass measurement
                and relative uncertainty due to uncertainty in P50X coincidence rate measurement and
                activity remaining in vial
                
            '''
            mlist = [x/1000. for x in range(mgmi,mgma,dmg)]
            M1,M2,M3,MTOT,CheckM,UAR,AL,Utot,Uest = [],[],[],[],[],[],[],[],[]
            for m2 in mlist:
                Q = m2 - AR/C 
                m1 = (MR - Q)*m2/Q
                m3 = m1 + m2 - MR
                mtot = m1 + m2
                sAR = (mtot + m3)/(mtot - m3) * sM / m3 # relative uncertainty in activity dispensed
                Aleft = C*m2/(m1+m2)*(MR) # activity remaining in vial in Bq
                if m1>0 and m3>0:
                    #print 'm1,m2,m3',m1,m2,m3,'sAR',sAR
                    M1.append(m1)
                    M2.append(m2)
                    M3.append(m3)
                    MTOT.append(mtot)
                    CheckM.append(m1+m2-m3)
                    UAR.append( 1000.*sAR )
                    Utot.append( 1000.*math.sqrt(sAR*sAR + coincUnc*coincUnc)) # total rel. unc.
                    Uest.append( 1000.*math.sqrt(sAR*sAR + dcoincRate*dcoincRate/Aleft/Aleft) )# est. rel. un.
                    AL.append(Aleft)
                    if GOAL is None: GOAL = [sAR, m1,m2,m3,Aleft,Uest[-1]/1000.]
                    if goalUnc>0:
                        if abs(GOAL[0]-goalUnc) > abs(sAR-goalUnc): GOAL = [sAR, m1,m2,m3,Aleft,Uest[-1]/1000.]
                    else:
                        if Uest[-1]/1000.<GOAL[5]: GOAL = [sAR, m1,m2,m3,Aleft,Uest[-1]/1000.]


            X,Y = numpy.array(M2),numpy.array(M3)    
            plt.plot(X,Y,'bo'+line,label='m3 mass removed from vial',markersize=markersize)

            Y = numpy.array(M1)
            plt.plot(X,Y,'gx'+line,label='m1 mass in vial before spike',markersize=markersize)

            Y = numpy.array(MTOT)
            plt.plot(X,Y,'ko'+line,label='Total mass in vial after spike',markersize=markersize)


            Y = numpy.array(CheckM)
            plt.plot(X,Y,'r'+line,label='mass remaining in vial after dispense')

            Y = numpy.array(UAR)
            plt.plot(X,Y,'cD'+line,label= permille+' uncertainty(mass msmt) in AD1 spike activity',markersize=markersize)
            #Y = numpy.array(Utot)
            #plt.plot(X,Y,'c'+line,label= permille+' total unc. in AD1 spike activity(est1)',markersize=markersize)

            Y = numpy.array(Uest)
            plt.plot(X,Y,'cd'+line,label= permille+' est. total unc. in AD1 spike activity',markersize=markersize)
            

            Y = numpy.array(AL)
            plt.plot(X,Y,'m'+line,label='Activity remaining in vial(Bq) ')
            
            #plt.plot(X,X,'b-',label='m3=m2')
        plt.plot([GOAL[2],GOAL[2]],[ymin,ymax],'r--',label='Operating point') # vertical dotted line at operating point
        plt.legend(prop={'size':13-3},loc=9) #loc=2=upper left
        
        # report calculated spike, dispense mass
        words = 'Spike mass {2:.3f} g, Dispensed mass {3:.3f} g, Dispensed activity rel. unc. {0:.3f}(mass only), Dispensed activity rel. unc. {5:.3f}(est.total), Activity left in vial {4:.1f} Bq, Initial mass in vial {1:.2f} g'.format(GOAL[0],GOAL[1],GOAL[2],GOAL[3],GOAL[4],GOAL[5])
        print words
        x1 = 1.3+.01
        y1 = 6.5
        dy = 0.9
        for word in words.split(','):
            plt.text(x1,y1,word,size='medium',fontname='serif')
            y1 -= dy


            
        if makePDF:
            if goalUnc>0:
                name += ' {0:.3f} goal mass unc'.format(goalUnc)
            else:
                name += ' minimize total unc'
            pdf = self.figdir + name.replace(' ','_').replace('.','x').replace(',','') + '.pdf'
            plt.savefig(pdf)
            print 'spike.main Wrote',pdf
        else:
            plt.show()

                
        return
  
if __name__ == '__main__' :
    pdf = False
    goalUnc = 0.005
    rateUnc = 0.05 
    if len(sys.argv)>1:
        w = sys.argv[1].lower()
        pdf = (w=='print' or w=='draw' or w=='pdf')
    if len(sys.argv)>2: goalUnc = float(sys.argv[2])
    if len(sys.argv)>3: rateUnc = float(sys.argv[3])
    A = spike()
    A.main(makePDF=pdf,goalUnc=goalUnc,rateUnc=rateUnc)
