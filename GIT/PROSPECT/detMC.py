#!/usr/bin/env python
'''
simulation to estimate expected spectrum from 227Ac
20161215
'''
import math
import sys
import random
import numpy
import scipy
from scipy.stats.mstats import chisquare
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import datetime,os

class detMC():
    def __init__(self,pe_MeV=10.):

        self.DautNames = ['227Ac','227Th','223Ra','219Rn','215Po','211Pb','211Bi','207Tl']
        self.colors    = ['b'    ,'g'    ,'r'    ,'c'    ,'m'    ,'y'    ,'k'    ,'p']

        self.DautAlpha = {'227Ac': [()],
                          '227Th': [(6038.,0.242), (6008.8,0.029), (5977.7,0.235), (5959.7,0.030),
                                   (5866.6,0.0127),(5756.87,0.204),(5713.2,0.0489),(5708.8,0.083),
                                   (5700.8,0.0363),(5694.,0.0150),(5668.0,0.0206)],
                          '223Ra': [(5871.3,0.010),(5747.,0.090),(5716.23,0.252),(5606.73,0.252),
                                    (5539.80,0.090),(5501.6,0.010),(5433.6,0.0222)],
                          '219Rn': [(6819.1,0.794),(6552.6,0.129),(6425.0,0.075)],
                          '215Po': [(7386.1,1.0)],
                          '211Pb': [()],
                          '211Bi': [(6622.9,0.8354),(6278.2,0.1619)],
                          '207Tl': [()]
                                   }
        for daut in self.DautNames:
            self.DautAlpha[daut] = self.read(daut)

        self.pePerMeVa = pe_MeV #
        print 'detMC.__init__ estimated number of photoelectrons per MeV of alpha energy=',pe_MeV
        
        self.figdir = 'Figures/detMC/'
        makeSubDir = False
        if makeSubDir:
            now = datetime.datetime.now()
            fmt = '%Y%m%d_%H%M%S_%f'
            self.start_time = cnow = now.strftime(fmt)
            self.figdir = self.figdir + cnow+'/'
            if os.path.isdir(self.figdir):
                pass
            else:
                try:
                    os.makedirs(self.figdir)
                except IOError,e:
                    print 'detMC__init__',e
                else:
                    print 'detMC__init__ created',self.figdir
        return
    def expt(self,events=10):
        '''
        perform a toy experiment 
        '''
        print 'detMC.expt Generate',events,'per isotope'
        decays = {}
        hists = {}
        Ndecays = events
        Nbin = 600
        emi,ema = 3000.,9000.
        d = (ema-emi)/float(Nbin)
        bins = [float(x) for x in range(int(emi),int(ema+d),int(d))]
        Nbin = bins
        Eall = numpy.array([])
        nrows,ncols=4,2
        fs = 10
        fst=  6
        fig,axs = plt.subplots(nrows=nrows,ncols=ncols)
        irow,icol = 0,0
        for isotope in self.DautNames:
            decays[isotope] = self.fill(Branch=isotope,Ndecays=Ndecays)
            print 'isotope,len(decays[isotope])',isotope,len(decays[isotope])
            if len(decays[isotope])>0:
                label = isotope#+' '+str(len(decays[isotope])) + ' entries'
                axs[irow,icol].set_xlim([emi,ema])
                axs[irow,icol].hist(decays[isotope],Nbin,label=label,facecolor='red',histtype='step',color='red',alpha=1.,fill='True')
                #axs[irow,icol].set_title(label,fontsize=fs,loc='left',weight='bold')
                axs[irow,icol].text(emi+0.05*(ema-emi),1000.,label,fontsize=fs,weight='bold')
                axs[irow,icol].set_xlabel('Alpha energy (keV)',fontsize=fs)
                irow += 1
                if irow==nrows:
                    icol += 1
                    irow = 0

                Eall = numpy.append(Eall,decays[isotope])
        # sum of all distributions
        axs[irow,icol].set_xlim([emi,ema])
        axs[irow,icol].hist(Eall,Nbin,label='Sum',facecolor='red',histtype='step',color='red',alpha=1.,fill='True')

        axs[irow,icol].text(emi+0.1*(ema-emi),1000.,'Sum',fontsize=fs,weight='bold')
        axs[irow,icol].set_xlabel('Alpha energy (keV)',fontsize=fs)
        ymi,yma = axs[irow,icol].get_ylim()
        for irow in range(nrows):
            for icol in range(ncols):
                axs[irow,icol].set_ylim([ymi,yma])
                axs[irow,icol].grid()
                axs[irow,icol].tick_params(axis='both',labelsize=fst)



        fig.subplots_adjust(hspace=0.4)
        fig.suptitle(str(self.pePerMeVa) +' pe/MeV of alpha energy')
        #plt.tight_layout()
        pdf = self.figdir + str(int(self.pePerMeVa)) + 'pePerMeValpha.pdf'
        plt.savefig(pdf)
        print 'detMC.expt Wrote',pdf
        #        plt.show()
        return
    def readAll(self):
        ''' for testing reading of files '''
        for isotope in self.DautNames:
            data = self.read(isotope)
        return
    def read(self,Branch='227Th'):
        '''
        read alpha decay NNDC data that has been extracted
        '''
        debug = False
        fn = 'Isotopes/'+Branch+'.data'
        data = []
        tot = 0.
        if os.path.isfile(fn):
            f = open(fn,'r')
            for l in f:
                s = l.split()
                if debug : print s[0],s[1],s[2]
                E,note,pct = s[0],s[1],s[2]
                E = float(E)
                if '%' in pct: pct = note
                pct = float(pct)
                tot += pct
                data.append( (E, pct/100.) )
            f.close()
        print 'detMC.read',fn,'found',len(data),'lines Total alpha Br',tot
        return data
    def fill(self,Branch='222Th',Ndecays=1000):
        '''
        return energies of Ndecays decays for input branch
        '''
        E = numpy.array([])
        alphaBranches = self.DautAlpha[Branch]
        #print 'Branch,len(alphaBranches)', Branch,len(alphaBranches)
        if len(alphaBranches)==0: return E
        bsum = 0.
        for pair in alphaBranches:
            if len(pair)>0:
                e,b = pair
                bsum += b
                s = self.Res(e)
                #print 'Branch,b,e,s',Branch,b,e,s
                E = numpy.append(E, numpy.random.normal(e,s,int(Ndecays*b)))
        E = numpy.asarray(E)
        #print 'Branch,bsum',Branch,bsum
        return E
    def Res(self,E):
        Emev = E/1000.
        Npe = self.pePerMeVa * Emev
        s = math.sqrt(Npe)
        s = s/Npe * Emev
        s = s * 1000.
        return s

if __name__ == '__main__' :
    pe_MeV = 10.
    Nev    = 100000
    if len(sys.argv)>1: pe_MeV = float(sys.argv[1])
    if len(sys.argv)>2: Nev    = int(sys.argv[2])
    dMC = detMC(pe_MeV=pe_MeV)

    
    #dMC.readAll()
    #sys.exit()
    
    dMC.expt(events=Nev)
    sys.exit()
