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
import cutsAndConstants

class detMC():
    def __init__(self,pe_MeV=10.):

        self.cAC = cutsAndConstants.cutsAndConstants()

        self.DautNames = ['227Ac','227Th','223Ra','219Rn','215Po','211Pb','211Bi','207Tl','219At']
        self.colors    = ['b'    ,'g'    ,'r'    ,'c'    ,'m'    ,'y'    ,'k'    ,'p'    , 'lime']

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

        # 227Ac and all significant (Br>1e-7) daughters
        self.parentBr = {'227Ac': [ ['227Th', 0.9862], ['223Fr', 0.0138] ],
                         '227Th': [ ['223Ra', 1.0] ],
                         '223Ra': [ ['219Rn', 1.0] ],
                         '219Rn': [ ['215Po', 1.0] ],
                         '215Po': [ ['211Pb', 1.0], ['215At', 2.3e-6] ],
                         '211Pb': [ ['211Bi', 1.0] ],
                         '211Bi': [ ['211Po', 0.0028], ['207Tl',0.9972] ],
                         '207Tl': [ ['207Pb', 1.0] ],
                         '211Po': [ ['207Pb', 1.0] ],
                         '223Fr': [ ['223Ra', 0.9999], ['219At', 6.e-5] ],
                         '215At': [ ['211Bi', 1.0] ],
                         '219At': [ ['215Bi', 0.97], ['219Rn', 0.03] ],
                         '215Bi': [ ['215Po', 1.0] ],
                         '207Pb': [ [None, 1.] ], # stable
                         }
        # normalize total br to unity for each parent
        parentBr = self.parentBr
        for parent in parentBr:
            sum = 0.
            for dautBr in parentBr[parent]:
                daut,Br = dautBr
                if daut is not None:
                    sum += Br
            for i,dautBr in enumerate(parentBr[parent]):
                daut,Br = dautBr
                if daut is not None:
                    parentBr[parent][i][1] = Br/sum
        self.parentBr = parentBr
        # total isotope branching fractions
        self.alphaBF = {}
                                    
        
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
    def getDauts(self,parent):
        '''
        return names of daughters, used by alphaAccounting
        '''
        if parent not in self.parentBr: print 'detMC.getDauts',parent,'is not a valid parent'
        dauts = []
        for d in self.parentBr[parent]: dauts.append(d[0])
        return dauts
    def getDautBr(self,parent,daut):
        '''
        return branching fraction of parent into daut, used by alphaAccounting
        '''
        for db in self.parentBr[parent]:
            if db[0]==daut : return db[1]
        return None
    def alphaAccounting(self):
        '''
        determine all possible decay chains from 227Ac to a stable isotope
        calculate the cumulative branching fraction for that chain
        and the number of alphas in each chain
        '''
        parent = '227Ac'
        Chains = [parent]
        while parent in Chains:
            dauts = self.getDauts(parent[-5:])
            if dauts[0] is None:
                break
            
            for d in dauts:
                Chains.append( parent + d )
            Chains.pop(0)
            parent = Chains[0]

        # turn string describing chain into a list of strings
        ParsedChains = []    
        for X in Chains:
            chain = [X[i:i+5] for i in range(0,len(X),5)]
            ParsedChains.append(chain)
        # error check
        for i,c1 in enumerate(ParsedChains):
            for j,c2 in enumerate(ParsedChains):
                if j>i and c1==c2:
                    print 'detMC.alphaAccounting DUPLICATE i,j',i,j
                    sys.exit('detMC.alphaAccounting ERROR Duplicate chain')
            
        totAlpha = 0.
        alphaBF = {}
        PoBF = 0.
        coinc =  '219Rn215Po211Pb'
        for isotope in self.DautNames: alphaBF[isotope] = 0.
        for chain in ParsedChains:
            print '%s' % '=>'.join(map(str, chain)),
            prodBr = 1.
            sumAlpha = 0.
            for i,parent in enumerate(chain):
                if i+1<len(chain):
                    daut = chain[i+1]
                    prodBr *= self.getDautBr(parent,daut)
                    if self.isAlpha(parent,daut): sumAlpha += 1.
                
            print '{0:.6f} {1:.1f}'.format(prodBr,sumAlpha),
            if coinc in '%s' % ''.join(map(str, chain)):
                PoBF += prodBr
                print 'has',coinc,'coincidence'
            else:
                print ''
            totAlpha += prodBr*sumAlpha
            for isotope in alphaBF:
                if isotope in chain: alphaBF[isotope] += prodBr
            
        print 'total alphas',totAlpha,'into',coinc,' coincidences {0:.6f}'.format(PoBF)
        print 'isotope: product branching fraction',
        for isotope in alphaBF:
            print '{0}: {1:.6f}'.format(isotope,alphaBF[isotope]),
        print ''
        self.alphaBF = alphaBF
        return
    def isAlpha(self,parent,daut):
        ''' return True if parent -> daut, alpha '''
        if daut is None: return False
        ip = int(parent[:3])
        id = int(daut[:3])
        return ip-id==4
    def expt(self,events=10,units='keV',cf=1.):
        '''
        perform a toy experiment
        cf = 0.521e-3 pC/keV
        '''
        self.alphaAccounting()
        
        print 'detMC.expt Generate',events,'per isotope'
        decays = {}
        hists = {}
        Ndecays = events
        Nbin = 600

        emi,ema = 3000.,9000.
        if units=='pC':
            emi,ema = cf*emi, cf*ema
            print 'detMC.expt Results in',units,'assuming',cf,units+'/MeV'
        
        d = (ema-emi)/float(Nbin)
        bins = []
        for i in range(Nbin+1):
            bins.append(emi+float(i)*d)
#        bins = [float(x) for x in range(int(emi),int(ema+d),int(d))]
        Nbin = bins
        Eall = numpy.array([])
        nrows,ncols=4,2
        fs = 10
        fst=  6
        fig,axs = plt.subplots(nrows=nrows,ncols=ncols)
        irow,icol = 0,0
        for isotope in self.DautNames:
            decays[isotope] = cf*self.fill(Branch=isotope,Ndecays=Ndecays)
            print 'isotope,len(decays[isotope])',isotope,len(decays[isotope])
            if isotope=='219Rn' and units=='pC':
                vmi,vma = self.cAC.promptChargeCut[0]*1e3, self.cAC.promptChargeCut[1]*1e3 # 30.,40. # pC
                effy = self.cutEffy(decays[isotope],vmi,vma)
                print isotope,'effy',effy,'cut range('+units+')',vmi,vma
            if len(decays[isotope])>0:
                label = isotope#+' '+str(len(decays[isotope])) + ' entries'
                axs[irow,icol].set_xlim([emi,ema])
                axs[irow,icol].hist(decays[isotope],Nbin,label=label,facecolor='red',histtype='step',color='red',alpha=1.,fill='True')
                #axs[irow,icol].set_title(label,fontsize=fs,loc='left',weight='bold')
                axs[irow,icol].text(emi+0.05*(ema-emi),1000.,label,fontsize=fs,weight='bold')
                axs[irow,icol].set_xlabel('Alpha energy ('+units+')',fontsize=fs)
                irow += 1
                if irow==nrows:
                    icol += 1
                    irow = 0

                Eall = numpy.append(Eall,decays[isotope])
        # effy of cut on sum of all distributions
        if units=='pC':
            vmi,vma = self.cAC.lowChargeCut*1.e3, self.cAC.Qmax*1.e3
            effy = self.cutEffy(Eall,vmi,vma)
            print 'low charge cut effy',effy,'cut range('+units+')',vmi,vma
        axs[irow,icol].set_xlim([emi,ema])
        axs[irow,icol].hist(Eall,Nbin,label='Sum',facecolor='red',histtype='step',color='red',alpha=1.,fill='True')

        axs[irow,icol].text(emi+0.1*(ema-emi),1000.,'Sum',fontsize=fs,weight='bold')
        axs[irow,icol].set_xlabel('Alpha energy ('+units+')',fontsize=fs)
        ymi,yma = axs[irow,icol].get_ylim()
        for irow in range(nrows):
            for icol in range(ncols):
                axs[irow,icol].set_ylim([ymi,yma])
                axs[irow,icol].grid()
                axs[irow,icol].tick_params(axis='both',labelsize=fst)



        fig.subplots_adjust(hspace=0.4)
        fig.suptitle(str(self.pePerMeVa) +' pe/MeV of alpha energy')
        #plt.tight_layout()
        pdf = self.figdir + str(int(self.pePerMeVa)) + 'pePerMeValpha_units_'+units+'.pdf'
        plt.savefig(pdf)
        print 'detMC.expt Wrote',pdf
        #        plt.show()
        return
    def cutEffy(self,values,vmi,vma):
        '''
        determine fraction of values in range [vmi,vma]
        '''
        tot = float(len(values))
        sum = 0.
        for v in values:
            if vmi<=v and v<=vma: sum+=1.
        effy = -1.
        if tot>0: effy = sum/tot
        return effy
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
        BF = self.alphaBF[Branch]
        #print 'Branch,len(alphaBranches)', Branch,len(alphaBranches)
        if len(alphaBranches)==0: return E
        bsum = 0.
        for pair in alphaBranches:
            if len(pair)>0:
                e,b = pair
                bsum += b
                s = self.Res(e)
                #print 'Branch,b,e,s',Branch,b,e,s
                E = numpy.append(E, numpy.random.normal(e,s,int(Ndecays*b*BF)))
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
    units  = 'keV'
    cf     = 1.0
    if len(sys.argv)>1: pe_MeV = float(sys.argv[1])
    if len(sys.argv)>2: Nev    = int(sys.argv[2])
    if len(sys.argv)>3:
        units = 'pC'
        cf    = 38.5/7386.1 # pC/keV from run98 analysis 215Po mean 38.5 pC for 7386.keV
    dMC = detMC(pe_MeV=pe_MeV)

    
    #dMC.alphaAccounting()
    #sys.exit()
    
    dMC.expt(events=Nev,units=units,cf=cf)
    sys.exit()
