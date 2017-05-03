#!/usr/bin/env python
'''
plot neutron production rate for some nuclei
Table 1 R.Heaton NIMA276 (1989) 529 
20170502
'''
import math
import sys
#import random
import numpy
import copy
#import scipy
#from scipy.stats.mstats import chisquare
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class alphan():
    def __init__(self,makeFigure=False):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/alphan/'

        self.dataFile = 'neutrons_per_alpha.txt'

        self.jendlFile= 'F019.dat'
        self.MTdescrip= {}

        self.makeFigure = makeFigure
        
        return
    def readJENDL(self,idebug=0):
        '''
        read JENDL file
         MF=3  (a,n) Reaction Cross Sections                               925 1451   18
   MT=4    Neutron Production Cross Section                        925 1451   19
   
   MT=50   (a,n) Cross Section to the Ground State                 925 1451   78
       Calculated by using the EGNASH-2 program.                   925 1451   79
                                                                   925 1451   80
   MT=51-77  (a,n) Cross Sections to the Discrete Excited Levels   925 1451   81
       Calculated by using the EGNASH-2 program.                   925 1451   82

        '''
        validMT = [4]
        for i in range(50,78): validMT.append(i)
        self.validMT = validMT
        self.MTdescrip[4] = 'n prod xsec'
        self.MTdescrip[50]= 'xsec to gs'
        for i in validMT:
            if i not in self.MTdescrip: self.MTdescrip[i] = 'xsec to es'
                
        iHEAD = 0
        iTAB1 = iHEAD + 1
        iINTR = iTAB1 + 1
        
        
        f = open(self.jendlFile,'r')
        print 'alphan.readJENDL Opened',self.jendlFile
        first = True
        Buffer = {}
        for l in f:
            if not first:
                MMM = l[67:75]
                if idebug>1: print 'l',l,'\nMMM',MMM,
                I1=67
                I2=I1+4
                MAT = int(l[I1:I2])
                I1=I2
                I2=I1+2
                MF = int(l[I1:I2])
                I1=I2
                I2=I1+3
                MT = int(l[I1:I2])
                if idebug>1: print 'MAT',MAT,'MT',MT,'MF',MF
                if MF==3 and MT in validMT: 
                    if MT not in Buffer:
                        if idebug>2: print 'initializing MT',MT
                        Buffer[MT] = []
                    Buffer[MT].append(l)
            first = False
        f.close()

        pData = {}
        for MT in validMT:

            if idebug>2: print 'Buffer[MT]',Buffer[MT]

            HEAD = Buffer[MT][iHEAD]
            ZA,QWR = self.floty(HEAD[:11]),self.floty(HEAD[11:22])
            A = int(ZA)%1000
            Z = (int(ZA)-A)/1000
            TAB1 = Buffer[MT][iTAB1]
            QM,QI = self.floty(TAB1[:11]),self.floty(TAB1[11:22])
            INTR = Buffer[MT][iINTR]
            Energy,Xsection = [],[]
            for i in range(iINTR+1,len(Buffer[MT])):
                DATA = Buffer[MT][i][:67]
                if idebug>1: print 'MT',MT,'i',i,'DATA',DATA
                E,S = self.parseDATA(DATA)
                Energy.extend(E)
                Xsection.extend(S)
            pData[MT] = [Z, A,QWR, QM,QI, Energy,Xsection]
            if idebug>0: print 'alphan.readJENDL MT',MT,'ZA',ZA,'Z,A',Z,A,'QWR,QM,QI',QWR,QM,QI,'# energies,xsections',len(Energy)
        return pData
    def parseDATA(self,DATA):
        '''
        DATA is a string of the form
        A1 B1 A2 B2 ...
        A1 = energy in eV
        B1 = cross-section in barn
        '''

        #print 'alphan.parseDATA DATA',DATA
        E,S = [],[]
        for i in range(0,len(DATA.strip())-1,11*2):
            #print 'alphan.parseDATA i',i,'DATA[i:i+11],DATA[i+11:i+11+11]',DATA[i:i+11],DATA[i+11:i+11+11]
            energy,xsec = self.floty(DATA[i:i+11]),self.floty(DATA[i+11:i+11+11])
            E.append(energy)
            S.append(xsec)
        return E,S
    def floty(self,s):
        #print 'alphan.floty s',s
        q,i = None,None
        if '+' in s:
            i = s.index('+')
            q = 1.
        elif '-' in s:
            i = s.index('-')
            q = -1.
            
        f = float(s[:i]) * math.pow(10,q*float(s[i+1:]))
        return f
    def read(self):
        ''' read file of neutrons/alpha from Heaton'''
        f = open(self.dataFile,'r')
        NperA = {}
        nucleus = None
        for l in f:
            line = l.replace('\n','')
            s = line.split(' ')
            if '#' in line:
                print 'alphan.read',line[:-1]
            elif '*' in line:
                nucleus = s[1]
                if nucleus in NperA: sys.exit('alphan.read ERROR duplicate nucleus ' + line[:-1])
                NperA[nucleus] = []
            elif len(s)>1:
                e,npera = float(s[0]),float(s[1])
                if nucleus is None: sys.exit('alphan.read ERROR no nucleus ' + line[:-1] )
                NperA[nucleus].append( [e,npera] )
        f.close()
        return NperA 
    def main(self,mode=2):
        '''
        main routine

       '''
        if mode==1: # heaton plots
            NperA = self.read()
            ip = 0
            fig,ax = plt.subplots()

            for ic,nucleus in enumerate(NperA):
                x = [pair[0] for pair in NperA[nucleus]]
                y = [pair[1] for pair in NperA[nucleus]]
                X,Y = numpy.array(x),numpy.array(y)    
                name = nucleus

                ax.plot(X,Y,self.colors[ic]+self.points[ip],label=nucleus)
                ax.set_yscale('log')

            ax.set_ylabel('Neutrons per incident alpha')
            ax.set_xlabel('Alpha energy (MeV)')
            plt.legend(loc='right')
            plt.grid()
            plt.title('neutron yield per alpha for thick targets R.Heaton etal NIMA276 (1989)529')
            
            if not self.makeFigure:
                plt.show()
            else:
                pdf = self.figdir + 'alphan.pdf'
                plt.savefig(pdf)
                print 'alphan.main created',pdf
        if mode==2: # JENDL
            ### JENDL[MT] = [Z, A,QWR, QM,QI, Energy,Xsection]
            JENDL = self.readJENDL()


            
            MT = 50 # excitation to gs
            gsE,gsXS = copy.deepcopy(JENDL[MT][5]),copy.deepcopy(JENDL[MT][6]) # total cross-section to  ground state
            exciteE,exciteXS = JENDL[MT][5],JENDL[MT][6] # total cross-section to excited states & ground state

            
            fig,ax = plt.subplots()
            for MT in self.validMT:
                E = JENDL[MT][5]
                XS= JENDL[MT][6]
                QM= JENDL[MT][3]/1.e3 # change to keV
                QI= JENDL[MT][4]/1.e3 # change to keV
                levelE = -(QI-QM)
                X,Y = numpy.array(E),numpy.array(XS)
                X = X/1.e3 # change to keV
                ic,ip = MT%len(self.colors),MT%len(self.points)
                label = self.MTdescrip[MT] + ' ' + str(int(levelE))
                ax.plot(X,Y,self.colors[ic]+self.points[ip]+'-',label=label)
                if MT!=4 and MT!=50:
                    e0 = E[1] # skip first point
                    for i,e in enumerate(exciteE):
                        #print 'i,e,e0',i,e,e0
                        if abs(e-e0)<10.: break
                    for j in range(len(E)-1):
                        #print 'MT',MT,'i',i,'j',j
                        exciteXS[i+j] += XS[j+1]

            # total xsec to gs and excited states
            X,Y = numpy.array(exciteE),numpy.array(exciteXS)
            X = X/1.e3 # change to keV
            ic,ip = 0,0
            label = 'xsec to gs + sum(es)'
            ax.plot(X,Y,self.colors[ic]+self.points[ip]+'-',label=label)
            print label,Y

            # total xsec to excited states
            X,Y = numpy.array(gsE),numpy.array(exciteXS)-numpy.array(gsXS)
            X = X/1.e3
            ic,ip = 1,1
            label = 'xsec to sum(es)'
            ax.plot(X,Y,self.colors[ic]+self.points[ip]+'-',label=label)
            print label,Y

            # ratio of excited states to total
            denom = numpy.array([max(1.e-5,x) for x in exciteXS])
            X,Y = numpy.array(gsE),(numpy.array(exciteXS)-numpy.array(gsXS))/denom
            X = X/1.e3
            ic,ip = 0,3
            label = 'sum(es)/total'
            ax.plot(X,Y,self.colors[ic]+self.points[ip]+'--',label=label)
            print label,Y
            
                        
            ax.set_yscale('log')
            ax.set_ylabel('Cross-section (b)')
            ax.set_xlabel('alpha energy (keV)')
            plt.legend(prop={'size':10-1}, numpoints=1,labelspacing=0.1,
                        columnspacing=0.1,handletextpad=0.1)
            plt.grid()
            if not self.makeFigure:
                plt.show()
            else:
                name = 'xsec_vs_alpha_energy_19F'
                pdf = self.figdir + name + '.pdf'
                plt.savefig(pdf)
                print 'alphan.main Wrote',pdf

                name += '_liny'
                ax.set_yscale('linear')
                pdf = self.figdir + name + '.pdf'
                plt.savefig(pdf)
                print 'alphan.main Wrote',pdf
 
        return

    
  
if __name__ == '__main__' :
    makeFigure = len(sys.argv)>1
    A = alphan(makeFigure=makeFigure)
    mode = 2 # 1=heaton,2=jendl
    A.main(mode=mode)
