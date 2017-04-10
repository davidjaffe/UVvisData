#!/usr/bin/env python
'''
define cuts and constants used for 227Ac setup and processing
20161224
'''
import math
import sys
import numpy
import time,os
import datetime


class cutsAndConstants():
    def __init__(self):
        self.Po215halflife = 0.001781 # from NNDC, units = seconds
        self.Po215lifetime = self.Po215halflife/math.log(2.)  
        self.tOffset = 10.*self.Po215lifetime

        self.Ac227halflife = 21.772 # from Nucl.Data Sheets 93, 763 (2001), units = years
        self.Ac227lifetime = self.Ac227halflife/math.log(2.)
        self.Ac227lifetime_sec = self.Ac227lifetime * 365. * 24. * 60. * 60. # lifetime in seconds
        
        self.Qmax = 0.05*2 # 20170316 increase Qmax range
        self.psdCut = 0.35
        self.lifeRange = [1,2,3,5]
        self.lowChargeCut = 0.012
        self.promptChargeCut = [0.03,0.04]
        nsig = self.nsigmaCutPo215Peak = 3.0
        mean,sigma = 0.03839, 0.00265 # these values from fit to run74 data
        self.delayChargeCut  = [mean-nsig*sigma,mean+nsig*sigma]

        # exclude high trigger rate at the start of runs after run607 when threshold lowered
        # apply this cut to all data just to make it easy
        self.startTimeCut = 0.01 
        
        self.maxTimeWindow = float(max(self.lifeRange))*self.Po215lifetime

        # measured activity and mass from Eckert & Ziegler see doc1482
        self.initialAc227Activity = 3.711e4 # Bq
        self.initialAc227ActivityRelUnc = 1.32e-2 # combined type A,B uncertainties in quad.
        self.initialAc227Date = '20160906' # at noon EDT, but that can't matter
        self.initialAc227mass = 10.22710 # grams
        self.dispensedAc227mass = 0.503 # grams, from notebook 20161215
        self.dispensedAc227massUnc = 0.010 # grams, estimated
        self.totalLiLSmass = 192. # grams, according to Richard

        self.Ac227MassPerLiLSGram = self.dispensedAc227mass/self.initialAc227mass / self.totalLiLSmass
        
        self.LiLS2mass = 10.030 # grams, according to notebook 20161215
        self.LiLSmassUnc = 0.010 # grams, estimated, common for all measurements
        self.LiLS2massUnc = self.LiLSmassUnc

        # this will be the database for sample masses
        self.SampleMass = {'LiLS2':self.LiLS2mass}

        self.SampleMass['LiLS3'] = self.LiLS3mass = 9.98 # g, https://elog-phy.bnl.gov/PROSPECT+Ac227/19
        self.SampleMass['LiLS4'] = self.LiLS4mass = 9.98
        self.SampleMass['LiLS5'] = self.LiLS5mass = 9.999
        self.SampleMass['LiLS6'] = self.LiLS6mass = 9.99
        self.SampleMass['LiLS7'] = self.LiLS7mass = 9.981
        self.SampleMass['LiLS8'] = self.LiLS8mass =10.011

        # materials 
        self.SampleMaterial = {'LiLS2':'Reference'}
        self.SampleMaterial['LiLS3'] = 'UVTacrylic'
        self.SampleMaterial['LiLS4'] = 'FEP'
        self.SampleMaterial['LiLS5'] = 'PLA'
        self.SampleMaterial['LiLS6'] = 'PEEK'
        self.SampleMaterial['LiLS7'] = 'RG188'
        self.SampleMaterial['LiLS8'] = 'Viton'

        # sample mass divided by reference sample mass
        if self.SampleMaterial['LiLS2']!='Reference': sys.exit('cutsAndConstants.__init__ ERROR LiLS2 is not the Reference sample')
        rmass = self.SampleMass['LiLS2']
        self.SampleMassByRef = {}
        for key in self.SampleMass:
            self.SampleMassByRef[key] = self.SampleMass[key]/rmass

        self.spikeLiLSmass = None
        
        # derived quantities
        self.totalAlphas = 5.
        self.lowChargeCutEffy =  0.999998600339 # cut range(pC) 12.0 50.0, from detMC
        self.promptChargeCutEffy =  0.947929001184 # 219Rn cut range(pC) 30.0 40.0 from detMC

        # bad runs
        # run 56 has event with dt = 6.7 second and log file "Comments: Accidently hit keys while running"
        # run 158 has incorrect sample name of LiLS#1, it should be LiLS#2. Tag as bad until software fix implemented.
        # runs 162, 163 have empty log files. Tag as bad until corrected files obtained
        # runs 258,259,260,261,262 vial raised in well for optics test
        # run 263,264, 265 black paper @ 0,90,180 degrees for optics test
        # run 685, 653, 663 high trigger rates
        # run 676, 266 no HV
        # run 544 wrong threshold
        # run 689 high trigger rate
        # Runs 287-295 are garbage test runs made during modification testing
        # run 658 LiLS6 very low apparent rate, PSD vs charge shows large peak at charge~.001, PSD~0.86
        self.badRuns = [56, 158, 258,259,260,261,262, 263,264,265,266, 676, 653, 663, 685, 544, 689, 658 ]
        self.badRuns.extend(  range(287,295+1) )
        
        
        print 'cutsAndConstants.__init__ Initialized'
        return
    def getSampleMass(self,inputSampleName=None):
        '''
        return best guess at sample mass given name
        '''
        if inputSampleName is None:
            sys.exit('cutsAndConstants.getSampleMass ERROR No sampleName given')
            
        if inputSampleName is list:
            sampleName = inputSampleName[0]
            if len(inputSampleName)>1 : print 'cutsAndConstants.getSampleMass WARNING inputSampleName is not a simple string. It is',inputSampleName
        else:
            sampleName = inputSampleName
                        
        sampleName = sampleName.replace(' ','') # remove blanks
        if sampleName in self.SampleMass: return self.SampleMass[sampleName]


        sN = sampleName.replace('#','') # remove # from string
        if sN in self.SampleMass: return self.SampleMass[sN]

        # try for match ignoring case and removing # from key
        for key in self.SampleMass:
            if key.lower()==sampleName.lower(): return self.SampleMass[key]
            if key.lower()==sN.lower()        : return self.SampleMass[key]
            if key.replace("#","").lower()==sN.lower() : return self.SampleMass[key]

        # Danielle initially used simple name for initial spiked sample
        # or bad typing run225
        # or deal with 'LiLS#2pedestal' and similar runs 260-265
        if sampleName=='LiLS' or sampleName=='LiLS#2s' or 'LiLS#2' in sampleName:
            key = 'LiLS2'
            return self.SampleMass[key]
                

        # can't do it
        print 'cutsAndConstants.getSampleMass FAILED. Input sampleName',sampleName,'not in existing names',self.SampleMass.keys()
        sys.exit('cutsAndConstants.getSampleMass FAILED. Input sampleName '+sampleName+' not in existing names')
        return None
    def setSpikeLiLSmass(self,mass):
        self.spikeLiLSmass = mass
        return
    def expectAc227Rate(self,inputDay=None,mass=None, sampleName=None, massUnc=None):
        '''
        return expected Ac227 activity in Bq given date and spiked LiLS mass or sample name
        Use today's date if no date is given.
        If spiked LiLS mass is not given, use set mass from cutsAndConstants.setSpikeLiLSmass
        Can accept inputDay as datetime.dateime object or as string, otherwise fails
        Abort if date prior to initial Ac227 activity measurement date.
        '''
        debug = False
        if debug : print 'cutsAndConstants.expectAc227Rate INPUTS inputDay,mass,sampleName',inputDay,mass,sampleName

        fmtYMD = '%Y%m%d' # 20161225, for example
        
        
        if mass is None:
            if sampleName is None:
                mass = self.spikeLiLSmass
            else:
                mass = self.getSampleMass(sampleName)
                
        if massUnc is None: massUnc = self.LiLSmassUnc
            
        day = None
        if inputDay is None: day = now = datetime.datetime.now()
        # fmt = '%Y%m%d_%H%M%S_%f'
        # cnow = now.strftime(fmt)

        if type(inputDay)==datetime.datetime: day = inputDay

        if type(inputDay)==str : day = datetime.datetime.strptime( inputDay, fmtYMD )

        if type(inputDay)==int : day = datetime.datetime.fromtimestamp( inputDay )
            

        if day is None:
            print 'cutsAndConstants.expectAc227Rate type(inputDay)',type(inputDay),'inputDay',inputDay
            sys.exit('cutsAndConstants.expectAc227Rate Cannot process inputDay')

        day0 = datetime.datetime.strptime( self.initialAc227Date, fmtYMD)

        deltaT = day - day0
        deltaTsec = deltaT.total_seconds()
        #deltaTyears = deltaTsec/60./60./24./365.
        if deltaTsec<0.:
            print 'cutsAndConstants.expectAc227Rate inputDay',inputDay
            sys.exit('cutsAndConstants.expectAc227Rate ERROR Input day is before initial date '+self.initialAc227Date)

        rate = self.initialAc227Activity * mass*self.Ac227MassPerLiLSGram * math.exp(-deltaTsec/self.Ac227lifetime_sec)


        
        drate = rate * math.sqrt(math.pow(self.initialAc227ActivityRelUnc,2) + math.pow(self.dispensedAc227massUnc/self.dispensedAc227mass,2) + math.pow(massUnc/mass,2) )
            
        return rate,drate
    def allSampleRates(self,inputDay=None):
        '''
        calculates all sample rates given inputDay.
        if no Day given, calculates rate today
        '''
        print 'cutsAndConstants.allSampleRates sample,material,mass(g), mass/mass(ref), area(cm2), rate(Hz), rate/rate(ref) on',inputDay
        sr,sr2, sm,sm2 = 0.,0.,0.,0.
        refrate,drefrate = self.expectAc227Rate(sampleName='LiLS2',inputDay=inputDay)
        for sN in sorted(self.SampleMass.keys()):
            area = self.sampleArea(sampleName=sN)
            rate,drate = self.expectAc227Rate(sampleName=sN,inputDay=inputDay)
            ratewrtRef = rate/refrate
            dratewrtRef= ratewrtRef * drate/rate
            sr += rate
            sr2+= rate*rate
            mat = self.SampleMaterial[sN]
            m = self.SampleMass[sN]
            mr = self.SampleMassByRef[sN]
            sm += m
            sm2+= m*m
            dm= self.LiLSmassUnc
            print ' {0} {1:10} {2:6.3f}({3:5.3f}) {9:.4f} {4:4.1f} {5:.2f}({6:.2f}) {7:.4f}({8:.4f})'.format(
                sN,mat,m,dm,area,rate,drate,ratewrtRef,dratewrtRef,mr)
        N = float(len(self.SampleMass))
        sr = sr/N
        sr2= sr2/N
        dr = math.sqrt((sr2-sr*sr)/(N-1.))
        sm = sm/N
        sm2= sm2/N
        dm = math.sqrt((sm2-sm*sm)/(N-1.))
        print 'mean mass {0:.3f}({1:.3f}) rate {2:.2f}({3:.2f}) Hz '.format(sm,dm,sr,dr)
        return
    def sampleArea(self,sampleName=None):
        '''
        calculate sample total surface area in square cm from measurements in logbook 20170224
        '''
        area = -1.
        if sampleName=='LiLS3': # UVT acrylic
            length = 1.0
            width  = 1.15
            thick  = 0.1
            area = 2*(length*width + length*thick + thick*width)
        if sampleName=='LiLS4': # FEP
            length = 1.5
            width  = 1.5
            thick  = 0.003*2.54
            area = 2*(length*width + length*thick + thick*width)
        if sampleName=='LiLS5': # PLA (10 disks)
            diam = 0.5
            thick= 0.1
            n    = 10.
            area = n*(2.*math.pi*diam*diam/4. + math.pi*diam*thick) # 2ends+ sides
        if sampleName=='LiLS6': # PEEK nut
            e = 1.1  # corner to corner
            d = 1.0  # side to side = 2 x apothem
            D = 0.5  # diameter of hole
            t = 0.5  # thickness
            c = math.sqrt(e*e-d*d)
            area = 2.*(3./2.*c*d - math.pi*D*D/4) # top + bottom area with hole subtracted
            area+= 6*c*t       # outer sides
            area+= math.pi*D*t # hole wall
        if sampleName=='LiLS7': # RG-188 cable
            h = 2.54 # height of liquid
            d = 2.54/2. # ID of vial (estimated, OD=28mm according to spec sheet)
            t = 0.2489 # diameter of cable from spec sheet
            length = 2.*math.sqrt(h*h + d*d/4.)
            area = math.pi*t*length
        if sampleName=='LiLS8': # viton (10 o-rings)
            OD = 0.6
            ID = 0.3
            OR = OD/2.
            IR = ID/2.
            n = 10.
            area = n*math.pi*math.pi*(OR-IR)*(OR+IR)
        return area
        
        
if __name__ == '__main__':
    cAC = cutsAndConstants()

    if 0:
        sn = 'LiLS#2'
        m = cAC.getSampleMass(sn)
        cAC.setSpikeLiLSmass(m)

        inputDay = '20170103'
        x,dx = cAC.expectAc227Rate(inputDay=inputDay)
        print sn,'inputDay',inputDay,'rate',x,'+/-',dx
    if 1:
        cnow = datetime.datetime.now().strftime('%Y%m%d')
        cAC.allSampleRates(inputDay=cnow)
