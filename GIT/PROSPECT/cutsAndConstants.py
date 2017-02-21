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
        
        self.Qmax = 0.05
        self.psdCut = 0.35
        self.lifeRange = [1,2,3,5]
        self.lowChargeCut = 0.012
        self.promptChargeCut = [0.03,0.04]
        nsig = self.nsigmaCutPo215Peak = 3.0
        mean,sigma = 0.03839, 0.00265 # these values from fit to run74 data
        self.delayChargeCut  = [mean-nsig*sigma,mean+nsig*sigma]

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
        self.LiLS2massUnc = 0.010 # grams, estimated

        # this will be the database for sample masses
        self.SampleMass = {'LiLS2':self.LiLS2mass}

        self.spikeLiLSmass = None
        
        # derived quantities
        self.totalAlphas = 5.
        self.lowChargeCutEffy =  0.999998600339 # cut range(pC) 12.0 50.0, from detMC
        self.promptChargeCutEffy =  0.947929001184 # 219Rn cut range(pC) 30.0 40.0 from detMC

        # bad runs
        # run 56 has event with dt = 6.7 second and log file "Comments: Accidently hit keys while running"
        # run 158 has incorrect sample name of LiLS#1, it should be LiLS#2. Tag as bad until software fix implemented.
        # runs 162, 163 have empty log files. Tag as bad until corrected files obtained
        # runs 259,260,261,262 vial raised in well for optics test
        # run 263,264 black paper @ 0,90,180 degrees for optics test
        self.badRuns = [56, 158, 259,260,261,262, 263,264,265 ]

        
        
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
    def expectAc227Rate(self,inputDay=None,mass=None, sampleName=None):
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


        
        drate = rate * math.sqrt(math.pow(self.initialAc227ActivityRelUnc,2) + math.pow(self.dispensedAc227massUnc/self.dispensedAc227mass,2) + math.pow(self.LiLS2massUnc/self.LiLS2mass,2) )
            
        return rate,drate
if __name__ == '__main__':
    cAC = cutsAndConstants()

    if 1:
        sn = 'LiLS#2'
        m = cAC.getSampleMass(sn)
        cAC.setSpikeLiLSmass(m)

        inputDay = '20170103'
        x,dx = cAC.expectAc227Rate(inputDay=inputDay)
        print sn,'inputDay',inputDay,'rate',x,'+/-',dx
