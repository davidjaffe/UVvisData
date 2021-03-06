#!/usr/bin/env python
'''
read log files 227Ac setup and get run information
20161220
'''
import math
import sys
import numpy
import datetime,os


class readLogFile():
    def __init__(self):
        print 'readLogFile Initialized'
        self.Comments = None
        return
    def readFile(self,fn=None):
        '''
        if logfile exists, then read log file
        get timestamp, run time, sample name, source(s), etc.
        If timestamp not in log file, try to take timestamp from filename.
        If actual run time not in log file, try to use requested run time
        20170117 Add treatment for log files with known typos
        20170829 pass comments in global
        '''
        if fn is None:
            sys.exit('readLogFile.readFile ERROR No input file specified')
        if not os.path.isfile(fn):
            sys.exit('readLogFile.readFile ERROR '+fn+' does not exist')
            
        f = open(fn,'r')
        bn = os.path.basename(fn)
        timestamp,sources,sample,runtime = None,None,None,None
        request,actual = None,None
        for l in f:
            if 'Timestamp' in l:
                try: 
                    timestamp = int(l.split()[-1])
                except ValueError:
                    pass
            if ' source ' in l:
                newl = self.cleanLine(l)
                sources = newl.split(':')[-1].split()
                if sources is None:
                    print 'readLogFile.readFile in',fn,'No sources in line:',newl
            if 'Sample being measured' in l:
                sample = l.split(':')[-1][:-1]
            if 'Comments:' in l:
                self.Comments = l.split(':')[1][:-1]
            if 'Request run time' in l:
                try: 
                    request = float(l.split()[-1])
                except ValueError:
                    pass
            if 'Actual run time' in l:
                try:
                    actual = float(l.split()[-1])
                except ValueError:
                    pass
        f.close()
        
        if request is not None: runtime = request
        if actual  is not None:
            if actual<1.e10: runtime = actual

        if timestamp is None: # get timestamp from filename
            bn = os.path.basename(fn)
            ts = (bn.split('_')[1]).split('.')[0]
            ts = ts.replace('ts','')
            try: 
                timestamp = int(ts)
            except ValueError:
                pass

        ### special fixes
        debug = False
        if debug: 
            print 'readLogFile.readFile bn',bn,'sample',sample,
            if sample is not None:
                print 'len(sample)',len(sample)
            else:
                print ''
        if sample is not None:
            slash = "\\"
            if slash in sample: 
                sample = sample.replace(slash,'')
            if bn=='run00158_ts1483996034.log' and sample=='LiLS#1': sample = 'LiLS#2'
            if bn=='run00224_ts1484615681.log' and (sample=='LiLS#2s' or 'LiLS#2s' in sample) :sample = 'LiLS#2'
            if  'LiLS#2 pedestal' in sample or 'LiLS#2 black' in sample: sample = 'LiLS#2'
            if sample==' ': sample = None
            
        return timestamp,sources,sample,runtime
    def getSampleNumber(self,sample):
        '''
        return sample number (int) given sample name (str)
        if no sample number, return -1
        '''
        #print 'readLogFile.getSampleNumber sample',sample
        if sample is None: return -1
        if '#' in sample:
            i = sample.index('#')
            sn = int( sample[i+1:] )
            return sn
        else:
            for i in range(len(sample)-1,0,-1):
                try:
                    sn = int(sample[i:])
                except ValueError:
                    if i+1>=len(sample): return -1
                    sn = int(sample[i+1:])
                    return sn
            return sn
    def cleanLine(self,line):
        '''
        return line with spurious characters removed
        '''
        bogus = [',','\\']
        newline = ''
        for c in line:
            if c not in bogus: newline += c
        return newline
    def numberNine(self,Number=9,debug=False):
        '''
        return list of run_timestamp and list of timestamp read from file for good runs for sample==Number
        

        AVOID GIGANTIC RUN 1052
        AVOID RUNS <1000
        '''
        fp = '/Users/djaffe/work/WaveDumpData/Copied_logfiles/'
        print 'readLogFile Check all files in',fp
        import get_filepaths
        gfp = get_filepaths.get_filepaths()
        file_paths = gfp.get_filepaths(fp)

        RTlist,TSlist = [],[]
        
        for fn in sorted( file_paths ):
            ts,s,sam,rt = self.readFile(fn=fn)
            #print fn
            bn = os.path.basename(fn)
            bnNoLog = bn.replace('.log','')
            samNum = self.getSampleNumber(sam)
            if samNum==Number:
                reject = ('BAD' in self.Comments) or ('short' in self.Comments) or ('~' in bn) or ('run01052' in bn) or ('run00' in bn)
                if not reject:
                    if debug: print sam,rt,bnNoLog,self.Comments
                    RTlist.append( bnNoLog )
                    TSlist.append( ts )
        return RTlist,TSlist
    def makeSCP(self):
        '''
        generate scp commands from Temple PC

        assumes ssh tunnel opened with
         ssh -L 1234:130.199.23.223:22 djaffe@gateway.phy.bnl.gov
        '''
        L = []
        for Number in [9,1]:
            LN,TS = self.numberNine(Number=Number)
            print 'readLogFile.makeSCP Found',len(LN),'sample',Number,'files'
            L.extend(LN)
            
        cnow = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')   #_%f')
        scriptName = 'scp_'+cnow+'.sh'
        fsn = open(scriptName,'w')
        

        origin = 'napolitanolab@localhost:~/BNL/RootTrees/'
        dest   = '/Users/djaffe/work/WaveDumpData/rootfiles/'
        suffix = '.root'
        NCMD = 0
        for rt in L:
            df = dest + rt + suffix
            if os.path.isfile(df):
                print 'readLogFile.makeSCP File exists:',df
            else:
                NCMD += 1
                cmd = 'scp -P 1234 '+ origin +rt+suffix+ ' '+ df + '\n'
                #print cmd
                fsn.write(cmd)
        fsn.close()
        print 'readLogFile.makeSCP Wrote',NCMD,'lines to',scriptName
        
        return
            
if __name__ == '__main__' :
    rLF = readLogFile()
    rLF.makeSCP()
    
    sys.exit('too boring to continue')
    
    simple = False
    if simple:
        fn = '/Users/djaffe/work/WaveDumpData/logfiles/run00078_ts1482161661.log'
        if len(sys.argv)>1 : fn = sys.argv[1]
            
        ts,s,sam,rt = rLF.readFile(fn=fn)
        sn = rLF.getSampleNumber(sam)
        print fn,'timestamp,sources,sample,runtime,sample#',ts,s,sam,rt,sn
    else:
        fp = '/Users/djaffe/work/WaveDumpData/logfiles/'
        print 'readLogFile Check all files in',fp
        import get_filepaths
        gfp = get_filepaths.get_filepaths()
        file_paths = gfp.get_filepaths(fp)

        
        RunTimeByWeek = {}
        for sn in range(0,10): RunTimeByWeek[sn] = []
        weekNumbers = []
        for fn in sorted( file_paths ):
            ts,s,sam,rt = rLF.readFile(fn=fn)
            timeStamp = ts
            if ts>2490712523.:timeStamp= ts/1000.
            weekNumber =  datetime.datetime.fromtimestamp(timeStamp).isocalendar()[1]
            if weekNumber not in weekNumbers: weekNumbers.append(weekNumber)
            sn = rLF.getSampleNumber(sam)
            if sn in RunTimeByWeek: RunTimeByWeek[sn].append( (weekNumber, rt) )
            #print fn,'timestamp,sources,sample,runtime,sample#',ts,s,sam,rt,sn
            Show = False
            Flag = 0
            if ts is None or s is None or sam is None or rt is None: Flag += 1
            if s is not None and 'Ac-227' not in s: Flag += 10
#            if sam is not None and 'LiLS#2'not in sam: Flag += 100
            if sam=='LiLS#2s': Flag += 1000
            if sam<1 : Flag += 10000
            if Flag: print fn,'timestamp,sources,sample,runtime',ts,s,sam,rt,sn,'Flag',Flag
        print 'Flag: 1=any=None, 10=bad src, 100=bad sample, 1000=LiLS#2s, 10000=bad sample#'

        weekNumbers.sort()
        for week in weekNumbers:
            print 'week#',week,'sample#,runTime(s)',
            for sn in range(0,10):
                runTime = 0.
                for pair in RunTimeByWeek[sn]:
                    if pair[0]==week: runTime += pair[1]
                print sn,runTime,
            print ''
            
