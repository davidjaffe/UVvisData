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
        return
    def readFile(self,fn=None):
        '''
        if logfile exists, then read log file
        get timestamp, run time, sample name, source(s), etc.
        If timestamp not in log file, try to take timestamp from filename.
        If actual run time not in log file, try to use requested run time
        '''
        if fn is None:
            sys.exit('readLogFile.readFile ERROR No input file specified')
        if not os.path.isfile(fn):
            sys.exit('readLogFile.readFile ERROR '+fn+' does not exist')
            
        f = open(fn,'r')
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
            if 'Sample being measured' in l:
                sample = l.split()[-1]
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
        if actual  is not None: runtime = actual

        if timestamp is None: # get timestamp from filename
            bn = os.path.basename(fn)
            ts = (bn.split('_')[1]).split('.')[0]
            ts = ts.replace('ts','')
            try: 
                timestamp = int(ts)
            except ValueError:
                pass
            
        return timestamp,sources,sample,runtime
    def cleanLine(self,line):
        '''
        return line with spurious characters removed
        '''
        bogus = [',','\\']
        newline = ''
        for c in line:
            if c not in bogus: newline += c
        return newline
if __name__ == '__main__' :
    rLF = readLogFile()
    fn = '/Users/djaffe/work/WaveDumpData/logfiles/run00078_ts1482161661.log'
    if len(sys.argv)>1 : fn = sys.argv[1]
    ts,s,sam,rt = rLF.readFile(fn=fn)
    print fn,'timestamp,sources,sample,runtime',ts,s,sam,rt
