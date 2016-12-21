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
        read log file
        get timestamp, run time, sample name, source(s), etc.
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
                sources = l.split(':')[-1].split()
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
        return timestamp,sources,sample,runtime
if __name__ == '__main__' :
    rLF = readLogFile()
    fn = '/Users/djaffe/work/WaveDumpData/logfiles/run00078_ts1482161661.log'
    if len(sys.argv)>1 : fn = sys.argv[1]
    ts,s,sam,rt = rLF.readFile(fn=fn)
    print fn,'timestamp,sources,sample,runtime',ts,s,sam,rt
