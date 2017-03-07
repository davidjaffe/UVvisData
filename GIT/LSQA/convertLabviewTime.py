#!/usr/bin/env python
'''
convert labview time as recorded by labview for PSD data
to python datetime object or string

secret knowledge from http://blogs.cae.tntech.edu/mwr/2008/04/14/converting-national-instruments-lvm-timestamps-to-excel/ "NI counts a real-valued number of seconds since January 1, 1904"
20170224
'''
import sys
# the ONETON dir contains wfanal.py 
#sys.path.append('/Users/djaffe/work/GIT/ONETON')
# the PROSPECT dir contains get_filepaths.py
#sys.path.append('/Users/djaffe/work/GIT/PROSPECT')


import datetime


class convertLabviewTime():
    def __init__(self):
        self.start = datetime.datetime(1904,1,1)
        return
    def convert(self,tsin,fmt='%Y%m%d %H:%M:%S'):
        '''
        convert timestamp recorded by labview and used to create the run number for the PSD analysis
        return date/time according to input format
        if fmt is None, return datetime object
        otherwise return string 
        '''
        ts = tsin
        if type(ts) is str: ts = int(tsin)
        ts = float(ts) 
        step = datetime.timedelta(seconds=tsin)
        dt = self.start + step
        if fmt is None: return dt
        return dt.strftime(fmt)
if __name__ == '__main__':
    cLT = convertLabviewTime()
    ist = 3570538570
    s = cLT.convert(ist)
    print 'ist',ist,'s',s,'type(s)',type(s) 
    dt = cLT.convert(ist,fmt=None)
    print 'ist',ist,'dt',dt,'type(dt)',type(dt)
    
