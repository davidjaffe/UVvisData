#!/usr/bin/env python
'''
read temperature data for 227Ac setup
20170309
'''
import math
import sys
import numpy
import time,os
import datetime
import get_filepaths
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import SU

class readTemperature():
    def __init__(self):
        self.tDir = '/Users/djaffe/work/WaveDumpData/Temperature/'
        self.gfp = get_filepaths.get_filepaths()
        return
    def readfile(self,fn):
        '''
        read csv file, return time and temperature lists
        '''
        hdr = '\xef\xbb\xbf"Plot Title: BNL_Start_03032017_18_13"\r\n'
        f = open(fn,'r')
        num,datet,temp = [],[],[]
        for line in f:
            if line[0]!=hdr[0] and 'Logged' not in line:
                s = line.split(',')
                if s[0]!='"#"':
                    #print s
                    n,t,T = int(s[0]),s[1],float(s[2])
                    tt = t.replace('/','')
                    dt = datetime.datetime.strptime(tt,'%m%d%y %H:%M:%S %p')
                    num.append(n)
                    datet.append(dt)
                    temp.append(T)
        f.close()
        print 'readTemperature read',len(num),'entries from',fn
        return num,datet,temp
    def readAll(self):
        '''
        read all files in self.tDir
        '''
        fp = self.gfp.get_filepaths(self.tDir)
        num,datet,temp = [],[],[]
        for fn in fp:
            if '.csv' in fn:
                N,D,T = self.readfile(fn)
                num.extend(N)
                datet.extend(D)
                temp.extend(T)
        return num,datet,temp
    def main(self):
        num,datet,temp = self.readAll()
        X,Y = numpy.array(datet),numpy.array(temp)
        fig,ax = plt.subplots()
        ax.plot(X,Y,marker='o')
        days = mdates.DayLocator()
        weeks= mdates.WeekdayLocator(SU)
        ax.xaxis.set_minor_locator(days)
        ax.xaxis.set_major_locator(weeks)
        ax.xaxis.set_major_formatter( mdates.DateFormatter('%y%m%d') )
        plt.grid()
        plt.show()
        return
if __name__=='__main__':
    rT = readTemperature()
    rT.main()
