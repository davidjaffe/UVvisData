#!/usr/bin/env python
'''
read log files 227Ac setup and get run information
20161220
'''
import math
import sys
import numpy
import time,os


class timing():
    def __init__(self):
        self.b = numpy.array([0,1,2,3,4])
        self.n = ['a','b','c','d','e']
        return
    def fillD(self,x):
        d = {}
        for i,k in enumerate(self.n):
            d[k] = self.b[i]
        return d
    def fillA(self,x):
        A = self.b
        return A
    def simple(self,N=10000):




        t0 = time.clock()
        for i in range(N):
            A = self.fillA(i)
            L = len(A)
            for j in range(L):
                x = A[j]
        t1 = time.clock()
        for i in range(N):
            D = self.fillD(i)
        t2 = time.clock()
        A = self.fillA(i)
        for i in range(N):
            L = len(A)
            for j in range(L):
                x = A[j]
        t3 = time.clock()
        print '\n# of loops',N
        print 'array dt',t1-t0,'s ',(t1-t0)*1e6/float(N),'us/loop'
        print ' dict dt',t2-t1,'s ',(t2-t1)*1e6/float(N),'us/loop'
        print 'loclA dt',t3-t2,'s ',(t3-t2)*1e6/float(N),'us/loop'
        return
if __name__ == '__main__':
    T = timing()
    for N in [1000, 10000, 100000]:
        T.simple(N=N)
    sys.exit()
