#!/usr/bin/env python
'''
calculate r_b, etc. from P.K.Resmi et al, JHEP01(2018)082 

20190320
'''
import sys
import math
import numpy


class D0decay():
    def __init__(self):
        self.bins = range(1,10)
        
        b = self.binKK = {} # contents of table 6
        b[1] = [ [0.2224, 0.0187], [0.1768, 0.0168] ] # bin1, [Ki,dKi], [Kbari,dKbari]
        b[2] = [ [0.3933, 0.9219], [0.1905, 0.0173] ]
        b[3] = [ [0.0886, 0.0128], [0.3176, 0.0205] ]
        b[4] = [ [0.0769, 0.0119], [0.0469, 0.0093] ]
        b[5] = [ [0.0576, 0.0105], [0.0659, 0.0109] ]
        b[6] = [ [0.0605, 0.0107], [0.0929, 0.0128] ]
        b[7] = [ [0.0454, 0.0094], [0.0450, 0.0091] ]
        b[8] = [ [0.0233, 0.0068], [0.0195, 0.0061] ]
        b[9] = [ [0.0319, 0.0079], [0.0447, 0.0091] ]

        c = self.bincs = {} # contents of table 9
        c[1] = [ [-1.11, 0.09, 0.02, 0.01], [0., 0., 0., 0.] ] # bin1, c,dc,dc+sys,dc-sys, s ibid
        c[2] = [ [-0.30, 0.05, 0.01, 0.01], [-0.03, 0.09, 0.01, 0.02] ]
        c[3] = [ [-0.41, 0.07, 0.02, 0.01], [ 0.04, 0.12, 0.01, 0.02] ]
        c[4] = [ [-0.79, 0.09, 0.05, 0.05], [-0.44, 0.18, 0.06, 0.06] ]
        c[5] = [ [-0.62, 0.12, 0.03, 0.02], [ 0.42, 0.20, 0.06, 0.06] ]
        c[6] = [ [-0.19, 0.11, 0.02, 0.02], [0., 0., 0., 0.] ]
        c[7] = [ [-0.82, 0.11, 0.03, 0.03], [-0.11, 0.19, 0.04, 0.03] ]
        c[8] = [ [-0.63, 0.18, 0.03, 0.03], [ 0.23, 0.41, 0.04, 0.03] ]
        c[9] = [ [-0.69, 0.15, 0.15, 0.12], [0., 0., 0., 0.] ]
        return
    def doR(self,KK):
        K,aK = KK
        a,da = K
        b,db = aK
        r = b/a
        dr = r*math.sqrt( da*da/a/a + db*db/b/b)
        return r,dr
    def main(self):
        print '{0:5} {1:5} {2:5} {3:5} {4:5}'.format('r_b','c_b','s_b','xpre','ypre')
        sumK,sumaK = 0.,0.
        for b in self.bins:
            KK = self.binKK[b]
            sumK += KK[0][0]
            sumaK+= KK[1][0]
            r,dr = self.doR(KK)
            cs = self.bincs[b]
            c = cs[0][0]
            s = cs[1][0]
            ypre = (1.-r)*c*math.sqrt(r)
            xpre = (1.+r)*s*math.sqrt(r)
            print '{0:5.2f} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f}'.format(r,c,s,xpre,ypre)
        print 'sumK',sumK,'sumaK',sumaK
        return 
        
if __name__ == '__main__':
    D0d = D0decay()
    D0d.main()
