#!/usr/bin/env python
'''
daya bay publication stats
20161207
'''
import math
import sys
import numpy
import matplotlib.pyplot as plt

class pubstats():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/'

        self.durNames = ['IntlRev', '1stCollRev', '2ndCollRev', 'Accepted', 'Published','Rejected']
        # name :  [2wk-IR, 1wk-2wk, sub-1wk, acc-sub, pub-acc]
        self.pubs   = {'LongReac'     : [54, 46, 20, 94, 0,0],
                       'LongOsc'      : [24, 29,125,  0, 0,0],
                       '68ADSter'     : [31, 32, 40, 43,51,0],
                       'DYBMINSteril' : [ 0, 43, 19, 54,40,0],
                       'Decoh'        : [29, 17, 30, 0,0,24+95], # not accepted
                       'FuelEvol'     : [84,  0,  0,  0,  0,0]
                       }
        self.order = ['LongReac', 
                       'LongOsc', 
                       '68ADSter', 
                       'DYBMINSteril', 
                       'Decoh', 
                       'FuelEvol'
                       ]
        self.orderLong = {'LongReac' : 'Long Reactor\nCPC',
                        'LongOsc'    : 'Long Osc.\nPRD', 
                       '68ADSter'    : '6+8AD\nsterile\nPRL', 
                       'DYBMINSteril': 'DYB/MINOS\nsterile\nPRL', 
                       'Decoh'      : 'Wave\n packet\n $PLB\!\!\!\!\!///$',
                       'FuelEvol'     : 'Fuel\nevolution\nPRL'
                       }
        
        return
    def plot(self):
        '''
        make horizontal,stacked bar chart
        '''
        plt.rcdefaults() # why?

        width = 0.5
        offset = 0.1
        
        left = numpy.array( range(len(self.order)) ) + offset
        bottom = [0 for x in range(len(self.order))]
        p = []
        for idur in range(len(self.durNames)):
            height = []
            for name in self.order:
                h = self.pubs[name][idur]
                height.append( h )
            p.append( plt.bar(left,height,color=self.colors[idur],bottom=bottom,width=width) )
            for i,x in enumerate(height):
                bottom[i] += x
        xlab = []
        for o in self.order: xlab.append(self.orderLong[o])
        plt.xticks(left+width/2.,xlab)
        plt.ylabel('Days')
        plt.xlabel('Paper')
        plt.grid()
        plt.title('Publication history and status since March 2016')
        plt.legend(p, self.durNames, bbox_to_anchor=(1.1, 1.05-.1))
        pdf = self.figdir + 'test.pdf'
        plt.savefig(pdf)
        print 'pubstats.plot Wrote',pdf
        return
      
if __name__ == '__main__' :
    P = pubstats()
    P.plot()
