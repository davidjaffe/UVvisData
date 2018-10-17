#!/usr/bin/env python
'''
daya bay publication stats
20161207
20171110 update
20181013 update
'''
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

class pubstats():
    def __init__(self):
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {0:'o', 1:'s',2:'D',3:'+'} # filled circle, square, diamond
        
        self.figdir = 'Figures/'

        self.durNames = ['IntlRev', '1stCollRev', '2ndCollRev', 'Accepted', 'Published','Rejected']
        # name :  [2wk-IR, 1wk-2wk, sub-1wk, acc-sub, pub-acc]
        self.pubs   = {'LongReac'     : [54, 46, 20, 94, 67,0],
                       'LongOsc'      : [24, 29,125,  105, 53,0],
                       '68ADSter'     : [31, 32, 40, 43,51,0],
                       'DYBMINSteril' : [ 0, 43, 19, 54,40,0],
                       'Decoh'        : [29, 17, 30, 301,105,0], # not accepted at PLB, acc by EJPC
                       'FuelEvol'     : [84,  138, 51,  24,  52,0],
                       'MuonMod'      : [132, 85, 175, 208, 35, 0],
                       'NeutProd'     : [82,  62,  39, 54, 68, 0],
                       'Sidereal'     : [45, 71, 45, 0,0,0],
                       'ImpEffReac'   : [155, 32, 60,0,0,0],
                       'Osc2018'      : [42, 14, 48,0,0,0]
                       }
        self.order = ['LongReac', 
                       'LongOsc', 
                       '68ADSter', 
                       'DYBMINSteril', 
                       'Decoh', 
                       'FuelEvol',
                       'MuonMod',
                       'NeutProd',                       
                       'Sidereal',
                       'ImpEffReac',
                       'Osc2018'    
                       ]
        self.orderLong = {'LongReac' : 'Long\nReac\nCPC',
                        'LongOsc'    : 'Long\nOsc.\nPRD', 
                       '68ADSter'    : '6+8AD\nsterile\nPRL', 
                       'DYBMINSteril': 'D/M\nsterile\nPRL', 
                       'Decoh'      : 'Wave\n packet\n EJPC',
                       'FuelEvol'     : 'Fuel\nevol\nPRL',
                       'MuonMod'      : 'Muon\nMod\nJCAP',
                       'NeutProd'    :  'Neut\nProd\nPRD',
                       'Sidereal'   : 'Side\nMod\nPRD',
                       'ImpEffReac' : 'Eff\nReac\nPRD',
                       'Osc2018'     : 'Osc\n2018\nPRL'
                       }
        
        return
    def timetosub(self):
        '''
        histogram of the time taken to submit to a journal
        '''
        plt.rcdefaults() # why?


        fig,ax = plt.subplots()

        x = []
        for name in self.pubs:
            a = sum(self.pubs[name][:3])
            x.append(a)
        nbin = 8
        hist,bins = np.histogram(x,bins=nbin)

        ax.set_yticks(np.arange(0,4))
        ax.set_ylim([0,max(hist)*1.1])
        ax.set_xlim([0,int(max(x)/nbin+1)*nbin])
        ax.set_xlabel('Days to submission to journal')
        ax.set_title(str(int(np.mean(x)))+' days to submission on average')

        width = 0.7/0.9 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.grid()

        pdf = self.figdir + 'average_days.pdf'
        plt.savefig(pdf)
        print 'pubstats.timetosub Wrote',pdf
        return

    def comments(self):
        '''
        trend of # of comments/manuscript vs time
        '''
        plt.rcdefaults() # why?
        fig,ax = plt.subplots()
        y = [4,9,7,3,12,6,6,7,4]
        x = np.arange(1,len(y)+1)
        meany = np.average(y)

        xmi,xma = 0.5,max(x)+0.5
        ax.set_xlim([xmi,xma])
        ax.set_ylim([0,max(y)*1.1])
        ax.set_xlabel('Chronological order of submission to journal')
        ax.set_ylabel('The number of collaborators who contributed comments')
        ax.set_title('{0:.2f}'.format(meany)+' comments/manuscript on average')
        plt.grid()


        l = plt.plot(x,y, 'ro')
        plt.setp(l, markersize=30/2)
        plt.setp(l, markerfacecolor='red')

        l = plt.plot([xmi,xma],[meany,meany],'b--')


        pdf = self.figdir + 'trend.pdf'
        plt.savefig(pdf)
        print 'pubstats.trend Wrote',pdf
        
        return

    def plot(self):
        '''
        make horizontal,stacked bar chart
        '''
        plt.rcdefaults() # why?

        width = 0.5 # width of bars in bar chart
        offset = 0.1 # horiz offset of first bar from ordinate axis
        
        left = np.array( range(len(self.order)) ) + offset
        bottom = [0 for x in range(len(self.order))]
        p = []
        maxheight = 0
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
        
        plt.xticks(left+width/2.,xlab,fontsize=10)
        plt.ylabel('Days')
        plt.xlabel('Paper')
        #plt.figure()
        plt.grid()
        plt.title('Publication history and status since March 2016')
        plt.legend(p, self.durNames, bbox_to_anchor=(0.40-.07,1.05-.1))#(1.1, 1.05-.1))
        pdf = self.figdir + 'test.pdf'
        plt.savefig(pdf)
        print 'pubstats.plot Wrote',pdf
        return
      
if __name__ == '__main__' :
    P = pubstats()
    P.plot()
    P.timetosub()
    P.comments()
