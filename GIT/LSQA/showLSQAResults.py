#!/usr/bin/env python
'''
read, plot lsqa results
20170303
'''
import math
import sys
import numpy
import datetime
import matplotlib.pyplot as plt
import convertLabviewTime

class showLSQAResults():
    def __init__(self,mf):
        self.makeFigure = mf
        self.cLT = convertLabviewTime.convertLabviewTime()
        
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {}
        for i,p in enumerate(['s','o','v','^','<','>','X','p','D']):
            self.points[i] = p
        self.lencolors = len(self.colors)
        self.lenpoints = len(self.points)

        self.resFile= 'Samples/resultsSummaryFile.txt'
        
        self.figdir = 'Samples/ResultsFigures/'
        return
    def colorAndPoint(self,i):
        '''
        return color and point given index
        rotate through colors, then points
        '''
        j = i%(self.lenpoints*self.lencolors)
        ip = j/self.lencolors
        ic = j%self.lencolors
        return ic,ip
    def reader(self):
        '''
        parse results file
        '''
        fn = self.resFile
        f = open(fn,'r')

        self.cols = cols = ['run','EperQ','FOM','Z','fast','total','gDate','gTime','nDate','nTime','gFile','nFile','jobTime']
        self.results = results = {}

        cjT = 'jobTime'
        N = 0
        
        for line in f:
            if line[0]!='*' : # not a comment line
                N += 1
                s = line.split()
                r = {}
                for i,key in enumerate(cols):
                    r[key] = s[i]
                if s[0] in results:
                    t1,t2 = results[s[0]][cjT],r[cjT]
                    print 'showLSQAResults.reader replacing results from',cjT,t1,'with',t2
                    if t1>t2:
                        sys.exit('showLSQAResults.reader ERROR Replacing data with later job time with data from earlier job time!')
                    
                results[ s[0] ] = r
        f.close()
        print 'showLSQAResults.reader Read',N,'lines. Found',len(results),'results from',fn
        return
    def translate(self):
        '''
        translate results from file into reals, integers, useful strings
        '''
        results = self.results
        useful = {}
        for rn in sorted( results.keys() ):
            r = results[rn]
            u = {}
            for col in self.cols:
                if col in ['EperQ']:
                    u[col] = float(r[col])
                elif col in ['fast','total']:
                    u[col] = int(r[col])
                elif col in ['FOM','Z']:
                    u[col] = [float(x) for x in r[col].replace('(',' ').replace(')',' ').split()]
                elif col=='run':
                    u[col] = int(r[col].replace('run',''))
                elif col=='nFile':
                    u['sampleName'] = r[col].split('/')[-2]
                    u[col] = r[col]
                else:
                    u[col] = r[col]
            useful[rn] = u
        #print 'showLSQAResults.translate useful',useful
        self.useful = useful
        return
    def compactSampleNames(self,samples):
        '''
        return dict of compact sample names given list of full sample names
        '''
        cSN = {}
        for s in samples:
            c = s.replace('LiLS_','').replace('batch','B').replace('sample','S')
            cSN[s] = c
        return cSN
    def main(self):
        '''
        main routine

        '''
        self.reader()
        self.translate()
        useful = self.useful


        samples = []
        for rn in useful:
            s = useful[rn]['sampleName']
            if s not in samples: samples.append(s)
        cSN = self.compactSampleNames(samples)
        
        # plot Quantity, unc for all runs


        Quantities = ['Z','FOM','EperQ']

        fig, axes = plt.subplots(nrows=len(Quantities))
        for iQ,Q in enumerate(Quantities):
            ax = axes[iQ]
            tmi,tma = None,None
            for s in samples:

#                ic = ip = 1 + samples.index(s)
                ic,ip = self.colorAndPoint(samples.index(s))
                x,y,dx,dy = [],[],[],[]
                for rn in sorted( useful.keys() ):
                    u = useful[rn]
                    if u['sampleName']==s:
                        Z = u[Q]
                        if Q in ['FOM','Z']:
                            z,dz = Z
                        elif Q=='EperQ':
                            z,dz = Z,0.
                        else:
                            print 'showLSQAResults.main Unknown quantity',Q
                        y.append(z)
                        dy.append(dz)
                        T = u['run'] # integer version of labview time
                        if tmi==None: tmi,tma = T,T
                        tmi = min(T,tmi)
                        tma = max(T,tma)
                        t = self.cLT.convert((u['run']),fmt=None)
                        x.append(t)
                        dx.append(0.)
                x,y,dx,dy = numpy.array(x),numpy.array(y),numpy.array(dx),numpy.array(dy)
                ax.errorbar(x, y, color=self.colors[ic],marker=self.points[ip], yerr=dy, linestyle='empty',label=cSN[s])
            #plt.plot(x,y,yerr=dy)
            # limits
            dt = (tma-tmi)/20
            t1 = self.cLT.convert(tmi-dt,fmt=None)
            t2 = self.cLT.convert(tma+dt,fmt=None)
            ax.set_xlim((t1,t2))
            ax.set_ylabel(Q)
            ax.grid()
            if iQ==0:
                #print 'ax.legend loc next'
                ax.legend(loc='upper center')
                #print 'ax.legend bbox next'
                ncol = len(samples)/2+1
                #print 'ncol',ncol
                ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,  ncol=ncol, borderaxespad=0., numpoints=1,labelspacing=0.1,columnspacing=0.1,handletextpad=0.1) #mode="expand"

                
        #plt.title('Z vs time all samples')
        fmt = '%Y%m%d_%H%M%S_%f'
        cnow = datetime.datetime.now().strftime(fmt)
        fig.suptitle('PSD QA jobtime '+cnow,y=0.05)
        Show = not self.makeFigure
        if Show: 
            plt.show()
        else:
            pdf = self.figdir + cnow + '.pdf'
            fig.set_size_inches(8.5,11)
            plt.savefig(pdf)
            print 'showLSQAResults.main Wrote',pdf
 
        return

    
  
if __name__ == '__main__' :
    '''
    usage to create pdf: 
    $ python showLSQAResults.py FIGURE
    to show only:
    $ python showLSQAResults.py
    '''
    makeFigure = len(sys.argv)>1
    A = showLSQAResults(makeFigure)
    if 0: # test points, colors
        for i in range(100):
            ic,ip = A.colorAndPoint(i)
            print 'i,ic,ip',i,ic,ip
        sys.exit('--------------')
    A.main()
