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
from matplotlib import dates
from matplotlib.ticker import AutoMinorLocator
import convertLabviewTime

class showLSQAResults():
    def __init__(self,mf):
        self.makeFigure = mf
        self.cLT = convertLabviewTime.convertLabviewTime()
        
        self.colors = {0:'k', 1:'b',2:'g',3:'r',4:'c',5:'m',6:'y'}
        self.points = {}
        for i,p in enumerate(['s','o','v','^','<','>','*','D','p','x','+']):
            self.points[i] = p
        self.lencolors = len(self.colors)
        self.lenpoints = len(self.points)

        self.resFile= 'Samples/resultsSummaryFile.txt'
        self.curFile= 'Samples/currentResultsSummaryFile.txt'
        
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
        #print 'showLSQAResults.colorAndPoint i,j,ip,lenpoints,ic,lencolors',i,j,ip,self.lenpoints,ic,self.lencolor
        return ic,ip
    def reader(self):
        '''
        parse results file
        '''
        fn = self.resFile
        outfn = self.curFile
        f = open(fn,'r')
        outf = open(outfn,'w')

        self.cols = cols = ['run','EperQ','FOM','Z','fast','total','gDate','gTime','nDate','nTime','gFile','nFile','jobTime']
        self.results = results = {}

        current = {'1A':None}
        line0 = None

        cjT = 'jobTime'
        N = 0
        
        for line in f:
            if current['1A'] is None: current['1A'] = line
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
                current[ s[0] ] = line
        f.close()
        print 'showLSQAResults.reader Read',N,'lines. Found',len(results),'results from',fn
        for key in sorted( current.keys()) :
            #print 'key',key,'current',current[key]
            outf.write(current[key])
        print 'showLSQAResults.reader Wrote',len(current),'lines to',outfn
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
            c = s.replace('LiLS_','').replace('batch','B').replace('sample','S').replace('_repeat','r').replace('drum','D')
            cSN[s] = c
        return cSN
    def MYindex(self,x):
        ''' used to sort sample names '''
        words = ['batch','drum']
        for w in words:
            if w in x: return x.index(w)+len(w)
        return None
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
        # make P50-1 the first sample and sort the rest of the samples
        cP50 =  'P50-1'
        samples.remove(cP50)
#        samples.sort(key=lambda x :int(x[x.index('batch')+len('batch'):x.index('_s')].replace('_repeat',''))*100+int(x[x.index('sample')+len('sample'):].replace('_repeat',''))  )
        samples.sort(key=lambda x :int(x[self.MYindex(x):x.index('_s')].replace('_repeat',''))*100+int(x[x.index('sample')+len('sample'):].replace('_repeat',''))  )
        samples.insert(0,cP50)
        cSN = self.compactSampleNames(samples)

        print 'showLSQAResults.main samples',samples
        
        # plot Quantity, unc for all runs


        Quantities = ['Z','FOM','EperQ']
        yAxLimits  = [ [94., 110.], [1.25, 1.60], [22.,30.] ]
        factors    = [1000., 1., 1.]
        histlabels = []
        for f,Q in zip(factors,Quantities):
            lab = Q
            if f!=1.:
                lab = str(int(f))
                lab += r'$\times$'+Q
            histlabels.append(lab)
#xlabel(r'\textbf{time} (s)')                

        showHist = True
        
        
        if showHist:
            fig = plt.figure()
            axes = []
            NCOL = len(Quantities)
            NROW = len(Quantities)+1
            for i in range(len(Quantities)):
                if i==0:
                    ax = plt.subplot2grid( (NROW,NCOL), (i,0), colspan=NCOL)
                else:
                    ax = plt.subplot2grid( (NROW,NCOL), (i,0), colspan=NCOL, sharex=axes[0])
                    plt.setp(axes[i-1].get_xticklabels(), visible=False)

                    
                axes.append(ax)
            Haxes = []
            for i in range(len(Quantities)):
                ax = plt.subplot2grid( (NROW,NCOL), (NROW-1,i))
                Haxes.append(ax)
            #Hfig,Haxes= plt.subplots(ncols=len(Quantities))
            #fig, axes = plt.subplots(nrows=len(Quantities))
        else:
            fig, axes = plt.subplots(nrows=len(Quantities))
        for iQ,Q in enumerate(Quantities):
            f = factors[iQ]
            ax = axes[iQ]
            if showHist:    Hax=Haxes[iQ]
            y50, ybs = [],[]
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
                            z,dz = f*z,f*dz
                        elif Q=='EperQ':
                            z,dz = f*Z,f*0.
                        else:
                            print 'showLSQAResults.main Unknown quantity',Q
                        y.append(z)
                        if s==cP50:
                            y50.append(z)
                        else:
                            ybs.append(z)
                        dy.append(dz)
                        T = u['run'] # integer version of labview time
                        if tmi==None: tmi,tma = T,T
                        tmi = min(T,tmi)
                        tma = max(T,tma)
                        t = self.cLT.convert((u['run']),fmt=None)
                        x.append(t)
                        dx.append(0.)
                x,y,dx,dy = numpy.array(x),numpy.array(y),numpy.array(dx),numpy.array(dy)
                yLimits = yAxLimits[iQ]
                if min(y-dy)<yLimits[0] or max(y+dy)>yLimits[1]:
                    print 'showLSQAResults.main ERROR',Q,'yLimits',yLimits,'ymin,ymax',min(y-dy),max(y+dy)
                    sys.exit('showLSQAResults.main ERROR Fix y limits for ' + Q)
                ax.errorbar(x, y, color=self.colors[ic],marker=self.points[ip], yerr=dy, linestyle='empty',label=cSN[s])
                ax.set_ylim(yLimits)

            # histograms
            if showHist:
                YY = numpy.concatenate((y50,ybs))
                ymi,yma = min(YY)-(max(YY)-min(YY))/10.,max(YY)+(max(YY)-min(YY))/10.
                bins = numpy.linspace(ymi,yma,20)

                Hax.hist([y50,ybs],bins,histtype='stepfilled',stacked=True,
                         label=[cP50,'batches'],color=['black','red'])

                Hax.set_xlabel(histlabels[iQ])
                Hax.grid()
                ylo,yhi = Hax.get_ylim()
                yhi = yhi + (yhi-ylo)/10.
                Hax.set_ylim( (ylo,yhi) )
                xlo,xhi = Hax.get_xlim()
            #Hax.xaxis.set_ticks(numpy.linspace(xlo,xhi,4))
            
            #plt.plot(x,y,yerr=dy)
            # limits
            dt = (tma-tmi)/20
            t1 = self.cLT.convert(tmi-dt,fmt=None)
            t2 = self.cLT.convert(tma+dt,fmt=None)
            ax.set_xlim((t1,t2))
            hfmt = dates.DateFormatter('%y%m%d')
            ax.xaxis.set_major_formatter(hfmt)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_ylabel(histlabels[iQ])
            ax.grid()
            if iQ==0:
                ax.legend(loc='upper center')
                nrow = 2
                if len(samples)>10: nrow +=1 
                if len(samples)>25: nrow +=1
                if len(samples)>40: nrow +=1
#                if len(samples)>50: nrow +=1
                ncol = len(samples)/nrow+1
                ax.legend(bbox_to_anchor=(0.-.05, 1.0002, 1.+.05, .102+.01), loc=3,  prop={'size':8}, #9
                          ncol=ncol, borderaxespad=0., numpoints=1,labelspacing=0.1,
                          columnspacing=0.1,handletextpad=0.1) #mode="expand"
                
            if iQ==len(Quantities)-1:
                if showHist: Hax.legend(loc='upper right')

        plt.tight_layout(h_pad=-0.5,w_pad=-0.2,pad=2.0)
#        plt.xticks(rotation=30)#'vertical')
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
