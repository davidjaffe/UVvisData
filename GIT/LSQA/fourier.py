#!/usr/bin/env python
'''
fourier analysis of exponentials and such.
to analysis LS decay time
20170127
'''
import sys
# the ONETON dir contains wfanal.py 
sys.path.append('/Users/djaffe/work/GIT/ONETON')
# the PROSPECT dir contains get_filepaths.py
from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt

class fourier():
    def __init__(self):
        return
    def test1(self):
    # Number of sample points
        N = 600*2
        # sample spacing
        T = 1.0 / 800.0
        x = np.linspace(0.0, N*T, N)
        y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
        y1 = np.exp(-x/2.)
        y2 =  0.5*np.exp(-x/8.)
        y = y1+y2
        plt.plot(x,y,'r.', x,y1,'bx', x,y2,'g+')
        plt.show()
        yf1 = fft(y1)
        yf2 = fft(y2)
        yf = fft(y)
        xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
        plt.semilogy(xf, 2.0/N * np.abs(yf[0:N/2]),'r.', xf, 2.0/N * np.abs(yf1[0:N/2]), 'bx',
                     xf, 2.0/N * np.abs(yf2[0:N/2]),'g+')
        plt.grid()
        plt.show()
        return
if __name__ == '__main__':
    F = fourier()
    F.test1()
