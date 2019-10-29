#!/usr/bin/env python
'''
investigate chi2 behavior for pmt calibration algorithm
20191026
'''
import math
import sys,os,shutil

import datetime
import numpy
from scipy.stats import poisson,norm,relfreq
#import copy

import random

import re
import glob # used in __init__

import matplotlib.pyplot as plt
import math
from scipy.stats import betaprime

a,b = 2.,3.
bins = range(100)
for b in [2.2]:
    fig, axs = plt.subplots(4)
    a = 10.
    for i in range(4):
        y = betaprime.rvs(a,b,size=10000)
        lab = 'a={0} b={1} mean={2:.2f} rms={3:.2f}'.format(a,b,y.mean(),math.sqrt(y.var()))
        x = axs[i].hist( betaprime.rvs(a,b,size=10000), bins,label=lab )
        a += 20./3.
        axs[i].legend(loc='best')
    plt.show()
