#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:39:30 2019

@author: Dartoon

To recover the plot for the omh fig1
"""
import numpy as np
import matplotlib.pyplot as plt
import operator as op
from functools import reduce

z = np.array([0.07 , 0.1  , 0.12 , 0.17 , 0.179, 0.199, 0.2  , 0.27 , 0.28 ,
       0.35 , 0.352, 0.4  , 0.44 , 0.48 , 0.57 , 0.593, 0.6  , 0.68 ,
       0.73 , 0.781, 0.875, 0.88 , 0.9  , 1.037, 1.3  , 1.43 , 1.53 ,
       1.75 , 2.34 ])
h = np.array([0.69 , 0.69 , 0.686, 0.83 , 0.75 , 0.75 , 0.729, 0.77 , 0.888,
       0.827, 0.83 , 0.95 , 0.826, 0.97 , 0.929, 1.04 , 0.879, 0.92 ,
       0.973, 1.05 , 1.25 , 0.9  , 1.17 , 1.54 , 1.68 , 1.77 , 1.4  ,
       2.02 , 2.22 ])


Omh_list = []
for i in range(len(h)):
    for j in range(i+1, len(h)):
        omh_ij = (h[i]**2 - h[j]**2)/((1+z[i])**3 - (1+z[j])**3)
        Omh_list.append(omh_ij)

Omh = np.asanyarray(Omh_list)
        
def nchoosek(n, r):
    '''
    See https://www.mathworks.com/help/matlab/ref/nchoosek.html
    '''
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom


#The PDF of each calcualted Omh, assuming as binomial distribution.
N = len(Omh)      
p_i = []
for i in range(1,N+1):
    p_i.append((2**(-N)*nchoosek(N,i)))
p_i = np.asarray(p_i)

plt.plot(Omh[Omh.argsort()],p_i)
plt.xlim(0.12,0.19)
plt.ylim(0.,0.045)
plt.show()