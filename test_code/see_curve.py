#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:55:30 2017

@author: Dartoon
"""

'''
To see how y change with x
'''
import numpy as np
import matplotlib.pyplot as plt

#x=np.linspace(0, 1)
##print x
##y=5.*x*(4.-x)**3/256
##print y
#def gaussian(m1, mu, sigma):
#    poss =  1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(m1 - mu)**2/(2*sigma**2))
#    return poss
#
#y = gaussian(x, 0.5,0.15)
#
#plt.plot(1/x,y)
#plt.xlabel("$1/x$")
#plt.ylabel("Gaussain")
#plt.xlim([0,5])
#plt.show()

x=np.linspace(0, 17)
y = 2.35 + 0.7 * x/(1+x)
plt.plot(x,y)
plt.show()


x=np.linspace(0, 17)
y = 2.35 + 0.04*x
plt.plot(x,y)
plt.show()