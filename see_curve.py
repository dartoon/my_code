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

x=np.linspace(0, 4)
print x
y=5.*x*(4.-x)**3/256
print y

plt.plot(x,y)
plt.xlabel("$\Theta$")
plt.ylabel("P($\Theta$)")
plt.show()
