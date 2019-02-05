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

plt.subplots(figsize=(12, 7))

def poss_gaussian(m1, mu, sigma):
    poss =  1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(m1 - mu)**2/(2*sigma**2))
    return poss
'''
alpha = -1
x=np.linspace(5, 80)
print x
y=x**(alpha)
print y
plt.plot(x,y, 'b')


x1, x2 = 5**(alpha), 80**(alpha)
y1, y2=np.linspace(0, x1), np.linspace(0,x2)
plt.plot(5 + y1*0,y1, 'b')
plt.plot(80 + y2*0,y2, 'b')

mu, sigma= 25, 5
x3 = np.linspace(mu-3*sigma, mu+3*sigma)
y3 = 0.5 * 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(x3 - mu)**2/(2*sigma**2))
y3 = 0.5 * poss_gaussian(x3,mu, sigma)
plt.plot(x3, y3, 'g')
import scipy
d = scipy.zeros(len(y3))

plt.xlabel("$m$",fontsize=15)
plt.ylabel("P($m$)",fontsize=15)
#plt.yticks([])
plt.tick_params(labelsize=15)
plt.show()
'''

def poss_ln_gaussian(m1, mu, sigma):
    poss =  1/(m1 * sigma * np.sqrt(2 * np.pi)) * np.exp(-(np.log(m1) - mu)**2/(2*sigma**2))
    return poss
#x_top = 40
mu, sigma = np.log(20), 0.4          #mu is the median point where CDF: Fx = 1/2. sigma 40% level,
x4 = np.linspace(0.001,50,500)
y4 = poss_ln_gaussian(x4, mu, sigma)
plt.plot(x4,y4, 'b')

expected_m1 = np.exp(mu+sigma**2/2.) #The mean point (expected value)
mode_m1 = np.exp(mu-sigma**2)   # The mode is the point of global maximum of the probability density function. In particular, it solves the equation (f)`=0
y_l=np.linspace(0, 2)

mu_star = np.exp(mu) 
sig_star = np.exp(sigma)

plt.plot(mu_star*sig_star + y_l*0,y_l, 'g--')
plt.plot(mu_star/sig_star + y_l*0,y_l, 'g--')
plt.plot(expected_m1 + y_l*0,y_l, 'b')     #Blue, the expected mean.
plt.plot(mu_star + y_l*0,y_l, 'g')         #Green, the np.exp(mu) is the MEDIAN point.
plt.plot(mode_m1 + y_l*0,y_l, 'r')       #Red, the mode
plt.ylim(0,y4.max()*1.2)
plt.xlabel("$m$",fontsize=15)
plt.ylabel("P($m$)",fontsize=15)
#plt.yticks([])
plt.tick_params(labelsize=15)
#plt.show()

samples = np.random.lognormal(np.log(20), 0.4,100000)
plt.hist(samples,bins=100,normed=True)
plt.show()

print mu_star/sig_star,mu_star,mu_star*sig_star
