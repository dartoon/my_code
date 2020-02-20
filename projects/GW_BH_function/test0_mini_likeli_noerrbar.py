#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Test if whether the likelihood would recover the para
if all the m1 are measureable.
Without error bar
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from BH_mass_function import BH_mass_function
from cal_likelihood import poss_mass_fun
import scipy.optimize as op
a, mbh_max, mbh_min = 2.35, 80., 5.
vol = 1000
BHBH_mass = BH_mass_function(a = a, mbh_max=mbh_max, mbh_min=mbh_min, vol= vol)
m1 = BHBH_mass.gen_dm()

#poss_mass_fun_test = poss_mass_fun(m1, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
def point_Chisq(para, m1):
    a = para
    mbh_max, mbh_min = [80., 5.]
    if 1.1 < a < 3:
        poss_m1_i = poss_mass_fun(m1, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
#        print a, mbh_max, mbh_min
        chisq = -np.sum(np.log(poss_m1_i))
        return chisq
    else:
        return np.inf
para_ini = 2.35
nll = lambda *args: point_Chisq(*args)
bnds = (0, None)
best_p = []
rounds = 1000
for loop in range(rounds):
    m1_loop = BHBH_mass.gen_dm()
    result = op.minimize(nll, para_ini, method='SLSQP', args=(m1_loop))
    best_p.append(result["x"][0])
    if loop/10 > (loop-1)/10 :
        print loop
l =np.percentile(best_p,16,axis=0)
m =np.percentile(best_p,50,axis=0)
h =np.percentile(best_p,84,axis=0)
y = np.linspace(0,rounds,100)
x_true = y*0+ 2.35
x_l = y*0+ l
x_m = y*0 + m
x_h = y*0 + h
plt.figure(figsize=(6,5))
plt.plot(x_true, y, 'red',linewidth=2.0)
plt.plot(x_l, y, 'k--',linewidth=2.0)
plt.plot(x_m, y, 'blue',linewidth=2.0)
plt.plot(x_h, y, 'k--',linewidth=2.0)
a,b,_ = plt.hist(best_p,bins=30)
plt.axis([2.2,2.5,0,a.max()*1.2])
plt.xlabel('$alpha$',fontsize=15)
plt.ylabel('PDF',fontsize=15)
plt.tick_params(labelsize=10)
plt.yticks([])
#plt.savefig('test0.pdf')
plt.show()

