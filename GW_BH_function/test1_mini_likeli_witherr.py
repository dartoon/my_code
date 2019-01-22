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
import matplotlib.pyplot as plt
from BH_mass_function import BH_mass_function
from cal_likelihood import likelihood as likelihood_i
likelihood = np.vectorize(likelihood_i)
from scipy.optimize import fmin
import time
#t1 = time.time()
#a, mbh_max, mbh_min = 2.35, 80., 5.
#vol = 1000
#BHBH_mass = BH_mass_function(a = a, mbh_max=mbh_max, mbh_min=mbh_min, vol= vol)
#m_noise_level = 0.2
#def point_Chisq(para, m1_obs,m1_sig):
#    a,mbh_max, mbh_min  = para
#    if 1.1 < a < 3 and 50 < mbh_max < 100 and 2 < mbh_min < 8:
#        poss_m1_i = likelihood(m1_obs,m1_sig, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
#        chisq = -np.sum(np.log(poss_m1_i))
##        print a, mbh_max, mbh_min
#        return chisq
#    else:
#        return np.inf
#para_ini = [2.35,80, 5.]
#best_p = []
#para_result =  open('test1_para.txt','r+') 
#rounds = 500
#for loop in range(rounds):
#    m1 = BHBH_mass.gen_dm()
#    m1_sig= m1 * m_noise_level
#    m1_obs = m1 + np.random.normal(0, m1_sig, size=m1_sig.shape)
#    mini=fmin(point_Chisq,para_ini,maxiter=1000, args=(m1_obs, m1_sig))
#    if loop > 0:
#        para_result = open('test1_para.txt','r+')
#        para_result.read()
#    para_result.write(repr(mini)+"\n")
#    para_result.close()
#    best_p.append(mini)
#    t2 = time.time()
#    time_sp = t2-t1
#    time_ave = (t2-t1)/(loop+1)
#    time_total = time_ave * rounds
#    t_left = time_total - time_sp
#    print mini
#    print "loop:", loop
#    print "Finish percent:",round(time_sp/time_total*100,2),"%" ,"total time needed :", round(time_total/60,2), "mins", "time_left", round(t_left/60,2), 'mins'

import corner
with open('test1_maxbin161_un20.txt') as f:
        content = f.readlines()
lines = [x.strip() for x in content] 
numbers = np.asarray([np.float_(lines[i].split(' ')) for i in range(len(lines))])
samples = numbers
fig = corner.corner(samples, labels=["$alpha$","$Mbh_{max}$", "$Mbh_{min}$"],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},truths=[2.35,80,5],
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f')
#fig.savefig("test1_norm_likeli.pdf")
#####
plt.show()  