#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Test if whether the likelihood would recover the para
if all the m1 are measureable.
With error bar
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
from BH_mass_function import BH_mass_function
#from cal_likelihood import likelihood_lognorm as likelihood_i
likelihood = np.vectorize(likelihood_i)
from scipy.optimize import fmin
import time
t1 = time.time()
a, mbh_max, mbh_min = 2.35, 80., 5.
vol = 1000
BHBH_mass = BH_mass_function(a = a, mbh_max=mbh_max, mbh_min=mbh_min, vol= vol)
m_noise_level = 0.20     #The actually sigma_star is np.exp(0.2)

use_method='med' #raw_input('Consider random.lognormal generate value is med or exp?:\n')

if use_method == 'med':
    filename = 'test1_conv_lognorm_{0}.txt'.format(int(m_noise_level*100))   #obsun means the std are from the obs data rather than true data
#elif use_method == 'exp':
#    filename = 'test2_1_trans_exp_obsun{0}.txt'.format(int(m_noise_level*100))

def point_Chisq(para, m1_obs,m1_sig):
    a,mbh_max,mbh_min  = para
    if 1.1 < a < 3 and 50 < mbh_max < 100 and 2 < mbh_min < 8:
        x = np.logspace(np.log10(m1_obs.min()/np.exp(m_noise_level)**5/1.1), np.log10(m1_obs.max()*np.exp(m_noise_level)**5*1.1),300)
        y = cov_twofun(x, a=a, mbh_max=mbh_max, mbh_min=mbh_min,sigma=m_noise_level)
        f = interpolate.interp1d(x, y)
#        poss_m1_i = likelihood(m1_obs,m1_sig, a=a, mbh_max=mbh_max, mbh_min=mbh_min, use_method='med',f=f)
        poss_m1_i = f(m1_obs)
        chisq = -np.sum(np.log(poss_m1_i))
#        print a, mbh_max, mbh_min, chisq
        return chisq
    else:
        return np.inf
#print point_Chisq([2.35,80, 5.], m1_obs, m1_sig_ratio)

def each_chisq(para, m1_obs,m1_sig):
    a,mbh_max,mbh_min  = para
    if 1.1 < a < 3 and 50 < mbh_max < 100 and 2 < mbh_min < 8:
        x = np.logspace(np.log10(m1_obs.min()/np.exp(m_noise_level)**5/1.1), np.log10(m1_obs.max()*np.exp(m_noise_level)**5*1.1),300)
        y = cov_twofun(x, a=a, mbh_max=mbh_max, mbh_min=mbh_min,sigma=m_noise_level)
        f = interpolate.interp1d(x, y)
#        poss_m1_i = likelihood(m1_obs,m1_sig, a=a, mbh_max=mbh_max, mbh_min=mbh_min, use_method='med', f=f)
        poss_m1_i = f(m1_obs)
#        print a, mbh_max, mbh_min, chisq
        return poss_m1_i
    else:
        return np.inf
    
from scipy import interpolate
from cal_likelihood import cov_twofun    

para_ini = [2.35,80.,5]
best_p = []
if_file = glob.glob(filename)
if if_file == []:
    para_result =  open(filename,'w') 
else:
    para_result =  open(filename,'r+') 
rounds = 1
for loop in range(rounds):
    m1 = BHBH_mass.gen_dm()   #The expected position
    if use_method == 'med':
        m1_mu = np.log(m1)   # 0_med np.log(mu_star) = mu 
    elif use_method == 'exp':
        m1_mu = np.log(m1) - m_noise_level**2/2.   # 1_trans_exp np.log(mu) = np.log(exp-sig**2/2)
    
    m1_sigstar= np.exp(m_noise_level)  #Not useful in the generate generation, but useful in understand the upper lower level.
    m1_obs = np.random.lognormal(m1_mu, m_noise_level, size=m1.shape)  #Generating for the mu as med, 
    m1_sig_fake = m1_obs * m_noise_level      #The fake "sigma", (m1_obs * m_noise_level) and (m1 * m_noise_level)
    print "m1_obs.min(), m1_obs.max():",m1_obs.min(), m1_obs.max()
    mini=fmin(point_Chisq,para_ini,maxiter=1000, args=(m1_obs, m1_sig_fake))
    print mini

m1_put = m1_obs
dis_likeli_true = each_chisq(para_ini, m1_put, m1_put * m_noise_level)
dis_likeli_mini = each_chisq(mini, m1_put, m1_put * m_noise_level)

plt.figure(figsize=(10, 8))
plt.plot(m1_put, dis_likeli_true, 'r.')
plt.plot(m1_put, dis_likeli_mini, '.')

plt.xlim([0,20])
plt.show()
print dis_likeli_true.sum()
print dis_likeli_mini.sum()
#print m1_put.min() * np.exp(m_noise_level)**2