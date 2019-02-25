#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 17:00:17 2019

@author: dxh

Test prior distribute
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from BH_mass_function import gene_BHBH, dl, solve_z
from cal_likelihood import fac_s_eff_v
import pickle
import glob
import random
import scipy.optimize as op
import time

a, mbh_max, mbh_min = 2.35, 80., 5.
filename = 'sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
if_file = glob.glob(filename)  
test = gene_BHBH(h0=70)
if if_file==[]:
    event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min)
    dl_zs = dl(zs_detected)
    sim_data = [event_rate, zs_detected, masses, rhos_detected, dl_zs]
    pickle.dump(sim_data, open(filename, 'wb'))
else:
    event_rate, zs_detected, masses, rhos_detected, dl_zs=pickle.load(open(filename,'rb'))
    
zs_all, chirp_mass_all, m1_all, m2_all, lumi_dis_all = zs_detected, masses[:,0], masses[:,1], masses[:,2], dl_zs                                                      

from scipy import interpolate
from cal_likelihood import cov_twofun
from scipy.optimize import fmin
from cal_likelihood import random_Theta  #The prior given dl and chirpmass, no errorbar

def select_effect(m1_obs, fname = 'select_effect_MBHmin5_cov_lognorm0.2'):
    x, y  = pickle.load(open(fname,'rb'))
    f = interpolate.interp1d(x, y)
    prior = f(m1_obs)
    return prior

#thetas = random_Theta()

#filename = 'inveres_prior_scatter.txt'
#if_file = glob.glob(filename)
#if if_file == []:
#    scatter_result =  open(filename,'w') 
#else:
#    scatter_result =  open(filename,'r+') 
#invprior_true, invprior_mean, invprior_median, m1_list, m1_obs_list, dl_list, massChirp_list = [], [], [], [], [], [], []
#index = np.arange(len(m1_all))
#for i in range(1):
#    idx = random.sample(index, 1)
#    m1 = m1_all[idx]
#    dl = lumi_dis_all[idx]
#    mass_Chirp = chirp_mass_all[idx]
#    prior_true = fac_s_eff_v(dl=dl, mass_Chirp=mass_Chirp, thetas=thetas)
#    m1_obs_i = np.random.lognormal(np.log(m1), 0.2)
##    prior_m1_obs = select_effect(m1_obs)
#    dl_noised =  np.random.lognormal(np.log(dl), 0.35, size=5000)
#    mass_Chirp_noised = np.random.lognormal(np.log(mass_Chirp), 0.17, size=5000)
#    prior_noised = fac_s_eff_v(dl=dl_noised,
#                               mass_Chirp=mass_Chirp_noised, thetas=thetas)
##    prior_noised = prior_noised * 1.3
#    hist = plt.hist(prior_noised, bins=50)
#    y = np.linspace(0, hist[0].max()*1.2) 
#    x = y*0 + prior_true
#    plt.plot(x,y,'r')
#    x_mean = y*0 + np.mean(prior_noised)
#    plt.plot(x_mean,y,'y')
#    #x1 = y*0 + prior_m1_obs
#    #plt.plot(x1,y,'b')
#    plt.show()
#    prior_noised = prior_noised[prior_noised!=0]
#    inves = 1/prior_noised
#    inves = inves[inves<100]
#    hist = plt.hist(inves, bins=50)
#    y = np.linspace(0, hist[0].max()*1.2) 
#    x0 = y*0 + 1/prior_true
#    plt.plot(x0,y,'r')
#    x_mean = y*0 + np.mean(inves)
#    plt.plot(x_mean,y,'y')
#    #x1 = y*0 + 1/prior_m1_obs
#    #plt.plot(x1,y,'b')
#    plt.show()
#    invprior_true.append(1/float(prior_true))
#    invprior_mean.append(np.mean(1/prior_noised))
#    invprior_median.append(np.median(1/prior_noised))
#    m1_list.append(float(m1))
#    m1_obs_list.append(float(m1_obs_i))
#    dl_list.append(float(dl))
#    massChirp_list.append(float(mass_Chirp))
##    if i > 0:
##        scatter_result = open(filename,'r+')
##        scatter_result.read()
##    scatter_result.write(repr(invprior_true[i])+' '+repr(invprior_mean[i])+' '+\
##                         repr(invprior_median[i])+' '+ repr(m1_list[i]) +' '+ repr(m1_obs_list[i]) +' '+\
##                         repr(dl_list[i]) +' ' + repr(massChirp_list[i]) +"\n")
##    scatter_result.close()
#    if i/5 > (i-1)/5: 
#        print "invprior_True, invprior_mean, invprior_median, m1_obs"
#        print i,':', invprior_true[i], invprior_median[i], invprior_mean[i], m1_obs_list[i]
#        print "sigma", np.sqrt(2*np.log(np.array(invprior_mean[i])/np.array(invprior_median[i])))
#    
# =============================================================================
# ### Sigma of the Lognorm  
# =============================================================================
lines = np.loadtxt('inveres_prior_scatter.txt')
invprior_true, invprior_mean, invprior_median, m1_list, m1_obs_list, dl_list, massChirp_list = [lines[:,i] for i in range(len(lines.T))]
#plt.plot(np.log(np.array(m1_obs_list)), np.sqrt(2*np.log(np.array(invprior_mean)/np.array(invprior))),'.')
#plt.show()
#plt.plot(np.log(1/np.array(dl_list)), np.sqrt(2*np.log(np.array(invprior_mean)/np.array(invprior))),'.')
#plt.show()
#plt.plot(np.log(np.array(massChirp_list)**(5/6.)), np.sqrt(2*np.log(np.array(invprior_mean)/np.array(invprior))),'.')
#plt.show()

sigma = np.sqrt(2*np.log(np.array(invprior_mean)/np.array(invprior_median)))

#plt.scatter(np.array(massChirp_list), sigma,c=np.array(dl_list))
#plt.show()
solve_z = np.vectorize(solve_z)
z = solve_z(np.array(dl_list))
plt.plot((np.array(dl_list)*(1/np.array(massChirp_list)**(5/6.))/(1+z)**(5/6.) ), sigma,'.')
plt.show()
#plt.plot(np.array(dl_list)*(1/np.array(massChirp_list)**(5/6.)), sigma,'.')
#plt.show()
