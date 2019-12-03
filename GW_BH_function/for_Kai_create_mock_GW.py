#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Test if whether the likelihood would recover the para
if all the m1 are measureable.
With error bar
"""


#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Generate the simulating data
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
test = gene_BHBH(h0=70)
test.rho0 = 0

event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min)
dl_zs = dl(zs_detected)
sim_data = [event_rate, zs_detected, masses, rhos_detected, dl_zs]

zs_all, chirp_mass_all, m1_all, m2_all, lumi_dis_all = zs_detected, masses[:,0], masses[:,1], masses[:,2], dl_zs                                                      
#==============================================================================
#Setting the noise level 
#==============================================================================
from scipy import interpolate
from cal_likelihood import cov_twofun
from scipy.optimize import fmin
from cal_likelihood import random_Theta  #The prior given dl and chirpmass, no errorbar
#def select_effect(m1_obs, fname = 'select_effect_MBHmin5_cov_lognorm0.2'):
#    x, y  = pickle.load(open(fname,'rb'))
#    f = interpolate.interp1d(x, y)
#    prior = f(m1_obs)
#    return prior

thetas = random_Theta()

rounds = 200
index = np.arange(len(m1_all))

# =============================================================================
# Infer the information of the sigma.
# =============================================================================
#lines = np.loadtxt('inveres_prior_scatter.txt')
#invprior_true, invprior_mean, invprior_median, m1_list, m1_obs_list, dl_list, massChirp_list = [lines[:,i] for i in range(len(lines.T))]
#solve_z = np.vectorize(solve_z)
#z = solve_z(np.array(dl_list))
#x = np.log(np.array(dl_list)*(1/np.array(massChirp_list)**(5/6.))/(1+z)**(5/6.) )
#sigma = np.sqrt(2*np.log(np.array(invprior_mean)/np.array(invprior_median)))
#y = sigma
#y = y[x!=x.min()] #Delete the mini point
#x = x[x!=x.min()]
#v_min, v_max = x.min(), x.max()
#fit1d = np.poly1d(np.polyfit(x, y, 30))

idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
dl = lumi_dis_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
dl_noised =  np.random.lognormal(np.log(dl), 0.35, size=dl.shape)

save_filename = 'for_kai_simulated_data.txt'
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR \n") 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close()  

