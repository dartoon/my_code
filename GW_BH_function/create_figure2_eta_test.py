#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:10:46 2019

@author: Dartoon

Test the distribution of 
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from BH_mass_function import gene_BHBH, dl, solve_z
import time
import pickle

from cal_likelihood import poss_ln_gaussian

solve_z = np.vectorize(solve_z)
#a0, a1, mbh_max, mbh_min = 2.35, 0.1, 80., 5.
a0, a1, mbh_max, mbh_min = 2.4, 0.7, 50., 5.   #mode1

seed = 2
#filename = 'test3_seed{0}_data'.format(seed)
filename = 'sim_a0_{0}_a1_{1}_max_{2}'.format(a0, a1, mbh_max)
if_file = glob.glob(filename)  
if if_file==[]:
    np.random.seed(seed)
    test = gene_BHBH(h0=70)
    event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=[a0, a1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode1')
    dl_zs = dl(zs_detected)
    sim_data = [event_rate, zs_detected, masses, rhos_detected, dl_zs]
    pickle.dump(sim_data, open(filename, 'wb'))
else:
    event_rate, zs_detected, masses, rhos_detected, dl_zs=pickle.load(open(filename,'rb'))

np.random.seed(seed)
#event_rate0, zs_detected0, masses0, rhos_detected0 = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'const')
#event_rate1, zs_detected1, masses1, rhos_detected1 = test.mc_year_rate(a=[2.35,0.1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode0')
#event_rate2, zs_detected2, masses2, rhos_detected2 = test.mc_year_rate(a=[a0, a1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode1')

zs_all, chirp_mass_all, m1_all, m2_all, lumi_dis_all = zs_detected, masses[:,0], masses[:,1], masses[:,2], dl_zs                                                      
#
##%% Test if the data slope follows the redshift trend.
#plt.hist(m1_all[zs_all==2.141], m1_all[zs_all==5.531],normed=True, log=True)
#plt.show()

from cal_likelihood import random_Theta  #The prior given dl and chirpmass, no errorbar
thetas = random_Theta()
#%%
from scipy import interpolate
from cal_likelihood import cov_twofun
from scipy.optimize import fmin
m_noise_level = 0.20

para_ini = [a0, a1, mbh_max, mbh_min]
t1 = time.time()

index = np.arange(len(m1_all))
from cal_likelihood import fac_s_eff_v
rounds = 250
count = 0

part = 0

np.random.seed(2)
idx = np.random.choice(index, 1)[0]
m1 = np.ones(1000)*m1_all[idx]
dl = np.ones(1000)*lumi_dis_all[idx]
mass_Chirp = chirp_mass_all[idx]
dl_noised =  np.random.lognormal(np.log(dl), 0.35, size=dl.shape)
mass_Chirp_noised = np.random.lognormal(np.log(mass_Chirp), 0.17, size=mass_Chirp.shape)
m1_mu = np.log(m1)   # 0_med np.log(mu_star) = mu 
m1_sigstar= np.exp(m_noise_level)  #Not useful in the generate generation, but useful in understand the upper lower level.
m1_obs = np.random.lognormal(m1_mu, m_noise_level, size=m1.shape)  #Generating for the mu as med, 
m1_sig_fake = m1_obs * m_noise_level      #The fake "sigma", (m1_obs * m_noise_level) and (m1 * m_noise_level)
prior = fac_s_eff_v(dl=dl, mass_Chirp=mass_Chirp, thetas=thetas)[0]
sf_factor = 1/prior
prior_noi = fac_s_eff_v(dl=dl_noised, mass_Chirp=mass_Chirp_noised, thetas=thetas)
prior_noi[prior_noi==0] = 0.001
sf_factor_noi = 1/prior_noi


#%%
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

fig, ax = plt.subplots(figsize=(12, 10))

a, b, c =plt.hist(sf_factor_noi,bins=60, histtype=u'step', linewidth = 2, normed=True, 
         label= r"1$/{\rm \eta}$ histogram.")

x = sf_factor
y = np.linspace(0,100,200)
#plt.plot(y*0+x , y, linewidth = 4,color='orange',linestyle=('dashed'), label = "True")

x_infer = np.median(sf_factor_noi)
x_line = np.linspace(x_infer/3,x_infer*4,200)
sf_sigma = np.log(x_infer)/3.
x_guess = x_infer/np.exp(sf_sigma**2/2)

x = x_guess
#plt.plot(y*0+x , y, linewidth = 4,color='blue',linestyle=('dashed'), label = "Corrected")

x = x_infer
#plt.plot(y*0+x , y, linewidth = 4,color='red',linestyle=('dashed'), label = "median")


y_line = poss_ln_gaussian(x_line, np.log(x_infer), sf_sigma)
plt.plot(x_line , y_line, linewidth = 4,color='orange',linestyle=('dashed'),
         label= r"The Log-Normal distribution with $\sigma = -\log(\eta_{median})/3$.")

plt.xlim(0.5,8)
plt.ylim(0,a.max()*1.2)
plt.xlabel(r"Probability distribution of 1$/{\rm \eta}$",fontsize=27)
plt.ylabel("PDF",fontsize=37)
plt.tick_params(labelsize=30)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], prop={'size':18})
#print x_guess, sf_factor
plt.savefig("eta_distribution.pdf")
plt.show()
#sf_sigma = np.log(sf_factor)/3