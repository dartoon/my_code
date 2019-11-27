#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 23:15:38 2019

@author: Dartoon

Test 3, considering the alpha is evolving with redshift.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from BH_mass_function import gene_BHBH, dl, solve_z
import random
import time
solve_z = np.vectorize(solve_z)
#a0, a1, mbh_max, mbh_min = 2.35, 0.1, 80., 5.
a0, a1, mbh_max, mbh_min = 2.35, 0.7, 80., 5.   #mode1

#filename = 'test3_sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
#if_file = glob.glob(filename)  

seed = 1
np.random.seed(seed)
test = gene_BHBH(h0=70)
#event_rate0, zs_detected0, masses0, rhos_detected0 = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'const')
#event_rate1, zs_detected1, masses1, rhos_detected1 = test.mc_year_rate(a=[2.35,0.1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode0')
event_rate2, zs_detected2, masses2, rhos_detected2 = test.mc_year_rate(a=[a0, a1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode1')

dl_zs = dl(zs_detected2)
zs_all, chirp_mass_all, m1_all, m2_all, lumi_dis_all = zs_detected2, masses2[:,0], masses2[:,1], masses2[:,2], dl_zs                                                      
#
##%% Test if the data slope follows the redshift trend.
#plt.hist(m1_all[zs_all==2.141], m1_all[zs_all==5.531],normed=True, log=True)
#plt.show()

#%%
from scipy import interpolate
from cal_likelihood import cov_twofun
from scipy.optimize import fmin
from cal_likelihood import random_Theta  #The prior given dl and chirpmass, no errorbar
   
step = 0
 
def posterior(para, m1_obs,m_noise_level,sf_factor,z):
    a0, a1, mbh_max,mbh_min  = para
#    a_z = a0 + a1 * z
#    if 1.1 < a0 < 3 and -0.5 < a1 < 0.5  and 50 < mbh_max < 100 and 2 < mbh_min < 8:
    a_z = a0 + a1 * z/(1+z)
    if 1.1 < a0 < 3 and -0.5 < a1 < 1.5 and 50 < mbh_max < 100 and 2 < mbh_min < 8:
        x = np.logspace(np.log10(m1_obs.min()/1.1), np.log10(m1_obs.max()*1.1),50)  #!!! 50 could be big enough?
        post, az_lib, y_lib = [], [], []
        for i in range(len(z)): 
            if a_z[i] not in az_lib:
                y = cov_twofun(x, a=a_z[i], mbh_max=mbh_max,
                               mbh_min=mbh_min,sigma=m_noise_level) #!!!
            else:
                j = np.where(np.asarray(az_lib) == a_z[i])[0][0]
                y = y_lib[j]
            az_lib.append(a_z[i])
            y_lib.append(y)
            f = interpolate.interp1d(x, y)
            poss_m1_i = f(m1_obs[i])
            post.append(poss_m1_i)
            if i/300 > (i-1)/300:
                print i
        post = np.array(post)
        chisq = -0.5*np.sum(np.log(post)*sf_factor)
        t2 = time.time()
        t_cost = (t2-t1)
        t_per_cost = t_cost/step
        t_remain = t_per_cost * (step_total-step)
        step=np.sum(sampler.flatchain!=0)/(ndim*nwalkers)
        if step != 0:
            print "step:",step,"percent",round(step/(step_total/100),2),"%","\r","para:",para, "time remain:", round(t_remain/60), 'mins'
        return chisq
    else:
        return np.inf       
m_noise_level = 0.20
thetas = random_Theta()

para_ini = [a0, a1, mbh_max, mbh_min]
#filename = 'test3_mode1_siglogdiv3_{0}.txt'.format(int(m_noise_level*100)) 
#if_file = glob.glob(filename)
#if if_file == []:
#    para_result =  open(filename,'w') 
#else:
#    para_result =  open(filename,'r+') 
t1 = time.time()

index = np.arange(len(m1_all))
from cal_likelihood import fac_s_eff_v

np.random.seed(seed)
idx = np.random.choice(index, 1000)
m1 = m1_all[idx]
dl = lumi_dis_all[idx]
dl_noised =  np.random.lognormal(np.log(dl), 0.35, size=dl.shape)
mass_Chirp = chirp_mass_all[idx]
mass_Chirp_noised = np.random.lognormal(np.log(mass_Chirp), 0.17, size=mass_Chirp.shape)
m1_mu = np.log(m1)   # 0_med np.log(mu_star) = mu 
m1_sigstar= np.exp(m_noise_level)  #Not useful in the generate generation, but useful in understand the upper lower level.
m1_obs = np.random.lognormal(m1_mu, m_noise_level, size=m1.shape)  #Generating for the mu as med, 
m1_sig_fake = m1_obs * m_noise_level      #The fake "sigma", (m1_obs * m_noise_level) and (m1 * m_noise_level)
#    prior = select_effect(m1_obs)
#    prior_true = fac_s_eff_v(dl=dl, mass_Chirp=mass_Chirp, thetas=thetas)
#    prior = prior_true
prior = fac_s_eff_v(dl=dl_noised, mass_Chirp=mass_Chirp_noised, thetas=thetas)
prior[prior==0] = 0.001
sf_factor = 1/prior
z_inf = solve_z(np.array(dl_noised))
sf_sigma = np.log(sf_factor)/3
sf_factor = sf_factor/np.exp(sf_sigma**2/2)
print "m1_obs.min(), m1_obs.max():",m1_obs.min(), m1_obs.max()
#mini=fmin(posterior,para_ini,maxiter=1000, args=(m1_obs, m_noise_level, sf_factor,z_inf))

ndim, nwalkers = 4, 70
pos = [[a0, a1, mbh_max, mbh_min] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

step_total = 500.
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior, args=(m1_obs, m_noise_level, sf_factor,z_inf))  ###the input args() should change
sampler.run_mcmc(pos, step_total)

###############save the data###############
samples = sampler.chain.reshape((-1, ndim))   #sampler.chain[:, 50:, :] to cut the data
burn=sampler.chain[:,:,:].T
import pickle
pickle.dump([sampler.chain,burn], open('test3_mode1_seed1.txt','wb'))
