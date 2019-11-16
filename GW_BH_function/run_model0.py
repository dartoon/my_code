#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Test if whether the likelihood would recover the para
if all the m1 are measureable.
With error bar

The script for run scenario0
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

#a, mbh_max, mbh_min = 2.35, 80., 5.
a, mbh_max, mbh_min = 0.8, 50., 5.
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
##averaged Chirp mass
#print 'averaged Chirp mass, m1, m2: ', np.average(masses[:,0]), np.average(masses[:,1]), np.average(masses[:,2])
#Generate the select biased m1.
zs_all, chirp_mass_all, m1_all, m2_all, lumi_dis_all = zs_detected, masses[:,0], masses[:,1], masses[:,2], dl_zs                                                      
#plt.hist(masses[:,:2])
#plt.hist(zs_detected)
#plt.hist(m1)
#plt.hist(rhos_detected)
#plt.plot(m1,m2,'.')
#plt.yscale('log', nonposy='clip')
#plt.show()
#import corner
#fig = corner.corner(np.column_stack((m1_all,m2_all)), labels=["$m1$", "$m2$"],title_kwargs={"fontsize": 12})
##                    quantiles=[0.16, 0.84],show_titles=True,
##                    title_kwargs={"fontsize": 12})
##                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
##                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2))
##                    range=[(0.15,0.45),(60,90)] )
#plt.show()
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
    
def posterior(para, m1_obs,m_noise_level,sf_factor):
    a,mbh_max,mbh_min  = para
    if 0.5 < a < 3 and 30 < mbh_max < 80 and 2 < mbh_min < 8:   #!!!
        x = np.logspace(np.log10(m1_obs.min()/1.1), np.log10(m1_obs.max()*1.1),300)
        y = cov_twofun(x, a=a, mbh_max=mbh_max, mbh_min=mbh_min,sigma=m_noise_level)
        f = interpolate.interp1d(x, y)
        poss_m1_i = f(m1_obs)
        post = poss_m1_i
#        if prior_true is not None:
#            prior[np.where(prior<prior_true.min())]  = prior_true.min()
#        sf_factor[np.where(sf_factor>75)]  = 75  # Do it in other place
        chisq = -0.5*np.sum(np.log(post)*sf_factor)
        return chisq
    else:
        return np.inf       
#    even

#def each_likeli(para, m1_obs,m1_sig,prior):
#    a,mbh_max,mbh_min  = para
#    if 1.1 < a < 3 and 30 < mbh_max < 80 and 2 < mbh_min < 8:
#        x = np.logspace(np.log10(m1_obs.min()/np.exp(m_noise_level)**5/1.1), np.log10(m1_obs.max()*np.exp(m_noise_level)**5*1.1),300)
#        y = cov_twofun(x, a=a, mbh_max=mbh_max, mbh_min=mbh_min,sigma=m_noise_level)
#        f = interpolate.interp1d(x, y)
##        poss_m1_i = likelihood(m1_obs,m1_sig, a=a, mbh_max=mbh_max, mbh_min=mbh_min, use_method='med', f=f)
#        poss_m1_i = f(m1_obs)
#        post = poss_m1_i #(1/ prior)
##        print a, mbh_max, mbh_min, chisq
#        return post
#    else:
#        return np.inf


m_noise_level = 0.20
thetas = random_Theta()

para_ini = [a, mbh_max, mbh_min]
#filename = 'test2_select-eff_correct_sigmalogdiv3_{0}.txt'.format(int(m_noise_level*100))   #First file
filename = '201911_newrun/model0_a-{1}_mbhmax-{2}_noizl-{0}.txt'.format(int(m_noise_level*100), a , mbh_max) 
if_file = glob.glob(filename)
if if_file == []:
    para_result =  open(filename,'w') 
else:
    para_result =  open(filename,'r+') 
t1 = time.time()
rounds = 500
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

for loop in range(rounds):
    idx = random.sample(index, 1000)
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
    prior_true = fac_s_eff_v(dl=dl, mass_Chirp=mass_Chirp, thetas=thetas)
    prior = fac_s_eff_v(dl=dl_noised, mass_Chirp=mass_Chirp_noised, thetas=thetas)
    prior[prior==0] = 0.001
    sf_factor = 1/prior
#    sf_factor[np.where(sf_factor>300)]  = 300
#    z_inf = solve_z(np.array(dl_noised))
#    x_value = np.log(np.array(dl_noised)*(1/np.array(mass_Chirp_noised)**(5/6.))/(1+z_inf)**(5/6.))
#    sf_sigma = fit1d(x_value)
#    sf_sigma[x_value<v_min] = fit1d(v_min)
#    sf_sigma[x_value>v_max] = fit1d(v_max)
#    sf_sigma = 0.434302399  # The mean value of the sigma
    sf_sigma = np.log(sf_factor)/3
    sf_factor = sf_factor/np.exp(sf_sigma**2/2)
    print "m1_obs.min(), m1_obs.max():",m1_obs.min(), m1_obs.max()
    mini=fmin(posterior,para_ini,maxiter=1000, args=(m1_obs, m_noise_level, sf_factor))
#    mini_true = fmin(posterior,para_ini,maxiter=1000, args=(m1_obs, m_noise_level, prior))
#    print "ini Chisq:", posterior(para_ini, m1_obs,m_noise_level,prior)
#    print "Minimazied Chisq:", posterior(mini, m1_obs,m_noise_level,prior)
    if loop > 0:
        para_result = open(filename,'r+')
        para_result.read()
    para_result.write(repr(mini)+"\n")
    para_result.close()
    t2 = time.time()
    time_sp = t2-t1
    time_ave = (t2-t1)/(loop+1)
    time_total = time_ave * rounds
    t_left = time_total - time_sp
    print mini
    print "loop:", loop, "m_noise_level", m_noise_level, "Test alpha:", a
    print "Finish percent:",round(time_sp/time_total*100,2),"%" ,"total time needed :", round(time_total/60,2), "mins", "time_left", round(t_left/60,2), 'mins'

