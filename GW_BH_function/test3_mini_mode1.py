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
import time
import pickle

solve_z = np.vectorize(solve_z)
#a0, a1, mbh_max, mbh_min = 2.35, 0.1, 80., 5.
a0, a1, mbh_max, mbh_min = 2.35, 0.7, 80., 5.   #mode1

seed = 2
filename = 'test3_seed{0}_data'.format(seed)
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

#%%
from scipy import interpolate
from cal_likelihood import cov_twofun
from scipy.optimize import fmin
from cal_likelihood import random_Theta  #The prior given dl and chirpmass, no errorbar
def posterior(para, m1_obs,m_noise_level,sf_factor,z):
    a0, a1, mbh_max,mbh_min  = para
#    a_z = a0 + a1 * z
#    if 1.1 < a0 < 3 and -0.5 < a1 < 0.5  and 50 < mbh_max < 100 and 2 < mbh_min < 8:
    a_z = a0 + a1 * z/(1+z)
    if 1.1 < a0 < 3.5 and -0.5 < a1 < 1.5 and 50 < mbh_max < 115 and 2 < mbh_min < 8:
        post = []
        for i in range(len(z)): 
            poss_m1_i = cov_twofun(m1_obs[i], a=a_z[i], mbh_max=mbh_max, mbh_min=mbh_min,sigma=m_noise_level)
            post.append(poss_m1_i)
        post = np.array(post)
        chisq = -0.5*np.sum(np.log(post)*sf_factor)
#        print para, chisq
    	global mini_count
        if mini_count/20 > (mini_count-1)/20:
            print "State of count {0}".format(mini_count), para, chisq
        mini_count = mini_count+1
        return chisq
    else:
        return np.inf       
m_noise_level = 0.20
thetas = random_Theta()

para_ini = [a0, a1, mbh_max, mbh_min]
t1 = time.time()

index = np.arange(len(m1_all))
from cal_likelihood import fac_s_eff_v
rounds = 250
count = 0

part = 0
for loop in range(part*rounds,(part+1)*rounds):
    mini_count = 0 
    print "Calculating loop:", loop
    seed_i = loop
    np.random.seed(seed_i)
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
    prior = fac_s_eff_v(dl=dl_noised, mass_Chirp=mass_Chirp_noised, thetas=thetas)
    prior[prior==0] = 0.001
    sf_factor = 1/prior
    z_inf = solve_z(np.array(dl_noised))
    sf_sigma = np.log(sf_factor)/3
    sf_factor = sf_factor/np.exp(sf_sigma**2/2)
    print "m1_obs.min(), m1_obs.max():",m1_obs.min(), m1_obs.max()
    mini=fmin(posterior,para_ini,maxiter=1000, args=(m1_obs, m_noise_level, sf_factor,z_inf))
    datafile = 'test3_mode1_take2_level{0}_p{1}.txt'.format(int(m_noise_level*100),part) 
    if count ==0:
        if_file = glob.glob(datafile)
        if if_file == []:
            para_result =  open(datafile,'w') 
        else:
#            print "HAHAHA"
            para_result =  open(datafile,'r+')
            para_result.read()
    if count > 0:
        para_result = open(datafile,'r+')
        para_result.read()
    para_result.write("seed = {0}, ".format(seed_i))    
    para_result.write(repr(mini)+"\n")
    para_result.close()    
    t2 = time.time()
    time_sp = t2-t1
    time_ave = (t2-t1)/(count+1)
    time_total = time_ave * rounds
    t_left = time_total - time_sp
    print "Finish percent:",round(time_sp/time_total*100,2),"%" ,"total time needed :", round(time_total/60,2), "mins", "time_left", round(t_left/60,2), 'mins'
    count = count+1
