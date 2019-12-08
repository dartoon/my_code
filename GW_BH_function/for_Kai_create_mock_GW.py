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

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHBH', scenario=2)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHBH', scenario=1)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHNS', scenario=2)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHNS', scenario=1)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'NSNS', scenario=2)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate))  
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate

#%%
test = gene_BHBH(h0=70, rho0=0, m_type = 'NSNS', scenario=1)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=1)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
save_file =  open(save_filename,'w') 
save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
for i in range(len(zs)):
    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
save_file.close() 
print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate
