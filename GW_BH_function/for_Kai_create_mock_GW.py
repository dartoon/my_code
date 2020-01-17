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

event_rate_list = []
zs_all_list = []
#%%
sce = 2
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHBH', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
sce = 1
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHBH', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_high_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
sce = 2
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHNS', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_low_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
sce = 1
test = gene_BHBH(h0=70, rho0=0, m_type = 'BHNS', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_high_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
sce = 2
test = gene_BHBH(h0=70, rho0=0, m_type = 'NSNS', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_low_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate))  
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_high_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
sce = 1
test = gene_BHBH(h0=70, rho0=0, m_type = 'NSNS', scenario=sce)
event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, itera=5)
zs_all, chirp_mass_all, m1_all, m2_all = zs_detected, masses[:,0], masses[:,1], masses[:,2]                                                      
index = np.arange(len(m1_all))
idx = random.sample(index, 5000)
m1 = m1_all[idx]
m2 = m2_all[idx]
rho = rhos_detected[idx]
zs = zs_all[idx]
#save_filename = 'simulated_{0}_high_end.txt'.format(test.m_type)
#save_file =  open(save_filename,'w') 
#save_file.write("#zs, m1, m2, SNR, total_rate = {0:.3f} \n".format(event_rate)) 
#for i in range(len(zs)):
#    save_file.write("{0:.3f} {1:.3f} {2:.3f} {3:.3f}\n".format(zs[i],m1[i],m2[i], rho[i])) 
#save_file.close() 
#print  'simulated_{0}_high_end.txt'.format(test.m_type), event_rate
zs_all_list.append(zs_all)
event_rate_list.append(event_rate)
#%%
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

event_rate_list = np.asarray(event_rate_list)
zs_all_list = np.asarray(zs_all_list)
plt.figure(figsize=(12,10))
classes = ['BH-BH', 'BH-NS', 'NS-NS']
color = ['orange', 'green', 'blue']
sce = ['low-end metallicity ' , 'high-end metallicity ']
line_type = ['-', '--']
zs_grid = test.z

for i in range(len(zs_all_list)):
    idx = random.sample(index, int(event_rate_list[i]/(event_rate_list.max())*40000) )
    zs_dis_g = zs_all_list[i][idx]
    zs_dis = []
    for j in zs_dis_g:
        zs_dis.append(np.random.uniform(zs_grid[zs_grid<j].max(),j))
    plt.hist(zs_dis, normed=False, histtype=u'step', bins = 30,
             label=(classes[i/2]+' system, '+sce[i%2]), linewidth = 2,linestyle=line_type[i%2], color=color[i/2])
plt.yticks([])
plt.xlabel("z$_s$",fontsize=27)
plt.ylabel("relative distribution",fontsize=27)
plt.legend(prop={'size':20})
plt.tick_params(labelsize=20)
plt.savefig('zs.pdf')
plt.show()
