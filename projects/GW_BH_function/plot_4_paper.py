#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 17:38:50 2019

@author: Dartoon

Plot the figures for the paper
"""
import numpy as np
import matplotlib.pyplot as plt
#%%
from BH_mass_function import gene_BHBH, dl, solve_z
import pickle
import glob
import corner

a, mbh_max, mbh_min = 2.35, 80., 5.
filename = 'sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
if_file = glob.glob(filename)  
if if_file==[]:
    test = gene_BHBH(h0=70)
    event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min)
    dl_zs = dl(zs_detected)
    sim_data = [event_rate, zs_detected, masses, rhos_detected, dl_zs]
    pickle.dump(sim_data, open(filename, 'wb'))
else:
    event_rate, zs_detected, masses, rhos_detected, dl_zs=pickle.load(open(filename,'rb'))
 
#fig = corner.corner(np.asarray([masses[:,0], rhos_detected]).T, labels=["Chirp mass", "rho"],
#                    #truths=[2.35, 0.7, 80., 5.],
#                    quantiles=[0.16, 0.84],show_titles=True,
#                    title_kwargs={"fontsize": 12},#truths=[2.35,80,5],
##                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
#                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
#                    title_fmt='.2f')    
#%%
from scipy.integrate import quad
def poss_ln_gaussian(m1, mu, sigma):
    if m1 >= 0:
        poss =  1/(m1 * sigma * np.sqrt(2 * np.pi)) * np.exp(-(np.log(m1) - mu)**2/(2*sigma**2))
    else:
        poss = 0
    return poss


def mass_fun_i(m1, a=2.35, mbh_max=80, mbh_min=5):
    if m1>=mbh_min and m1<=mbh_max:
        N = (mbh_max)**(-a+1)/(1-a) - (mbh_min)**(-a+1)/(1-a)             # The intergral function of m**a for normalize the power law function
        return m1**(-a)/N
    else:
        return 0.
poss_mass_fun=np.vectorize(mass_fun_i)    

def joint_twofun(tau, t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    return mass_fun_i(tau, a=a, mbh_max=mbh_max, mbh_min=mbh_min) * poss_ln_gaussian(t, mu=np.log(tau), sigma = sigma)
def cov_twofun(t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    inter = quad(joint_twofun, 0, 100, args=(t, a, mbh_max, mbh_min,sigma))[0]
    return inter
cov_twofun=np.vectorize(cov_twofun)   

x = np.logspace(np.log10(0.2), np.log10(85),300)
y = poss_mass_fun(x)
plt.plot(x, y)
x = np.logspace(np.log10(0.2), np.log10(85),300)
y = cov_twofun(x, sigma=0.04)
plt.plot(x, y)
plt.xlim(3,20)
plt.show()


#%% The corner plot:
filename = 'test2_select-eff_correct_sigmalogdiv3_20.txt'
with open(filename) as f:
        content = f.readlines()
lines = [x.strip() for x in content] 
import re
lines = [re.findall("\d+\.\d+", lines[i]) for i in range(len(lines))]
lines = [lines[i] for i in range(len(lines)) if len(lines[i]) ==3]
numbers = np.asarray(lines).astype(float)
samples = numbers#[numbers[:,1]<3]
fig = corner.corner(samples, labels=["a", "$MBH_{max}$", "$MBH_{min}$"],
                    truths=[2.35, 80., 5.],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},#truths=[2.35,80,5],
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f')
fig.savefig("fig_results_3para.pdf")
#####
plt.show()  

filename = 'test3_result/test3_mode1_level20.txt'
numbers = np.loadtxt(filename)
samples = numbers#[numbers[:,1]<3]
fig = corner.corner(samples, labels=["a0", "a1", "mbh_max", "mbh_min"],
                    truths=[2.35, 0.7, 80., 5.],
                    quantiles=[0.16, 0.84],show_titles=True,
                    title_kwargs={"fontsize": 12},#truths=[2.35,80,5],
#                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2),
                    title_fmt='.2f')
#fig.savefig("fig_results_4para.pdf")
#####
plt.show()  

