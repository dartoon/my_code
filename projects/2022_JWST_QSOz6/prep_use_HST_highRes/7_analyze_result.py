#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 10:42:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import pickle
import glob

target_info = pickle.load(open('target_info.pkl','rb'))
res_files = glob.glob('/Users/Dartoon/Downloads/sim_results/qsoID*filt_f35*pkl')
# res_files = glob.glob('sim_results/*seed2_small*_result.pkl')
res_files.sort()

res_scatter_diffPSF, res_scatter_samePSF = [], []
res_reff = []
y = 'mag (inf$-$true)'
# y = 'reff ratio (inf/true)'
ratio = []
for file in res_files:
    res = pickle.load(open(file,'rb'))
    # print(round(res['true_host_flux'],1), round(res['inferred_host_flux'],1))
    # print(res['PSF_id_true'], res['PSF_id_model'])
    # print(res['host_Reff_kpc'], res['host_flux_ratio'])
    res_reff.append(res['assumed_host_Re_kpc'])
    if y == 'mag (inf$-$true)':
        res_scatter_diffPSF.append( res['inferred_magnitude_diff_psf'] - res['true_host_mag'] )
        res_scatter_samePSF.append( res['inferred_magnitude_same_psf'] - res['true_host_mag'] )
        # fit_run = pickle.load(open(file[:-11] + 'same_psf_fit.pkl','rb'))
        # res_scatter_samePSF.append( -2.5*np.log10(fit_run.final_result_galaxy[0]['flux_sersic_model']) + 28.98244840385909- res['true_host_mag'] )
    elif y == 'reff ratio (inf/true)':
        res_scatter_diffPSF.append( res['inferred_R_sersic_diff_psf'] / res['galfit_Re'] )
        res_scatter_samePSF.append( res['inferred_R_sersic_same_psf'] / res['galfit_Re'] )
    ratio.append( res['true_host_flux_ratio'] )
#%%

plt.figure(figsize=(11,11))
plt.scatter(ratio, res_scatter_samePSF, label = 'use True PSF')
plt.scatter(ratio, res_scatter_diffPSF, label = 'use diff PSF')
if y == 'mag (inf$-$true)':
    plt.ylim([-2.5, 2.5])
else:
    plt.ylim([0, 3])
plt.tick_params(labelsize=20)
plt.xlabel('flux ratio',fontsize=27)
plt.ylabel(y, fontsize=27)
plt.legend(prop={'size':28})
plt.show()

#%%
plt.figure(figsize=(11,11))
plt.scatter(res_reff, res_scatter_samePSF, label = 'use True PSF')
plt.scatter(res_reff, res_scatter_diffPSF, label = 'use diff PSF')
if y == 'mag (inf$-$true)':
    plt.ylim([-2.5, 2.5])
else:
    plt.ylim([0, 3])
plt.tick_params(labelsize=20)
plt.xlabel('Host Reff (kpc)',fontsize=27)
plt.ylabel(y, fontsize=27)
plt.legend(prop={'size':28})
plt.show()

#%%
if y == 'mag (inf$-$true)':
    res_scatter_samePSF = np.array(res_scatter_samePSF)
    res_scatter_samePSF = res_scatter_samePSF[abs(res_scatter_samePSF)<3]
    res_scatter_diffPSF = np.array(res_scatter_diffPSF)
    res_scatter_diffPSF = res_scatter_diffPSF[abs(res_scatter_diffPSF)<3]
plt.figure(figsize=(8,6))
if y == 'mag (inf$-$true)':
    dis0, dis1 = res_scatter_samePSF, res_scatter_diffPSF
else:
    y = 'log( '+y +')'
    dis0, dis1 = np.log10(res_scatter_samePSF), np.log10(res_scatter_diffPSF)

high0, x0, _ = plt.hist(dis0, density=True, histtype=u'step',
         label=('use True PSF'), linewidth = 2)
high1, x1, _ = plt.hist(dis1, density=True, histtype=u'step',
          label=('use diff PSF'), linewidth = 2)
plt.xlabel(y,fontsize=27)
# plt.ylabel("Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])
plt.xlim([-2.5, 2.5])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()
print(round(np.mean(dis0),2), round(np.std(dis0),2))
print(round(np.mean(dis1),2), round(np.std(dis1),2))

#%%
# Check Sersic index
n_diff_sP, n_diff_dP = [], []
for file in res_files:
    res = pickle.load(open(file,'rb'))
    inf_n_samePSF = res['inferred_n_sersic_same_psf']
    inf_n_diffPSF = res['inferred_n_sersic_diff_psf']
    # inf_n_noaddPSF = res['inferred_n_sersic_diff_psf']
    source_id = res['source_id']
    hostfile = glob.glob('test_int/*ID{0}*seed1*pkl'.format(source_id))[0]
    host_inf  = pickle.load(open(hostfile,'rb'))
    inf_n_noaddPSF = host_inf['inferred_n_sersic_same_psf']
    n_diff_sP.append( inf_n_samePSF - inf_n_noaddPSF )
    n_diff_dP.append( inf_n_diffPSF - inf_n_noaddPSF )
    
    
    

plt.figure(figsize=(11,11))
plt.scatter(ratio, n_diff_sP, label = 'use True PSF')
plt.scatter(ratio, n_diff_dP, label = 'use diff PSF')
plt.ylim([-2.5, 2.5])
# else:
#     plt.ylim([0, 3])
plt.tick_params(labelsize=20)
plt.xlabel('flux ratio',fontsize=27)
plt.ylabel('Sersic n (inf - true)', fontsize=27)
plt.legend(prop={'size':28})
plt.show()


plt.figure(figsize=(8,6))
high0, x0, _ = plt.hist(n_diff_sP, density=True, histtype=u'step',
         label=('use True PSF'), linewidth = 2)
high1, x1, _ = plt.hist(n_diff_dP, density=True, histtype=u'step',
          label=('use diff PSF'), linewidth = 2)
plt.xlabel('Sersic n (inf - true)',fontsize=27)
# plt.ylabel("Density",fontsize=27)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])
# plt.xlim([-2.5, 2.5])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()
print(round(np.mean(dis0),2), round(np.std(dis0),2))
print(round(np.mean(dis1),2), round(np.std(dis1),2))