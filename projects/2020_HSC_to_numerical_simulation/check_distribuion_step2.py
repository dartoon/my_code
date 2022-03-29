#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 10:18:00 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from sum_prep import imf_dict

zs = 0.3
consider_noise = False
# Imag = [19.5, 20.0, 20.5, 21.0, 21.5, 22.0][2]

# zs 0.3 Imag 19.5; zs 0.5 Imag 20.5; zs 0.7 Imag 21.5
Imag_list = {0.3: 19.5, 0.5: 20.5, 0.7: 21.5}
Imag= Imag_list[zs]

files = glob.glob("offset_result/*{0}.txt".format(zs))
if zs == 0.5 or zs == 0.7:
    files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
files.sort()
m_ml, b_ml = (0.981139684856507, -2.545890295477823)

print(zs, Imag)
n_flag = ''
if consider_noise == False:
    n_flag = '_nonoise'
dis = {'MBII': 142.7, 'Illustris': 106.5, 'TNG100': 111, 'TNG300': 302, 'Horizon':142, 'SAM':100}
hists = []
labels = []
for i in [3,2,1,4,5,0]:
    label = files[i].split("/")[1].split('_')[0] 
    
    imf = imf_dict[label]
    detlaM  = 0
    if imf == 'Sal': 
        detlaM = 0.23
        
    # if Imag == 22.0:
    #     obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}{1}.pkl'.format(zs,n_flag))[0],'rb'))
    # else:
    obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{2}{1}.pkl'.format(zs,n_flag,Imag))[0],'rb'))
        
    if zs == 0.5 and Imag!=22.0:
        _, sim_dict_ = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{2}{1}.pkl'.format(zs,n_flag,21.0))[0],'rb'))
        sim_dict['MBII'] = sim_dict_['MBII']
    
    sm_sim, bh_sim = sim_dict[label]['Stellar_Mass_nois_sl']- detlaM, sim_dict[label]['BH_Mass_nois_sl']
    # print(label)
    # sm_sim, bh_sim = sim_dict[label]['Stellar_Mass_nois_sl'][:500]- detlaM, sim_dict[label]['BH_Mass_nois_sl'][:500]
    off_sim = bh_sim - (m_ml*sm_sim+b_ml)
    print(label, len(sim_dict[label]['Stellar_Mass_nois_sl']), '({0:.2f}, {1:.2f})'.format(np.mean(off_sim), np.std(off_sim)) ), #/(dis[label])**3)
    # print('({0:.2f}, {1:.2f})'.format(np.mean(off_sim), np.std(off_sim)) )
    
    dis0 = sim_dict[label]['logLbol_nois_sl']
    dis0 = sim_dict[label]['BH_Mass_nois_sl']
    # dis0 = sim_dict[label]['Stellar_Mass_nois_sl']
    if zs > 1.5:
        dis1 = obs_dict[imf]['logLbol']
        dis1 = obs_dict[imf]['MBHs']
        # dis1 = obs_dict[imf]['Mstar']
    if zs <=1:         
        dis1 = obs_dict[imf]['HSC_Lbol']
        dis1 = obs_dict[imf]['HSC_MBHs']
    # dis1 = obs_dict[imf]['HSC_Mstar']
    dis0 = dis0[dis0>0]
    hists.append(dis0)
    labels.append(label)
    # plt.figure(figsize=(8,6))
    # high0, x0, _ = plt.hist(dis0,density=True, histtype=u'step',
    #           label=('Simulation sample'), linewidth = 2, color='orange')
    # high1, x1, _ = plt.hist(dis1,density=True, histtype=u'step',
    #           label=('HSC sample'), linewidth = 2, color='green')
    # plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':12},ncol=2,handletextpad=0)
    # plt.close()
# print('obs', len(obs_dict[imf]['HSC_MBHs']) )
plt.figure(figsize=(8,6))
high1, x1, _ = plt.hist(dis1,density=True, histtype=u'step',
          label=('Obs sample'), linewidth = 5)
for i in range(len(labels)):
    dis0 = hists[i]
    high0, x0, _ = plt.hist(dis0,density=True, histtype=u'step',
              label=(labels[i]), linewidth = 2)
plt.legend(scatterpoints=1,numpoints=1,prop={'size':12},ncol=2,handletextpad=0)
plt.show()
print('Imag break', Imag)