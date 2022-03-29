#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 18:16:30 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mpl
from sum_prep import load_HSC_comp, load_HST_comp, imf_dict
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

import glob
import pickle
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns

m_ml, b_ml = (0.981139684856507, -2.545890295477823)
#%%
eps = 1.5
I_mag_break=21.5
fig = plt.figure(figsize=(16*eps, 20*eps))
gs = fig.add_gridspec(5, 4, hspace=0, wspace=0)
axs = gs.subplots()
i, j = 0, 0
obsname = ['HSC', 'HST']
colors_sim = ['deepskyblue', 'steelblue', 'c', 'deeppink', 'm']
for zs in [0.3, 0.5, 0.7, 1.5]:
    if zs == 0.3:
        I_mag_break = 19.5
    if zs == 0.5 :
        I_mag_break = 20.5
    if zs == 0.7 :
        I_mag_break = 21.5
    if glob.glob('comb_zs{0}_Imag_{1}.pkl'.format(zs, I_mag_break)) != []:
        obs_dict, sim_dict = pickle.load(open(glob.glob('comb_zs{0}_Imag_{1}.pkl'.format(zs,I_mag_break))[0],'rb'))
    # else:
    #     if zs<1:
    #         obs_dict, sim_dict = load_HSC_comp(zs =zs, I_mag_break=I_mag_break)
    #     elif zs>1: 
    #         obs_dict, sim_dict = load_HST_comp()
    # pickle.dump([obs_dict, sim_dict], open('comb_zs{0}_Imag_{1}.pkl'.format(zs,I_mag_break), 'wb'))    
        
    for sname in ['SAM', 'MBII', 'Illustris', 'TNG100', 'Horizon']:
        sim = sim_dict[sname]
        imf = imf_dict[sname]
        obs = obs_dict[imf]
        
        if sname == 'MBII':
            if zs==0.5: 
                obs = obs_dict['zs06']   
                zs_ = 0.6
            if zs == 0.7:
                axs[i,j].axes.yaxis.set_ticklabels([])
                axs[i,j].axes.xaxis.set_ticklabels([])
                axs[i,j].xaxis.set_visible(False) 
                axs[i,j].yaxis.set_visible(False) 
                i = i+1
                continue
        else:
            zs_ = zs
        
        panel=axs[i,j].hist2d(sim['Stellar_Mass'], sim['BH_Mass'],
                          norm=mpl.colors.LogNorm(), density = True, cmap='summer',bins=50,zorder=-1,
                          alpha=0.5, cmin = 0.001)
        if i == 0:
            sns.kdeplot(sim['Stellar_Mass'], sim['BH_Mass'], linewidths = 2*eps, color = 'seagreen', 
                        fill=True, alpha=0.5, zorder = -10, ax = axs[i,j])

        s, alpha = 150, 0.7
        if sname == 'Horizon':
            sname = 'Horizon-AGN'
        axs[i,j].scatter(sim['Stellar_Mass_nois_sl'][0], sim['BH_Mass_nois_sl'][0],c=colors_sim[i],
                    s=s*eps, marker=".",zorder=0, 
                    edgecolors='black', alpha = alpha, label='{1} z={0}'.format(zs_, sname))
        axs[i,j].scatter(obs['Mstar'][0], #[obs['HSC_ps_mag']<22],
                    obs['MBHs'][0], #[obs['HSC_ps_mag']<22],
                    c='orange',
                    s=s*eps, marker=".",zorder=0.5, edgecolors='black', 
                    alpha = alpha, label=obsname[zs>1])
        #Random place the zorder for zs 0.5 and zs 0.7
        scatters_x = np.concatenate((sim['Stellar_Mass_nois_sl'][1:500], obs['Mstar'][1:])) 
        scatters_y = np.concatenate((sim['BH_Mass_nois_sl'][1:500], obs['MBHs'][1:])) 
        colors = [colors_sim[i]]*len(sim['Stellar_Mass_nois_sl'][1:500]) + ['orange']*len(obs['MBHs'][1:])
        colors = np.asanyarray(colors)
        orders = np.arange(len(scatters_x))
        if zs ==0.5 or zs == 0.7:
            np.random.shuffle(orders)
        axs[i,j].scatter(scatters_x[orders], 
                    scatters_y[orders], 
                    c=colors[orders],
                    s=s*eps, marker=".",zorder=0.5, edgecolors='black', 
                    alpha = alpha)
        xl = np.linspace(5, 13, 100)
        detlaM  = 0
        if imf == 'Sal': 
            detlaM = 0.23
        axs[i,j].plot(xl+detlaM, m_ml*xl+b_ml, color="k", linewidth=1.5*eps,zorder=-0.5)
        # plt.tick_params(which='both', width=2, top=True, right=True,direction='in')
        axs[i,j].grid(linestyle='--')
        axs[i,j].tick_params(labelsize=15*eps)
        axs[i,j].tick_params(which='major', length=7*eps,direction='in')
        axs[i,j].tick_params(which='minor', length=5*eps,direction='in')
        axs[i,j].legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':17*eps},ncol=1, 
                        handletextpad=0.5, handlelength=0, frameon=False )
        if zs<1:
            axs[i,j].set_xlim(9,12.5)
            axs[i,j].set_ylim(6.0,10.3)
        elif zs>1:
            axs[i,j].set_xlim(9.7,11.9)  #
            axs[i,j].set_ylim(7.2, 9.4)  #
        axs[i,j].xaxis.set_minor_locator(AutoMinorLocator())
        axs[i,j].yaxis.set_minor_locator(AutoMinorLocator())
        if i == 4:
            axs[i,j].set_xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=18*eps)
        if j == 0:
            axs[i,j].set_ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=18*eps)
        if i != 4:
            axs[i,j].axes.xaxis.set_ticklabels([])
        if j ==1 or j==2:
            axs[i,j].axes.yaxis.set_ticklabels([])
        if j == 3:
            box = axs[i,3].get_position()
            box.x0 = box.x0 + 0.03
            box.x1 = box.x1 + 0.03
            axs[i,3].set_position(box)
        i = i+1
    j = j+1
    i = 0
plt.show()    


