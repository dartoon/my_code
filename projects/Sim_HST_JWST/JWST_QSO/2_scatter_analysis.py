#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:58:39 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

import sys
from matplotlib.colors import LogNorm
import pickle

#fig, ax = plt.subplots(figsize=(11,8))
#for ID in range(id_range[0], id_range[1]):
#    key = folder_type + '{0}'.format(ID)
#    gamma_bias = result_dic[key][1]['kwargs_lens'][0]['gamma'] - result_dic[key][0]['kwargs_lens'][0]['gamma']
#    plt.scatter(ID, gamma_bias,
#                c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#    ax.set_xticks(range(id_range[0]-1, id_range[1]+1,3)) 
#    plt.plot(np.linspace(id_range[0]-1, id_range[1]+1), np.linspace(id_range[0]-1, id_range[1]+1)*0)
#    plt.xlabel("ID",fontsize=27)
#    plt.ylabel("$\gamma$ bias (inferred - truth)",fontsize=27)
#    plt.ylim(-0.23,0.23)
#    plt.tick_params(labelsize=20)
#plt.show()
#
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
color = {'F150W': 'cyan', 'F200W':'lime', 'F356W':'tomato', 'F444W':'firebrick'}
seed = 0
use_true_PSF = False
if use_true_PSF == True:
    save_name = 'fit_result_truePSF'
else:
    save_name = 'fit_result'
for ID in [1, 2, 3, 4, 5]:
    fig, ax = plt.subplots(figsize=(11,8))
    for filt_i in range(4):    
        filt  = ['F444W', 'F356W', 'F200W', 'F150W'][filt_i]
        zp = [28., 27.9, 26.7, 27.75][filt_i]
        folder = 'sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed)    
#        print("The truth:")
        f = open(folder+"/sim_info.txt","r")
        string = f.read()
        lines = string.split('\n')   # Split in to \n
#        true_host_flux_text = float([lines[i] for i in range(len(lines)) if 'host_flux' in lines[i]][0].split('\t')[1])
        true_host = pyfits.getdata(folder+'/Drz_HOSTclean_image.fits')
        framesize = [101, 101, 81, 81][filt_i]
        half_r = int(framesize/2)
        peak = np.where(true_host==true_host.max())
        peak = [peak[0][0], peak[1][0]]
        true_host = true_host[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]        
        true_host_flux = true_host.sum()
        true_host_mag = -2.5*np.log10(true_host_flux) + zp
        true_host_ratio = [lines[i] for i in range(len(lines)) if 'host_flux_ratio' in lines[i]][0].split('\t')[1]
        true_total_flux = true_host_flux/float(true_host_ratio[:-1])*100
        true_total_mag = -2.5*np.log10(true_total_flux) + zp
        result = pickle.load(open(folder+'/{0}.pkl'.format(save_name),'rb'))
        fit_host_mag = result['host_mag'] 
        plt.scatter(filt_i, fit_host_mag - true_host_mag,
                c=color[filt],s=580, marker=".",zorder=0, edgecolors='k',alpha=0.7)
        plt.text(filt_i-0.1, fit_host_mag - true_host_mag+ 0.05, '{0}'.format(round(true_total_mag,2)), color='g',fontsize=14)
        plt.text(filt_i-0.1, fit_host_mag - true_host_mag- 0.05, '{0}%'.format(round(float(true_host_ratio[:-1]),1)), color='crimson',fontsize=14)
        plt.plot(np.linspace(-1, 5), np.linspace(-1, 5)*0,'k')   
    plt.title('Inferred mag bias for ID '+repr(ID), fontsize=27)
    #legend:
    plt.scatter(2.2, -0.4, c='dimgrey',s=580, marker=".",zorder=0, edgecolors='k',alpha=0.7)
    plt.text(2.2-0.1, -0.4+ 0.05, 'value: True AGN total magnitude', color='g',fontsize=14)
    plt.text(2.2-0.1, -0.4- 0.05, 'value: True host flux ratio', color='crimson',fontsize=14)
    labels = ['F444W', 'F356W', 'F200W', 'F150W']
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels)
    plt.ylabel("$\Delta$mag (inferred - truth)",fontsize=27)
    plt.xlim(-0.2, 3.2)
    plt.ylim(-0.6, 0.6)    
    plt.tick_params(labelsize=20)
#    plt.savefig("host_mag_bias_ID{0}_seed{1}.pdf".format(ID, seed))
    plt.show()