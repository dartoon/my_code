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
#seed = 0
#use_true_PSF =False
use_true_PSF ='same_drizzle'
zp_dic = {'F444W':27.3012, 'F356W':27.1841, 'F200W':27.0383, 'F150W':26.8627} #Using mirage

#save_name = 'fit_result_samedrizzle_biggerframe'
#save_name = 'fit_result_samedrizzle'
#save_name = 'fit_result_diffPSF_biggerframe'
#save_name = 'fit_result_samedrizzle_biggerframecalrms'
#folder_suf = 'sim_seed201_5000s/'
# folder_suf = 'sim_seed301_2500s_reff12n2/'
# folder_suf = 'sim_seed301_1000s_reff13n14/'
# folder_suf = 'sim_seed301_5000s_reff13n14_ID246_noboost/'
folder_suf = 'sim_seed301_3000s_ID3/'

save_name = 'fit_result_diffPSF'
#if use_true_PSF == True:
#    save_name = 'fit_result_truePSF'
#elif use_true_PSF == False:
#    save_name = 'fit_result_diffPSF'
#elif use_true_PSF == 'same_drizzle':
#    save_name = 'fit_result_samedrizzle'

ID_range = [2,3]   
seed_range = [301,351]
#%%Mag bias

for ID in [3]:
    fig, ax = plt.subplots(figsize=(11,8))
    for filt_i in [1,2]:    
        filt  = ['F150W', 'F200W', 'F356W','F444W'][filt_i]
        zp = zp_dic[filt]
        bias_list = []
        for seed in range(seed_range[0],seed_range[1]):
            folder = folder_suf + 'sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed)    
            f = open(folder+"/sim_info.txt","r")
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            true_host = pyfits.getdata(folder+'/Drz_HOSTclean_image.fits')
            result, framesize = pickle.load(open(folder+'/{0}.pkl'.format(save_name),'rb'))
            half_r = int(framesize/2)
            peak = np.where(true_host==true_host.max())
            peak = [peak[0][0], peak[1][0]]
            true_host = true_host[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]        
            true_host_flux = true_host.sum()
            true_host_mag = -2.5*np.log10(true_host_flux) + zp
            # print(true_host_mag)
            true_host_ratio = [lines[i] for i in range(len(lines)) if 'host_flux_ratio' in lines[i]][0].split('\t')[1]
            true_total_flux = true_host_flux/float(true_host_ratio[:-1])*100  #!!!
            true_total_mag = -2.5*np.log10(true_total_flux) + zp
            fit_host_mag = result['host_mag'] 
            bias= fit_host_mag - true_host_mag
            bias_list.append(bias)
            plt.scatter(filt_i, bias,
                    c=color[filt],s=200, marker=".",zorder=0, edgecolors='w',alpha=0.5)
        if filt == 'F356W':
            seed = np.where(bias_list == np.min(bias_list))[0][0]
            print("seed:", np.where(bias_list == np.min(bias_list))[0][0])
            print("sim_ID{0}_{1}_seed{2}".format(ID, filt, seed))
#        plt.scatter(filt_i+0.1, np.mean(bias_list), c=color[filt], s=380, marker=".",zorder=0, edgecolors='black')
        if np.mean(bias_list)<1.0:
            y_loc = np.mean(bias_list)
        else:
            y_loc = 0.85
        plt.errorbar(filt_i-0.1, np.mean(bias_list), yerr=np.std(bias_list), color=color[filt],ecolor='black', fmt='o',zorder=-500,markersize=10)
        plt.text(filt_i-0.55, y_loc-0.05, '{0:.2f}\n$\pm${1:.2f}'.format(np.mean(bias_list),np.std(bias_list)), color='black',fontsize=24)
        plt.text(filt_i, -1.4+0.15, '{0}'.format(round(true_total_mag,2)), color='g',fontsize=24)
        plt.text(filt_i, -1.4-0.05, '{0}%'.format(round(float(true_host_ratio[:-1]),1)), color='crimson',fontsize=24)
    plt.plot(np.linspace(-1, 5), np.linspace(-1, 5)*0,'k')   
    plt.title('Tested mag bias for J085907.19+002255.9 ', fontsize=32)
    #legend:
    plt.text(-0.5, -1.4+0.15, 'AGN total magnitude', color='g',fontsize=24)
    plt.text(-0.5, -1.4-0.05, 'Host flux ratio', color='crimson',fontsize=24)
    labels = ['F200W', 'F356W']
    ax.set_xticks([1,2])
    plt.tick_params(which='both', width=3, length = 7)
    ax.set_xticklabels(labels)
    plt.ylabel("$\Delta$mag (inferred - truth)",fontsize=37)
    plt.xlim(-0.5, 3.2)
    plt.ylim(-1.5, 1.)    
    plt.tick_params(labelsize=40)
    plt.savefig("ID3_host_mag_bias_ID{0}_3000s.pdf".format(ID), bbox_inches = "tight")
    plt.show()
    
# #%%Color bias
# compare = [2, 0]    
# filt_list = ['F150W', 'F200W', 'F356W','F444W']
# for ID in range(ID_range[0],ID_range[1]):
#     fig, ax = plt.subplots(figsize=(8,8))
#     true_filt_mag_list = []
#     infe_filt_mag_list = []
#     for filt_i in range(4):    
#         filt  = filt_list[filt_i]
#         zp = zp_dic[filt]
#         true_mag_list = []
#         infe_mag_list = []
#         for seed in range(seed_range[0],seed_range[1]):
#             folder = folder = folder_suf + 'sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed)
#             f = open(folder+"/sim_info.txt","r")
#             string = f.read()
#             lines = string.split('\n')   # Split in to \n
#             true_host = pyfits.getdata(folder+'/Drz_HOSTclean_image.fits')
#             result, framesize = pickle.load(open(folder+'/{0}.pkl'.format(save_name),'rb'))
#             half_r = int(framesize/2)
#             peak = np.where(true_host==true_host.max())
#             peak = [peak[0][0], peak[1][0]]
#             true_host = true_host[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]        
#             true_host_flux = true_host.sum()
#             true_host_mag = -2.5*np.log10(true_host_flux) + zp
#             fit_host_mag = result['host_mag'] 
#             true_mag_list.append(true_host_mag)
#             infe_mag_list.append(fit_host_mag)
#         true_filt_mag_list.append(np.array(true_mag_list))
#         infe_filt_mag_list.append(np.array(infe_mag_list))
#     plt.scatter(true_filt_mag_list[compare[0]]-true_filt_mag_list[compare[1]], 
#                 infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]],
#                 s=200, marker=".",zorder=0, edgecolors='w',alpha=0.5)    
#     plt.errorbar(np.mean(true_filt_mag_list[compare[0]]-true_filt_mag_list[compare[1]])-0.05, 
#                   np.mean(infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]]), 
#                   yerr=np.std(infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]]), 
#                   color=color[filt],ecolor='black', fmt='o',zorder=-500,markersize=10)    
#     plt.text(np.mean(true_filt_mag_list[compare[0]]-true_filt_mag_list[compare[1]])-0.1, 
#               np.mean(infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]]), 
#               '{0}\n$\pm${1}'.format(round(np.mean(infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]]),2),round(np.std(infe_filt_mag_list[compare[0]]-infe_filt_mag_list[compare[1]]),2)), color='black',fontsize=14)
     
#     plt.plot(np.linspace(-2, 0), np.linspace(-2,0)*1,'k')   
#     plt.title('Inferred color bias for ID '+repr(ID), fontsize=27)
#     #legend:
#     plt.xlabel("True $\Delta$mag ({0} - {1})".format(filt_list[compare[0]], filt_list[compare[1]]),fontsize=20)
#     plt.ylabel("Inferred $\Delta$mag ({0} - {1})".format(filt_list[compare[0]], filt_list[compare[1]]),fontsize=20)
#     plt.tick_params(labelsize=20)
#     plt.xlim(-2.0, -0.5)
#     plt.ylim(-2.0, -0.5)
# #    plt.savefig("host_color_bias_ID{0}.pdf".format(ID))
#     plt.show()    
    
# #%%Other parameter bias
# #para = 'n'   
# para = 'Reff' 
# for ID in [3]:
#     fig, ax = plt.subplots(figsize=(11,8))
#     for filt_i in [1,2]:    
#         filt  = ['F150W', 'F200W', 'F356W','F444W'][filt_i]
#         bias_list = []
#         for seed in range(seed_range[0],seed_range[1]):
#             folder = folder_suf + 'sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed) 
#             f = open(folder+"/sim_info.txt","r")
#             string = f.read()
#             lines = string.split('\n')   # Split in to \n
#             true_para = [lines[i] for i in range(len(lines)) if 'host_'+para+':' in lines[i]][0].split('\t')[1]
#             true_q = [lines[i] for i in range(len(lines)) if '(phi, q):' in lines[i]][0].split('\t')[1]
#             true_q = float(true_q.split(', ')[1].split(')')[0])
#             true_para = float(true_para[:5]) / np.sqrt(true_q)
#             result, framesize = pickle.load(open(folder+'/{0}.pkl'.format(save_name),'rb'))
#             result_key_dic = {'Reff': 'R_sersic', 'n': 'n_sersic'}
#             fit_para = result[result_key_dic[para]] 
# #            print(fit_para)
#             bias= fit_para/true_para
#             bias_list.append(bias)
#             plt.scatter(filt_i, bias,
#                     c=color[filt],s=100, marker=".",zorder=0, edgecolors='w',alpha=0.5)
# #        if filt == 'F200W':
# #            seed = np.where(bias_list == np.min(bias_list))[0][0]
# #            print("seed:", np.where(bias_list == np.min(bias_list))[0][0])
# #            print("sim_ID{0}_{1}_seed{2}".format(ID, filt, seed))
# #        plt.scatter(filt_i+0.1, np.mean(bias_list), c=color[filt], s=380, marker=".",zorder=0, edgecolors='black')
#         plt.errorbar(filt_i-0.1, np.mean(bias_list), yerr=np.std(bias_list), color=color[filt],ecolor='black', fmt='o',zorder=-500,markersize=10)
#         plt.text(filt_i-0.35, np.mean(bias_list)-0.05, '{0}\n$\pm${1}'.format(round(np.mean(bias_list),2),round(np.std(bias_list),2)), color='black',fontsize=14)
#     if para == 'n':
#         plt.ylim(-2, 2)       
#     elif para == 'Reff':
#         plt.ylim(0.7, 1.3)    
# #    for filt_i in range(4):                   
# #        y_min, y_max=ax.get_ylim()
# #        plt.text(filt_i, y_min+(y_max-y_min)*0.1, '{0}'.format(round(true_total_mag,2)), color='g',fontsize=14)
# #        plt.text(filt_i, y_min+(y_max-y_min)*0.05, '{0}%'.format(round(float(true_host_ratio[:-1]),1)), color='crimson',fontsize=14)
#     plt.plot(np.linspace(-1, 5), np.linspace(-1, 5)*0+1,'k')   
#     plt.title('Inferred {0} bias for ID '.format(para)+repr(ID), fontsize=27)
#     #legend:
# #    plt.text(-0.5,y_min+(y_max-y_min)*0.25, 'Total magnitude', color='g',fontsize=14)
# #    plt.text(-0.5,y_min+(y_max-y_min)*0.2, 'Host flux ratio', color='crimson',fontsize=14)
#     labels = ['F200W', 'F356W']
#     ax.set_xticks([1,2])
#     ax.set_xticklabels(labels)
#     plt.ylabel("{0} ratio (inferred/truth)".format(para),fontsize=27)
#     plt.xlim(-0.5, 3.2)
#     plt.tick_params(labelsize=20)
#     # plt.savefig("host_Reff_bias_ID{0}_5000.pdf".format(ID))
#     plt.show()
    
