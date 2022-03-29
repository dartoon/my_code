#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 18:30:42 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from matplotlib.ticker import AutoMinorLocator
from scipy import stats
import pickle
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

zs = 1.5
print_p = True
print_dis = False


consider_noise = True
n_flag = ''
if consider_noise == False:
    n_flag = '_nonoise'

Imag_list = {0.3: 19.5, 0.5: 20.5, 0.7: 21.5}

files = glob.glob("offset_result/*{0}.txt".format(zs))
if zs == 0.5 or zs == 0.7:
    files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
files.sort()

ct = 0
m_ml, b_ml = (0.981139684856507, -2.545890295477823)
# fig, ax = plt.subplots(figsize=(8,7))
# colors = ['green','steelblue','c','deeppink','plum','m']
# colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
colors = ['orange','deepskyblue', 'steelblue', 'c', 'deeppink', 'deeppink', 'm']
fig, axs = plt.subplots(len(files)+1, figsize=(5,8), sharex = True)

from sum_prep import imf_dict
for i in [3,2,1,4,5,0]:
    label = files[i].split("/")[1].split('_')[0] 
    imf = imf_dict[label]
    detlaM  = 0
    if imf == 'Sal': 
        detlaM = 0.23
    # data = np.loadtxt(files[i]) #nosl, sl, obs, sl_M*, sl_BH, obs_M*, obs_BH
    # if zs == 1.5:
    #     data[:,2] = -data[:,2]
    #     data[:,1] = -data[:,1]
    
    # if Imag == 22.0:
    #     obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}{1}.pkl'.format(zs,n_flag))[0],'rb'))
    # else:
    if zs<1:
        Imag= Imag_list[zs]
        obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}{1}_Imag_{2}.pkl'.format(zs,n_flag,Imag))[0],'rb'))
        _, sim_int_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}_nonoise.pkl'.format(zs,Imag))[0],'rb'))
    elif zs>1:
        obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}.pkl'.format(zs))[0],'rb'))
        _, sim_int_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_nonoise.pkl'.format(zs))[0],'rb'))
        
    if zs == 0.5:
        _, sim_dict_ = pickle.load(open(glob.glob('pickles/comb_zs{0}{1}_Imag_{2}.pkl'.format(zs,n_flag,21.0))[0],'rb'))
        sim_dict['MBII'] = sim_dict_['MBII']

    if zs >1:
        # dis1 = obs_dict[imf]['logLbol']
        # dis1 = obs_dict[imf]['MBHs']
        off_obs = obs_dict['Cha']['MBHs'] - (m_ml*obs_dict['Cha']['Mstar']+b_ml),
        # dis1 = obs_dict[imf]['Mstar']
    elif zs <=1:         
        # dis1 = obs_dict[imf]['HSC_Lbol']
        # dis1 = obs_dict[imf]['HSC_MBHs']
        off_obs = obs_dict['Cha']['HSC_MBHs'] - (m_ml*obs_dict['Cha']['HSC_Mstar']+b_ml),
    off_sim = sim_dict[label]['BH_Mass_nois_sl']- (m_ml*sim_dict[label]['Stellar_Mass_nois_sl']+b_ml)
    # dis1 = obs_dict[imf]['HSC_Mstar']
    
    if ct == 0:
        his_xy_ = axs[0].hist(off_obs,histtype= u'barstacked', 
                  density=True, color = colors[0],
                  linewidth = 3, alpha = 0.8)
        obs_mean = np.mean(off_obs)
        axs[0].plot([obs_mean,obs_mean], [0, 2], linewidth = 3,color = colors[0])
    axs[0].set_yticks([])
    axs[0].set_ylim([0,np.max(his_xy_[0])*1.1])
    # axs[0].set_xlim([his_xy_[1][0]-0.2 ,his_xy_[1][-1]+0.2])
    axs[0].set_xlim([-1.5, 2])
    axs[0].tick_params(which='major', width=2, length=4, direction='in')
    
    # obs_dict, sim_dict = pickle.load(open(glob.glob('pickles/comb_zs{0}_Imag_{1}.pkl'.format(zs,Imag))[0],'rb'))
    
    if zs == 1.5:
        axs[0].set_xlim([-1., 2])
    sm_sim_int, bh_sim_int = sim_int_dict[label]['Stellar_Mass_nois_sl'][:]- detlaM, sim_int_dict[label]['BH_Mass_nois_sl'][:]
    off_int_sim = bh_sim_int - (m_ml*sm_sim_int+b_ml),
    sm_sim, bh_sim = sim_dict[label]['Stellar_Mass_nois_sl'][:]- detlaM, sim_dict[label]['BH_Mass_nois_sl'][:]
    off_sim = bh_sim - (m_ml*sm_sim+b_ml),
    
    bins = None
    if zs == 1.5 and label =='TNG300' or zs == 1.5 and label =='SAM' :
        bins = 30
    his_xy = axs[ct+1].hist(off_int_sim, 
              density=True,histtype=u'step',
              linewidth = 2,alpha = 0.8,
                color = 'k', zorder = 100,bins = bins)
    #Print no-noise simulation result
    if print_dis == True:        
        print(label, len(sim_dict[label]['Stellar_Mass_nois_sl']), 
              '({0:.2f}, {1:.2f})'.format(np.mean(off_int_sim), np.std(off_int_sim)),
              '({0:.2f}, {1:.2f})'.format(np.mean(off_sim), np.std(off_sim))) #/(dis[label])**3)

    
    if label == 'Horizon':
        label = label + '\n-AGN'
    axs[ct+1].hist(off_sim, histtype=u'barstacked',
             density=True,
              linewidth = 14,alpha = 0.8,
               color = colors[ct+1])
    axs[ct+1].set_yticks([])
    if zs == 0.3: # or zs == 0.7:
        axs[ct+1].set_ylabel(label, fontsize = 20)#, rotation=0)
        axs[0].set_ylabel('Obs.', fontsize = 20)
    axs[ct+1].plot([obs_mean,obs_mean], [0, 10], linewidth = 3,color = colors[0])
    axs[ct+1].plot([np.mean(off_sim)]*2, [0, 10],'--' , linewidth = 2 , color = colors[ct+1], alpha= 0.8)
    
    axs[ct+1].set_ylim([0,np.max(his_xy[0])*1.1])
    # axs[ct+1].set_xlim([his_xy_[1][0]-0.2 ,his_xy_[1][-1]+0.2])
    axs[ct+1].tick_params(which='major', width=2, length=4, direction='in')
    ct += 1
    # n = len(data)
    # for j in range(1, len(data)):
    #     if data[j-1,1] == -99 and data[j,1] != -99:
    #             n = j-1
    #             break
    # data_c = data[:n+1]        
    # pvalue = stats.ks_2samp(off_obs[:n+1] , off_sim[:n+1]).pvalue
    
    if print_p == True:    
        pvalue = stats.ks_2samp(off_obs[0], off_sim[0]).pvalue
        if pvalue>1.e-10:
            print(zs, files[i].split("/")[1].split('_')[0], "%e"%pvalue )
        else:
            print(zs, files[i].split("/")[1].split('_')[0], "<1e-10")

        
# axs[0].set_title(r"Offset compared with the local, z={0}".format(zs), fontsize = 25)
# fig.suptitle(r"Offset compared with the local, z={0}".format(zs), fontsize = 20)
fig.suptitle(r"Offset compared with the local, z={0}".format(zs), fontsize = 20)
plt.tight_layout()
plt.tick_params(labelsize=20)
# plt.legend(prop={'size':20})
# if zs == 1.5:
    # plt.legend(scatterpoints=1,numpoints=1,loc=1,prop={'size':20},ncol=1,handletextpad=0)
    # plt.xlim([-1.,2.5])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# plt.tick_params(which='major', width=2, length=10, direction='in')
# plt.tick_params(which='major', length=10)
# plt.tick_params(which='minor', length=6)#, color='r’)
plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=20)
plt.savefig('offset_dis_z{0}.pdf'.format(zs), bbox_inches = "tight")
plt.show()

# #%%
# zs = 1.5
# files = glob.glob("offset_result/*{0}.txt".format(zs))
# if zs == 0.5 or zs == 0.7:
#     files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
# files.sort()
# ct = 0
# fig, ax = plt.subplots(figsize=(8,7))
# for i in [3,2,1,4,5,0]:
#     data = np.loadtxt(files[i]) #nosl, sl, obs, sl_M*, sl_BH, obs_M*, obs_BH
#     if zs == 1.5:
#         data[:,2] = -data[:,2]
#     if ct == 0:
#         plt.hist(data[:,6][abs(data[:,6])!=99],histtype=u'step',density=True,
#                   label=' Observation', linewidth = 4)
        
#     plt.hist(data[:,4][abs(data[:,4])!=99], histtype=u'step',density=True,
#               label=' '+files[i].split("/")[1].split('_')[0], linewidth = 2,alpha = 0.8)
#               # color = colors[i])
#     ct += 1
#     pvalue = stats.ks_2samp(data[:,2][abs(data[:,2])!=99], data[:,1][abs(data[:,1])!=99] ).pvalue
#     if pvalue>1.e-10:
#         print(zs, files[i].split("/")[1].split('_')[0], "%e"%pvalue )
#     else:
#         print(zs, files[i].split("/")[1].split('_')[0], "<1e-10")

# plt.title(r"Observed MBH function z={0}".format(zs), fontsize = 25)
# plt.tick_params(labelsize=20)
# # plt.legend(prop={'size':20})
# if zs == 1.5:
#     plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=1,handletextpad=0)
# plt.yticks([])
# # plt.xlim([-2.5,3.5])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# plt.tick_params(which='both', width=2, top=True,direction='in')
# plt.tick_params(which='major', length=10)
# plt.tick_params(which='minor', length=6)#, color='r’)
# plt.xlabel(r'M$_{\rm BH}$',fontsize=30)
# plt.savefig('MBH_dis_z{0}.pdf'.format(zs))
# plt.show()

# #%%
# zs = 1.5
# files = glob.glob("offset_result/*{0}.txt".format(zs))
# if zs == 0.5 or zs == 0.7:
#     files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
# files.sort()
# ct = 0
# sal_list = [0.2266, 0.2266, 0, 0, 0, 0.2266]
# fig, ax = plt.subplots(figsize=(8,7))
# for i in [3,2,1,4,5,0]:
#     data = np.loadtxt(files[i]) #nosl, sl, obs, sl_M*, sl_BH, obs_M*, obs_BH
#     if zs == 1.5:
#         data[:,2] = -data[:,2]
#     if ct == 0: #This is the SAM, sal
#         plt.hist(data[:,5][abs(data[:,5])!=99] - sal_list[ct],histtype=u'step',density=True,
#                   label=' Observation', linewidth = 4)
        
#     plt.hist(data[:,3][abs(data[:,3])!=99]- sal_list[ct], histtype=u'step',density=True,
#               label=' '+files[i].split("/")[1].split('_')[0], linewidth = 2,alpha = 0.8)
#               # color = colors[i])
#     ct += 1
#     pvalue = stats.ks_2samp(data[:,2][abs(data[:,2])!=99], data[:,1][abs(data[:,1])!=99] ).pvalue
#     if pvalue>1.e-10:
#         print(zs, files[i].split("/")[1].split('_')[0], "%e"%pvalue )
#     else:
#         print(zs, files[i].split("/")[1].split('_')[0], "<1e-10")

# plt.title(r"Observed stellar mass distribution z={0}".format(zs), fontsize = 25)
# plt.tick_params(labelsize=20)
# # plt.legend(prop={'size':20})
# if zs == 1.5:
#     plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=1,handletextpad=0)
# plt.yticks([])
# plt.xlim([9.0,12])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
# plt.tick_params(which='both', width=2, top=True,direction='in')
# plt.tick_params(which='major', length=10)
# plt.tick_params(which='minor', length=6)#, color='r’)
# plt.xlabel(r'M$_*$',fontsize=30)
# plt.savefig('Mstar_dis_z{0}.pdf'.format(zs))
# plt.show()