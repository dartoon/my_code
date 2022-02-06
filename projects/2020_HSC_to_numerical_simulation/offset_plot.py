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

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

zs = 0.7
files = glob.glob("offset_result/*{0}.txt".format(zs))
if zs == 0.5 or zs == 0.7:
    files = files+glob.glob("offset_result/*{0}.txt".format(0.6))
files.sort()
ct = 0
# fig, ax = plt.subplots(figsize=(8,7))
# colors = ['green','steelblue','c','deeppink','plum','m']
# colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
colors = ['orange','deepskyblue', 'steelblue', 'c', 'deeppink', 'deeppink', 'm']
fig, axs = plt.subplots(len(files)+1, figsize=(5,8), sharex = True)
for i in [3,2,1,4,5,0]:
    data = np.loadtxt(files[i]) #nosl, sl, obs, sl_M*, sl_BH, obs_M*, obs_BH
    if zs == 1.5:
        data[:,2] = -data[:,2]
        data[:,1] = -data[:,1]
    
    if ct == 0:
        his_xy_ = axs[0].hist(data[:,2][abs(data[:,2])!=99],histtype= u'barstacked', 
                 density=True, color = colors[0],
                 linewidth = 3, alpha = 0.8)
        obs_mean = np.mean(data[:,2][abs(data[:,2])!=99])
        obs_ks_sample = data[:,2][abs(data[:,2])!=99]
        axs[0].plot([obs_mean,obs_mean], [0, 2], linewidth = 3,color = colors[0])
    axs[0].set_yticks([])
    axs[0].set_ylim([0,np.max(his_xy_[0])*1.1])
    # axs[0].set_xlim([his_xy_[1][0]-0.2 ,his_xy_[1][-1]+0.2])
    axs[0].set_xlim([-1.5, 2])
    axs[0].tick_params(which='major', width=2, length=4, direction='in')
    if zs == 1.5:
        axs[0].set_xlim([-1., 2])
    label = files[i].split("/")[1].split('_')[0] 
    if label == 'Horizon':
        label = label + '\n-AGN'
    his_xy = axs[ct+1].hist(data[:,1][abs(data[:,1])!=99], histtype=u'barstacked',
             density=True,
              linewidth = 14,alpha = 0.8,
               color = colors[ct+1])
    axs[ct+1].set_yticks([])
    if zs == 0.3 or zs == 0.7:
        axs[ct+1].set_ylabel(label, fontsize = 20)#, rotation=0)
        axs[0].set_ylabel('Obs.', fontsize = 20)
    axs[ct+1].plot([obs_mean,obs_mean], [0, 2], linewidth = 3,color = colors[0])
    # axs[ct+1].plot([np.mean(data[:,1][abs(data[:,1])!=99])]*2, [0, 2], linewidth = 3, color = colors[ct+1])
    axs[ct+1].plot([np.mean(data[:,1][abs(data[:,1])!=99])]*2, [0, 2],'--' , linewidth = 2 , color = 'black', alpha= 0.8)
    axs[ct+1].set_ylim([0,np.max(his_xy[0])*1.2])
    # axs[ct+1].set_xlim([his_xy_[1][0]-0.2 ,his_xy_[1][-1]+0.2])
    axs[ct+1].tick_params(which='major', width=2, length=4, direction='in')
    ct += 1
    n = len(data)
    for j in range(1, len(data)):
        if data[j-1,1] == -99 and data[j,1] != -99:
                n = j-1
                break
    data_c = data[:n+1]        
    # pvalue = stats.ks_2samp(data[:,2][:n+1] , data[:,1][:n+1]).pvalue
    pvalue = stats.ks_2samp(obs_ks_sample, data_c[:,1][abs(data_c[:,1])!=99] ).pvalue
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