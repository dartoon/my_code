#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:05:43 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import scipy.stats as st
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

#%%
from prep_comparison import HSC_set, comp_plot
# from prep_comparison import TNG_set as Illustris_set

filenames = glob.glob('EAGLE/*.txt') 
filenames.sort()

#Note that Spurious=1 should be excluded. The MassType_Star is the stellar mass within the whole 
#galaxy, in unites of Msun. You can use BlackHoleMassAccretionRate(Msun/yr) to calculate the L_bol.

#0GalaxyID,1ApertureSize,2Mass_Star,3Mass_BH:
texts0 = np.loadtxt(filenames[0], delimiter=',')
#0GalaxyID,1Redshift,2Snapnum,3Spurious,4BlackHoleMass,5BlackHoleMassAccretionRate,6MassType_Star:
texts1 = np.loadtxt(filenames[1], delimiter=',')
texts1 = texts1[texts1[:,3]!=1]
# zs = 0.7
# zs_range = np.array([zs-0.1, zs+0.1])
zs = 1.5
zs_range = np.array([zs-0.2, zs+0.2])
print(zs)
zs_bool = (texts1[:,1]>zs_range[0])*(texts1[:,1]<zs_range[1])
texts1 = texts1[zs_bool]
#%%
# IDs = texts1[:,0]
# idxs = []
# mass0,mass1,mass2,mass3 = [], [], [], []
# for i in range(len(IDs)):
#     try:
#         idx = np.where(IDs[i] == texts0[:,0])[0][0]
#         mass0.append([texts1[i, 6], texts1[i, 4]])  #All subhalo
#         mass1.append([texts1[i, 6], texts0[idx, 3]])  #
#         mass2.append([texts0[idx, 2], texts1[i, 4]])  #
#         mass3.append([texts0[idx, 2], texts0[idx, 3]])  #
#     except:
#         pass
# mass0 = np.array(mass0)
# mass1 = np.array(mass1)
# mass2 = np.array(mass2)
# mass3 = np.array(mass3)
# plt.scatter(np.log10(mass2[:,0]), np.log10(mass2[:,1]))
# plt.show()
#Mass 3 符合最好？
#%%

# 0.63 * 10**46 * BlackHoleMassAccretionRate = Lbol
IDs = texts1[:,0]
# write_file = open('EAGLE/simdata_{0}.txt'.format(zs),'w') 
# write_file.write("#ID, Stellar_Mass, BH_Mass, Eddington_ratio, logLbol \n")
Eddington_ratios, BH_Masses, logLbols = [], [], []
for i in range(len(IDs)):
    try:
        idx = np.where(IDs[i] == texts0[:,0])[0][0]
        smass = np.log10(texts0[idx, 2])
        smass_ = np.log10(texts1[i, 6])
        BH_Mass = np.log10(texts1[i, 4])
        logLbol = np.log10(0.63 * 10**46 * texts1[i, 5] * 0.9  )
        logLedd = 38. + np.log10(1.2) + BH_Mass
        Eddington_ratio = logLbol - logLedd
        # logLbol = logLedd + Eddington_ratio
        if texts1[i, 4] !=0:
            Eddington_ratios.append(Eddington_ratio)
            BH_Masses.append(BH_Mass)
            logLbols.append(logLbol)
            # write_file.write("{0} {1} {2} {3} {4} {5} {6}\n".format(int(IDs[i]), 
            #                                                     smass, BH_Mass,Eddington_ratio, logLbol, 
            #                                                     smass_, texts1[i, 1] ) )
    except:
        pass
# write_file.close()

#%%
# plt.scatter(logLbols, 38. + np.log10(1.2) +BH_Masses,c='orange',alpha=0.2,zorder = 0.5, label = 'EAGLE QSO distribution')
# plt.xlabel('AGN logLbol', fontsize=25)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.xlim(35, 60)
# plt.ylim(5.8,10)
plt.figure(figsize=(9,7))     
BH_Masses = np.array(BH_Masses)
Eddington_ratios = np.array(Eddington_ratios)
cal_M_range = np.concatenate([np.arange(5, 8, 0.3), np.arange(8.1,9.5,0.5)])
sim_stat = []
for i in range(len(cal_M_range)-1):
    s_bool = (BH_Masses>cal_M_range[i])*(BH_Masses<cal_M_range[i+1])
    Edd_s = 10**Eddington_ratios[s_bool]
    sim_stat.append( [np.mean(Edd_s), np.median(Edd_s)] )
sim_stat = np.array(sim_stat)
plt.plot(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2, np.log10(sim_stat[:,0]) , color = 'red', 
      zorder = 10, linewidth = 3, linestyle= '--')
plt.plot(cal_M_range[:-1]+ (cal_M_range[1]-cal_M_range[0])/2, np.log10(sim_stat[:,1]) , color = 'blue', 
      zorder = 10, linewidth = 3)

xspace = np.linspace(5,10)
plt.plot(xspace, 0*xspace,'k--',linewidth=1)
plt.plot(xspace, 0*xspace-1.5,'k--',linewidth=1)
y_line3 = -1.1*(xspace-7.5) -0.5
plt.plot(xspace, y_line3,'k--',linewidth=1)
yspace = np.linspace(-5,2)
plt.plot(yspace*0+7.7, yspace,'k--',linewidth=1)
plt.plot(xspace*0+8.6, yspace,'k--',linewidth=1)
xfill = np.linspace(7.7, 8.6)
yfill_sline = -1.1*(xfill-7.5) -0.5
y_sline1 = xfill*0
y_sline2 = xfill*0-1.5
y4 = np.maximum(yfill_sline, y_sline2)
plt.fill_between(xfill, y4, y2=0, color='steelblue', alpha=0.5, zorder=10)

plt.scatter(BH_Masses, Eddington_ratios,c='orange',alpha=0.2,zorder = 0.5, label = 'EAGLE QSO distribution')
plt.tick_params(labelsize=25)
plt.grid(linestyle='--')
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)', fontsize=25)
plt.ylabel(r"log(L$_{\rm bol}$/L$_{\rm Edd}$)", fontsize=25)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.xlim(6, 10)
plt.ylim(-5,0)
# plt.legend(scatterpoints=1,numpoints=1,prop={'size':20},ncol=2,handletextpad=0)
plt.show()

#%%

texts = np.loadtxt('Horizon/outt00266_halo_gal_centralBHs_2reff')
Stellar_Mass = np.log10(texts[:, 1])+ 11
BH_Mass = np.log10(texts[:, 2]) + 8
Eddington_ratio = np.log10(texts[:, 3])
# logLedd = 38. + np.log10(1.2) + BH_Mass
# logLbol = logLedd + Eddington_ratio
logLbol = texts[:, 4]
plt.scatter(BH_Mass, Eddington_ratio,c='blue',alpha=0.2,zorder = 0.5, label = 'Horizion QSO distribution')
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)', fontsize=25)
plt.ylabel(r"log(L$_{\rm bol}$/L$_{\rm Edd}$)", fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(6, 10)
plt.ylim(-5,0)
# plt.legend(scatterpoints=1,numpoints=1,prop={'size':20},ncol=2,handletextpad=0)
plt.show()
