#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:21:57 2021

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import scipy.stats as st

#%%
from prep_comparison import HSC_set, Horizon_set, comp_plot

filenames = glob.glob('Horizon/*') 
filenames.sort()
filename = filenames[2]
zs = 0.5   #Actually from z = 1.0

# filename = filenames[3]
# zs = 0.5

# filename = filenames[4]
# zs = 0.3    #Actually from z = 0, not work, too few sample


if zs <= 0.5:
    # HSC_Mstar = HSC_Mstar_overall[HSC_z<0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z<0.5]
    I_mag_break = 20.5  #z~0.3
if zs > 0.5:    
    # HSC_Mstar = HSC_Mstar_overall[HSC_z>0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z>0.5]
    I_mag_break = 22.0    #z~0.7

detlaM  = 0
imf = 'Sal'
if imf == 'Sal':
    detlaM = 0.23

HSC = HSC_set(zs, core = False,imf = imf)

Horizon = Horizon_set(filename, HSC_Lbol_overall=HSC['HSC_Lbol_overall'], HSC_MBHs_overall=HSC['HSC_MBHs_overall'],
              zs = zs, I_mag_break = I_mag_break, imf =  imf)


#%%
comp_plot(Horizon['BH_Mass'], Horizon['Eddington_ratio'], 'BH_Mass', 'Eddington_ratio')
# plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.ylabel('AGN Eddington ratio', fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(5,10)
plt.show()

# comp_plot(Horizon['logLbol'], Horizon['BH_Mass'], 'logLbol', 'BH_Mass')
comp_plot(Horizon['logLbol_nois'], Horizon['BH_Mass_nois'], 'logLbol', 'BH_Mass')
plt.scatter(HSC['HSC_Lbol_overall'], HSC['HSC_MBHs_overall'], alpha = 0.1, color = 'orange')
# plt.xlabel('AGN logLbol', fontsize=25)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.show()

#%% 
comp_plot(Horizon['logLbol_nois'], Horizon['BH_Mass_nois'], 'logLbol_nois', 'BH_Mass_nois', alpha = 0.2, label = 'Horizon')
plt.scatter(Horizon['logLbol_nois_sl'], Horizon['BH_Mass_nois_sl'], color = 'green',alpha=0.2, zorder = 1, label = 'selected Horizon')
plt.scatter(HSC['HSC_Lbol_overall'], HSC['HSC_MBHs_overall'],c='orange',alpha=0.2,zorder = 0.5, label = 'HSC QSO distribution')
plt.xlabel('AGN logLbol', fontsize=25)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=2,handletextpad=0)
plt.show()


#%%
plt.figure(figsize=(11.5,12))      
plt.scatter(Horizon['Stellar_Mass_nois'], Horizon['BH_Mass_nois'],c='gray',
            s=220, marker=".",zorder=-10, edgecolors='k', alpha = 0.2)
plt.scatter(Horizon['Stellar_Mass_nois_sl'], Horizon['BH_Mass_nois_sl'],c='steelblue',
            s=220, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='Horizon sample z={0}'.format(zs))
# plt.scatter(HSC['HSC_Mstar'],HSC['HSC_MBHs'],c='orange',
#             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
plt.scatter(HSC['HSC_Mstar'][HSC['HSC_ps_mag']<I_mag_break],HSC['HSC_MBHs'][HSC['HSC_ps_mag']<I_mag_break],c='orange',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')

m_ml, b_ml = (0.981139684856507, -2.545890295477823)
xl = np.linspace(5, 13, 100)
plt.plot(xl+ detlaM, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.xlim(9,12.5)
plt.ylim(6.0,10.3)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()
#%%
Horizon_scatter = (Horizon['BH_Mass_nois_sl'] - ( m_ml* (Horizon['Stellar_Mass_nois_sl'] - detlaM)+b_ml ) )
Horizon_scatter_nosl = (Horizon['BH_Mass_nois'] - ( m_ml* (Horizon['Stellar_Mass_nois'] -detlaM) +b_ml ) )
HSC_scatter = (HSC['HSC_MBHs'] - ( m_ml*(HSC['HSC_Mstar'] - detlaM)+b_ml ) )
#Plot the 1-D scatter for MM.
fig, ax = plt.subplots(figsize=(8,7))
plt.hist(Horizon_scatter_nosl, histtype=u'step',density=True,
          label=('Horizon sample scatter nosl'), linewidth = 2, color='gray')
plt.hist(Horizon_scatter,histtype=u'step',density=True,
          label=('Horizon sample scatter'), linewidth = 2, color='steelblue')
plt.hist(HSC_scatter, histtype=u'step',density=True,
          label=('HSC sample scatter'), linewidth = 2, color='orange')
# plt.hist(Horizon_scatter_noselect,histtype=u'step',density=True,
#           label=('Horizon sample scatter, no selection'), linewidth = 2, color='gray')
plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 25)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':10})
plt.yticks([])
# ax.xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(which='both', width=2, top=True,direction='in')
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=6)#, color='r’)
plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM_MBIIonly.pdf')
plt.show()
from scipy import stats
sim_scatter_std = np.std(Horizon_scatter)
obs_scatter_std = np.std(HSC_scatter)
print("obs scatter:", obs_scatter_std)
print("sim scatter:", sim_scatter_std)
print("KS p-value:", stats.ks_2samp(Horizon_scatter, HSC_scatter).pvalue)
print(np.mean(Horizon_scatter_nosl), np.mean(Horizon_scatter) )


print("for paper Observation", 'zs=', zs)
print('{0:.2f}, {1:.2f}'.format(np.mean(HSC_scatter), np.std(HSC_scatter)))

print("for paper Horizon", 'zs=', zs)
print('{0:.2f}, {1:.2f}'.format(np.mean(Horizon_scatter), np.std(Horizon_scatter)))
