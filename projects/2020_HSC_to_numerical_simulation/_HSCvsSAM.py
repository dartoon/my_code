#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:27:56 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob

#Load HSC sample:
line_means = ['id','z','mbh','mbh_err','stellar_mass','lbol','spectra','bit','ps_gmag','ps_rmag','ps_imag','ps_rmag','ps_zmag','ps_ymag','host_gmag','host_rmag','host_imag','host_zmag','host_ymag']
infers  = np.loadtxt('HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,4].astype(np.float)
HSC_MBHs_ = infers[:,2].astype(np.float)
HSC_MBHs_err_ = infers[:,3].astype(np.float)
HSC_label_ = infers[:,-1]
HSC_Lbol = infers[:,5].astype(np.float)

HSC_galaxy_i_mag = infers[:,16].astype(np.float)
HSC_agn_i_mag = infers[:,10].astype(np.float)

HSC_z = HSC_z_
HSC_Mstar = HSC_Mstar_
HSC_MBHs = HSC_MBHs_

# plt.figure(figsize=(12.5,12))      
# plt.scatter(HSC_Mstar,HSC_MBHs,c=HSC_z,
#             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')
# clb = plt.colorbar()
m_ml, b_ml = (0.981139684856507, -2.545890295477823)
# xl = np.linspace(5, 13, 100)
# plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
# plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
# plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
# plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
# plt.xlim(9,12.5)
# plt.ylim(6.0,10.3)
# plt.grid(linestyle='--')
# plt.tick_params(labelsize=25)
# plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
# plt.show()


#%% load_SAM 
filename = 'SAM/catalogue_1st.dat'
text  = np.loadtxt(filename)
text = text[text[:,5] != 0]
# print(text.shape)
redshift = text[:, 0]
logM_mass =  text[:, 1]
logMBH =  text[:, 2]
g_galaxy =  text[:, 3]
r_galaxy =  text[:, 4]
AGN_bol_10_45 =  text[:, 5]
stat_weight =  text[:, 6]
 
L_5100_10_45 = AGN_bol_10_45 / 9.26
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dl = cosmo.luminosity_distance(redshift).value  # Mpc
dis = dl * 3.085677581e+24 #from Mpc to cm
F_lam_5100 = L_5100_10_45*10**45 / (4*np.pi * dis**2) / 5100

wave_lam_g = 4700 * (1+redshift)  #rest-frame in A  #vanden berk 2001, UV-optical SED slope as -0.44
F_lam_g = F_lam_5100 / (5100*(1+redshift) )**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A  #!!! 5100 rest-frame?
F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
g_pointsource = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
# g_totall = -2.5*np.log10(10**(-0.4*g_galaxy) + 10**(-0.4*g_pointsource))
g_totall = g_pointsource
    
wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
r_pointsource = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
# r_totall = -2.5*np.log10(10**(-0.4*r_galaxy) + 10**(-0.4*r_pointsource))
r_totall = r_pointsource

#%%
zs = 0.7  # for z < 0.5  #!!!
HSC_Mstar = HSC_Mstar[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]
HSC_MBHs = HSC_MBHs[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]
HSC_galaxy_i_mag = HSC_galaxy_i_mag[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]
HSC_agn_i_mag = HSC_agn_i_mag[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]
HSC_Lbol = HSC_Lbol[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]
HSC_z = HSC_z[(HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))]

import copy
#Transfer the magnitude to absolute ones:
I_mag_break = np.ones_like(redshift) * 20.5
I_mag_break[redshift>0.5] = 21.5
agn_mag = copy.deepcopy(r_totall)
agn_mag[redshift>0.5] = copy.deepcopy(g_totall[redshift>0.5])  #agn mag is r when z<0.5 and g when z>0.5

abs_Mags_break = I_mag_break - 5*(np.log10(dl * 10**6)-1)   #Compare with and with r when z<0.5,  g when z >0.5
dMBH, dmag, dMstar= 0.4, 0.3, 0.2  #dmag is for host magnitude. 
# dMBH, dmag, dMstar= 0.01, 0.01, 0.01  #dmag is for host magnitude. 

redshift_bool = (redshift<(zs+0.1)) * (redshift>(zs-0.1))
weight = stat_weight[redshift_bool]
logMBH = logMBH[redshift_bool]
AGN_bol_10_45 = AGN_bol_10_45[redshift_bool]
logM_mass = logM_mass[redshift_bool]
g_galaxy = g_galaxy[redshift_bool]
r_galaxy = r_galaxy[redshift_bool]
abs_Mags_break = abs_Mags_break[redshift_bool]
agn_mag = agn_mag[redshift_bool]
redshift = redshift[redshift_bool]

c_weight = np.ones_like(weight)
for i in range(len(weight)):
    c_weight[i] = np.sum(weight[:i+1])
#%% Before adding noise and selection effect:
plt.figure(figsize=(11.5,12))      
import matplotlib
cmap_r = matplotlib.cm.get_cmap('RdBu_r')
z_range = np.arange(0.2, 1.0, 0.05)
mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
mstar_cut = mstar_cut_range[zs > z_range][-1]
    
Stellar_Mass_nois_nosl, BH_Mass_nois_nosl, agn_mag_reali, abs_Mags_break_reali, redshift_reali,AGN_bol_reali  = [], [], [], [], [],[]
for i in range(len(HSC_Mstar)):
    seed = np.random.uniform(0, c_weight[-1])
    j = np.sum(seed > c_weight)  #Random based on the weighting
    Stellar_Mass_nois_nosl.append(logM_mass[j] + np.random.normal(0, dMstar) )
    BH_Mass_nois_nosl.append(logMBH[j] + np.random.normal(0, dMBH))
    agn_mag_reali.append(agn_mag[j])
    abs_Mags_break_reali.append(abs_Mags_break[j])
    AGN_bol_reali.append(np.log10(AGN_bol_10_45[j]) + 45 ) 
    redshift_reali.append(redshift[j])
Stellar_Mass_nois_nosl = np.array(Stellar_Mass_nois_nosl)
BH_Mass_nois_nosl = np.array(BH_Mass_nois_nosl)
AGN_bol_reali = np.array(AGN_bol_reali)
agn_mag_reali = np.array(agn_mag_reali)
abs_Mags_break_reali = np.array(abs_Mags_break_reali)

select_bool = (Stellar_Mass_nois_nosl>mstar_cut) * (agn_mag_reali<abs_Mags_break_reali) 
Lbol_boll = (AGN_bol_reali>43.6) * (AGN_bol_reali<46.2)
# select_bool = select_bool*Lbol_boll
BH_Mass_nois = BH_Mass_nois_nosl[select_bool]
Stellar_Mass_nois = Stellar_Mass_nois_nosl[select_bool]

SAM_scatter = (BH_Mass_nois - ( m_ml*Stellar_Mass_nois+b_ml ) )
SAM_scatter_nosl = (BH_Mass_nois_nosl - ( m_ml*Stellar_Mass_nois_nosl+b_ml ) )

plt.scatter(Stellar_Mass_nois_nosl, BH_Mass_nois_nosl,c='gray',
            s=220, marker=".",zorder=1, edgecolors='k', alpha = 0.2, label='SAM sample no sl')

plt.scatter(Stellar_Mass_nois, BH_Mass_nois,c='green',
            s=220, marker=".",zorder=1, edgecolors='k', alpha = 0.7, cmap='autumn', label='SAM sample z={0}'.format(zs))

plt.scatter(HSC_Mstar,HSC_MBHs,c='orange',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.7, label='HSC sample')

m_ml, b_ml = (0.981139684856507, -2.545890295477823)
xl = np.linspace(5, 13, 100)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)
plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.xlim(9,12.5)
plt.ylim(6.0,10.3)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=1,handletextpad=0)
plt.show()

print(np.mean(SAM_scatter_nosl), np.mean(SAM_scatter))

#%%
HSC_scatter = (HSC_MBHs - ( m_ml*HSC_Mstar+b_ml ) )

#Plot the 1-D scatter for MM.
fig, ax = plt.subplots(figsize=(8,7))
plt.hist(HSC_scatter, histtype=u'step',density=True,
          label=('HSC sample scatter'), linewidth = 2, color='orange')
plt.hist(SAM_scatter,histtype=u'step',density=True,
          label=('SAM sample scatter'), linewidth = 2, color='green')
plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])

# ax.xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(which='both', width=2, top=True,direction='in')
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=6)#, color='r’)

plt.xlabel(r'$\Delta$log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM_MBIIonly.pdf')
plt.show()

from scipy import stats
sim_scatter_std = np.std(SAM_scatter)
obs_scatter_std = np.std(HSC_scatter)
print("obs scatter:", obs_scatter_std)
print("sim scatter:", sim_scatter_std)
print("KS p-value:", stats.ks_2samp(SAM_scatter, HSC_scatter).pvalue)


#%%
logM_mass_plot, g_galaxy_plot, BH_mass_plot, agn_mag_plot, AGN_bol_plot  = [], [], [], [], []
for i in range(5000):
    seed = np.random.uniform(0, c_weight[-1])
    j = np.sum(seed > c_weight)  #Random based on the weighting    
    # if np.log10(AGN_bol_10_45[j])+45>43.6 and np.log10(AGN_bol_10_45[j])+45<46.2:
    if 1>0:
        logM_mass_plot.append(logM_mass[j])# + np.random.normal(0, 0.01) )
        g_galaxy_plot.append(g_galaxy[j])# + np.random.normal(0, 0.01))
        BH_mass_plot.append(logMBH[j])# + np.random.normal(0, 0.01))
        agn_mag_plot.append(agn_mag[j])
        AGN_bol_plot.append(AGN_bol_10_45[j])
    
logM_mass_plot = np.array(logM_mass_plot)
g_galaxy_plot = np.array(g_galaxy_plot)
AGN_bol_plot = np.array(AGN_bol_plot)
AGN_bol_plot = np.log10(AGN_bol_plot)+45
# plt.scatter(logM_mass, g_galaxy) #Before consider weighting.

dl_HSC = cosmo.luminosity_distance(HSC_z).value  # Mpc
HSC_galaxy_i_abs_mag = HSC_galaxy_i_mag - 5*(np.log10(dl_HSC * 10**6)-1)   #Compare HSC i with and with r when z<0.5,  g when z >0.5
HSC_agn_i_mag = HSC_agn_i_mag - 5*(np.log10(dl_HSC * 10**6)-1)  
plt.scatter(HSC_Mstar, HSC_galaxy_i_abs_mag, c = 'orange', alpha=0.2)
plt.scatter(logM_mass_plot, g_galaxy_plot, c = 'green',alpha=0.1)
plt.xlim(9,11.8)
plt.ylim(-25, -17.5)
plt.xlabel('M*')
plt.ylabel('host galaxy g magnitude')
plt.title("Stellar mass VS galaxy magnitude")
plt.close()

plt.figure(figsize=(11.5,10))      
plt.scatter(BH_mass_plot, agn_mag_plot, c='green',alpha=0.1)
plt.scatter(HSC_MBHs, HSC_agn_i_mag, c='orange',alpha=0.1)
plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + np.mean(abs_Mags_break_reali) )
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
plt.title("SAM simulation",fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(5,10)
plt.ylim(-28, -14)
plt.close()

#%%

plt.figure(figsize=(11.5,10))      
plt.scatter(AGN_bol_plot, BH_mass_plot, c='green',alpha=0.2)
plt.scatter(HSC_Lbol, HSC_MBHs,c='orange',alpha=0.2)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.xlabel('Lbol', fontsize=25)
plt.title("SAM simulation", fontsize=25)
plt.tick_params(labelsize=25)
# plt.xlim(5,10)
# plt.ylim(-28, -14)
plt.show()

plt.figure(figsize=(11.5,10))      
plt.scatter(AGN_bol_plot, logM_mass_plot, c='green',alpha=0.2)
plt.scatter(HSC_Lbol, HSC_Mstar,c='orange',alpha=0.2)
plt.ylabel(r'log(M$_{*}$/M$_{\odot}$)',fontsize=30)
plt.xlabel('Lbol', fontsize=25)
plt.title("SAM simulation", fontsize=25)
plt.tick_params(labelsize=25)
# plt.xlim(5,10)
# plt.ylim(-28, -14)
plt.show()