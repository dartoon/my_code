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

#%%Load HSC sample:
line_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
         'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
         'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
         'sed_chisq', 'logMBH', 'logMBH_err']
infers  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,20].astype(np.float)
HSC_MBHs_ = infers[:,22].astype(np.float)
HSC_MBHs_err_ = infers[:,23].astype(np.float)

HSC_RA_ = infers[:,2].astype(np.float)
HSC_DEC_ = infers[:,3].astype(np.float)

flags_  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1_catalog_flag.txt', dtype=str)
flags = flags_[:,0]

IDs, HSC_z, HSC_Mstar, HSC_MBHs, HSC_MBHs_err =[], [], [], [], []
HSC_RA, HSC_DEC  = [], []
for i in range(len(IDs_)):
    idx = np.where(IDs_[i] == flags)[0][0]
    if flags_[idx][1] == 'y':
        IDs.append(IDs_[i])
        HSC_z.append(HSC_z_[i])
        HSC_Mstar.append(HSC_Mstar_[i])  #Uncertainty as 0.2 dex
        HSC_MBHs.append(HSC_MBHs_[i])
        HSC_MBHs_err.append(HSC_MBHs_err_[i])
        HSC_RA.append(HSC_RA_[i])
        HSC_DEC.append(HSC_DEC_[i])
HSC_Mstar = np.array(HSC_Mstar)
HSC_MBHs = np.array(HSC_MBHs)
HSC_z = np.array(HSC_z)
HSC_Mstar = HSC_Mstar[HSC_z>0.5]
HSC_MBHs = HSC_MBHs[HSC_z>0.5]

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
wave_lam_g = 4700 * (1+redshift)  #rest-frame in A
F_lam_g = F_lam_5100 / (5100*(1+redshift) )**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A  #!!! 5100 rest-frame?
F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
g_pointsource = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
g_totall = -2.5*np.log10(10**(-0.4*g_galaxy) + 10**(-0.4*g_pointsource))
    

wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
r_pointsource = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
r_totall = -2.5*np.log10(10**(-0.4*r_galaxy) + 10**(-0.4*r_pointsource))

#%% Before adding noise and selection effect:
plt.figure(figsize=(11.5,12))      
plt.scatter(logM_mass, logMBH,c='red',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.24, label='SAM sample')

# plt.scatter(HSC_Mstar,HSC_MBHs,c='blue',
#             s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.24, label='HSC sample')

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
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.close()
#%%
I_mag_break = 21.7
#Transfer the magnitude to absolute ones:
abs_Mags_break = I_mag_break - 5*(np.log10(dl * 10**6)-1)   #Compare to g when z >0.5 and r when z<0.5
dMBH, dmag, dMstar= 0.4, 0.3, 0.2  #dmag is for host magnitude. 

weight_s05 = stat_weight[(redshift<0.5)*(r_totall<abs_Mags_break)]
weight_l05 = stat_weight[(redshift>0.5)*(g_totall<abs_Mags_break)]

logMBH_s05 = logMBH[(redshift<0.5)*(r_totall<abs_Mags_break)]
logMBH_l05 = logMBH[(redshift>0.5)*(g_totall<abs_Mags_break)]
logM_mass_s05 = logM_mass[(redshift<0.5)*(r_totall<abs_Mags_break)]
logM_mass_l05 = logM_mass[(redshift>0.5)*(g_totall<abs_Mags_break)]

Stellar_Mass_s05_nois = logM_mass_s05 + np.random.normal(0, dMstar, size=logM_mass_s05.shape)
BH_Mass_s05_nois = logMBH_s05 + np.random.normal(0, dMBH, size=logMBH_s05.shape)
Stellar_Mass_l05_nois = logM_mass_l05 + np.random.normal(0, dMstar, size=logM_mass_l05.shape)
BH_Mass_l05_nois = logMBH_l05 + np.random.normal(0, dMBH, size=logMBH_l05.shape)

#%% Before adding noise and selection effect:
plt.figure(figsize=(11.5,12))      
import matplotlib
cmap_r = matplotlib.cm.get_cmap('RdBu_r')

plt.scatter(Stellar_Mass_l05_nois, BH_Mass_l05_nois, c= weight_l05/weight_l05.max(), # c='red',
            s=220, marker=".",zorder=1, edgecolors='k', alpha = 0.4, cmap=cmap_r, label='SAM sample z>0.5')

plt.scatter(Stellar_Mass_s05_nois, BH_Mass_s05_nois,c= weight_s05/weight_s05.max(), #c='yellow',
            s=220, marker="*",zorder=1, edgecolors='k', alpha = 0.4, cmap=cmap_r, label='SAM sample z<0.5')

plt.scatter(HSC_Mstar,HSC_MBHs,c='red',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.24, label='HSC sample')

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

#%%
TNG_scatter_0 = (BH_Mass_s05_nois - ( m_ml*Stellar_Mass_s05_nois+b_ml ) )
TNG_scatter_1 = (BH_Mass_s05_nois - ( m_ml*Stellar_Mass_s05_nois+b_ml ) )
TNG_scatter = np.concatenate((TNG_scatter_0, TNG_scatter_1))

HSC_scatter = (HSC_MBHs - ( m_ml*HSC_Mstar+b_ml ) )

#Plot the 1-D scatter for MM.
fig, ax = plt.subplots(figsize=(8,7))
plt.hist(HSC_scatter, histtype=u'step',density=True,
          label=('HSC sample scatter'), linewidth = 2, color='orange')
plt.hist(TNG_scatter,histtype=u'step',density=True,
          label=('SAM sample scatter'), linewidth = 2, color='green')
plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
plt.tick_params(labelsize=20)
plt.legend(prop={'size':20})
plt.yticks([])

# ax.xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(which='both', width=2, top=True,direction='in')
plt.tick_params(which='major', length=10)
plt.tick_params(which='minor', length=6)#, color='r’)

plt.xlabel('$\Delta$log(M$_{*}$/M$_{\odot}$)',fontsize=30)
#plt.savefig('comp_scatter_MM_MBIIonly.pdf')
plt.show()

from scipy import stats
sim_scatter_std = np.std(TNG_scatter)
obs_scatter_std = np.std(HSC_scatter)
print("obs scatter:", obs_scatter_std)
print("sim scatter:", sim_scatter_std)
print("KS scatter:", stats.ks_2samp(TNG_scatter, HSC_scatter).pvalue)


#%%
s = stat_weight/stat_weight.max()
plt.scatter(logM_mass[s>0.2], g_galaxy[s>0.2], c= s[s>0.2])
plt.xlabel('M*')
plt.ylabel('host galaxy g magnitude')
plt.title("SAM simulation")
plt.show()