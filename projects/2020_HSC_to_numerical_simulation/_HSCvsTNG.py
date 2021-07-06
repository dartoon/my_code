#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 16:46:42 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import scipy.stats as st
#%%Load HSC sample:
# line_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
#           'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
#           'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
#           'sed_chisq', 'logMBH', 'logMBH_err']
# infers  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1.txt', dtype=str)
# IDs_ = infers[:, 0]
# HSC_z_ = infers[:,1].astype(np.float)
# HSC_Mstar_ = infers[:,20].astype(np.float)
# HSC_gmag_galaxy_ = infers[:,9].astype(np.float)
# HSC_MBHs_ = infers[:,22].astype(np.float)
# HSC_MBHs_err_ = infers[:,23].astype(np.float)

# HSC_RA_ = infers[:,2].astype(np.float)
# HSC_DEC_ = infers[:,3].astype(np.float)

# flags_  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1_catalog_flag.txt', dtype=str)
# flags = flags_[:,0]

# IDs, HSC_z, HSC_Mstar, HSC_MBHs, HSC_MBHs_err, HSC_g_mag_galaxy =[], [], [], [], [], []
# HSC_RA, HSC_DEC  = [], []
# for i in range(len(IDs_)):
#     idx = np.where(IDs_[i] == flags)[0][0]
#     if flags_[idx][1] == 'y':
#         IDs.append(IDs_[i])
#         HSC_z.append(HSC_z_[i])
#         HSC_Mstar.append(HSC_Mstar_[i])  #Uncertainty as 0.2 dex
#         HSC_g_mag_galaxy.append(HSC_gmag_galaxy_[i])
#         HSC_MBHs.append(HSC_MBHs_[i])
#         HSC_MBHs_err.append(HSC_MBHs_err_[i])
#         HSC_RA.append(HSC_RA_[i])
#         HSC_DEC.append(HSC_DEC_[i])
# plt.scatter(HSC_Mstar, HSC_g_mag_galaxy) #The correlation between M* and g_band_mag
# plt.xlabel('M*')
# plt.ylabel('host galaxy g magnitude')
# plt.title("TNG simulation")
# plt.xlim(9,11.8)
# plt.ylim(-23.5, -17.5)
# plt.show()

#%%

# HSC_Mstar = np.array(HSC_Mstar)
# HSC_MBHs = np.array(HSC_MBHs)
# HSC_z = np.array(HSC_z)


# line_means = ['id', 'z', 'mbh', 'mbh_err', 'mgal', 'ps_mag', 'spectra', 'bit']
line_means = ['id','z','mbh','mbh_err','stellar_mass','lbol','spectra','bit','ps_gmag','ps_rmag','ps_imag','ps_rmag','ps_zmag','ps_ymag','host_gmag','host_rmag','host_imag','host_zmag','host_ymag']
infers  = np.loadtxt('HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_Mstar_ = infers[:,4].astype(np.float)
HSC_MBHs_ = infers[:,2].astype(np.float)
HSC_ps_mag_overall = infers[:,10].astype(np.float)
HSC_MBHs_err_ = infers[:,3].astype(np.float)
HSC_label_ = infers[:,-1]
HSC_Lbol_overall = infers[:,5].astype(np.float)


HSC_z = HSC_z_
HSC_Mstar_overall = HSC_Mstar_
HSC_MBHs_overall = HSC_MBHs_
#If use the core:
# HSC_z, HSC_Mstar, HSC_MBHs = [], [], []
# for i in range(len(IDs_)):
#     if HSC_label_[i] in ['eboss_core', 'boss_core', 'ugri']:
#         HSC_z.append(HSC_z_[i])
#         HSC_Mstar_overall.append(HSC_Mstar_[i])
#         HSC_MBHs_overall.append(HSC_MBHs_[i])
# HSC_z = np.array(HSC_z)
# HSC_Mstar_overall = np.array(HSC_Mstar)
# HSC_MBHs_overall = np.array(HSC_MBHs)

#%%Load TNG sample:
filenames = glob.glob('TNG100/*') 
filenames.sort()
# for i in range(len(filenames)):
# for i in [1]:     #!!! 
for i in [3]:    
    text  = np.load(filenames[i])
    zs = float(filenames[i].split("_z")[1][:4])
    filename  = filenames[i]
    print(text.shape)
    # Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource = np.load('TNG100/data_z0.70.npy')
    Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource, Eddington_ratio = np.load(filename)

BH_Mass = np.log10(BH_Mass)
Stellar_Mass = np.log10(Stellar_Mass)

#%%
plt.figure(figsize=(11.5,10))      
plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(5,10)
plt.ylim(-26, -14)
plt.close()

plt.figure(figsize=(11.5,10))      
plt.scatter(BH_Mass, np.log10(Eddington_ratio))
# plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.ylabel('AGN Eddington ratio', fontsize=25)
plt.tick_params(labelsize=25)
plt.xlim(5,10)
plt.close()

logLedd_overall = 38. + np.log10(1.2) + BH_Mass
logLbol = logLedd_overall + np.log10(Eddington_ratio)

plt.figure(figsize=(11.5,10))      
plt.scatter(logLbol + np.random.normal(0, 0.2, size=logLbol.shape) , BH_Mass + np.random.normal(0, 0.4, size=BH_Mass.shape) ,alpha=0.2)
plt.scatter(HSC_Lbol_overall, HSC_MBHs_overall,c='orange',alpha=0.2,zorder = 1)
# plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
plt.xlabel('AGN logLbol', fontsize=25)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.close()

#%%
if zs <= 0.5:
    # HSC_Mstar = HSC_Mstar_overall[HSC_z<0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z<0.5]
    I_mag_break = 20.5  #z~0.3
if zs > 0.5:    
    # HSC_Mstar = HSC_Mstar_overall[HSC_z>0.5]
    # HSC_MBHs = HSC_MBHs_overall[HSC_z>0.5]
    I_mag_break = 22.0    #z~0.7
    
redshift_bool = (HSC_z>(zs-0.1))*(HSC_z<(zs+0.1))
HSC_Mstar = HSC_Mstar_overall[redshift_bool]
HSC_MBHs = HSC_MBHs_overall[redshift_bool]
HSC_ps_mag = HSC_ps_mag_overall[redshift_bool]

import numpy as np
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)               #Perform the function EE in array style
#Calculate the luminosity distance:
dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #zs is a list of redshift of the sources (input as 1 D numpy array)

#%%Transfer the TNG magnitude to absolute ones:
   
i_wavelen = 7496
g_wavelen = 4700
r_wavelen = 6177

obs_sdss_i_pointsource = sdss_i_pointsource + 5*(np.log10(dl)-1)
F_mu_i = 10**(0.4*(-48.60 - obs_sdss_i_pointsource))   
F_lam_i = F_mu_i /(i_wavelen*(1+zs)) **2 * (2.9979 * 10**18)

# F_lam_g = F_lam_i / (i_wavelen*(1+zs))**(-(-0.44+2)) * (g_wavelen*(1+zs))**(-(-0.44+2))
# F_mu_g = F_lam_g * (g_wavelen*(1+zs)) **2 / (2.9979 * 10**18)
# obs_sdss_g_pointsource = (-2.5 * np.log10(F_mu_g)) - 48.60
# abs_sdss_g_pointsource = obs_sdss_g_pointsource - 5*(np.log10(dl)-1)
F_lam_r = F_lam_i / (i_wavelen*(1+zs))**(-(-0.44+2)) * (r_wavelen*(1+zs))**(-(-0.44+2))
F_mu_r = F_lam_r * (r_wavelen*(1+zs)) **2 / (2.9979 * 10**18)
obs_sdss_r_pointsource = (-2.5 * np.log10(F_mu_r)) - 48.60
abs_sdss_r_pointsource = obs_sdss_r_pointsource - 5*(np.log10(dl)-1)
               
# =============================================================================
# F_lam_5100 = L_5100_10_45*10**45 / (4*np.pi * dis**2) / 5100
# wave_lam_g = 4700 * (1+redshift)  #rest-frame in A
# F_lam_g = F_lam_5100 / (5100*(1+redshift))**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A
# F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
# obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
# g_pointsource = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
# =============================================================================
abs_Mags = I_mag_break -5*(np.log10(dl)-1)   #dl is the luminosity distance which is a function of redshift:
# abs_Mags = 0    

#!!! Following changes, not total Mag, but AGN mag to select.    
if zs >= 0.5:    
    sdss_mag_totall = sdss_g_pointsource #-2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))
elif zs<0.5:
    sdss_mag_totall = abs_sdss_r_pointsource #-2.5*np.log10(10**(-0.4*sdss_r_galaxy) + 10**(-0.4*abs_sdss_r_pointsource))

add_noise_select = 1
if add_noise_select == 0:
    dMBH, dmag, dMstar= 0.0001, 0.0001, 0.0001  #dmag is for host magnitude. 
elif add_noise_select == 1:
    dMBH, dmag, dMstar= 0.4, 0.3, 0.2  #dmag is for host magnitude. 


z_range = np.arange(0.2, 1.0, 0.05)
mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
mstar_cut = mstar_cut_range[zs > z_range][-1]


Stellar_Mass_nois_nosl = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
BH_Mass_nois_nosl = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)
logLbol_nois_nosl = logLbol + np.random.normal(0, 0.1, size=logLbol.shape)
selet_bool = (sdss_mag_totall<abs_Mags) * (Stellar_Mass_nois_nosl>mstar_cut)

def quasar_filter(group_list):
    xmin, xmax = int(HSC_Lbol_overall.min()-2), int(HSC_Lbol_overall.max()+3)
    ymin, ymax = int(HSC_MBHs_overall.min()-2), int(HSC_MBHs_overall.max()+3)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    # positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([HSC_Lbol_overall, HSC_MBHs_overall])
    kernel = st.gaussian_kde(values)
    # f = np.reshape(kernel(positions).T, xx.shape)
    # from matplotlib import cm 
    # plt.contourf(xx, yy, f, cmap=cm.Blues, alpha=0.5)
    t = [kernel.pdf([HSC_Lbol_overall[i] , HSC_MBHs_overall[i]]) for i in range(len(HSC_Lbol_overall))]
    min_pdf = np.min(t)
    bools = [ (kernel.pdf([group_list[0][i] , group_list[1][i]])>min_pdf)[0] for i in range(len(group_list[0]))]  
    return np.array(bools)
type1_bools = quasar_filter([logLbol_nois_nosl, BH_Mass_nois_nosl])
plt.figure(figsize=(11.5,10))      
plt.scatter(logLbol_nois_nosl, BH_Mass_nois_nosl,alpha=0.2)
plt.scatter(logLbol_nois_nosl[type1_bools], BH_Mass_nois_nosl[type1_bools],alpha=0.2)
plt.scatter(HSC_Lbol_overall, HSC_MBHs_overall,c='orange',alpha=0.2,zorder = 1)
# plt.scatter(BH_Mass, sdss_g_pointsource)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
plt.xlabel('AGN logLbol', fontsize=25)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
plt.tick_params(labelsize=25)
plt.xlim(30, 46.5)
plt.ylim(5.8,10)
plt.show()


selet_bool = selet_bool * type1_bools

Stellar_Mass_nois = Stellar_Mass_nois_nosl[selet_bool]
BH_Mass_nois = BH_Mass_nois_nosl[selet_bool]
sdss_g_galaxy_ = sdss_g_galaxy[selet_bool]
sdss_mag_totall_ = sdss_mag_totall[selet_bool]
        
plt.figure(figsize=(11.5,12))      
plt.scatter(Stellar_Mass_nois_nosl, BH_Mass_nois_nosl,c='gray',
            s=220, marker=".",zorder=-10, edgecolors='k', alpha = 0.2, label='TNG sample overall z={0}'.format(zs))
plt.scatter(Stellar_Mass_nois, BH_Mass_nois,c='steelblue',
            s=220, marker=".",zorder=0, edgecolors='k', alpha = 0.7, label='TNG sample z={0}'.format(zs))
plt.scatter(HSC_Mstar[HSC_ps_mag<I_mag_break],HSC_MBHs[HSC_ps_mag<I_mag_break],c='orange',
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
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()


#%%

TNG_scatter = (BH_Mass_nois - ( m_ml*Stellar_Mass_nois+b_ml ) )
TNG_scatter_nosl = (BH_Mass_nois_nosl - ( m_ml*Stellar_Mass_nois_nosl+b_ml ) )
HSC_scatter = (HSC_MBHs - ( m_ml*HSC_Mstar+b_ml ) )

#Plot the 1-D scatter for MM.
fig, ax = plt.subplots(figsize=(8,7))
plt.hist(TNG_scatter_nosl, histtype=u'step',density=True,
          label=('HSC sample scatter nosl'), linewidth = 2, color='gray')
plt.hist(TNG_scatter,histtype=u'step',density=True,
          label=('TNG sample scatter'), linewidth = 2, color='steelblue')
plt.hist(HSC_scatter, histtype=u'step',density=True,
          label=('HSC sample scatter'), linewidth = 2, color='orange')
# plt.hist(TNG_scatter_noselect,histtype=u'step',density=True,
#           label=('TNG sample scatter, no selection'), linewidth = 2, color='gray')
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
sim_scatter_std = np.std(TNG_scatter)
obs_scatter_std = np.std(HSC_scatter)
print("obs scatter:", obs_scatter_std)
print("sim scatter:", sim_scatter_std)
print("KS p-value:", stats.ks_2samp(TNG_scatter, HSC_scatter).pvalue)

print(np.mean(TNG_scatter_nosl), np.mean(TNG_scatter) )

      
# """
# #%% simulation
# plt.scatter(Stellar_Mass_nois, sdss_g_galaxy_) #The correlation between M* and g_band_mag
# plt.xlabel('M*')
# plt.ylabel('host galaxy g magnitude')
# plt.title("TNG simulation")
# plt.xlim(9,11.8)
# plt.ylim(-23.5, -17.5)
# plt.show()

# #%%
# line_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
#           'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
#           'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
#           'sed_chisq', 'logMBH', 'logMBH_err']
# infers  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1.txt', dtype=str)
# IDs_ = infers[:, 0]
# HSC_z_ = infers[:,1].astype(np.float)
# HSC_Mstar_ = infers[:,20].astype(np.float)
# HSC_gmag_galaxy_ = infers[:,9].astype(np.float)
# HSC_MBHs_ = infers[:,22].astype(np.float)
# HSC_MBHs_err_ = infers[:,23].astype(np.float)

# HSC_RA_ = infers[:,2].astype(np.float)
# HSC_DEC_ = infers[:,3].astype(np.float)

# flags_  = np.loadtxt('HSC_fitting/sdss_quasar_decomposition_v1_catalog_flag.txt', dtype=str)
# flags = flags_[:,0]

# IDs, HSC_z, HSC_Mstar, HSC_MBHs, HSC_MBHs_err, HSC_g_mag_galaxy =[], [], [], [], [], []
# HSC_RA, HSC_DEC  = [], []
# for i in range(len(IDs_)):
#     idx = np.where(IDs_[i] == flags)[0][0]
#     if flags_[idx][1] == 'y':
#         IDs.append(IDs_[i])
#         HSC_z.append(HSC_z_[i])
#         HSC_Mstar.append(HSC_Mstar_[i])  #Uncertainty as 0.2 dex
#         HSC_g_mag_galaxy.append(HSC_gmag_galaxy_[i])
#         HSC_MBHs.append(HSC_MBHs_[i])
#         HSC_MBHs_err.append(HSC_MBHs_err_[i])
#         HSC_RA.append(HSC_RA_[i])
#         HSC_DEC.append(HSC_DEC_[i])
# HSC_MBHs = np.array(HSC_MBHs)
# HSC_z = np.array(HSC_z)
# HSC_Mstar = np.array(HSC_Mstar)
# HSC_g_mag_galaxy = np.array(HSC_g_mag_galaxy)

# dl_HSC=(1+HSC_z)*c*vec_EE(HSC_z)/h0 *10**6   #zs is a list of redshift of the sources (input as 1 D numpy array)
# HSC_g_abs_Mags = HSC_g_mag_galaxy -5*(np.log10(dl_HSC)-1)   #dl is the luminosity distance which is a function of redshift:
    
    
# plt.scatter(HSC_Mstar[([HSC_z>0.5][0]*[HSC_z<0.8][0])], HSC_g_abs_Mags[([HSC_z>0.5][0]*[HSC_z<0.8][0])]) #The correlation between M* and g_band_mag
# plt.xlabel('M*')
# plt.ylabel('host galaxy g magnitude')
# plt.xlim(9,11.8)
# plt.ylim(-23.5, -17.5)
# plt.show()

# # plt.scatter(HSC_Mstar[HSC_z>0.5], HSC_MBHs[HSC_z>0.5]) #The correlation between M* and g_band_mag
# # plt.show()
# """
# #%%
# plt.figure(figsize=(8,7))      
# plt.scatter(BH_Mass, sdss_g_pointsource, alpha = 0.1) #Overall sample
# plt.scatter(BH_Mass_nois, sdss_mag_totall_)
# plt.plot(np.linspace(5,10), np.linspace(5,10)*0 + abs_Mags)
# plt.xlabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=30)
# plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
# plt.tick_params(labelsize=25)
# plt.xlim(5,10)
# plt.ylim(-26, -14)
# plt.show()

# plt.figure(figsize=(8,7))      
# plt.scatter(Stellar_Mass, sdss_g_pointsource, alpha = 0.1) #Overall sample
# plt.scatter(Stellar_Mass_nois, sdss_mag_totall_)
# plt.plot(np.linspace(5,15), np.linspace(5,10)*0 + abs_Mags)
# plt.xlabel(r'log(M$_{*}$/M$_{\odot}$)',fontsize=30)
# plt.ylabel('AGN abs magnitude rest-frame g band', fontsize=25)
# plt.tick_params(labelsize=25)
# plt.xlim(8,12)
# plt.ylim(-26, -14)
# plt.show()

