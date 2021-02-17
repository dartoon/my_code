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

#%%Load TNG sample:
filenames = glob.glob('TNG100/*') 
filenames.sort()
# for i in range(len(filenames)):
for i in [1]:    
    text  = np.load(filenames[i])
    zs = float(filenames[i].split("_z")[1][:4])
    filename  = filenames[i]
    print(text.shape)
    # Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource = np.load('TNG100/data_z0.70.npy')
    Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource = np.load(filename)

I_mag_break = 21.7
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
#Transfer the magnitude to absolute ones:
   
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

if zs > 0.5:    
    sdss_mag_totall = -2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))
elif zs<0.5:
    sdss_mag_totall = -2.5*np.log10(10**(-0.4*sdss_r_galaxy) + 10**(-0.4*abs_sdss_r_pointsource))

add_noise_select = 1
if add_noise_select == 0:
    dMBH, dmag, dMstar= 0.0001, 0.0001, 0.0001  #dmag is for host magnitude. 
elif add_noise_select == 1:
    dMBH, dmag, dMstar= 0.4, 0.3, 0.2  #dmag is for host magnitude. 


BH_Mass = np.log10(BH_Mass)
Stellar_Mass = np.log10(Stellar_Mass)


Stellar_Mass_nois = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
BH_Mass_nois = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)
Stellar_Mass_nois = Stellar_Mass_nois[sdss_mag_totall<abs_Mags]
BH_Mass_nois = BH_Mass_nois[sdss_mag_totall<abs_Mags]

        
plt.figure(figsize=(11.5,12))      
plt.scatter(Stellar_Mass_nois, BH_Mass_nois,c='red',
            s=220, marker=".",zorder=-1, edgecolors='k', alpha = 0.24, label='TNG sample z={0}'.format(zs))

plt.scatter(HSC_Mstar,HSC_MBHs,c='blue',
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
plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()


# #%%

# TNG_scatter = (BH_Mass_nois - ( m_ml*Stellar_Mass_nois+b_ml ) )
# HSC_scatter = (HSC_MBHs - ( m_ml*HSC_Mstar+b_ml ) )

# #Plot the 1-D scatter for MM.
# fig, ax = plt.subplots(figsize=(8,7))
# plt.hist(HSC_scatter, histtype=u'step',density=True,
#          label=('HSC sample scatter'), linewidth = 2, color='orange')
# plt.hist(TNG_scatter,histtype=u'step',density=True,
#          label=('TNG sample scatter'), linewidth = 2, color='green')
# plt.title(r"The offset comparison for the M$_{\rm BH}$-M$_{*}$ relation", fontsize = 20)
# plt.tick_params(labelsize=20)
# plt.legend(prop={'size':20})
# plt.yticks([])

# # ax.xaxis.set_minor_locator(AutoMinorLocator())
# plt.tick_params(which='both', width=2, top=True,direction='in')
# plt.tick_params(which='major', length=10)
# plt.tick_params(which='minor', length=6)#, color='r’)

# plt.xlabel('$\Delta$log(M$_{*}$/M$_{\odot}$)',fontsize=30)
# #plt.savefig('comp_scatter_MM_MBIIonly.pdf')
# plt.show()

# from scipy import stats
# sim_scatter_std = np.std(TNG_scatter)
# obs_scatter_std = np.std(HSC_scatter)
# print("obs scatter:", obs_scatter_std)
# print("sim scatter:", sim_scatter_std)
# print("KS scatter:", stats.ks_2samp(TNG_scatter, HSC_scatter).pvalue)

# #%%
# plt.scatter(Stellar_Mass, sdss_g_galaxy) #The correlation between M* and g_band_mag
# plt.xlabel('M*')
# plt.ylabel('host galaxy g magnitude')
# plt.title("TNG simulation")
# plt.show()