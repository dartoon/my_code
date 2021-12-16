#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:33:43 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import scipy.stats as st
import glob
import copy
from astropy.cosmology import FlatLambdaCDM

def HSC_set(zs = 0.5, core = False, imf = 'Cha'):
    HSC = {}
    line_means = ['id','z','mbh','mbh_err','stellar_mass','lbol','spectra','bit','ps_gmag','ps_rmag','ps_imag','ps_rmag','ps_zmag','ps_ymag','host_gmag','host_rmag','host_imag','host_zmag','host_ymag']
    infers  = np.loadtxt('HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
    IDs_ = infers[:, 0]
    HSC_z_overall = infers[:,1].astype(np.float)
    HSC_Mstar_overall = infers[:,4].astype(np.float)
    if imf == 'Sal':
        HSC_Mstar_overall = HSC_Mstar_overall+ 0.23
    HSC_MBHs_overall = infers[:,2].astype(np.float)
    HSC_ps_mag_overall = infers[:,10].astype(np.float)  #'ps_imag'
    HSC_MBHs_err_overall  = infers[:,3].astype(np.float)
    HSC_Lbol_overall = infers[:,5].astype(np.float)
    HSC_i_mag_galaxy_overall = infers[:,16].astype(np.float)
    if core == True:
        HSC_label_ = infers[:,7]
        HSC_z, HSC_Mstar, HSC_MBHs, HSC_ps_mag, HSC_MBHs_err, HSC_Lbol, HSC_i_mag_galaxy, HSC_label= [], [], [], [], [], [], [], []
        for i in range(len(IDs_)):
            if HSC_label_[i] in ['eboss_core', 'boss_core', 'ugri']:
                HSC_z.append(HSC_z_overall[i])
                HSC_Mstar.append(HSC_Mstar_overall[i])
                HSC_MBHs.append(HSC_MBHs_overall[i])
                HSC_ps_mag.append(HSC_ps_mag_overall[i])
                HSC_MBHs_err.append(HSC_MBHs_err_overall[i])
                HSC_Lbol.append(HSC_Lbol_overall[i])
                HSC_i_mag_galaxy.append(HSC_i_mag_galaxy_overall[i])
                HSC_label.append([HSC_label_[i]])
        HSC_z_overall = np.array(HSC_z)
        HSC_Mstar_overall = np.array(HSC_Mstar)
        HSC_MBHs_overall = np.array(HSC_MBHs)   
        HSC_ps_mag_overall = np.array(HSC_ps_mag)   
        HSC_MBHs_err_overall  = np.array(HSC_MBHs_err)   
        HSC_Lbol_overall = np.array(HSC_Lbol)   
        HSC_i_mag_galaxy_overall = HSC_i_mag_galaxy
        HSC['label'] = np.array(HSC_label)
    HSC['HSC_z_overall'] = HSC_z_overall
    HSC['HSC_Mstar_overall'] = HSC_Mstar_overall
    HSC['HSC_MBHs_overall'] = HSC_MBHs_overall
    HSC['HSC_ps_mag_overall'] = HSC_ps_mag_overall
    HSC['HSC_MBHs_err_overall']  = HSC_MBHs_err_overall
    HSC['HSC_Lbol_overall'] = HSC_Lbol_overall
    
    i_mag_cut = np.zeros_like(HSC['HSC_z_overall'])
    redshift_bool = (HSC_z_overall>(zs-0.1))*(HSC_z_overall<(zs+0.1))
    if core == True:
        for i in range(len(i_mag_cut)):
            if HSC['HSC_z_overall'][i]<0.5:
                i_mag_cut[i] = 20.5 
            else:
                i_mag_cut[i] = 22.0
            if HSC['label'][i] == 'ugri':
                i_mag_cut[i] = 19.1
            if HSC['label'][i] == 'eboss_core' or HSC['label'][i] == 'boss_core':
                i_mag_cut[i] = 22.0
        i_mag_cut_bool = (HSC['HSC_ps_mag_overall'] <i_mag_cut)
        select_bool = i_mag_cut_bool * redshift_bool
    else:
        select_bool = redshift_bool
    
    HSC['HSC_Mstar'] = HSC_Mstar_overall[select_bool ]
    HSC['HSC_MBHs'] = HSC_MBHs_overall[select_bool]
    HSC['HSC_ps_mag'] = HSC_ps_mag_overall[select_bool]  #'ps_imag'
    HSC['HSC_Lbol'] = HSC_Lbol_overall[select_bool]
    HSC['HSC_MBHs_err'] = HSC_MBHs_err_overall[select_bool]
    HSC['HSC_z'] = HSC_z_overall[select_bool]
    
    dl_HSC=(1+HSC_z_overall)*c*vec_EE(HSC_z_overall)/h0 *10**6   #zs is a list of redshift of the sources (input as 1 D numpy array)
    HSC_galaxy_abs_iMags_overall = HSC_i_mag_galaxy_overall -5*(np.log10(dl_HSC)-1)   #dl is the luminosity distance which is a function of redshift:
    HSC['HSC_galaxy_abs_iMags_overall'] = HSC_galaxy_abs_iMags_overall
    HSC['HSC_galaxy_abs_iMags'] = HSC_galaxy_abs_iMags_overall[select_bool]
    HSC['HSC_ps_abs_iMags_overall'] = HSC['HSC_ps_mag_overall'] - 5*(np.log10(dl_HSC)-1)
    HSC['HSC_ps_abs_iMags'] = HSC['HSC_ps_abs_iMags_overall'][select_bool]
    return HSC

h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
from scipy.integrate import quad
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)               #Perform the function EE in array style


def Horizon_set(filename, HSC_Lbol_overall, HSC_MBHs_overall, zs, I_mag_break = 20.5, consider_type1= True, imf ='Sal'):
    Horizon = {}
    # Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource, Eddington_ratio = np.load(filename)
    texts = np.loadtxt(filename)
    Stellar_Mass = np.log10(texts[:, 1])+ 11
    BH_Mass = np.log10(texts[:, 2]) + 8
    Eddington_ratio = np.log10(texts[:, 3])
    # logLedd = 38. + np.log10(1.2) + BH_Mass
    # logLbol = logLedd + Eddington_ratio
    logLbol = texts[:, 4]  #44

    Horizon['Stellar_Mass'] = Stellar_Mass
    Horizon['BH_Mass'] = BH_Mass
    # Horizon['sdss_g_pointsource'] = sdss_g_pointsource
    Horizon['Eddington_ratio'] = Eddington_ratio
    Horizon['logLbol'] = logLbol
    
    L_5100 = 10**logLbol / 9.26
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    redshift = zs  #!!!
    dl = cosmo.luminosity_distance(redshift).value  # Mpc
    dis = dl * 3.085677581e+24 #from Mpc to cm
    F_lam_5100 = L_5100 / (4*np.pi * dis**2) / 5100
    wave_lam_g = 4700 * (1+redshift)  #rest-frame in A  #vanden berk 2001, UV-optical SED slope as -0.44
    F_lam_g = F_lam_5100 / (5100*(1+redshift) )**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A  #!!! 5100 rest-frame?
    F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
    obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
    g_pointsource = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    # g_totall = -2.5*np.log10(10**(-0.4*g_galaxy) + 10**(-0.4*g_pointsource))
    # g_totall = g_pointsource
    wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
    F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
    F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
    obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
    r_pointsource = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    
    abs_Mags = I_mag_break -5*(np.log10(dl * 10**6 )-1)   #dl is the luminosity distance which is a function of redshift:
    if zs >= 0.5:    
        sdss_mag = g_pointsource #-2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))
    elif zs<0.5:
        sdss_mag= r_pointsource #-2.5*np.log10(10**(-0.4*sdss_r_galaxy) + 10**(-0.4*abs_sdss_r_pointsource))
    dMBH, dmag, dMstar, dLbol= 0.4, 0.3, 0.2, 0.1  #dmag is for host magnitude. 
    # dMBH, dmag, dMstar, dLbol= 0.001, 0.001, 0.001, 0.001  #dmag is for host magnitude. 
    
    z_range = np.arange(0.2, 1.0, 0.05)
    mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
    if imf =='Sal':
        mstar_cut_range = mstar_cut_range+0.23
    mstar_cut = mstar_cut_range[zs > z_range][-1]
    Stellar_Mass_nois = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
    BH_Mass_nois = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)
    logLbol_nois = logLbol + np.random.normal(0, dLbol, size=logLbol.shape)
    select_bool = (Stellar_Mass_nois >mstar_cut) * (sdss_mag <abs_Mags) #!!!Assuming sdss_mag_totall are accurate.
    type1_bools = quasar_filter([logLbol_nois, BH_Mass_nois], HSC_Lbol_overall, HSC_MBHs_overall)
    if consider_type1 == True:
        select_bool = select_bool * type1_bools
    Stellar_Mass_nois_sl = Stellar_Mass_nois[select_bool]
    BH_Mass_nois_sl = BH_Mass_nois[select_bool]
    logLbol_nois_sl = logLbol_nois[select_bool]
    # sdss_g_galaxy_sl = sdss_g_galaxy[select_bool]
    Horizon['Stellar_Mass_nois'] = Stellar_Mass_nois
    Horizon['Stellar_Mass_nois_sl'] = Stellar_Mass_nois_sl
    Horizon['BH_Mass_nois'] = BH_Mass_nois
    Horizon['BH_Mass_nois_sl'] = BH_Mass_nois_sl
    Horizon['logLbol_nois'] = logLbol_nois
    Horizon['logLbol_nois_sl'] = logLbol_nois_sl
    # Horizon['sdss_g_galaxy'] = sdss_g_galaxy  #!!! #TODO: Change to r for z<0.5
    # Horizon['sdss_g_galaxy_sl'] = sdss_g_galaxy_sl
    Horizon['sdss_g_pointsource'] = g_pointsource
    Horizon['sdss_g_pointsource_sl'] = g_pointsource[select_bool]
    # Horizon['select_abs_Mags'] = abs_Mags
    Horizon['select_bool'] = select_bool
    return Horizon

def EAGLE_set(filename, HSC_Lbol_overall, HSC_MBHs_overall, zs, I_mag_break = 20.5, consider_type1= True, imf ='Sal'):
    EAGLE = {}
    # Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource, Eddington_ratio = np.load(filename)
    texts = np.loadtxt(filename)
    Stellar_Mass = texts[:, 1]
    BH_Mass = texts[:, 2]
    Eddington_ratio = texts[:, 3]  #In log
    # logLedd = 38. + np.log10(1.2) + BH_Mass
    # logLbol = logLedd + Eddington_ratio
    logLbol = texts[:, 4]

    EAGLE['Stellar_Mass'] = Stellar_Mass
    EAGLE['BH_Mass'] = BH_Mass
    # EAGLE['sdss_g_pointsource'] = sdss_g_pointsource
    EAGLE['Eddington_ratio'] = Eddington_ratio
    EAGLE['logLbol'] = logLbol
    
    L_5100 = 10**logLbol / 9.26
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    redshift = texts[:, -1]  #!!!
    dl = cosmo.luminosity_distance(redshift).value  # Mpc
    dis = dl * 3.085677581e+24 #from Mpc to cm
    F_lam_5100 = L_5100 / (4*np.pi * dis**2) / 5100
    wave_lam_g = 4700 * (1+redshift)  #rest-frame in A  #vanden berk 2001, UV-optical SED slope as -0.44
    F_lam_g = F_lam_5100 / (5100*(1+redshift) )**(-(-0.44+2)) * wave_lam_g **(-(-0.44+2)) #erg/s/cm^2/A  #!!! 5100 rest-frame?
    F_mu_g = F_lam_g *  (wave_lam_g) **2 / (2.9979 * 10**18)
    obs_mag_g = [(-2.5 * np.log10(F_mu_g[i]) - 48.60) for i in range(len(F_mu_g))]
    g_pointsource = obs_mag_g - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    # g_totall = -2.5*np.log10(10**(-0.4*g_galaxy) + 10**(-0.4*g_pointsource))
    # g_totall = g_pointsource
    wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
    F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
    F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
    obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
    r_pointsource = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    
    abs_Mags = I_mag_break -5*(np.log10(dl * 10**6 )-1)   #dl is the luminosity distance which is a function of redshift:
    if zs >= 0.5:    
        sdss_mag = g_pointsource #-2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))
    elif zs<0.5:
        sdss_mag= r_pointsource #-2.5*np.log10(10**(-0.4*sdss_r_galaxy) + 10**(-0.4*abs_sdss_r_pointsource))
    dMBH, dmag, dMstar, dLbol= 0.4, 0.3, 0.2, 0.1  #dmag is for host magnitude. 
    # dMBH, dmag, dMstar, dLbol= 0.001, 0.001, 0.001, 0.001  #dmag is for host magnitude. 
    
    z_range = np.arange(0.2, 1.0, 0.05)
    mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
    if imf =='Sal':
        mstar_cut_range = mstar_cut_range+0.23
    mstar_cut = mstar_cut_range[zs > z_range][-1]
    Stellar_Mass_nois = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
    BH_Mass_nois = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)
    logLbol_nois = logLbol + np.random.normal(0, dLbol, size=logLbol.shape)
    select_bool = (Stellar_Mass_nois >mstar_cut) * (sdss_mag <abs_Mags) #!!!Assuming sdss_mag_totall are accurate.
    type1_bools = quasar_filter([logLbol_nois, BH_Mass_nois], HSC_Lbol_overall, HSC_MBHs_overall)
    if consider_type1 == True:
        select_bool = select_bool * type1_bools
    Stellar_Mass_nois_sl = Stellar_Mass_nois[select_bool]
    BH_Mass_nois_sl = BH_Mass_nois[select_bool]
    logLbol_nois_sl = logLbol_nois[select_bool]
    # sdss_g_galaxy_sl = sdss_g_galaxy[select_bool]
    EAGLE['Stellar_Mass_nois'] = Stellar_Mass_nois
    EAGLE['Stellar_Mass_nois_sl'] = Stellar_Mass_nois_sl
    EAGLE['BH_Mass_nois'] = BH_Mass_nois
    EAGLE['BH_Mass_nois_sl'] = BH_Mass_nois_sl
    EAGLE['logLbol_nois'] = logLbol_nois
    EAGLE['logLbol_nois_sl'] = logLbol_nois_sl
    # EAGLE['sdss_g_galaxy'] = sdss_g_galaxy  #!!! #TODO: Change to r for z<0.5
    # EAGLE['sdss_g_galaxy_sl'] = sdss_g_galaxy_sl
    EAGLE['sdss_g_pointsource'] = g_pointsource
    EAGLE['sdss_g_pointsource_sl'] = g_pointsource[select_bool]
    # EAGLE['select_abs_Mags'] = abs_Mags
    EAGLE['select_bool'] = select_bool
    return EAGLE

def TNG_set(filename, HSC_Lbol_overall, HSC_MBHs_overall, I_mag_break = 20.5, consider_type1= True, imf ='Cha'):
    TNG = {}
    zs = float(filename.split("_z")[1][:4])
    TNG['zs'] = zs
    # Stellar_Mass, BH_Mass, sdss_i_galaxy, sdss_g_galaxy, sdss_r_galaxy, sdss_i_pointsource, sdss_g_pointsource, Eddington_ratio = np.load(filename)
    BH_Mass, Stellar_Mass, StellarMass_30kpc, sdss_i_stellar, sdss_g_stellar, sdss_r_stellar, sdss_i_pointsource, sdss_g_pointsource, BH_Lbol, Eddington_ratio = np.load(filename)
    BH_Mass = np.log10(BH_Mass)
    Stellar_Mass = np.log10(Stellar_Mass)
    Eddington_ratio = np.log10(Eddington_ratio)
    
    TNG['Stellar_Mass'] = Stellar_Mass
    TNG['BH_Mass'] = BH_Mass
    TNG['sdss_g_pointsource'] = sdss_g_pointsource
    TNG['Eddington_ratio'] = Eddington_ratio
    dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #zs is a list of redshift of the sources (input as 1 D numpy array)
    #Transfer the TNG magnitude to absolute ones:
    i_wavelen = 7496
    # g_wavelen = 4700
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
    abs_Mags = I_mag_break -5*(np.log10(dl)-1)   #dl is the luminosity distance which is a function of redshift:
    if zs >= 0.5:    
        sdss_mag = sdss_g_pointsource #-2.5*np.log10(10**(-0.4*sdss_g_galaxy) + 10**(-0.4*sdss_g_pointsource))
    elif zs<0.5:
        sdss_mag= abs_sdss_r_pointsource #-2.5*np.log10(10**(-0.4*sdss_r_galaxy) + 10**(-0.4*abs_sdss_r_pointsource))
    dMBH, dmag, dMstar, dLbol= 0.4, 0.3, 0.2, 0.1  #dmag is for host magnitude. 
    # dMBH, dmag, dMstar, dLbol= 0.001, 0.001, 0.001, 0.001  #dmag is for host magnitude. 
    logLedd = 38. + np.log10(1.2) + BH_Mass
    logLbol = logLedd + Eddington_ratio
    TNG['logLbol'] = logLbol
    
    z_range = np.arange(0.2, 1.0, 0.05)
    mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
    if imf == 'Sal':
        mstar_cut_range = mstar_cut_range+0.23
    mstar_cut = mstar_cut_range[zs > z_range][-1]
    Stellar_Mass_nois = Stellar_Mass + np.random.normal(0, dMstar, size=Stellar_Mass.shape)
    BH_Mass_nois = BH_Mass + np.random.normal(0, dMBH, size=BH_Mass.shape)
    logLbol_nois = logLbol + np.random.normal(0, dLbol, size=logLbol.shape)
    select_bool = (sdss_mag <abs_Mags) * (Stellar_Mass_nois >mstar_cut) #!!!Assuming sdss_mag_totall are accurate.
    type1_bools = quasar_filter([logLbol_nois, BH_Mass_nois], HSC_Lbol_overall, HSC_MBHs_overall)
    if consider_type1 == True:
        select_bool = select_bool * type1_bools
    Stellar_Mass_nois_sl = Stellar_Mass_nois[select_bool]
    BH_Mass_nois_sl = BH_Mass_nois[select_bool]
    logLbol_nois_sl = logLbol_nois[select_bool]
    sdss_g_galaxy_sl = sdss_g_stellar[select_bool]
    TNG['Stellar_Mass_nois'] = Stellar_Mass_nois
    TNG['Stellar_Mass_nois_sl'] = Stellar_Mass_nois_sl
    TNG['BH_Mass_nois'] = BH_Mass_nois
    TNG['BH_Mass_nois_sl'] = BH_Mass_nois_sl
    TNG['logLbol_nois'] = logLbol_nois
    TNG['logLbol_nois_sl'] = logLbol_nois_sl
    TNG['sdss_g_galaxy'] = sdss_g_stellar  #!!! #TODO: Change to r for z<0.5
    TNG['sdss_g_galaxy_sl'] = sdss_g_galaxy_sl
    TNG['sdss_g_pointsource'] = sdss_g_pointsource
    TNG['sdss_g_pointsource_sl'] = sdss_g_pointsource[select_bool]
    TNG['select_abs_Mags'] = abs_Mags
    TNG['select_bool'] = select_bool
    return TNG

def SAM_set(filename, zs, HSC_Lbol_overall, HSC_MBHs_overall, I_mag_break = [20.5,22], 
            consider_type1= True, imf = 'Sal'):
    SAM = {}
    text  = np.loadtxt(filename)
    text = text[text[:,5] != 0]
    redshift = text[:, 0]
    logMBH =  text[:, 1]
    logM_mass =  text[:, 2]
    g_galaxy =  text[:, 3]
    r_galaxy =  text[:, 4]
    AGN_bol_10_45 =  text[:, 5]
    stat_weight =  text[:, 6]
    L_5100_10_45 = AGN_bol_10_45 / 9.26
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
    # g_totall = g_pointsource
    wave_lam_r = 6100 * (1+redshift)  #rest-frame in A  
    F_lam_r = F_lam_5100 / 5100**(-(-0.44+2)) * wave_lam_r **(-(-0.44+2)) #erg/s/cm^2/A
    F_mu_r = F_lam_r *  (wave_lam_r) **2 / (2.9979 * 10**18)
    obs_mag_r = [(-2.5 * np.log10(F_mu_r[i]) - 48.60) for i in range(len(F_mu_r))]
    r_pointsource = obs_mag_r - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    # r_totall = -2.5*np.log10(10**(-0.4*r_galaxy) + 10**(-0.4*r_pointsource))
    # r_totall = r_pointsource    
    #Transfer the magnitude to absolute ones:
    I_mag_break = np.ones_like(redshift) * I_mag_break[0]
    I_mag_break[redshift>0.5] = I_mag_break[1]
    agn_mag = copy.deepcopy(r_pointsource)
    agn_mag[redshift>0.5] = copy.deepcopy(g_pointsource[redshift>0.5])  #agn mag is r when z<0.5 and g when z>0.5
    galaxy_mag = copy.deepcopy(r_galaxy)
    galaxy_mag[redshift>0.5] = copy.deepcopy(g_galaxy[redshift>0.5])  #agn mag is r when z<0.5 and g when z>0.5
    
    abs_Mags_break = I_mag_break - 5*(np.log10(dl * 10**6)-1)   #Compare with and with r when z<0.5,  g when z >0.5
    if zs<1:
        redshift_bool = (redshift<(zs+0.1)) * (redshift>(zs-0.1))
    elif zs>1:
        redshift_bool = (redshift<(zs+0.3)) * (redshift>(zs-0.3))
    weight = stat_weight[redshift_bool]
    logMBH = logMBH[redshift_bool]
    AGN_bol_10_45 = AGN_bol_10_45[redshift_bool]
    logM_mass = logM_mass[redshift_bool]
    galaxy_mag = galaxy_mag[redshift_bool]
    abs_Mags_break = abs_Mags_break[redshift_bool]
    agn_mag = agn_mag[redshift_bool]
    redshift = redshift[redshift_bool]
    dMBH, dmag, dMstar,dLbol= 0.4, 0.3, 0.2, 0.1 #dmag is for host magnitude. 
    # dMBH, dmag, dMstar, dLbol= 0.001, 0.001, 0.001, 0.001  #dmag is for host magnitude. 

    z_range = np.arange(0.2, 1.0, 0.05)
    mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
    if imf == 'Sal':
        mstar_cut_range = mstar_cut_range+0.23
    mstar_cut = mstar_cut_range[zs > z_range][-1]
      
    c_weight = np.ones_like(weight)
    for i in range(len(weight)):
        c_weight[i] = np.sum(weight[:i+1])    
    Stellar_Mass_reali, BH_Mass_reali, agn_mag_reali,galaxy_mag_reali, abs_Mags_break_reali, redshift_reali, AGN_bol_reali  = [], [], [], [], [], [],[]
    for i in range(3000):
        seed = np.random.uniform(0, c_weight[-1])
        j = np.sum(seed > c_weight)  #Random based on the weighting
        Stellar_Mass_reali.append(logM_mass[j] )
        BH_Mass_reali.append(logMBH[j])
        agn_mag_reali.append(agn_mag[j])
        galaxy_mag_reali.append(galaxy_mag[j])
        abs_Mags_break_reali.append(abs_Mags_break[j])
        AGN_bol_reali.append(np.log10(AGN_bol_10_45[j]) + 45 ) 
        redshift_reali.append(redshift[j])
    Stellar_Mass_reali= np.array(Stellar_Mass_reali)
    BH_Mass_reali= np.array(BH_Mass_reali)
    AGN_bol_reali = np.array(AGN_bol_reali)
    galaxy_mag_reali = np.array(galaxy_mag_reali)
    agn_mag_reali = np.array(agn_mag_reali)
    abs_Mags_break_reali = np.array(abs_Mags_break_reali)
    
    Stellar_Mass_nois = Stellar_Mass_reali + np.random.normal(0, dMstar, size=Stellar_Mass_reali.shape)
    BH_Mass_nois = BH_Mass_reali + np.random.normal(0, dMBH, size =BH_Mass_reali.shape)
    AGN_bol_nois  = AGN_bol_reali + np.random.normal(0, dLbol, size = AGN_bol_reali.shape)
    
    select_bool = (Stellar_Mass_nois>mstar_cut) * (agn_mag_reali<abs_Mags_break_reali) 
    # Lbol_boll = (AGN_bol_reali>43.6) * (AGN_bol_reali<46.2)  #!!! 
    type1_bools = quasar_filter([AGN_bol_nois, BH_Mass_nois], HSC_Lbol_overall, HSC_MBHs_overall)
    if consider_type1 == True:
        select_bool = select_bool*type1_bools
    print(np.sum(select_bool ), select_bool.shape)
    
    BH_Mass_nois_sl = BH_Mass_nois[select_bool]
    Stellar_Mass_nois_sl = Stellar_Mass_nois[select_bool]
    AGN_bol_nois_sl = AGN_bol_nois[select_bool]

    SAM['Stellar_Mass_reali'] = Stellar_Mass_reali
    SAM['BH_Mass_reali'] = BH_Mass_reali    
    SAM['logLbol_reali'] = AGN_bol_reali
    
    SAM['Stellar_Mass_nois'] = Stellar_Mass_nois
    SAM['Stellar_Mass_nois_sl'] = Stellar_Mass_nois_sl
    SAM['BH_Mass_nois'] = BH_Mass_nois
    SAM['BH_Mass_nois_sl'] = BH_Mass_nois_sl
    SAM['logLbol_nois'] = AGN_bol_nois
    SAM['logLbol_nois_sl'] = AGN_bol_nois_sl
    SAM['sdss_mag_galaxy_reali'] = galaxy_mag_reali
    SAM['sdss_mag_galaxy_sl'] = agn_mag_reali[select_bool]
    SAM['sdss_mag_pointsource_reali'] = agn_mag_reali
    SAM['sdss_mag_pointsource_sl'] = agn_mag_reali[select_bool]
    SAM['select_bool'] = select_bool
    SAM['abs_Mags_break_reali'] = abs_Mags_break_reali
    return SAM

def quasar_filter(group_list, HSC_Lbol_overall, HSC_MBHs_overall):
    xmin, xmax = int(HSC_Lbol_overall.min()-2), int(HSC_Lbol_overall.max()+3)
    ymin, ymax = int(HSC_MBHs_overall.min()-2), int(HSC_MBHs_overall.max()+3)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    values = np.vstack([HSC_Lbol_overall, HSC_MBHs_overall])
    kernel = st.gaussian_kde(values)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    f = np.reshape(kernel(positions).T, xx.shape)
    from matplotlib import cm 
    t = [kernel.pdf([HSC_Lbol_overall[i] , HSC_MBHs_overall[i]]) for i in range(len(HSC_Lbol_overall))]
    min_pdf = np.min(t)
    # print(np.mean(min_pdf))
    bools = [ (kernel.pdf([group_list[0][i] , group_list[1][i]]) > min_pdf/2)[0] for i in range(len(group_list[0]))]  
    # plt.scatter(HSC_Lbol_overall, HSC_MBHs_overall,c='blue')
    # plt.contourf(xx, yy, f, cmap=cm.Blues, alpha=0.5)
    # plt.scatter(group_list[0], group_list[1],c='gray',alpha=0.4)
    # plt.scatter(group_list[0][bools], group_list[1][bools],c='green',alpha=0.4)
    return np.array(bools)

def comp_plot(x, y, label_x='', label_y='',  alpha=1, label = '', c = 'blue', edgecolors=None):
    plt.figure(figsize=(11.5,10))      
    plt.scatter(x, y, alpha=alpha,label=label,c=c, edgecolors=edgecolors)
    plt.xlabel(label_x,fontsize=30)
    plt.ylabel(label_y, fontsize=25)
    plt.tick_params(labelsize=25)
