#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:44:21 2021

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pickle
# import scipy.stats as st
# import matplotlib as mat
# mat.rcParams['font.family'] = 'STIXGeneral'

def lfit(x,m = 1/0.9811,c = 2.545/0.9811):
    return m*x+c

import copy
ifplot = False
# zs = 0.3
imf_dict = {'MBII': 'Sal', 'SAM': 'Sal', 'TNG100': 'Cha', 'TNG300': 'Cha', 'Illustris': 'Cha', 'Horizon': 'Sal'}
#%% Prep HSC data with dict as Sal and Cha
from prep_comparison import HSC_set
from prep_comparison import TNG_set, SAM_set, Horizon_set   #MBII and Illustris are loaded using TNG_set function
def load_HSC_comp(zs, I_mag_break=None, no_noise = False):
    HSC_dict = {}
    HSC_dict['Sal'] = HSC_set(zs, core = True, imf = 'Sal')
    HSC_dict['Cha'] = HSC_set(zs, core = True, imf = 'Cha')
    HSC_dict['Sal']['Mstar'] = HSC_dict['Sal']['HSC_Mstar']
    HSC_dict['Sal']['MBHs'] = HSC_dict['Sal']['HSC_MBHs']
    HSC_dict['Cha']['Mstar'] = HSC_dict['Cha']['HSC_Mstar']
    HSC_dict['Cha']['MBHs'] = HSC_dict['Cha']['HSC_MBHs']
    I_mag_break_SAM = [20.5,22.0]
    if I_mag_break is None:
        if zs < 0.5:
            I_mag_break = 20.5  #z~0.3
        if zs >= 0.5:    
            I_mag_break = 22.0    #z~0.7
    else:
        I_mag_break_SAM = [I_mag_break, I_mag_break]
    filename_dict = {}
    zs_MBII = 0.3
    if zs!= 0.3:
        zs_MBII = 0.6
        HSC_dict['zs06'] = HSC_set(zs=zs_MBII, core = True, imf = 'Sal')
        HSC_dict['zs06']['Mstar'] = HSC_dict['zs06']['HSC_Mstar']
        HSC_dict['zs06']['MBHs'] = HSC_dict['zs06']['HSC_MBHs']
    filename_dict['MBII'] = glob.glob('MBII_data/*{0}*.npy'.format(zs_MBII))[0]
    filename_dict['SAM'] = 'SAM/catalogue.dat'
    if zs == 0.7:
        filename_dict['SAM'] = 'SAM/catalogue2.dat'
    filename_dict['TNG100'] = glob.glob('TNG100/*{0}*.npy'.format(zs))[0]
    filename_dict['TNG300'] = glob.glob('TNG300_data/*{0}*.npy'.format(zs))[0]
    filename_dict['Illustris'] = glob.glob('Illustris_data/*{0}*.npy'.format(zs))[0]
    Horizon_name_dict= {'0.3': 'Horizon/outt00638_halo_gal_centralBHs_2reff', 
                        '0.5': 'Horizon/outt00552_halo_gal_centralBHs_2reff', 
                        '0.7': 'Horizon/outt00439_halo_gal_centralBHs_2reff'} 
    filename_dict['Horizon'] = Horizon_name_dict['{0}'.format(zs)]
    sim_dict = {}
    sim_dict['MBII'] = TNG_set(filename_dict['MBII'], HSC_Lbol_overall=HSC_dict[imf_dict['MBII']]['HSC_Lbol_overall'], 
                   HSC_MBHs_overall=HSC_dict[imf_dict['MBII']]['HSC_MBHs_overall'],
                  I_mag_break = I_mag_break, imf = imf_dict['MBII'], no_noise=no_noise)
    sim_dict['SAM'] = SAM_set(filename_dict['SAM'], zs=zs, HSC_Lbol_overall=HSC_dict[imf_dict['SAM']]['HSC_Lbol_overall'], 
                  HSC_MBHs_overall=HSC_dict[imf_dict['SAM']]['HSC_MBHs_overall'],
                  I_mag_break = I_mag_break_SAM, imf =  imf_dict['SAM'], no_noise=no_noise)
    
    sim_dict['SAM']['BH_Mass'] = sim_dict['SAM']['BH_Mass_reali']
    sim_dict['SAM']['Stellar_Mass'] = sim_dict['SAM']['Stellar_Mass_reali']
    sim_dict['TNG100'] = TNG_set(filename_dict['TNG100'], HSC_Lbol_overall=HSC_dict[imf_dict['TNG100']]['HSC_Lbol_overall'], 
                  HSC_MBHs_overall=HSC_dict[imf_dict['TNG100']]['HSC_MBHs_overall'],
                  I_mag_break = I_mag_break, imf =  imf_dict['TNG100'], no_noise=no_noise)
    sim_dict['TNG300'] = TNG_set(filename_dict['TNG300'], HSC_Lbol_overall=HSC_dict[imf_dict['TNG300']]['HSC_Lbol_overall'], 
                  HSC_MBHs_overall=HSC_dict[imf_dict['TNG300']]['HSC_MBHs_overall'],
                  I_mag_break = I_mag_break, imf =  imf_dict['TNG300'], no_noise=no_noise)
    sim_dict['Illustris'] = TNG_set(filename_dict['Illustris'], HSC_Lbol_overall=HSC_dict[imf_dict['Illustris']]['HSC_Lbol_overall'], 
                        HSC_MBHs_overall=HSC_dict[imf_dict['Illustris']]['HSC_MBHs_overall'],
                        I_mag_break = I_mag_break, imf =  imf_dict['Illustris'], no_noise=no_noise)
    sim_dict['Horizon'] = Horizon_set(filename_dict['Horizon'], HSC_Lbol_overall=HSC_dict[imf_dict['Horizon']]['HSC_Lbol_overall'], 
                          HSC_MBHs_overall=HSC_dict[imf_dict['Horizon']]['HSC_MBHs_overall'],
                  zs = zs, I_mag_break = I_mag_break, imf =  imf_dict['Horizon'], no_noise=no_noise)
    return HSC_dict, sim_dict
#%%
def prep_sim_hst(bhmass_overall,mstar_overall,Eddr_overall, imf, no_noise=False):
    logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
    Lbol_overall = logLedd_overall + Eddr_overall
    if no_noise == False:
        dMBH, dMstar, dLbol= 0.4, 0.17, 0.1
    elif no_noise == True:  
        dMBH, dMstar, dLbol= 0.00004, 0.000017, 0.000003
    bhmass_overall_noi = bhmass_overall + np.random.normal(0, dMBH, size=bhmass_overall.shape)
    mstar_overall_noi = mstar_overall + np.random.normal(0, dMstar, size=mstar_overall.shape)
    Lbol_overall_noi = Lbol_overall + np.random.normal(0, dLbol, size= Lbol_overall.shape)
    logLedd_overall_noi = 38. + np.log10(1.2) + bhmass_overall_noi
    Eddr_overall_noi = Lbol_overall_noi - logLedd_overall_noi
    select_window = (bhmass_overall_noi>7.7)*(bhmass_overall_noi<8.6)*(Eddr_overall_noi<0.0)*\
    (Eddr_overall_noi > -1.1*(bhmass_overall_noi-7.5)-0.5 )* (Eddr_overall_noi>-1.5)
    bhmass_select_noi = bhmass_overall_noi[select_window]
    mstar_select_noi = mstar_overall_noi[select_window]
    Lbol_select_noi = Lbol_overall_noi[select_window]
    # return bhmass_select_noi, mstar_select_noi, Lbol_select_noi #, bhmass_overall_noi, mstar_overall_noi, Lbol_overall_noi
    sim = {}
    sim['Stellar_Mass'], sim['BH_Mass'], sim['logLbol'] =mstar_overall, bhmass_overall, Lbol_overall
    sim['Stellar_Mass_nois_sl'], sim['BH_Mass_nois_sl'], sim['logLbol_nois_sl'] =mstar_select_noi, bhmass_select_noi, Lbol_select_noi
    return sim

def load_HST_comp(no_noise = False):
    zs = 1.5
    HST_dict={}
    #Import HST sample
    HST_results = pickle.load(open('pickles/HST_saves.pkl','rb'))
    
    Lr, M_star, M_r, M_r_obs_err, bh_mass_obs, stellar_mass_obs_err, logLbol_obs = HST_results
    
    HST_dict = {}
    HST_dict['Cha'] = {}
    HST_dict['Cha']['Mstar'] = M_star
    HST_dict['Cha']['MBHs'] = bh_mass_obs
    HST_dict['Cha']['logLbol'] = logLbol_obs
    HST_dict['Sal'] = copy.deepcopy(HST_dict['Cha'])
    HST_dict['Sal']['Mstar'] = np.log10(10**M_star / 0.54 * 0.91)
    #Import sim sample
    sim_dict={}
    load = ['MBII', 'TNG100', 'TNG300', 'Illustris']
    for i in range(len(load)):
        filename = glob.glob('{1}*/*{0}*npy'.format(zs, load[i]))[0]
        BH_Mass, _, StellarMass_30kpc, _, _, _, _, _, _, Eddington_ratio = np.load(filename)
        
        
        bhmass_overall = np.log10(BH_Mass) #in log
        mstar_overall = np.log10(StellarMass_30kpc) #in log
        Eddr_overall = np.log10(Eddington_ratio) #in log
        if load[i] == 'Illustris':
            print('Illustris working')
            outl_bool = (lfit(bhmass_overall) - 0.5 < mstar_overall)
            print(len(outl_bool), np.sum(outl_bool))
            bhmass_overall, mstar_overall, Eddr_overall = bhmass_overall[outl_bool], mstar_overall[outl_bool],  Eddr_overall[outl_bool]
        sim_dict[load[i]] = prep_sim_hst(bhmass_overall,mstar_overall,Eddr_overall, imf=imf_dict[load[i]],no_noise = no_noise)
    filename = 'Horizon/outt00266_halo_gal_centralBHs_2reff'
    texts = np.loadtxt(filename)
    mstar_overall = np.log10(texts[:, 1])+ 11
    bhmass_overall = np.log10(texts[:, 2]) + 8
    Eddr_overall = np.log10(texts[:, 3])
    sim_dict['Horizon'] = prep_sim_hst(bhmass_overall,mstar_overall,Eddr_overall, imf=imf_dict['Horizon'],no_noise = no_noise)
    _HSC = HSC_set(zs)
    
    SAM = SAM_set('SAM/catalogue.dat', zs=zs, HSC_Lbol_overall=_HSC['HSC_Lbol_overall'], HSC_MBHs_overall=_HSC['HSC_MBHs_overall'],
                  I_mag_break = [20.5 +5 ,22.0 + 5 ])  #Will take the overall, just quickly load SAM sample. Selection will do later.
    bhmass_overall = SAM['BH_Mass_reali']
    mstar_overall = SAM['Stellar_Mass_reali']
    logLedd_overall = 38. + np.log10(1.2) + bhmass_overall
    Lbol_overall = SAM['logLbol_reali']
    Eddr_overall =  (Lbol_overall - logLedd_overall)
    sim_dict['SAM'] = prep_sim_hst(bhmass_overall,mstar_overall,Eddr_overall, imf=imf_dict['SAM'],no_noise = no_noise)
    return HST_dict, sim_dict

