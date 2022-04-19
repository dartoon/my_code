#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 10:12:49 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy

# # ---- RA ---  --- Dec ---  - Magnitudes -  - Discovery Spectroscopy -      M_dyn       M_bh
# # hh mm ss.ss   dd mm ss.s   zAB     yAB    ID  z/SpT    M1450   Paper  # (1e9 Msun)  (1e9 Msun)
# 22 16 44.47  -00 16 50.1  22.775  22.942   Q  6.10  -23.82 0.04   I   #     82        0.70
# 11 52 21.27  +00 55 36.6  21.827  21.614   Q  6.37  -25.31 0.04   I   #     14        0.63
# 08 59 07.19  +00 22 55.9  22.768  23.649   Q  6.39  -24.09 0.07   I   #     56        0.038
# 12 08 59.23  -02 00 34.8  22.121  22.045   Q  6.2   -24.73 0.02  II   #     13        0.71
# 22 39 47.47  +02 07 47.5  22.402  22.326   Q  6.26  -24.69 0.04  II   #    290        1.1


target_info = {'J223644.58+003256.9': {'z':6.4, 'M1450': -23.8, 'yABmag': 23.2, 'zABmag': 23.777,'T_NIRSpec':1.5, 'T_NIRCam': 1.8 },
'J225538.04+025126.6': {'z':6.34, 'M1450': -23.9, 'yABmag': 23.0, 'zABmag': 23.435,'T_NIRSpec':1.4, 'T_NIRCam': 2.3 },
'J152555.79+430324.0': {'z':6.27, 'M1450': -23.9, 'yABmag': 23.5, 'zABmag': 23.375,'T_NIRSpec':2.2, 'T_NIRCam': 1.8 },
'J114648.42+012420.1': {'z':6.27, 'M1450': -23.7, 'yABmag': 23.8, 'zABmag': 22.984,'T_NIRSpec':1.8, 'T_NIRCam': 1.9 },
'J084431.60-005254.6': {'z':6.25, 'M1450': -23.7, 'yABmag': 23.1, 'zABmag': 23.201,'T_NIRSpec':1.5, 'T_NIRCam': 1.8 },
'J021721.59-020852.6': {'z':6.2, 'M1450': -23.2, 'yABmag': 23.5, 'zABmag': 23.882,'T_NIRSpec':1.8, 'T_NIRCam': 2.3 },
'J091833.17+013923.4': {'z':6.19, 'M1450': -23.7, 'yABmag': 23.2, 'zABmag': 23.143,'T_NIRSpec':1.4, 'T_NIRCam': 1.9 },
'J151248.71+442217.5': {'z':6.18, 'M1450': -23.1, 'yABmag': 24.2, 'zABmag': 23.648,'T_NIRSpec':3.3, 'T_NIRCam': 1.9 },
'J084408.61-013216.5': {'z':6.18, 'M1450': -24.0, 'yABmag': 23.7, 'zABmag': 23.310,'T_NIRSpec':1.8, 'T_NIRCam': 2.3 },
'J142517.72-001540.8': {'z':6.18, 'M1450': -23.4, 'yABmag': 23.4, 'zABmag': 22.833,'T_NIRSpec':1.4, 'T_NIRCam': 2.3 },
'J114658.89-000537.7': {'z':6.3, 'M1450': -21.5, 'yABmag': 24.8, 'zABmag': 23.741,'T_NIRSpec':4.1, 'T_NIRCam': 1.8 },
'J091114.27+015219.4': {'z':6.07, 'M1450': -22.1, 'yABmag': 24.6, 'zABmag': 24.222,'T_NIRSpec':3.8, 'T_NIRCam': 1.8 }}


keys = []
for key in target_info.keys():
    keys.append(key)
    
#%% Use redshift and M1450 to estimate the AGN magnitude in F356w and F150w.
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
def cal_filt_flam(array_spec, fil):
    '''
    Calculate a filter f_lambda
    Parameter
    --------
    Return
    --------
        Filter flux
    '''
    fil[:,1] /= fil[:,1].max()
    f_filt = interp1d(fil[:,0], fil[:,1], kind='cubic')
    int_spec = array_spec[(array_spec[:,0]>fil[0,0]) * (array_spec[:,0]<fil[-1,0])]
    int_flux = 0
    int_filt = 0
    for i in range(len(int_spec)-1):
        int_flux += int_spec[i,1] * f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
        int_filt += f_filt(int_spec[i,0]) * (int_spec[i+1,0]-int_spec[i,0])
    filt_flam = int_flux/ int_filt   
    return filt_flam

#%%Load AGN typtical template
filename_p  = '/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/sdss_temp/spDR2-029.fit'  #QSO ID: 30
# filename_p  = 'sdss_temp/spDR2-030.fit'  #QSO with some BAL activity (QSO)
# filename_p  = 'sdss_temp/spDR2-031.fit'  #QSO with some BAL activity (QSO)
# filename_p  = 'sdss_temp/spDR2-032.fit'  #High-luminosity QSO (QSO) !!! Not work, no short Wavelen information
hdul = pyfits.open(filename_p)
spec_header = hdul[0].header
spec_data = hdul[0].data
wavelen_first_pi = spec_header['COEFF0']
delta_wave = spec_header['COEFF1']
array_spec = np.zeros((len(spec_data.T),2))
array_spec[0,0]  = wavelen_first_pi
for i in range(len(spec_data.T)-1):
    array_spec[i+1,0] = array_spec[i, 0] + delta_wave
array_spec[:,0] = 10** array_spec[:,0]      
array_spec[:,1] = np.sum(spec_data, axis=0)
array_spec = array_spec[100:]
array_spec_temp = copy.deepcopy(array_spec)
#Plot the Spec at the rest-frame
array_spec = copy.deepcopy(array_spec_temp)
plt.figure(figsize=(10, 6))    
plt.plot(array_spec[:,0], array_spec[:,1])
idx = np.where( abs(array_spec[:,0] - 1450.) == np.min(abs(array_spec[:,0] - 1450.)) )[0][0]
plt.scatter(1450, array_spec[idx,1], c='r')
plt.tick_params(labelsize=15)
plt.xlabel("A",fontsize=27)
plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/$\AA$",fontsize=27)
plt.show()

#%%Redshift and norm the spec.
for key in keys:
    z = target_info[key]['z']
    M1450 = target_info[key]['M1450']
    array_spec = copy.deepcopy(array_spec_temp)
    array_spec[:,0] *= (1+z)
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dl = cosmo.luminosity_distance(z=z).value  
    #m1450 = M1450 + 2.5*np.log10((dl*10**6/10)**2)   # Mv = m - 2.5 log[ (d/10)**2 ],  d in pc scale
    m1450 = M1450 + 5*(np.log10(dl*10**6)-1) - 2.5*np.log10(1+z)
    lam = 1450 * (1+z)
    f_lambda_1450 = 10**(-0.4*(m1450+2.402+5.0*np.log10(lam)))
    idx = np.where( abs(array_spec[:,0] - lam) == np.min(abs(array_spec[:,0] - lam)) )[0][0]
    norm = f_lambda_1450/array_spec[idx, 1]
    array_spec[:,1] *= norm
    
    z_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/317.fil')
    z_fil = z_fil[:,1:]
    z_fil[:,1] = z_fil[:,1]/z_fil[:,1].max()* array_spec[:,1].max()/6
    y_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/318.fil')
    y_fil = y_fil[:,1:]
    y_fil[:,1] = y_fil[:,1]/y_fil[:,1].max()* array_spec[:,1].max()/6

    plt.figure(figsize=(10, 6))    
    plt.plot(array_spec[:,0], array_spec[:,1])
    plt.plot(z_fil[:,0], z_fil[:,1], label='HSC Z filter response', c='green')
    plt.plot(y_fil[:,0], y_fil[:,1], label='HSC Y filter response', c='orange')
    plt.scatter(lam, f_lambda_1450, c='r', marker="*", s =200, zorder = 10, label='observed M1450')

    HSC_z_lam = 8908.2
    HSC_y_lam = 9775.1
    z_filt_flam = cal_filt_flam(array_spec, z_fil)
    plt.scatter(HSC_z_lam, z_filt_flam, c='green',  edgecolors='k', s =100, zorder = 100)
    y_filt_flam = cal_filt_flam(array_spec, y_fil)
    plt.scatter(HSC_y_lam, y_filt_flam, c='orange', edgecolors='k', s =100, zorder = 100)

    m_zAB = target_info[key]['zABmag']
    m_yAB = target_info[key]['yABmag']
    obs_z_filt_flam = 10**(-0.4*(m_zAB+2.402+5.0*np.log10(HSC_z_lam)))
    obs_y_filt_flam = 10**(-0.4*(m_yAB+2.402+5.0*np.log10(HSC_y_lam)))
    plt.scatter(HSC_z_lam, obs_z_filt_flam, c='green', marker="*", s =200, zorder = 10, label='observed z_mag')
    plt.scatter(HSC_y_lam, obs_y_filt_flam, c='orange', marker="*", s =200, zorder = 10, label='observed y_mag')
    
    lambda_mean = {'f150w':15104.23, 'f200w':20028.15, 'f356w':35934.49, 'f444w':44393.50}
    f444w_fil = np.loadtxt('/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F444W.dat')
    f444w_fil[:,1] = f444w_fil[:,1]/f444w_fil[:,1].max()* array_spec[:,1].max()/6
    plt.plot(f444w_fil[:,0], f444w_fil[:,1], label='JWST F444W filter response', c='firebrick')
    f444w_lam = lambda_mean['f444w'] # np.median(f444w_fil[:,0])
    f444w_filt_flam = cal_filt_flam(array_spec, f444w_fil) 
    f444w_mag = -2.5 * np.log10(f444w_filt_flam ) - 2.402 - 5.0 * np.log10(f444w_lam)

    f356w_fil = np.loadtxt('/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F356W.dat')
    f356w_fil[:,1] = f356w_fil[:,1]/f356w_fil[:,1].max()* array_spec[:,1].max()/6
    plt.plot(f356w_fil[:,0], f356w_fil[:,1], label='JWST F356W filter response', c='tomato')
    f356w_lam = lambda_mean['f356w'] # np.median(f356w_fil[:,0])
    f356w_filt_flam = cal_filt_flam(array_spec, f356w_fil) 
    f356w_mag = -2.5 * np.log10(f356w_filt_flam ) - 2.402 - 5.0 * np.log10(f356w_lam)

    f200w_fil = np.loadtxt('/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F200W.dat')
    f200w_fil[:,1] = f200w_fil[:,1]/f200w_fil[:,1].max()* array_spec[:,1].max()/6
    plt.plot(f200w_fil[:,0], f200w_fil[:,1], label='JWST F200W filter response', c='lime')
    f200w_lam = lambda_mean['f200w'] #np.median(f200w_fil[:,0])
    f200w_filt_flam = cal_filt_flam(array_spec, f200w_fil) 
    f200w_mag = -2.5 * np.log10(f200w_filt_flam ) - 2.402 - 5.0 * np.log10(f200w_lam)

    f150w_fil = np.loadtxt('/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F150W.dat')
    f150w_fil[:,1] = f150w_fil[:,1]/f150w_fil[:,1].max()* array_spec[:,1].max()/6
    plt.plot(f150w_fil[:,0], f150w_fil[:,1], label='JWST F150W filter response', c='cyan')
    f150w_lam = lambda_mean['f150w'] # np.median(f150w_fil[:,0])
    f150w_filt_flam = cal_filt_flam(array_spec, f150w_fil) 
    f150w_mag = -2.5 * np.log10(f150w_filt_flam ) - 2.402 - 5.0 * np.log10(f150w_lam)

    plt.tick_params(labelsize=15)
    plt.xlabel("$\AA$",fontsize=27)
    plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/A",fontsize=27)
    plt.xlim(800, 60000)
    plt.legend(prop={'size':15})
    plt.show()
    z_mag = -2.5 * np.log10(z_filt_flam ) - 2.402 - 5.0 * np.log10(HSC_z_lam)
    y_mag = -2.5 * np.log10(y_filt_flam ) - 2.402 - 5.0 * np.log10(HSC_y_lam)
    # print("\n\n================")
    # print("ID", key, "Input:")
    # print("M_1450:", round(M1450,3), "redshift:", z)
    # print("================")
    # print("Estimate AGN total mag")
    # print("z_mag", round(z_mag,3))
    # print("y_mag", round(y_mag,3))
    # print('AGN F444W mag', round(f444w_mag,3))
    # print('AGN F356W mag', round(f356w_mag,3))
    # print('AGN F200W mag', round(f200w_mag,3))
    # print('AGN F150W mag', round(f150w_mag,3))
    target_info[key]['mag_f356w'] = f356w_mag
    target_info[key]['mag_f150w'] = f150w_mag
    target_info[key]['mag_HSCz'] = z_mag
    target_info[key]['mag_HSCy'] = y_mag
# import pickle
# pickle.dump(target_info, open('target_info.pkl', 'wb')) 

#%%
plt.figure(figsize=(6,5))
for key in keys:
    plt.scatter(target_info[key]['T_NIRCam'], target_info[key]['M1450'])

plt.figure(figsize=(6,5))
    
for key in target_info.keys():
    # print(key,'redshift:', target_info[key]['z'], ', mag:' , round(target_info[key]['mag_'+filt],2) )
    print(key,'redshift:', target_info[key]['z'], ', HSC_obs_ymag:' , round(target_info[key]['yABmag'],2),'predict' ,round(target_info[key]['mag_HSCy'],2)  )
    plt.scatter(target_info[key]['yABmag'], target_info[key]['mag_HSCy'])
    plt.plot(np.linspace(22.5,25), np.linspace(22.5,25))
plt.xlabel('HSC obs y mag')
plt.ylabel('HSC predict y mag')
plt.show()

plt.figure(figsize=(6,5))
for key in target_info.keys():
    # print(key,'redshift:', target_info[key]['z'], ', mag:' , round(target_info[key]['mag_'+filt],2) )
    # print(key,'redshift:', target_info[key]['z'], ', HSC_obs_ymag:' , round(target_info[key]['yABmag'],2),'predict' ,round(target_info[key]['mag_HSCy'],2)  )
    plt.scatter(target_info[key]['zABmag'], target_info[key]['mag_HSCz'])
    plt.plot(np.linspace(22.5,25), np.linspace(22.5,25))
plt.xlabel('HSC obs z mag')
plt.ylabel('HSC predict z mag')
plt.show()
  