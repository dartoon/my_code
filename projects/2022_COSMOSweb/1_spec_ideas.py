#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 10:20:50 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM


# # 24      "Early-type" galaxy (Galaxy)    Template #24    
# # 25      Galaxy  Template #25    
# # 26      Galaxy  Template #26    
# # 27      Galaxy  Template #27    
# # 28      "Late-type" galaxy (Galaxy)
# # 30	QSO
# i = 30
# #Download from SDSS: http://classic.sdss.org/dr5/algorithms/spectemplates/
# filename_p  = '../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/sdss_temp/spDR2-0{0}.fit'.format(i-1)  #QSO ID: 30
# hdul = pyfits.open(filename_p)
# spec_header = hdul[0].header
# spec_data = hdul[0].data
# wavelen_first_pi = spec_header['COEFF0']
# delta_wave = spec_header['COEFF1']
# array_spec = np.zeros((len(spec_data.T),2))
# array_spec[0,0]  = wavelen_first_pi
# for i in range(len(spec_data.T)-1):
#     array_spec[i+1,0] = array_spec[i, 0] + delta_wave
# array_spec[:,0] = 10** array_spec[:,0]      
# array_spec[:,1] = np.sum(spec_data, axis=0)
# array_spec = array_spec[100:]
# array_spec_temp = copy.deepcopy(array_spec)


import glob
filt = 444
# filename  = '../Sim_HST_JWST/JWST_QSO/esti_galaxy_magnitide/f{0}w/summary_{0}{1}_PA00.fits'.format(filt,3)
# filename = glob.glob(filename)[0]
# hdul_sum = pyfits.open(filename)
# name_sum = hdul_sum[1].columns
# table_sum = hdul_sum[1].data
# spec1d = hdul_sum = pyfits.open("/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/templates/gsf_spec_{0}{1}.fits".format(filt,3))  
# name_spec = spec1d[1].columns
# table_spec = spec1d[1].data
# spec = table_spec['f_model_50']
# # wavel = table_spec['wave_model']/10000.
array_spec= []
z_org = 6.39
array_spec.append( table_spec['wave_model']/(1+z_org)  )
array_spec.append(spec)
array_spec = np.array(array_spec)
array_spec = array_spec.T
array_spec_temp = copy.deepcopy(array_spec)

#%%Plot the Spec at the rest-frame
array_spec = copy.deepcopy(array_spec_temp)
plt.figure(figsize=(10, 6))    
plt.plot(array_spec[:,0], array_spec[:,1])
idx = np.where( abs(array_spec[:,0] - 1450.) == np.min(abs(array_spec[:,0] - 1450.)) )[0][0]
plt.scatter(1450, array_spec[idx,1], c='r')
plt.tick_params(labelsize=15)
plt.xlabel("$\AA$",fontsize=27)
plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/$\AA$",fontsize=27)
plt.show()

#%%
# #%%Shift to redshift z = 6 and norm based on M1450 value
# ID = 3 #1, 2, 3, 4, 5
# #The data information of the five system
# if ID ==1:
#     M1450, m_zAB, m_yAB, z, M_dyn = -23.82, 22.775, 22.942, 6.10, 82.
# elif ID ==2:    
#     M1450, m_zAB, m_yAB, z, M_dyn = -25.31, 21.827, 21.614, 6.37, 14.
# elif ID ==3:
#     M1450, m_zAB, m_yAB, z, M_dyn = -24.09, 22.768, 23.649, 6.39, 56.
# elif ID ==4:    
#     M1450, m_zAB, m_yAB, z, M_dyn = -24.73, 22.121, 22.045, 6.2, 13.
# elif ID ==5:
#     M1450, m_zAB, m_yAB, z, M_dyn = -24.69, 22.402, 22.326, 6.26, 290.
# elif ID ==6:
#     M1450, m_zAB, m_yAB, z, M_dyn = -26.5, 22.402, 22.326, 6.26, 200.  #(M_dyn as 2*10^11, i.e. 200*10^9 )

M1450, z= -24.09, 6

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

# z_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/317.fil')
# z_fil = z_fil[:,1:]
# z_fil[:,1] = z_fil[:,1]/z_fil[:,1].max()* array_spec[:,1].max()/6
# y_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/318.fil')
# y_fil = y_fil[:,1:]
# y_fil[:,1] = y_fil[:,1]/y_fil[:,1].max()* array_spec[:,1].max()/6
plt.figure(figsize=(10, 6))    
plt.plot(array_spec[:,0], array_spec[:,1])
# plt.plot(z_fil[:,0], z_fil[:,1], label='HSC Z filter response', c='green')
# plt.plot(y_fil[:,0], y_fil[:,1], label='HSC Y filter response', c='orange')
# plt.scatter(lam, f_lambda_1450, c='r', marker="*", s =200, zorder = 10, label='observed M1450')

HSC_z_lam = 8908.2
HSC_y_lam = 9775.1

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


# z_filt_flam = cal_filt_flam(array_spec, z_fil)
# plt.scatter(HSC_z_lam, z_filt_flam, c='green')
# y_filt_flam = cal_filt_flam(array_spec, y_fil)
# plt.scatter(HSC_y_lam, y_filt_flam, c='orange')
# obs_z_filt_flam = 10**(-0.4*(m_zAB+2.402+5.0*np.log10(HSC_z_lam)))
# obs_y_filt_flam = 10**(-0.4*(m_yAB+2.402+5.0*np.log10(HSC_y_lam)))

lambda_mean = {'f150w':15104.23, 'f200w':20028.15, 'f356w':35934.49, 'f444w':44393.50,
               'f277w': 27844.64, 'f770w': 77111.39,'f115w':11623.89}
f770w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_MIRI.F770W.dat')
f770w_fil[:,1] = f770w_fil[:,1]/f770w_fil[:,1].max()* array_spec[:,1].max()/6
plt.plot(f770w_fil[:,0], f770w_fil[:,1], label='JWST F770W filter response', c='firebrick')
f770w_lam = lambda_mean['f770w'] # np.median(f770w_fil[:,0])
# f770w_filt_flam = cal_filt_flam(array_spec, f770w_fil) 
# f770w_mag = -2.5 * np.log10(f770w_filt_flam ) - 2.402 - 5.0 * np.log10(f770w_lam)

f444w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F444W.dat')
f444w_fil[:,1] = f444w_fil[:,1]/f444w_fil[:,1].max()* array_spec[:,1].max()/6
plt.plot(f444w_fil[:,0], f444w_fil[:,1], label='JWST F444W filter response', c='tomato')
f444w_lam = lambda_mean['f444w'] # np.median(f444w_fil[:,0])
f444w_filt_flam = cal_filt_flam(array_spec, f444w_fil) 
f444w_mag = -2.5 * np.log10(f444w_filt_flam ) - 2.402 - 5.0 * np.log10(f444w_lam)

# f356w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F356W.dat')
# f356w_fil[:,1] = f356w_fil[:,1]/f356w_fil[:,1].max()* array_spec[:,1].max()/6
# plt.plot(f356w_fil[:,0], f356w_fil[:,1], label='JWST F356W filter response', c='tomato')
# f356w_lam = lambda_mean['f356w'] # np.median(f356w_fil[:,0])
# f356w_filt_flam = cal_filt_flam(array_spec, f356w_fil) 
# f356w_mag = -2.5 * np.log10(f356w_filt_flam ) - 2.402 - 5.0 * np.log10(f356w_lam)
f277w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F277W.dat')
f277w_fil[:,1] = f277w_fil[:,1]/f277w_fil[:,1].max()* array_spec[:,1].max()/6
plt.plot(f277w_fil[:,0], f277w_fil[:,1], label='JWST F277W filter response', c='y')
f277w_lam = lambda_mean['f277w'] # np.median(f277w_fil[:,0])
f277w_filt_flam = cal_filt_flam(array_spec, f277w_fil) 
f277w_mag = -2.5 * np.log10(f277w_filt_flam ) - 2.402 - 5.0 * np.log10(f277w_lam)

f150w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F150W.dat')
f150w_fil[:,1] = f150w_fil[:,1]/f150w_fil[:,1].max()* array_spec[:,1].max()/6
plt.plot(f150w_fil[:,0], f150w_fil[:,1], label='JWST F150W filter response', c='lime')
f150w_lam = lambda_mean['f150w'] # np.median(f150w_fil[:,0])
f150w_filt_flam = cal_filt_flam(array_spec, f150w_fil) 
f150w_mag = -2.5 * np.log10(f150w_filt_flam ) - 2.402 - 5.0 * np.log10(f150w_lam)

f115w_fil = np.loadtxt('../Sim_HST_JWST/JWST_QSO/esti_filter_magnitude/JWST_NIRCam.F115W.dat')
f115w_fil[:,1] = f115w_fil[:,1]/f115w_fil[:,1].max()* array_spec[:,1].max()/6
plt.plot(f115w_fil[:,0], f115w_fil[:,1], label='JWST F115W filter response', c='cyan')
f115w_lam = lambda_mean['f115w'] #np.median(f115w_fil[:,0])
f115w_filt_flam = cal_filt_flam(array_spec, f115w_fil) 
f115w_mag = -2.5 * np.log10(f115w_filt_flam ) - 2.402 - 5.0 * np.log10(f115w_lam)

plt.text(80000, 4*10**-19, s='z='+str(z), fontsize=18)
plt.tick_params(labelsize=15)
plt.xlabel("$\AA$",fontsize=27)
plt.ylabel(r"f$_\lambda$ erg/s/cm$^2$/A",fontsize=27)
plt.xlim(800, 100000)
#plt.ylim(10**-19.5,10**-17.5)
#plt.ylim(10**-19.5,10**-17.5)
#plt.yscale('log',basey=10) 
plt.legend(prop={'size':15})
plt.show()

# z_mag = -2.5 * np.log10(z_filt_flam ) - 2.402 - 5.0 * np.log10(HSC_z_lam)
# y_mag = -2.5 * np.log10(y_filt_flam ) - 2.402 - 5.0 * np.log10(HSC_y_lam)


# print("ID", ID, "Input:")
# print("m_zAB", m_zAB)
# print("m_yAB", m_yAB)
print("M_1450", round(M1450,3))
print("================")
print("Estimate AGN total mag")
# print("z_mag", round(z_mag,3))
# print("y_mag", round(y_mag,3))
print('AGN F444W mag', round(f444w_mag,3))
print('AGN F356W mag', round(f277w_mag,3))
print('AGN F115W mag', round(f115w_mag,3))
print('AGN F150W mag', round(f150w_mag,3))