#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 18:24:10 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info, name_list
import glob
import pickle
remove_id = [24, 55]

HST_filt = 'F160W'
AGN_mag_list = []
host_mag_list = []
# for idx in range(59):  #CEERS
for idx in [35]:  #CEERS
    fit_run_F150_list = []
    fit_run_HST_list = []
    
    if idx in remove_id:
        continue
    target_id, z = load_info(idx)
    print(idx)
    # root_folder = '../*/*'  #Include HST
    files_F150 = glob.glob('fit_material/fit_run_idx{0}_F150W_psf*.pkl'.format(idx)) 
    files_HST = glob.glob('../*/fit_material/fit_run_idx{0}_{1}_psf*.pkl'.format(idx, HST_filt)) 
    if files_F150 == [] or files_HST == []:
        continue
    
    for i in range(len(files_F150)):
        fit_run_F150_list.append(pickle.load(open(files_F150[i],'rb')))    
    chisqs_f150 = np.array([fit_run_F150_list[i].reduced_Chisq for i in range(len(fit_run_F150_list))])
    sort_Chisq_F150 = chisqs_f150.argsort()  
    fit_run_F150 = fit_run_F150_list[sort_Chisq_F150[0]]
        
    for i in range(len(files_HST)):
        fit_run_HST_list.append(pickle.load(open(files_HST[i],'rb')))    
    chisqs_HST = np.array([fit_run_HST_list[i].reduced_Chisq for i in range(len(fit_run_HST_list))])
    sort_Chisq_HST = chisqs_HST.argsort()  
    fit_run_HST = fit_run_HST_list[sort_Chisq_HST[0]]
    AGN_mag_list.append([fit_run_F150.final_result_ps[0]['magnitude'], fit_run_HST.final_result_ps[0]['magnitude']]) 
    host_mag_list.append([fit_run_F150.final_result_galaxy[0]['magnitude'], fit_run_HST.final_result_galaxy[0]['magnitude']]) 
AGN_mag_list = np.array(AGN_mag_list)
plt.figure(figsize=(11, 11))
plt.scatter(AGN_mag_list[:,0], AGN_mag_list[:,1],
            c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("AGN mag F150W",fontsize=27)
plt.ylabel("AGN mag "+ HST_filt,fontsize=27)
plt.tick_params(labelsize=20)
plt.plot(np.linspace(20,35), np.linspace(20,35))
plt.xlim([20,32])
plt.ylim([20,32])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()

host_mag_list = np.array(host_mag_list)
plt.figure(figsize=(11, 11))
plt.scatter(host_mag_list[:,0], host_mag_list[:,1],
            c='green',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("host mag F150W",fontsize=27)
plt.ylabel("host mag "+HST_filt,fontsize=27)
plt.tick_params(labelsize=20)
plt.plot(np.linspace(20,35), np.linspace(20,35))
plt.xlim([20,32])
plt.ylim([20,32])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()

#%%
from galight.data_process import DataProcess
print("After remove candidates")
import glob
from galight.tools.cutout_tools import cutout
from galight.tools.measure_tools import measure_bkg
wht_max = {'F105W':11948740000, 'F125W':45840350000, 'F140W':858202000, 'F160W':45840350000.0}
zp_list = {'F606W':26.489, 'F814W':25.937, 'F105W':26.264, 
            'F125W':26.232, 'F140W':26.450, 'F160W':25.936,
            'F150W':29.50477, 'F115W': 29.50477}

#Load JWST
JWST_filt = 'F150W'
HST_filt = 'F140W'
PSF_lib_files = glob.glob('material/*'+JWST_filt[:-1]+'*_PSF_Library.pkl')[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))

#Load HST
HST_folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
HST_all_files= glob.glob(HST_folder+'/egs_all_wfc3_ir_{0}_030mas_v1.9_drz.fits'.format(HST_filt.lower()))  #For HST
HST_fitsFile = pyfits.open(HST_all_files[0])
fov_image_HST = HST_fitsFile[0].data 
header_HST = HST_fitsFile[0].header 

HST_PSF_list = []
for i in range(len(PSF_RA_DEC_list)):
    print(i)
    RA, Dec = PSF_RA_DEC_list[i]
    zp = zp_list[HST_filt]
    data_process = DataProcess(fov_image = fov_image_HST, target_pos = [RA, Dec],
                               pos_type = 'wcs', header = header_HST,
                               rm_bkglight = False, if_plot=False, #exptime= wht, 
                               zp = zp)
    fov_cutout = cutout(image=data_process.fov_image, center= data_process.target_pos, radius=200)
    try:
        bkglight = measure_bkg(fov_cutout, if_plot=False) # Remove bkg light
        # data_process.target_pos = data_process.target_pos + np.array(shift_list[idx][0][-1]) * exppix # Match to JWST
        data_process.generate_target_materials(radius=100,cut_kernel='center_bright', #bkg_std = bkg_std,
                                               create_mask = False, skip = True)
        ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
        data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]
        data_process.target_stamp = np.flip(data_process.target_stamp)
        HST_PSF_list.append(data_process.target_stamp)
    except:
        HST_PSF_list.append([1.e-9])
from galight.tools.measure_tools import flux_profile
F150W_AGN_mag, HST_AGN_mag = [], []
#%%
for i in range(len(PSF_list)):
    if np.mean(HST_PSF_list[i]) != 1.e-9:
        radius = len(PSF_list[i])/2 -2
        seeding_num = np.min([int(radius*2), 100])
        center = [len(PSF_list[i])/2, len(PSF_list[i])/2]
        r_flux, r_grids, _  =  flux_profile(PSF_list[i], center = center, radius = radius,
                                            x_gridspace = 'log', if_plot=False, fits_plot = False, grids=seeding_num )
        r_flux_s = r_flux[:-1]/r_flux[1:]
        r_flux_s = [r_flux_s>0.997]
        t_flux = 0
        for j in range(4, len(r_flux_s)):
            if r_flux_s[j]==True and r_flux_s[j-1]==True and r_flux_s[j-2]==True:
                t_flux = r_flux[j]
        if t_flux == 0:
            t_flux = r_flux[-1]
        F150W_AGN_mag.append(-2.5*np.log10(np.sum(t_flux)) + zp_list['F150W'])
        
        radius = len(HST_PSF_list[i])/2 -2
        seeding_num = np.min([int(radius*2), 100])
        center = [len(HST_PSF_list[i])/2, len(HST_PSF_list[i])/2]
        r_flux, r_grids, _  =  flux_profile(HST_PSF_list[i], center = center, radius = radius,
                                            x_gridspace = 'log', if_plot=False, fits_plot = False, grids=seeding_num )
        r_flux_s = r_flux[:-1]/r_flux[1:]
        r_flux_s = [r_flux_s>0.997]
        t_flux_hst = 0
        for j in range(4,len(r_flux_s)):
            if r_flux_s[j]==True and r_flux_s[j-1]==True and r_flux_s[j-2]==True:
                t_flux_hst = r_flux[j]
        if t_flux_hst == 0:
            t_flux_hst = r_flux[-1]
        HST_AGN_mag.append(-2.5*np.log10(np.sum(t_flux_hst)) + zp_list[HST_filt])

F150W_AGN_mag = np.array(F150W_AGN_mag)
HST_AGN_mag = np.array(HST_AGN_mag)

bools = (HST_AGN_mag<23) * abs(F150W_AGN_mag - HST_AGN_mag)<0.5
    
plt.figure(figsize=(11, 11))
plt.scatter(F150W_AGN_mag[bools], HST_AGN_mag[bools],
            c='darkred',s=280,marker=".",zorder=0, vmin=1.2, vmax=1.8, edgecolors='white',alpha=0.7)
#plt.title('', fontsize=27)
plt.xlabel("stars mag "+JWST_filt,fontsize=27)
plt.ylabel("stars mag "+HST_filt,fontsize=27)
plt.tick_params(labelsize=20)
plt.plot(np.linspace(15,35), np.linspace(15,35))
plt.xlim([18,23])
plt.ylim([18,23])
#plt.legend(scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.show()
np.mean( (F150W_AGN_mag[bools] - HST_AGN_mag[bools]))
    
