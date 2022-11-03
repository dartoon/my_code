#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 23:13:44 2022

@author: Dartoon

Similar to JWST's 1_generate_data_process.py for HST IR
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle
from galight.data_process import DataProcess
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.measure_tools import measure_bkg
from galight.tools.cutout_tools import cutout
import warnings
warnings.filterwarnings("ignore")

# shift_list = pickle.load(open('../model_JWST_stage3_Masafusa/material/jwst_shift_list.pkl','rb')) #For fov shift
# ignore_id = [10, 21, 30, 31,  41, 46, 47, 52]
# remove_id = [24, 55]

# f = open("../model_JWST_stage3_Masafusa/material/target_info.txt","r")
# string = f.read()
# lines = string.split('\n')   # Split in to \n
# result_folder = 'fit_result/'
# lines = lines[1:]


HST_folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_HST_data/'
# for idx in range(len(lines)):
wht_max = {'F105W':11948740000, 'F125W':45840350000, 'F140W':858202000, 'F160W':45840350000.0}

# for idx in range(48, len(lines)):
# for idx in range(10, 11):
    # if idx in remove_id:
    #     continue
    # line = lines[idx]
HST_all_files= glob.glob(HST_folder+'/egs_all_wfc3_ir_*_030mas_v1.9_drz.fits')  #For NIRCam
# target_id, _RA, _Dec, spec_z, photo_z = line.split(' ')
target_id, RA, Dec = 'SDSS', 214.82340490630813, 52.830309411297804
idx = 101
# i = -1
# cut_kernel = 'nearest_obj_center'
data_process_list = []
com_aper_l = []
for i in range(len(HST_all_files)):
    cut_kernel = None
    HST_fitsFile = pyfits.open(HST_all_files[i])
    print("Loading...,", 'idx', idx, HST_all_files[i].split('/')[-1])
    fov_image_HST = HST_fitsFile[0].data 
    header_HST = HST_fitsFile[0].header 
    data_process = DataProcess(fov_image = fov_image_HST, target_pos = [RA, Dec],
                               pos_type = 'wcs', header = header_HST,
                               rm_bkglight = False, if_plot=False, #exptime= wht, 
                               zp = 27)
    # del data_process.fov_noise_map
    wht = pyfits.open(HST_all_files[i].replace('drz','wht'))[0].data
    exp = header_HST['EXPTIME']
    print(header_HST['DATE-OBS'])
    exp_map = exp * wht[int(data_process.target_pos[1]), int(data_process.target_pos[0])] / wht_max[header_HST['filter']]
    print(exp_map)
    data_process.exptime = exp_map
    
    #estimate local bkg and remove:
    if fov_image_HST[int(data_process.target_pos[1]), int(data_process.target_pos[0])] == 0 :
        continue
    # data_process.generate_target_materials(radius=radius, create_mask = False, 
    #                                        cut_kernel = cut_kernel,
    #                                        npixels = 80)
    
    # data_process_file = pickle.load(open('fit_material/data_process_idx{0}_{1}_psf0.pkl'.format(idx, header_HST['filter']),'rb')) #For fov shift
    # if com_aper_l != []:
    #     shift_x = data_process.tbl[data_process.tbl['label']==0]['xcentroid']  - com_aper_l[0].positions[0]
    #     shift_y = data_process.tbl[data_process.tbl['label']==0]['ycentroid']  - com_aper_l[0].positions[1]
    # shift_x = 0
    # shift_y = 5
    # bkg_std = data_process_file.bkg_std
    # data_process.target_pos = data_process_file.target_pos
    
    fov_cutout = cutout(image=data_process.fov_image, center= data_process.target_pos, radius=200)
    bkglight = measure_bkg(fov_cutout, if_plot=False) # Remove bkg light
    
    # data_process.target_pos = data_process.target_pos + np.array(shift_list[idx][0][-1]) * exppix # Match to JWST
    data_process.generate_target_materials(radius=40,cut_kernel=None, #bkg_std = bkg_std,
                                           create_mask = False, skip = True)
    ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
    data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]
    
    #flip image to match with JWST
    data_process.target_stamp = np.flip(data_process.target_stamp)
    data_process.noise_map = np.flip(data_process.noise_map)
    
    del data_process.fov_image
    del data_process.exptime
    # del data_process.exp_map
    data_process.filt = header_HST['filter']
    if np.sum(data_process.target_stamp) != 0:
        data_process_list.append(data_process)
            
    filters =[data_process_list[i].filt for i in range(len(data_process_list))]
    if com_aper_l == []:
        com_aper_l = common_data_class_aperture(data_process_list, l_idx=0, return_idx=0)
print("Common aperture:")
for i in range(len(data_process_list)):
    print(idx, target_id, filters[i],':')
    plot_data_apertures_point(data_process_list[i].target_stamp * data_process_list[i].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                              com_aper_l, figsize=(4,3))
print("Above are for the", 'idx:', idx, target_id, 'filts:', filters)
# hold = input('Hold ... OK?\n')
pickle.dump([[data_process_list, com_aper_l]], open('material/'+'data_process+apertures_{0}_IR.pkl'.format(idx), 'wb'))

#%%
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
zp_list = {'F606W':26.489, 'F814W':25.937, 'F105W':26.264, 
            'F125W':26.232, 'F140W':26.450, 'F160W':25.936}

# print("work on", count, 'idx', idx, filt, "Total PSF NO.", len(idx_counts))
fit_run_all = []
aper = com_aper_l[0]
for i, filt in enumerate(filters):
    
    #Take the PSF from the otherband fitting:
    fit_files = glob.glob('fit_material/fit_run_idx48_{0}*.pkl'.format(filt))
    fit_run_list = []
    for k in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[k],'rb')))
    chisqs = np.array([fit_run_list[k].reduced_Chisq for k in range(len(fit_run_list))])
    sort_Chisq = chisqs.argsort()  
    # idx_counts = chisqs.argsort()  
    _fit_run = fit_run_list[sort_Chisq[0]]
    psf = _fit_run.fitting_specify_class.data_process_class.PSF_list[0]
    
    #Use PSF to run:
    data_process = data_process_list[i]
    data_process.zp = zp_list[filt]
    aper.positions = np.array([32, 50])
    aper.a, aper.b = 3, 3
    data_process.apertures = []
    data_process.PSF_list = [psf]
    psf[psf<0] = 0
    fit_sepc = FittingSpecify(data_process)
    target_stamp = data_process.target_stamp
    ps_pos = np.where(target_stamp == np.max(target_stamp))
    ps_pos = (ps_pos[0][0] - data_process.radius, ps_pos[1][0] - data_process.radius)
    ps_pos = [ps_pos[1], ps_pos[0]]    
    fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 1,
                                  ps_pix_center_list = [ ps_pos ]  ) #, fix_n_list= [[0,4],[1,1]])
    # fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
    fit_sepc.build_fitting_seq()
    print(filt)
    fit_sepc.plot_fitting_sets()
    fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep','deep'])
    fit_run.run(algorithm_list = ['PSO','PSO','PSO'])
    fit_run.plot_final_qso_fit(target_ID =target_id)
    fit_run_all.append(fit_run)
    # filt = data_process.filt
    # pickle.dump(fit_run , open('fit_material/'+'fit_run_idx{2}_{0}_psf{1}.pkl'.format(filt, i, idx), 'wb'))
pickle.dump([filters, fit_run_all] , open('run_z5_target.pkl', 'wb'))

#%%
from galight.tools.measure_tools import esti_bgkstd
from galight.tools.measure_tools import detect_obj, mask_obj
run_z5_target = pickle.load(open('run_z5_target.pkl','rb'))
filters, fit_run_all = run_z5_target
from galight.tools.measure_tools import measure_FWHM
for i in range(len(filters)):
    data_process = fit_run_all[i].fitting_specify_class.data_process_class
    apertures, segm_deblend, mask_apertures, tbl  = detect_obj(data_process.target_stamp)
    flux = tbl['kron_flux'][0]
    print(-2.5*np.log10(flux) + zp_list[filters[i]])
    # data_process = fit_run_all
    # print(filters[i], fit_run_all[i].final_result_ps[0]['magnitude'])
    # print(filters[i], np.mean(measure_FWHM(fit_run_all[i].fitting_specify_class.data_process_class.PSF_list[0])))
    # print(filters[i], fit_run_all[i].final_result_galaxy[0]['magnitude'])
    


    # import shutil
    # shutil.copyfile('1_generate_data_process.py', 'backup_files/1_generate_data_process_idx{0}.py'.format(idx))