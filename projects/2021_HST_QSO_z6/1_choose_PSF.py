#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 10:20:59 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import galight.tools.astro_tools as astro_tools
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
    
f = open("fit_files_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
info = [lines[i].split(' ') for i in range(len(lines))]

ra_dec_file = 'HST_download_sh/ALMAz6qso_radec_for_HST_archive.list'
f = open(ra_dec_file,"r")
string = f.read()
ra_dec_info = string.split('\n')   # Split in to \n

# i, psf_i, psf_pos
# [0, x, [793., 367.]]
# [1, x, [1001.,  218.]]
# [2, x, [1001.,  218.]]  #ID 1 and 2 are the same sources
# [3, x, [475., 185.]]
# [4, x, [641., 199.]]
# [5, x, [465., 513.]]
# [6, x, [156., 190.]]
# [7, x, [502., 222.]]
# [8, 34, [486., 778.]]
# [9, 3, [996., 140.]]
# [10, x, [532., 540.]]
# [11, x, [524., 442.]]
# [12, x, [657., 514.]]
# [13, x, [310., 855.]]
# [14, x, [751., 690.]]
# [15, x, [550., 146.]]
# [16, x, [310., 855.]]

# [17, x, [xxx, xxx]]
#%%
PSF_loc_dic = {'0': [793., 367.],
'1': [1001.,  218.],
'2': [1001.,  218.],
'3': [475., 185.],
'4': [641., 199.],
'5': [465., 513.],
'6': [156., 190.],
'7': [502., 222.],
'8': [486., 778.],
'9': [996., 140.],
'10': [532., 540.],
'11': [524., 442.],
'12': [657., 514.],
'13': [310., 855.],
'14': [751., 690.],
'15': [550., 146.],
'16': [310., 855.]}

for i in range(len(info)):
# for i in [10]:
    ID, filename, _ = info[i]
    if ID == 'J0842+1218':
        continue
    
    fitsFile = pyfits.open(filename)
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    
    wht = fitsFile[2].data # The WHT map
    exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
    pixel_scale = astro_tools.read_pixel_scale(fitsFile[1].header)  #Read pixel scale
    mean_wht = exp * (pixel_scale/0.135)**2
    exp_map = exp * wht/mean_wht
    # print(pixel_scale, exp, exp_map.max())
        
    idx = [i for i in range(len(ra_dec_info)) if ID in ra_dec_info[i]][0]
    # print(ID, idx)
    RA, Dec = ra_dec_info[idx].split(' ')[:2]
    RA, Dec = np.float(RA), np.float(Dec)
    if ID == 'J2100-1715':
        RA = 315.2279713
        Dec = -17.25608967

    fov_image = fitsFile[1].data # check the back grounp
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                          rm_bkglight = True, exptime = exp_map, if_plot=False, zp = 25.9463)  #!!! zp use F160W for now
    data_process.generate_target_materials(radius=20, cut_kernel = 'center_bright', create_mask = False, 
                                           # detect_tool = 'sep', if_select_obj= True, nsigma=2.5, thresh = 2.5, exp_sz= 1.2, npixels = 15, 
                                           if_plot=False)
    
    #Lines used to find PSF and set the PSF_loc_dic.
    data_process.find_PSF(radius = 30, user_option = True, if_filter=True, psf_edge =30)
    # data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
    #               prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )
    # PSF_loc = PSF_loc_dic[str(i)]
    data_process.plot_overview(label = ID, target_label = None)
    # data_process.find_PSF(radius = 30, PSF_pos_list = [PSF_loc])
    data_process.checkout()
    #Start to produce the class and params for lens fitting.
    # data_process.apertures = []
    print(ID, i)
    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
    # psf_error_map = np.ones_like(data_process.PSF_list[data_process.psf_id_for_fitting]) *0.01 # It is in the variance unit (std^2).
    # fit_sepc.prepare_fitting_seq(point_source_num = 1, psf_error_map = psf_error_map)
    fit_sepc.build_fitting_seq()
    #Plot the initial settings for fittings. 
    fit_sepc.plot_fitting_sets()

    # fit_run = FittingProcess(fit_sepc, savename = 'pkl_files/'+ ID+'_'+str(i), fitting_level='deep')
    # fit_run.run(algorithm_list = ['PSO'], setting_list=[None])
    #             # setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}, {'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
    # # fit_run.plot_all()
    # fit_run.plot_final_qso_fit()
    # fit_run.dump_result()
    # print(fit_run.final_result_galaxy[0])

#%% Co-add residual
import glob
import pickle
from galight.tools.astro_tools import plt_fits
tot_residual = None
for i in range(len(info)-1):
    ID, filename, _ = info[i]
    picklename = 'first_pkl_files/'+ID+'_'+str(i)+'.pkl'
    picklename = glob.glob(picklename)
    if picklename  != []:
        fitting_run_class = pickle.load(open(picklename[0],'rb'))
        data = fitting_run_class.fitting_specify_class.kwargs_data['image_data']
        psf = fitting_run_class.image_ps_list[0]
        nearby = fitting_run_class.image_host_list[1:]
        host_resi = data - psf
        # if i == 10:
        #     fitting_run_class.plot_final_qso_fit()
        for j in range(len(nearby)):
            host_resi = host_resi - nearby[j]
        # plt_fits(host_resi)
        if tot_residual is None:
            tot_residual = host_resi
        else:
            tot_residual += host_resi

plt_fits(tot_residual, colorbar = True)
