#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:50:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source')
from source_info import  source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise
import copy
import pickle
for source_id in range(1,25):
# source_id = 4
    source_name=source_list[0][source_id]
    Re = source_list[2][source_id]
    print(source_name, source_id)
    hdu = pyfits.open('/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source/fix/{0}_fix.fits'.format(source_name))
    hd_gal_img = hdu[0].data
    hdu.close()
    print("source name:", source_name)
    # from galight.tools.astro_tools import plt_fits
    # plt_fits(hd_gal_img)
    
    # from skimage.transform import resize
    # new_image_0 = resize(hd_gal_img, [25,25])
    # plt_fits(new_image_0)
    
    # from scipy.ndimage import zoom
    # new_image_1 = zoom(hd_gal_img, 0.1)
    # plt_fits(new_image_1)
    
    # project_gal_img = zoom(hd_gal_img, 100/len(hd_gal_img))
    # from galight.data_process import DataProcess
    # data_process = DataProcess(fov_image = project_gal_img, target_pos = [len(project_gal_img)/2]*2, pos_type = 'pixel',
    #                           rm_bkglight = False, exptime = 2400, if_plot=False, zp = 27.0)
    # data_process.generate_target_materials(radius=len(project_gal_img)/2, create_mask = False, nsigma=2.8, if_select_obj=False,
    #                                       exp_sz= 1.2, npixels = 30, if_plot=True, bkg_std=0.1)
    # data_process.apertures = [data_process.apertures[0]] 
    # PSF = copy.deepcopy(data_process.target_stamp[:11,:11]*0)
    # PSF[5,5]=1
    # data_process.PSF_list = [PSF]
    # data_process.deltaPix = 1
    # # data_process.noise_map[np.isnan(data_process.noise_map)] = bkg_std
    
    # from galight.fitting_specify import FittingSpecify
    # fit_sepc = FittingSpecify(data_process)
    # fit_sepc.prepare_fitting_seq(point_source_num = 0) #, fix_n_list= [[0,4],[1,1]])
    # fit_sepc.build_fitting_seq()
    # fit_sepc.plot_fitting_sets()  #The apertures shows how the images will be modelled.
    # from galight.fitting_process import FittingProcess
    # fit_run = FittingProcess(fit_sepc, savename = 'galight_fit_HD_galaxy/'+ source_name, fitting_level='norm')
    # fit_run.run(algorithm_list = ['PSO'], setting_list=[None])
    # fit_run.plot_final_galaxy_fit()
    # fit_run.dump_result()
    fit_run = pickle.load(open('galight_fit_HD_galaxy/'+ source_name+'.pkl','rb'))
    # print(round(fit_run.final_result_galaxy[0]['R_sersic']/100*len(hd_gal_img),2), Re)
    print(fit_run.final_result_galaxy[0]['n_sersic'])