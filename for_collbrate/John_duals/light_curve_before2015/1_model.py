#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 20:45:19 2022

@author: Dartoon
"""

import sys
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.astro_tools import plt_fits
import glob
#RA, DEC information of the QSO:
# QSO_RA, QSO_DEC = 35.2735048,  -4.6837584
#Calculate the zeropoint for HSC filters:
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle
from astropy.wcs import WCS
from _0_scp_data import data_list
# folders = glob.glob("HSC_data/*fits")
# folders.sort()
i = 0

band ='i'
# name = 'PSPS+sersic_fixpos_nearbyPSF_result'
# for file in folders: 
QSO_RA, QSO_DEC = 35.2735048,  -4.6837584
for i in range(len(data_list)):
# for i in [29]:
    if i !=4 and i != 16 and i != 22 and i!=30:
        b ,c = data_list[i][2:4]
        file = 'HSC_data/CORR-*{0}-*{1}.fits'.format(b, c)
        name = 'HSC_{0}_{1}'.format(b, c)
        file = glob.glob(file)[0]
        fitsFile = pyfits.open(file)
        file_header0 = fitsFile[0].header
        FLUXMAG0 = file_header0['FLUXMAG0']
        zp = 2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
        print(file, zp, data_list[i][4])
        #Load the fov image data:
        fov_image = fitsFile[1].data # check the back grounp
        #Derive the header informaion, might be used to obtain the pixel scale and the exposure time.
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        #Derive the fov noise level map:
        err_data= fitsFile[3].data ** 0.5
        #Load the PSF data:
        # PSF_files = glob.glob('psf*{0}*.fits'.format(HSCband))
        # PSF = pyfits.getdata(PSF_files[0])
        QSO_x, QSO_y = float(data_list[i][-2]), float(data_list[i][-1])
        wcs = WCS(header)
        cal_QSO_x, cal_QSO_y = wcs.all_world2pix([[QSO_RA, QSO_DEC]], 1)[0]
        data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_x, QSO_y],
                                    pos_type = 'pixel', header = header,
                                  rm_bkglight = False, if_plot=False, zp = zp)
    
        # PSF_loc = [35.29689501,-4.685896878]
        PSF_loc = [35.2807986,-4.682323596]
        cal_psf_x, cal_psf_y = wcs.all_world2pix([[PSF_loc[0], PSF_loc[1]]], 1)[0]
        if i ==6:
            cal_psf_x, cal_psf_y = [80, 314]
        # PSF_loc = [35.28084894, -4.6823214]  #The one nearby
        # try:
        cal_psf_x = cal_psf_x + (QSO_x - cal_QSO_x)
        cal_psf_y = cal_psf_y + (QSO_y - cal_QSO_y)
        data_process.find_PSF(radius = 15, PSF_pos_list = [[cal_psf_x, cal_psf_y]], pos_type = 'pixel',)
    
    #     #Generate the fitting materials
        data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8,
                                              exp_sz= 1.5, npixels = 40, if_plot=False)
        # print(name)
        plt_fits(data_process.PSF_list[0])
        # except:
        #     None
        # data_process.apertures = [] #Assuming there is no host (i.e., as constant.) #!!!
    
        # Manually input the PSF:
        # data_process.PSF_list = [PSF]
    # #%%
        # # Check if all the materials is given, if so to pass to the next step.
        data_process.checkout() #Check if all the materials is known.
        fit_sepc = FittingSpecify(data_process)
        # fit_sepc.prepare_fitting_seq(point_source_num = 2)
        fit_sepc.prepare_fitting_seq(point_source_num = 2)
        # if fit_sepc.kwargs_params['point_source_model'][0][0] == fit_sepc.kwargs_params['point_source_model'][0][1]:
        #     fit_sepc.prepare_fitting_seq(point_source_num = 2, ps_pix_center_list= [[-2.0, -1.0], [6.0, -1.0]])
        #     print(file)
        fit_sepc.build_fitting_seq()
        fit_sepc.plot_fitting_sets()
        # # if rerun != 's19a':
        # #     fit_sepc.kwargs_params['point_source_model'] = s21_res.fitting_specify_class.kwargs_params['point_source_model']
        fit_run = FittingProcess(fit_sepc, savename = 'results/'+name, fitting_level='deep')
        fit_run.run(algorithm_list = ['PSO','PSO'], setting_list=[None,None])
        fit_run.plot_final_qso_fit(save_plot=True, target_ID= name +'-'+ band)
        fit_run.cal_astrometry()
        fit_run.dump_result()