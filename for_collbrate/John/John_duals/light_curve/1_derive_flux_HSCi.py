#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:45:59 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from galight.tools.astro_tools import plt_fits
import glob
#RA, DEC information of the QSO:
QSO_RA, QSO_DEC = 35.2735048,  -4.6837584
#Calculate the zeropoint for HSC filters:
# file_header0 = fitsFile[0].header
# FLUXMAG0 = file_header0['FLUXMAG0']
zp = 27 # 2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle
folders = glob.glob("20*-*")
folders.sort()

name = 'results/PSPS+sersic_fixpos_nearbyPSF_result'
# folders = ['2021-10-31']
# for folder in folders:
# for folder in folders[:20]:
# for folder in folders[20:35]:
for folder in folders:    
    files = glob.glob(folder+"*/*/*/*/*/*fits")
    for file in files:
        HSCband = file.split('calexp-')[1][:5]
        band = HSCband[-1]
        # pklfiles = glob.glob(name+'-{0}-band{1}*pkl'.format(folder,band))
        # pklfile =  pklfiles[0]
        # fit_run = pickle.load(open(pklfile,'rb'))
        # if len(fit_run.final_result_galaxy) == 1:
        #     continue
        # else:
        #     print("fitting:", '{0}-band{1}'.format(folder,band))
        fitsFile = pyfits.open(file)
        file_header0 = fitsFile[0].header
        FLUXMAG0 = file_header0['FLUXMAG0']
        zp = 2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
        print(file, zp)
        #Load the fov image data:
        fov_image = fitsFile[1].data # check the back grounp
        #Derive the header informaion, might be used to obtain the pixel scale and the exposure time.
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        #Derive the fov noise level map:
        err_data= fitsFile[3].data ** 0.5
        #Load the PSF data:
        # PSF_files = glob.glob('psf*{0}*.fits'.format(HSCband))
        # PSF = pyfits.getdata(PSF_files[0])
        data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                                    pos_type = 'wcs', header = header,
                                  rm_bkglight = True, if_plot=False, zp = zp)
        
        PSF_loc = [35.2807986,-4.682323596]
        data_process.find_PSF(radius = 15, PSF_pos_list = [PSF_loc], pos_type = 'wcs',)
        
        #Generate the fitting materials
        data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8,
                                              exp_sz= 1.5, npixels = 40, if_plot=False)

        plt_fits(data_process.PSF_list[0])
        # data_process.apertures = [] #Assuming there is no host (i.e., as constant.) #!!!

        # Manually input the PSF:
        # data_process.PSF_list = [PSF]
        
        # Check if all the materials is given, if so to pass to the next step.
        data_process.checkout() #Check if all the materials is known.
        fit_sepc = FittingSpecify(data_process)
        # fit_sepc.prepare_fitting_seq(point_source_num = 2)
        fit_sepc.prepare_fitting_seq(point_source_num = 2, ps_pix_center_list= [[-2.0, -1.0], [6.0, -1.0]])
        # if fit_sepc.kwargs_params['point_source_model'][0][0] == fit_sepc.kwargs_params['point_source_model'][0][1]:
        #     fit_sepc.prepare_fitting_seq(point_source_num = 2, ps_pix_center_list= [[-2.0, -1.0], [6.0, -1.0]])
        #     print(file)
        fit_sepc.build_fitting_seq()
        fit_sepc.plot_fitting_sets()
        # if rerun != 's19a':
        #     fit_sepc.kwargs_params['point_source_model'] = s21_res.fitting_specify_class.kwargs_params['point_source_model']
        fit_run = FittingProcess(fit_sepc, savename = name+'-{0}-band{1}'.format(folder,band), fitting_level='deep')
        fit_run.run(algorithm_list = ['PSO','PSO'], setting_list=[None,None])
        fit_run.plot_final_qso_fit(save_plot=True, target_ID= folder +'-'+ band)
        fit_run.cal_astrometry()
        fit_run.dump_result()
