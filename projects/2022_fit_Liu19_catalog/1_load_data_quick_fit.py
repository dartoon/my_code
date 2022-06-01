#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:47:31 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#Select which data you plan to download
dr='dr4'
# rerun='s21a_dud'  #Deep or UltraDeep
rerun='s21a_wide'  #Wide 
#Select how many bands you plan to Download (G R I Z Y)
bands = 'GRIZY'  #Band that will be download
import os
f = open("catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
from galight.hsc_utils import hsc_image, hsc_psf
dr='dr4'
rerun='s21a_wide'  #Wide 
ct = 0
ctt = 0
from galight.data_process import DataProcess

for _, line in enumerate(lines[18:19]):
    ID, Ra, Dec = line.split(' ')
    # if i <10:
    #     hsc_psf.get_psfs(ID,Ra,Dec,'./',dr=dr,rerun=rerun,filters='i')
    for b in bands:
        img_globname = glob.glob('gfarm_data_download/{0}_HSC-{1}.fits'.format(ID,b)) + glob.glob('online_data_download/{0}/*cutout*-{1}-*.fits'.format(ID,b) )
        psf_globname = glob.glob('gfarm_data_download/{0}_HSC-{1}_psf.fits'.format(ID,b)) + glob.glob('online_data_download/{0}/*psf*-{1}-*.fits'.format(ID,b) )
        ctt = ctt +1
        if len(img_globname+psf_globname) == 2:
            ct = ct+1
        fitsFile = pyfits.open(img_globname[0])

        #Load the fov image data:
        fov_image = fitsFile[1].data # check the back grounp
        
        #Derive the header informaion, might be used to obtain the pixel scale and the exposure time.
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        
        #Derive the fov noise level map:
        err_data= fitsFile[3].data ** 0.5
        
        zp =  27.0   # This is something Xuheng can't make sure.
        
        #Load the PSF data:
        PSF = pyfits.getdata(psf_globname[0])
        
    
        #RA, DEC information of the QSO:
        QSO_RA, QSO_DEC = float(Ra), float(Dec)
        data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                                   pos_type = 'wcs', header = header,
                                  rm_bkglight = True, if_plot=False, zp = zp)
        
        #Generate the fitting materials
        data_process.generate_target_materials(radius=None, create_mask = False, nsigma=2.5, #if_select_obj=True,
                                              exp_sz= 1.5, npixels = 25, if_plot=True)        
        #Manually input the PSF:
        data_process.PSF_list = [PSF]
        data_process.clean_aperture()
        # Compare the 1D profile of all the components.
        data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
                          prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )
        
        #Check if all the materials is given, if so to pass to the next step.
        data_process.checkout() #Check if all the materials is known.
        
        # #Start to produce the class and params for lens fitting.
        # #For more details, see notebook galight_HST_QSO.ipynb
        # from galight.fitting_specify import FittingSpecify
        # fit_sepc = FittingSpecify(data_process)
        # fit_sepc.prepare_fitting_seq(point_source_num = 1)
        # #Using following line: want to fix Sersic_n as 4 for the source_id = 0, and if want to fix the QSO and host center:
        # # fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
        
        # #Plot the initial settings for fittings. 
        # fit_sepc.plot_fitting_sets()
        # fit_sepc.build_fitting_seq()
        
        # #Setting the fitting method and run.
        # from galight.fitting_process import FittingProcess
        
        # #Pass fit_sepc to FittingProcess,
        # # savename: The name of the saved files.    
        # fit_run = FittingProcess(fit_sepc, savename = 'Test-'+b, fitting_level='shallow')
        # fit_run.run(algorithm_list = ['PSO', 'PSO']) 
        # #Try also setting_list = [{'sigma_scale': 0.8, 'n_particles': 50, 'n_iterations': 50}, {'n_burn': 50, 'n_run': 100, 'walkerRatio': 10, 'sigma_scale': .1}]
        
        # # Plot all the fitting results, including:
        # #         run_diag() : The convergence of the chains.
        # #         model_plot(): The model plot (by lenstronomy)
        # #         plot_params_corner(): The mcmc corner for all the chains (MCMC should be peformed) 
        # #         plot_flux_corner(): The flux corner for all the component (MCMC should be peformed)
        # #         plot_final_qso_fit() or plot_final_galaxy_fit(): Plot the overall plot (data, model, data-ps, resudal, 1D profile)
        
        # #Save the fitting class as pickle format:
        # #     Note, if you use python3 (or 2), load with python3 (or 2)
        # fit_run.plot_final_qso_fit()
        # # fit_run.dump_result(savedata=False) 
            
        