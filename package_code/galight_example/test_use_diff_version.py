#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:01:17 2020

@author: Xuheng Ding

You can skip this step if the QSO stamp, noise level and the PSF is ready.
"""
#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
savepickle = input("savepickle?:\n")
if int(savepickle) == 1:
    fitsFile = pyfits.open('../example_data/HST/QSO/1147_final_drz.fits')
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    
    #%%
    wht = fitsFile[2].data # The WHT map
    import galight.tools.astro_tools as astro_tools
    exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
    mean_wht = exp * (0.0642/0.135)**2
    exp_map = exp * wht/mean_wht
    
    #%%Start to use galight
    from galight.data_process import DataProcess
    data_process = DataProcess(fov_image = fov_image, target_pos = [1142., 637.], pos_type = 'pixel', header = header,
                              rm_bkglight = False, exptime = exp_map, if_plot=False, zp = 27.0)
    data_process.generate_target_materials(radius=45, create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 15, if_plot=True)
    
    
    #%%PSF works.
    data_process.find_PSF(radius = 30, user_option = False, if_filter=True, nearyby_obj_filter=True)
    # data_process.find_PSF(radius = 50, PSF_pos_list = [[ 350., 1055.], [2078., 1910.]], user_option = False)
    data_process.plot_overview(label = 'Example', target_label = None)
    
    
    data_process.checkout() #Check if all the materials is known.
    
    #%%Start to produce the class and params for lens fitting.
    from galight.fitting_specify import FittingSpecify
    
    # data_process.apertures = []
    
    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor = 2, point_source_supersampling_factor = 2) #, fix_n_list= [[0,4],[1,1]])
    # psf_error_map = np.ones_like(data_process.PSF_list[data_process.psf_id_for_fitting]) *0.01 # It is in the variance unit (std^2).
    # fit_sepc.prepare_fitting_seq(point_source_num = 1, psf_error_map = psf_error_map)
    
    fit_sepc.build_fitting_seq()
    
    #Plot the initial settings for fittings. 
    fit_sepc.plot_fitting_sets()
    
    #%%Setting the fitting method and run.
    from galight.fitting_process import FittingProcess
    fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level='norm')
    fit_run.run(algorithm_list = ['PSO'])
                # setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}, {'n_burn': 100, 'n_run': 100, 'walkerRatio': 10,'sigma_scale': .1}])
    # fit_run.plot_all()
    fit_run.plot_final_qso_fit()
    
    # fit_run.dump_result()
    fit_run.dump_result_version_free()
    # # print(fit_run.final_result_galaxy[0])

#%%
# Test load pkl
import pickle
picklename = 'savename_version_free.pkl'
fitting_run_class = pickle.load(open(picklename,'rb'))
# fitting_run_class.plot_final_qso_fit()
# fitting_run_class.run_diag()

