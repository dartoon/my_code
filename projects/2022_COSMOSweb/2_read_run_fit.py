#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:32:37 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle,glob
from galight.tools.astro_tools import plt_fits
# import sys
# count = int(sys.argv[1]) - 1 # 1 - 12586
# idx = 10

filt_i = 0
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
# cata_list = pickle.load(open('material/cata_list.pkl','rb'))
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

# for i, item in enumerate(cata_list):
#     print(i, item)

#%%
# Quick check:
# fit_run.run(algorithm_list = ['PSO'], fitting_level=['norm'])
# fit_run.plot_final_qso_fit(target_ID =None, show_plot=True)
# fit_run.cal_astrometry()
# savename = fit_files[i].replace('_notrunyet_', '_run_')[:-4]+'_{0}.pkl'.format(i)
# pickle.dump(fit_run , open(savename, 'wb'))
top_psf_id = 0
# top_psf_id = 3 For idx 32 F444
# top_psf_id = 2 For idx 40 F444

# fit_file_folder ='/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
fit_file_folder ='./'
save_plot = False
save_psf = False

# check_name= 'cid_473'  #29
# check_name= 'cid_1210' #8
# check_name= 'cid_1245' #10

# for idx in range(50):
# [i for i in range(len(cata_list)) if cata_list[i][-1] == 'cid_473']
for idx in [29]:
    fit_run_list = []
    fit_files = glob.glob(fit_file_folder+'fit_result/fit2_run_{0}*idx{1}_psf*.pkl'.format(filt,idx))
    fit_files.sort()
    z = cata_list[idx][6]
    if z >0:
        zinfo = 'Zspec'+str(z)
    elif z <0:
        zinfo = 'Zphot'+str(cata_list[idx][5])
    
    if fit_files != []:
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        fit_run.cal_astrometry()
        target_ID = 'idx{0}_{1}_{2}_{3}'.format(idx,cata_list[idx][-1],filt,zinfo)
        
        savename = 'figures/idx{0}_{1}_{2}_{3}'.format(idx,cata_list[idx][-1],filt,zinfo)
        fit_run.fitting_specify_class.plot_fitting_sets(savename=savename+'_fitting_config.pdf')
        fit_run.savename = savename
        fit_run.plot_final_qso_fit(target_ID =target_ID, show_plot=True,save_plot=save_plot)
    if save_psf == True:
        pyfits.PrimaryHDU(fit_run.fitting_specify_class.data_process_class.PSF_list[0]).writeto('PSFs_library/PSF_{0}_id{1}.fits'.format(
            filt,sort_Chisq[top_psf_id]),overwrite=True)
    print(fit_run.reduced_Chisq)
    # plt_fits(fit_run.flux_2d_out['data'] - fit_run.flux_2d_out['data-point source'])
    
    
    # count_n = 5
    # Chisq_best = chisqs[sort_Chisq[top_psf_id]]
    # Chisq_last= chisqs[sort_Chisq[count_n-1]]
    # inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
    # weight = np.zeros(len(chisqs))
    # for i in sort_Chisq[:count_n]:
    #     weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
    # prop_name = 'magnitude'
    # # all_values = [fit_run_list[i].final_result_ps[0][prop_name] for i in range(len(fit_run_list))]
    # all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
    # weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
    # rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
    # print(cata_list[idx][-1], filt,': ' , round(weighted_value,2), '+-', round(rms_value,2))
    # print(prop_name, round(weighted_value,2), '+-', round(rms_value,2))
    
        
    # if fit_files == []:
    #     pos = cata_list[idx][3:5]
    #     filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
    #     # filename = '1727_cosmos_mosaic_miri_exptime_scale1.0.fits'
    #     fitsFile = pyfits.open(fit_file_folder+filename)
    #     img = fitsFile[1].data[int(pos[1])-100:int(pos[1])+100, int(pos[0])-100:int(pos[0])+100] #
    #     if save_plot == False:
    #         plt_fits(img)
    #     elif save_plot == True:
    #         plt_fits(img, savename='figures/idx{0}_{1}_{2}_{3}_notfit.pdf'.format(idx,cata_list[idx][-1],filt,zinfo))

#%% Read PSF FWHM
use_psf = fit_run.fitting_specify_class.data_process_class.PSF_list[0]
from galight.tools.measure_tools import measure_FWHM
FWHM = np.mean(measure_FWHM(use_psf))
delta_pixel = fit_run.fitting_specify_class.deltaPix
print(filt, round(FWHM*delta_pixel,3), 'arcsec')

#%%

# #%% Refit image but not adding point source
# from galight.fitting_specify import FittingSpecify
# from galight.fitting_process import FittingProcess
# fit_sepc = FittingSpecify(fit_run.fitting_specify_class.data_process_class)
# import copy
# # ps_pix_center_list = [copy.deepcopy(ps_pos)]
# # if idx == 15:
#     # ps_pix_center_list = None
# fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor = 3, apertures_center_focus=True )
# # fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
# # fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
# fit_sepc.plot_fitting_sets()
# # print(idx, filt, 'apertures', len(data_process.apertures) )
# fit_run = FittingProcess(fit_sepc, savename = cata_list[idx][-1])
# fit_run.run(algorithm_list = ['PSO','PSO','PSO'], fitting_level=['norm', 'norm','norm'])

# check_name = ''
# if idx == 8:
#     check_name = 'cid_1210'
# elif idx == 10:
#     check_name = 'cid_1245'

# fit_run.savename =  '/Users/Dartoon/Downloads/no_ps_fit_'+ check_name
# fit_run.plot_final_galaxy_fit(target_ID=' ', save_plot=True )
# print(fit_run.reduced_Chisq)

