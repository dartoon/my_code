#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:10:46 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle,glob
from galight.tools.astro_tools import plt_fits
import warnings
warnings.filterwarnings("ignore")
# import sys
# count = int(sys.argv[1]) - 1 # 1 - 12586
# idx = 10

# filt_i = 3
# filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
# cata_list = pickle.load(open('material/cata_list.pkl','rb'))
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

filename = 'figures/table_summay.txt'
if_file = glob.glob(filename)
if if_file == []:
    write_file =  open(filename,'w')
    write_file.write("idx name redshift\n")
    write_file.write("(0) filter (1) host magnitude (2) host Reff(arcsec) (3) host sersic_n (4) host q (5) QSO magnitude (6) poistional offset (arcsec) \n")
    write_file.close()

# for i, item in enumerate(cata_list):
#     print(i, item)

#%%
import shutil
# Quick check:
# fit_run.run(algorithm_list = ['PSO'], fitting_level=['norm'])
# fit_run.plot_final_qso_fit(target_ID =None, show_plot=True)
# fit_run.cal_astrometry()
# savename = fit_files[i].replace('_notrunyet_', '_run_')[:-4]+'_{0}.pkl'.format(i)
# pickle.dump(fit_run , open(savename, 'wb'))
top_psf_id = 0   
fit_file_folder ='/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
# fit_file_folder ='./'
save_plot = True
save_psf = False
# for idx in [0]:
for idx in range(50):
    for filt in ['F115W', 'F150W','F277W', 'F444W']:
        z = cata_list[idx][6]
        if z >0:
            zinfo = 'Zspec'+str(z)
        elif z <0:
            zinfo = 'Zphot'+str(cata_list[idx][5])
        fit_run_list = []
        fit_files = glob.glob('fit_result/fit2_run_{0}*idx{1}_psf*.pkl'.format(filt,idx))
        if fit_files == []:
            fit_files = glob.glob(fit_file_folder+'fit_result/fit2_run_{0}*idx{1}_psf*.pkl'.format(filt,idx))
        fit_files.sort()
        write_file =  open(filename,'r+') 
        write_file.read()
        if filt == 'F115W':
            write_file.write("idx{0} {1} {2}\n".format(idx, cata_list[idx][-1], zinfo))
        if fit_files != []:
            for i in range(len(fit_files)):
                fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
            chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
            fit_run = fit_run_list[sort_Chisq[top_psf_id]]
            fit_run.cal_astrometry()
            offset = [np.sqrt(np.sum((np.array(fit_run.final_result_galaxy[i]['position_xy']) - 
                                      np.array(fit_run.final_result_ps[0]['position_xy']))**2)) for i in range(len(fit_run.final_result_galaxy))]
            host_id = np.where(offset==np.min(offset))[0][0]
            write_file.write("{0} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f}\n".format(filt, 
                                                                    fit_run.final_result_galaxy[host_id]['magnitude'],
                                                                    fit_run.final_result_galaxy[host_id]['R_sersic'],
                                                                    fit_run.final_result_galaxy[host_id]['n_sersic'],
                                                                    fit_run.final_result_galaxy[host_id]['q'],
                                                                    fit_run.final_result_ps[0]['magnitude'],
                                                                    offset[host_id]))
            if fit_file_folder in fit_files[sort_Chisq[0]]:
                shutil.copy(fit_files[sort_Chisq[0]],'fit_result/')
                
            target_ID = 'idx{0}_{1}_{2}_{3}'.format(idx,cata_list[idx][-1],filt,zinfo)
            savename = 'figures/idx{0}_{1}_{2}_{3}'.format(idx,cata_list[idx][-1],filt,zinfo)
            if save_plot == True:
                fit_run.fitting_specify_class.plot_fitting_sets(savename=savename+'_fitting_config.pdf')
            fit_run.savename = savename
            fit_run.plot_final_qso_fit(target_ID =target_ID, show_plot=True,save_plot=save_plot)
        write_file.close()
        # if save_psf == True:
        #     pyfits.PrimaryHDU(fit_run.fitting_specify_class.data_process_class.PSF_list[0]).writeto('PSFs_library/PSF_{0}_id{1}.fits'.format(
        #         filt,sort_Chisq[top_psf_id]),overwrite=True)
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