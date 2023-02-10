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

for i, item in enumerate(cata_list):
    print(i, item)

#%%
# Quick check:
# fit_run.run(algorithm_list = ['PSO'], fitting_level=['norm'])
# fit_run.plot_final_qso_fit(target_ID =None, show_plot=True)
# fit_run.cal_astrometry()
# savename = fit_files[i].replace('_notrunyet_', '_run_')[:-4]+'_{0}.pkl'.format(i)
# pickle.dump(fit_run , open(savename, 'wb'))
top_psf_id = 0
fit_file_folder ='/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
for idx in range(50):
    fit_run_list = []
    fit_files = glob.glob(fit_file_folder+'fit_result/fit_run_{0}*idx{1}_psf*.pkl'.format(filt,idx))
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
        fit_run.savename = 'figures/idx{0}_{1}_{2}_{3}'.format(idx,cata_list[idx][-1],filt,zinfo)
        fit_run.plot_final_qso_fit(target_ID =target_ID, show_plot=True,save_plot=True)
    if fit_files == []:
        pos = cata_list[idx][3:5]
        filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
        # filename = '1727_cosmos_mosaic_miri_exptime_scale1.0.fits'
        fitsFile = pyfits.open(fit_file_folder+filename)
        img = fitsFile[1].data[int(pos[1])-100:int(pos[1])+100, int(pos[0])-100:int(pos[0])+100] #
        plt_fits(img, savename='figures/idx{0}_{1}_{2}_{3}_notfit.pdf'.format(idx,cata_list[idx][-1],filt,zinfo))