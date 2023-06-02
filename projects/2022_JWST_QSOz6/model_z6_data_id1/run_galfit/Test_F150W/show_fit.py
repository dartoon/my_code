#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 06:07:08 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

zp = 29.5385
#%%
folder = 'psf2'

PSF = pyfits.getdata('PSF_F150W_top{0}.fits'.format(folder[-1]))

#%%
fitsFile = pyfits.open(folder+'/imgblock.fits')
# fitsFile_host = pyfits.open('imgblock_host.fits')
qso_data = fitsFile[1].data
PS_image = fitsFile[2].data
# host_sersic_image = fitsFile_host[2].data

from galight.tools.astro_tools import plt_fits
plt_fits(qso_data-PS_image)

mag = -2.5*np.log10(np.sum(qso_data-PS_image)) + 29.5385
print(mag)


#%%
galfit_PSF_list, galight_PSF_list = [], []
for i in range(1,6):
    fitsFile = pyfits.open('psf{0}/imgblock.fits'.format(i))
    PS_image = fitsFile[2].data
    galfit_PSF_list.append(PS_image)
    

import glob, pickle
filt = 'F150W'
idx = 1
run_folder = '/Users/Dartoon/Astro/Projects/my_code/projects/2022_JWST_QSOz6/model_z6_data_id1/stage3_all/'
for top_psf_id in [0, 1, 2, 3, 4]:
        fit_run_list = []
        fit_files = glob.glob(run_folder+'*fit_material_super2/fit_run_fixn1__idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        galight_PSF_list.append(fit_run.image_ps_list[0])


from galight.tools.measure_tools import profiles_compare    
i = int(folder[-1])-1
profiles_compare([PSF]+galfit_PSF_list[i:i+1]+galight_PSF_list[i:i+1], if_annuli=True )
plt.show()

