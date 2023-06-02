#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:45:01 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle, copy, matplotlib
from matplotlib.colors import LogNorm

cmap = 'gist_heat'
my_cmap_0 = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
my_cmap_0.set_bad('black')
cmap = 'inferno'
my_cmap_1 = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
my_cmap_1.set_bad('black')

plt.rcParams["font.family"] = "sans-serif"

# Generate some random data
data = np.random.rand(4, 10, 10)

# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2,  figsize=(10, 8))

run_folder = '../model_z6_data_id{0}/stage3_all/'.format(1) 

filts = ['F356W', 'F150W']
usecodes = ['galight', 'galfit']
zp_list = [27.98081241339078, 29.538598992275794]

# Plot the data in each subplot
for i in range(2):
    if i == 0:
        my_cmap = my_cmap_0
    if i == 1:
        my_cmap = my_cmap_1
    for j in range(2):
        if j == 0:
            fit_run_list = []
            filt = filts[i]
            fit_files = glob.glob('run_galfit/Test_{0}/{0}_fit_psf?.pkl'.format(filt))#+\
            # fit_files = glob.glob('run_galfit/Test_{0}/{0}_fit_psf?_ssf1.pkl'.format(filt))#+\
            # fit_files = glob.glob(run_folder+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(1, filt))#+\
            fit_files.sort()
            for ii in range(len(fit_files)):
                fit_run_list.append(pickle.load(open(fit_files[ii],'rb')))
            chisqs = np.array([fit_run_list[ii].reduced_Chisq for ii in range(len(fit_run_list))])
            sort_Chisq = chisqs.argsort()  
            # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
            sort_Chisq = [0]
            fit_run = fit_run_list[sort_Chisq[0]]
            image = fit_run.flux_2d_out['data-point source']
            im = axs[i, j].imshow(image, norm = LogNorm(vmax = 1.1, vmin = 1.e-4), cmap=my_cmap, origin = 'lower')
            
            #If use host residual
            # t= axs[i, j].text(2,80, 'host mag {0:.2f}'.format(fit_run.final_result_galaxy[0]['magnitude']), color='white', fontsize=23)
            
            model_mag = -2.5*np.log10(fit_run.final_result_galaxy[0]['flux_sersic_model']) + fit_run.zp
            t= axs[i, j].text(2,80, 'host mag {0:.2f}'.format(model_mag), color='white', fontsize=23)
            best_chisq = fit_run.reduced_Chisq            
            # print(filt,'zp:', fit_run.zp)
        if j == 1:    
            zp = zp_list[i]
            all_log = glob.glob('run_galfit/Test_{0}/psf?/fit.log'.format(filt))
            all_log.sort()
            chisq_list = []
            mag_list = []
            for ii in range(len(all_log)):
                filename = all_log[ii]
                f = open(filename,"r")
                string = f.read()
                lines = string.split('\n')   # Split in to \n
                line = [lines[l] for l in range(len(lines)) if 'Chi^2/nu' in lines[l]][-1]
                chisq = float(line.split('=')[-1])
                chisq_list.append(chisq)
                # line_sersic = [lines[l] for l in range(len(lines)) if 'sersic' in lines[l]][-1]
                # mag = float(line_sersic.split(')')[-1][:10])
                # mag_list.append(mag)
                fitsFile = pyfits.open('run_galfit/Test_{0}/psf{1}/imgblock.fits'.format(filt,ii+1))
                sersic = fitsFile[2].data
                mag_list.append(-2.5*np.log10(np.sum(sersic))+zp)
            print(filt, mag_list)
            print(np.mean(mag_list), np.std(mag_list))
            # psf_i = np.where(chisq_list==np.min(chisq_list))[0][0]
            psf_i = sort_Chisq[0]
            best_chisq = chisq_list[psf_i]
            fitsFile = pyfits.open('run_galfit/Test_{0}/psf{1}/imgblock_ps.fits'.format(filt,psf_i+1))
            image = fitsFile[3].data
            fitsFile = pyfits.open('run_galfit/Test_{0}/psf{1}/imgblock.fits'.format(filt,psf_i+1))
            sersic = fitsFile[2].data
            # qso_data = fitsFile[1].data
            # PS_image = fitsFile[2].data
            # image = qso_data-PS_image
            im = axs[i, j].imshow(image, norm = LogNorm(vmax = 1.1, vmin = 1.e-4), 
                                  cmap=my_cmap, origin = 'lower')
            #If use host residual
            t= axs[i, j].text(2,80, 'host mag {0:.2f}'.format(-2.5*np.log10(np.sum(sersic))+zp), color='white', fontsize=23)
            
            # t= axs[i, j].text(2,80, 'host mag {0:.2f}'.format(mag_list[psf_i]), color='white', fontsize=23)
        t.set_bbox(dict(facecolor='black', alpha=0.8))    
        # axs[i, j].text(2,5, r'fitting $\chi^2$ = {0:.2f}'.format(best_chisq), color='white', fontsize=23)
        cb_ij = fig.colorbar(im, ax=axs[i, j], pad=0.01)            
        axs[i, j].set_title('J2236 fit by {1}, {0}'.format(filt, usecodes[j]), fontsize = 20)
        axs[i, j].tick_params(axis='both', labelsize = 15)
        cb_ij.ax.tick_params(labelsize=15) 
        
        
# Adjust the spacing between subplots
plt.tight_layout()
plt.savefig('../model_z6_data_id0/figures/galight_vs_galfit.pdf')# 
# Display the figure
plt.show()


#%%Pinrt the result for J2255
mag_list = []
filt = 'F356W'
zp = 27.98081241339078
for ii in range(5):
    fitsFile = pyfits.open('../model_z6_data_id0/run_galfit/Test_{0}/psf{1}/imgblock.fits'.format(filt,ii+1))
    sersic = fitsFile[2].data
    mag_list.append(-2.5*np.log10(np.sum(sersic))+zp)
print(np.mean(mag_list),np.std(mag_list))
mag_list = []

filt = 'F150W'
zp = 29.538598992275794
for ii in range(5):
    fitsFile = pyfits.open('../model_z6_data_id0/run_galfit/Test_{0}/psf{1}/imgblock.fits'.format(filt,ii+1))
    sersic = fitsFile[2].data
    mag_list.append(-2.5*np.log10(np.sum(sersic))+zp)
print(np.mean(mag_list),np.std(mag_list))
