#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 14:44:04 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle

import sys
sys.path.insert(0,'../model_z6_data_id0/')

run_folder = 'stage3_all/' #!!!
idx = 1
# filt = 'F356W'
filt = 'F150W'


#Load PSF information:
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))


#Load top 5 PSF idx
fit_files = glob.glob(run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
fit_files.sort()
fit_run_list = []
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
sort_Chisq = chisqs.argsort() 

prop_name = 'magnitude'
all_values = np.array([fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))])
print(all_values[sort_Chisq[:5]])


folder = '../NIRCam_data/Nov14/bkg_removed'
from target_info import target_info
info = target_info[str(idx)]

file = glob.glob(folder+'/*{1}*{0}*.fits'.format(filt, info['target_id'][:5]))[0]

#%%

from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess

run_test = False

if run_test == True:
    for i in range(5):
        for j in range(i+1,5):
            use_psf_id = int(fit_files[sort_Chisq[i]].split('psf')[1].split('_')[0])
            fit_psf_id = int(fit_files[sort_Chisq[j]].split('psf')[1].split('_')[0])
            print(i, j)
            print(use_psf_id, fit_psf_id)
            _fitsFile = pyfits.open(file)
            fov_image = _fitsFile[1].data # check the back grounp
            header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            flux_mjsr = header['PHOTMJSR']
            from galight.tools.astro_tools import read_pixel_scale
            pixscale = read_pixel_scale(header)
            zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
            wht = _fitsFile[4].data # The WHT map
            exp = _fitsFile[0].header['EFFEXPTM']
            print("Exp time:", exp)
            if _fitsFile[0].header['CHANNEL'] == 'LONG':
                gain_value = 2
                expsize = 1
            else:
                gain_value = 1.8
                expsize = 1.4
            exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
            print("Processing data...")
            psf = PSF_list_clean[use_psf_id]
            
            data_process = DataProcess(fov_image = fov_image, target_pos = [PSF_RA_DEC_list[fit_psf_id][0], PSF_RA_DEC_list[fit_psf_id][1]], #The final cut center is dete
                                            pos_type = 'wcs', header = header, rm_bkglight = False, 
                                            if_plot=False, zp = zp, exptime= exp_map, 
                                            fov_noise_map = None)
            if filt == 'F356W':
                radius = 40
            elif filt == 'F150W':
                radius = 20
            
            data_process.generate_target_materials(radius=radius, create_mask = False,
                                                    cut_kernel = None, if_select_obj=False,
                                                    if_plot=False)
            ct = int(len(psf)/2) - radius
            data_process.target_stamp = PSF_list_clean[fit_psf_id][ct:-ct,ct:-ct]
            
            data_process.apertures = []
            del data_process.fov_image
            del data_process.exptime
            data_process.filt = filt
            data_process.file = file
            data_process.plot_aperture()
            psf[psf<0] = 0.
            data_process.PSF_list = [psf]
            fit_sepc = FittingSpecify(data_process)
            fit_sepc.prepare_fitting_seq(point_source_num = 1)
            fit_run = FittingProcess(fit_sepc, savename = 'PSFfitPSF_use{0}fit{1}'.format(i,j))
            fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
            fit_run.plot_final_qso_fit(target_ID ='PSF fit the other PSF', save_plot=False)
            savename = 'PSFvPSF_checks/' + 'PSF_fit_PSF'+'_'+filt +'_use{0}fit{1}'.format(i,j)+'.pkl' 
            pickle.dump(fit_run , open(savename, 'wb'))
        
#%%

PSF_test_files = glob.glob('PSFvPSF_checks/*'+filt[:-1]+'*.pkl'.format(idx))
PSF_test_files.sort()
images = []
for i in range(len(PSF_test_files)):
    fit_run = pickle.load(open(PSF_test_files[i],'rb'))
    images.append([fit_run.flux_2d_out['data-point source'],fit_run.flux_2d_out['data']])
    
# from galight.tools.astro_tools import plt_many_fits
# plt_many_fits(images)

from matplotlib.colors import LogNorm

fig, (axs) = plt.subplots(2, 5, figsize=(15, 6))
for i in range(len(images)):
    _i = int(i / 5)
    _j = int(i % 5)
    im_i = axs[_i][_j].imshow(images[i][0]/np.sum(images[i][1]) * 280, origin='lower',vmin = -0.2, vmax = 0.2  , cmap='bwr')
    frame_size = len(images[i][0])
    # cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])
    fig.colorbar(im_i, ax=axs[_i][_j], pad=0.01, shrink=0.8)
    
    # label = labels[i]
    use_i = PSF_test_files[i].split('use')[1][0]
    fit_j = PSF_test_files[i].split('fit')[-1][0]
    label = 'use{0}fit{1}'.format(use_i,fit_j)
    plttext = axs[_i][_j].text(frame_size*0.05, frame_size*0.87, label,
              fontsize=15, weight='bold', color='black')
    plttext.set_bbox(dict(facecolor='white', alpha=0.5))
    # if texts is not None:
    #     plttext = axs[_i][_j].text(frame_size*0.05, frame_size*0.05, "{1} = {0}".format(round(texts[i],3), prop ),
    #              fontsize=label_size, weight='bold', color='black')
    #     plttext.set_bbox(dict(facecolor='white', alpha=0.5))
    axs[_i][_j].axes.xaxis.set_visible(False)
    axs[_i][_j].axes.yaxis.set_visible(False)
    fig.tight_layout()

# plt.savefig(savename,bbox_inches='tight')
plt.show()    

        
    