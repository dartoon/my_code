#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 16:40:31 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings, pickle
warnings.filterwarnings("ignore")
from astropy.wcs import WCS
from galight.tools.astro_tools import plt_many_fits
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle, copy, glob

filt_i = 0
filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
fitsFile = pyfits.open(filefolder+filename)
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
img = fitsFile[1].data #
print(img.shape)
wcs = WCS(header)
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

#%%
PSF_lib_files = glob.glob('material/'+filt+'_PSF_Library.pkl')[0]
PSF_org_list, PSF_clean_list, all_PSF_pos_list, PSF_RA_DEC_list = pickle.load(open(PSF_lib_files,'rb'))
# plt_many_fits(PSF_clean_list)
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
fov_noise_map = fitsFile[2].data
# for i in [0]:
#%%    
from galight.tools.astro_tools import plt_fits

# run_i = 46
# for i in range(run_i,len(cata_list)):
for i in [10]:
    pos = cata_list[i][3:5]
    plt_fits(img[int(pos[1])-100:int(pos[1])+100, int(pos[0])-100:int(pos[0])+100])
    radius_list = [60,70,80]
    npixels = 50
    if i in [6,11,16,34,36,37]:
        radius_list = [60,70]
    if i in [18,27,38]:
        radius_list = [80,90]
    if i in [17]:
        radius_list = [35]
    
    data_process = DataProcess(fov_image = img, target_pos = pos, pos_type = 'pixel', header = header,
                              rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = fov_noise_map)
    data_process.generate_target_materials(radius=80,radius_list=radius_list, 
                                           create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = npixels, if_plot=True)
    data_process.noise_map[data_process.noise_map == 0] = np.max(data_process.noise_map[data_process.noise_map!=0])
    if i == 27:
        del data_process.apertures[1]
    if i == 46:
        del data_process.apertures[2]
        # data_process.apertures[2] = data_process.apertures[3]
        # data_process.apertures[3] = ap
    data_process.filt = filt
    del data_process.fov_image
    del data_process.fov_noise_map
    # data_process.plot_aperture()
    # pickle.dump(data_process, open('material/'+'data_process_idx{0}_{1}.pkl'.format(i, filt), 'wb'))
    
    sz = 30
    ct = int((len(data_process.target_stamp) - 30)/2)
    test_img = data_process.target_stamp[ct:-ct,ct:-ct]
    ps_pos = np.array(np.where(test_img == test_img.max())) - len(test_img)/2
    ps_pos = ps_pos[::-1]
    target_id = cata_list[i][7]
    for j in range(len(PSF_clean_list)):
        savename = 'fit_material/fit_notrunyet_{2}_idx{0}_psf{1}.pkl'.format(i,j,filt)
        psf = PSF_clean_list[j]
        psf = psf[1:,1:]
        if filt == 'F150W':
            ct = int((len(psf) - len(data_process.target_stamp))/2)
            if ct >10:
                psf = psf[ct:-ct, ct:-ct]
        psf[psf<0] = 0.
        data_process.PSF_list = [psf]
        fit_sepc = FittingSpecify(data_process)
        
        fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3, apertures_center_focus=True,
                                      ps_pix_center_list = [copy.deepcopy(ps_pos)]  )
        # fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
        # fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
        # fit_sepc.kwargs_constraints['linear_solver'] = False
        if j == 0:
            fit_sepc.plot_fitting_sets()
            print(i, filt, 'apertures', len(data_process.apertures) )
            hold = input('Hold ... OK?\n')
        fit_run = FittingProcess(fit_sepc, savename = target_id)
        pickle.dump(fit_run , open(savename, 'wb'))
        # savename = savename.replace('_notrunyet_', '_run_')[:-4]+'_{0}.pkl'.format(i)
        # fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
        # fit_run.plot_final_qso_fit(target_ID =target_id)
        