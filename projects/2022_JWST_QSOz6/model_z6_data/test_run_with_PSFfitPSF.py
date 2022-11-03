#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:19:52 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt




import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
import pickle
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.measure_tools import measure_bkg
from astropy.wcs import WCS
from galight.tools.cutout_tools import common_data_class_aperture
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import cutout
import warnings
warnings.filterwarnings("ignore")


data_type = 'all' 
filt = 'F356W'
file_NO = 0

run_folder = 'stage3_all*/' #!!!
idx = 0
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))


idx = 0
folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_J2255_stage3_{0}/bkg_removed'.format(data_type)
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
jwst_all_filenames = glob.glob(folder+'/*{0}*.fits'.format(filt))
jwst_all_filenames.sort()
file = jwst_all_filenames[file_NO]
if data_type == 'all':
    run_folder = 'stage3_{0}/'.format(data_type)
elif data_type == 'half':
    if file_NO == 0:
        run_folder = 'stage3_first_half/'
    if file_NO == 1:
        run_folder = 'stage3_second_half/'

result_folder = run_folder + 'fit_result/'
#%%
cut_kernel = None #After pos correct then, do the nearest_obj_center
_fitsFile = pyfits.open(file)
fov_image = _fitsFile[1].data # check the back grounp
header = _fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
flux_mjsr = header['PHOTMJSR']
pixscale = read_pixel_scale(header)
zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
wht = _fitsFile[4].data # The WHT map
exp = _fitsFile[0].header['EFFEXPTM']
print("Exp time:", exp)
if _fitsFile[0].header['CHANNEL'] == 'LONG':
    gain_value = 2
    expsize = 1
    exppix = 1
else:
    gain_value = 1.8
    expsize = 1.4
    exppix = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
print("Processing data...")
data_process = DataProcess(fov_image = fov_image, target_pos = [PSF_RA_DEC_list[2][0], PSF_RA_DEC_list[2][1]], #The final cut center is dete
                               pos_type = 'wcs', header = header, rm_bkglight = False, 
                               if_plot=False, zp = zp, exptime= exp_map, 
                               fov_noise_map = None)

#%%

data_process.generate_target_materials(radius=40 * expsize, create_mask = False, nsigma=2.5, 
                                        cut_kernel = None, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 60 * expsize, if_plot=False)

data_process.apertures = []
# data_process.apertures[0].theta = 0
# data_process.apertures[0].b = data_process.apertures[0].b/2
# data_process.apertures[0].positions[0] = data_process.apertures[0].positions[0] - 1.5

del data_process.fov_image
del data_process.exptime
data_process.filt = filt
data_process.file = file
data_process.plot_aperture()

psf = PSF_list_clean[0]
psf[psf<0] = 0.
data_process.PSF_list = [psf]

from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3)
                              # ps_pix_center_list = [ps_pos]  ) #, fix_n_list= [[0,4],[1,1]])
# fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
# fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 0.6
# fit_sepc.kwargs_constraints['linear_solver'] = False
# fit_sepc.plot_fitting_sets()
fit_run = FittingProcess(fit_sepc, savename = 'PSFfitPSF')
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.plot_final_qso_fit(target_ID ='PSFfitPSF')
