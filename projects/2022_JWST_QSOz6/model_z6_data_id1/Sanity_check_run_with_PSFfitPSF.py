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
import sys
sys.path.insert(0,'../model_z6_data_id0/')

data_type = 'all' 
filt = 'F150W'
file_NO = 0

run_folder = 'stage3_all/' #!!!
idx = 1
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))


folder = '../NIRCam_data/Nov14/'
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# jwst_all_filenames = glob.glob(folder+'/*{0}*{1}*.fits'.format(target_id[:5], filts[0]))
jwst_all_filenames = glob.glob(folder+'/*{0}*.fits'.format(filt))
jwst_all_filenames.sort()
file = jwst_all_filenames[file_NO]
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
data_process = DataProcess(fov_image = fov_image, target_pos = PSF_RA_DEC_list[10], #array([10,  2, 13])
                               pos_type = 'wcs', header = header, rm_bkglight = False, 
                               if_plot=False, zp = zp, exptime= exp_map, 
                               fov_noise_map = None)

#%%

data_process.generate_target_materials(radius=30 * expsize, create_mask = False, nsigma=1.2, 
                                        cut_kernel = None, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 60 * expsize, if_plot=False)

data_process.apertures = []

del data_process.fov_image
del data_process.exptime
data_process.filt = filt
data_process.file = file
data_process.plot_aperture()

psf = PSF_list_clean[2]
psf[psf<0] = 0.
data_process.PSF_list = [psf]

from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3)
fit_run = FittingProcess(fit_sepc, savename = 'PSFfitPSF')
fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
fit_run.savename = 'figures/' + 'PSF_fit_PSF'+'_'+filt 
fit_run.plot_final_qso_fit(target_ID ='PSF fit the other PSF', save_plot=True)
