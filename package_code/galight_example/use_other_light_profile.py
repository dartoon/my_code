#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:45:16 2020

@author: Xuheng Ding
"""

#photutils in version > = 0.7.2
#astropy in version astropy-4.0.1

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

#%%
# object_id,ra,dec= '40158992189641871', 218.874090265527, -2.17422463704376
fitsFile = pyfits.open('./000017.88+002612.6/2-cutout-HSC-I-9469-s21a_wide.fits')

fov_image= fitsFile[1].data
header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
err_data= fitsFile[3].data ** 0.5

file_header0 = fitsFile[0].header
FLUXMAG0 = file_header0['FLUXMAG0']
zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
PSF = pyfits.getdata('./000017.88+002612.6/2-psf-calexp-s21a_wide-HSC-I-9469-4,2-0.07453-0.43684.fits')

#%%Start to use galight
from galight.data_process import DataProcess
QSO_RA, QSO_DEC = 0.07452999800443649, 0.4368380010128021
data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                           pos_type = 'wcs', header = header,
                          rm_bkglight = True, if_plot=False, zp = zp)

data_process.deltaPix = 1

data_process.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)

data_process.PSF_list = [PSF]

data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 0, supersampling_factor=3,
                               extend_source_model = ['LINEAR_ELLIPSE', 'LINEAR_ELLIPSE', 'LINEAR_ELLIPSE'],
                              apertures_center_focus = False)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
# # fit_sepc.plot_fitting_sets()
fit_sepc.build_fitting_seq()
fit_sepc.kwargs_params['lens_light_model'][0][-1]['k'] = fit_sepc.kwargs_params['lens_light_model'][0][-1].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][1][-1]['k'] = fit_sepc.kwargs_params['lens_light_model'][1][-1].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][3][-1]['k'] = fit_sepc.kwargs_params['lens_light_model'][3][-1].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][4][-1]['k'] = fit_sepc.kwargs_params['lens_light_model'][4][-1].pop('R_sersic')

fit_sepc.kwargs_params['lens_light_model'][0][-2]['k'] = fit_sepc.kwargs_params['lens_light_model'][0][-2].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][1][-2]['k'] = fit_sepc.kwargs_params['lens_light_model'][1][-2].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][3][-2]['k'] = fit_sepc.kwargs_params['lens_light_model'][3][-2].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][4][-2]['k'] = fit_sepc.kwargs_params['lens_light_model'][4][-2].pop('R_sersic')

fit_sepc.kwargs_params['lens_light_model'][0][-3]['k'] = fit_sepc.kwargs_params['lens_light_model'][0][-3].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][1][-3]['k'] = fit_sepc.kwargs_params['lens_light_model'][1][-3].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][3][-3]['k'] = fit_sepc.kwargs_params['lens_light_model'][3][-3].pop('R_sersic')
fit_sepc.kwargs_params['lens_light_model'][4][-3]['k'] = fit_sepc.kwargs_params['lens_light_model'][4][-3].pop('R_sersic')


#%%Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc)
fit_run.run(algorithm_list = ['PSO','PSO','PSO'], fitting_level=['norm','deep','deep'])
fit_run.plot_final_galaxy_fit()
# fit_run.dump_result()
# print(fit_run.final_result_galaxy[0])

#%%
plt.imshow(fit_run.image_host_list[1], origin='lower')
plt.show()
# #%%
# import pickle
# #links of file https://drive.google.com/file/d/1jE_6pZeDTHgXwmd2GW28fCRuPaQo8I61/view?usp=sharing
# fit_run_pkl = pickle.load(open('./HSC_QSO.pkl','rb'))
