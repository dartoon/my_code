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
                          rm_bkglight = False, if_plot=False, zp = zp)

data_process.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)

data_process.PSF_list = [PSF]

data_process.checkout() #Check if all the materials is known.

#%%Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor=3)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
fit_sepc.plot_fitting_sets()
fit_sepc.build_fitting_seq()


#%%Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'HSC_QSO', fitting_level='norm')
fit_run.run(algorithm_list = ['PSO'], setting_list=[None])
# fit_run.plot_all()
fit_run.dump_result()
# print(fit_run.final_result_galaxy[0])

#%%
import pickle
#links of file https://drive.google.com/file/d/1jE_6pZeDTHgXwmd2GW28fCRuPaQo8I61/view?usp=sharing
fit_run_pkl = pickle.load(open('./HSC_QSO.pkl','rb'))
from galight.tools.asymmetry_tools import CAS
CAS_class = CAS(fit_run_pkl, seg_cal_reg = 'or', obj_id=0)
# CAS_class.asy_segm(mask_type='aper')
# result = CAS_class.find_pos()
# asy = CAS_class.cal_asymmetry(rotate_pix = result["x"], if_remeasure_bkg=False ,if_plot=False, if_plot_bkg=False)
# print(asy)
from galight.tools.plot_tools import plt_fits
plt_fits(CAS_class.img,colorbar=True)
cas = CAS_class.cal_CAS(mask_type='segm', if_plot=True)
print(cas)

#%%
morph = fit_run_pkl.cal_statmorph(obj_id=0, segm=fit_run_pkl.fitting_specify_class.segm_deblend, if_plot = True)

from statmorph.utils.image_diagnostics import make_figure
fig = make_figure(morph)
print('xc_asymmetry =', morph.xc_asymmetry)
print('yc_asymmetry =', morph.yc_asymmetry)
print('ellipticity_asymmetry =', morph.ellipticity_asymmetry)
print('elongation_asymmetry =', morph.elongation_asymmetry)
print('orientation_asymmetry =', morph.orientation_asymmetry)
print('C =', morph.concentration)
print('A =', morph.asymmetry)
print('S =', morph.smoothness)