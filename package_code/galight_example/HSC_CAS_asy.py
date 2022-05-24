#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:08:21 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from galight.tools.asymmetry_tools import Measure_asy
from galight.tools.plot_tools import plt_fits
import pickle
# fit_run = pickle.load(open('./HSC_QSO.pkl','rb'))
fit_run = pickle.load(open('40968511920540721/result-band-I.pkl','rb'))
# fit_run.fitting_specify_class.plot_fitting_sets()
# data_process = fit_run.fitting_specify_class.data_process_class

# asy_class = Measure_asy(fit_run, seg_cal_reg = 'or', obj_id=0)
# asy_class.asy_segm(mask_type='aper')
# result = asy_class.find_pos()
# print(result["x"])
# plt_fits(asy_class.img,colorbar=True)
# asy_class.make_bkg(rotate_pix = result["x"], if_remeasure_bkg=False ,)
# asy = asy_class.cal_asymmetry(rotate_pix = result["x"], if_plot=True, if_plot_bkg=True)
# print('asymmetry :', asy)

from galight.tools.asymmetry_tools import CAS
CAS_class = CAS(fit_run, seg_cal_reg = 'or', obj_id=0, rm_ps=False)
# CAS_class.asy_segm(mask_type='aper')
# result = CAS_class.find_pos()
# asy = CAS_class.cal_asymmetry(rotate_pix = result["x"], if_remeasure_bkg=False ,if_plot=False, if_plot_bkg=False)
# print(asy)
# plt_fits(CAS_class.img,colorbar=True)
cas = CAS_class.cal_CAS(mask_type='segm',extend=1, if_plot=False)
print(CAS_class.find_pos_res["x"])
print(cas)

#%%
#Calculate the smoothness for the residual image:
from galight.tools.asymmetry_tools import CAS
CAS_class = CAS(fit_run, seg_cal_reg = 'or', obj_id=0, rm_model=True)  #Seeting rm_model = True for residual.
#Note that in the following equation, the image_org should be used to calculate the flux value with 0.25 and 1.25 of the petrosian radius.
img = fit_run.fitting_specify_class.kwargs_data['image_data']
cas = CAS_class.cal_CAS(mask_type='segm',extend=1, if_plot=False, if_residual=True, image_org = img)
print(CAS_class.find_pos_res["x"])
print(cas)
# print('smoothness (for abs diff and pos diff):', cas[1])