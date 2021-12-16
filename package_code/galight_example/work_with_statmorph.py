#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:00:49 2021

@author: Dartoon
"""


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import copy
import glob
import lenstronomy.Util.param_util as param_util
import pickle
from galight.tools.plot_tools import profile_plots

fitsFile_ = glob.glob('/Users/Dartoon/Astro/Projects/my_code/projects/2020_HSC_bulge_fit/SDSS_0.2-0.3/*_HSC-I.fits')
fitsFile_.sort()
# NO = 1
# for NO in range(21):
# for NO in range(21,42):
# for NO in range(0,63):
# for NO in [8]: 
# for NO in range(11,36):
#60, 34, 20    
fitsFile  = fitsFile_[34]
ID = fitsFile.split('/')[-1].split('_')[0]
PSF_filename = fitsFile.split('.fits')[0]+ '_psf.fits'
save_name = 'test_folder/'
f = open("/Users/Dartoon/Astro/Projects/my_code/projects/2020_HSC_bulge_fit/SDSS_0.2-0.3/SDSS_0.2-0.3.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
line = [i for i in range(len(lines)) if lines[i].split(' ')[0] == ID][0]
QSO_RA, QSO_DEC = lines[line].split(' ')[1:]

fits = pyfits.open(fitsFile)
fov_image= fits[1].data
header = fits[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
err_data= fits[3].data ** 0.5

file_header0 = fits[0].header
zp = 27.0# This is something Xuheng can't make sure.
PSF = pyfits.getdata(PSF_filename)

 # Start to use galight
picklename = save_name+ID+'_single_Sersic.pkl'
picklename_l = glob.glob(picklename)
if picklename_l != []:
    fit_run_0 = pickle.load(open(picklename_l[0],'rb'))
else:
    data_process_0 = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [float(QSO_RA), float(QSO_DEC)],
                                pos_type = 'wcs', header = header,
                              rm_bkglight = True, if_plot=False, zp = zp)
    
    data_process_0.noise_map = err_data
    
    data_process_0.generate_target_materials(radius=None, detect_tool='sep', create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 15, if_plot=True)
    data_process_0.PSF_list = [PSF]
    data_process_0.checkout() #Check if all the materials is known.
    
    fit_sepc_0 = FittingSpecify(data_process_0)
    fit_sepc_0.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    fit_sepc_0.plot_fitting_sets()
    fit_sepc_0.build_fitting_seq()
    # Setting the fitting method and run.
    fit_run_0 = FittingProcess(fit_sepc_0, savename = save_name+ID+'_single_Sersic')
    fit_run_0.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 50, 'n_iterations': 50}]*2)
    fit_run_0.plot_final_qso_fit()
    bic_0 = fit_run_0.fitting_seq.bic
    fit_run_0.dump_result()

#%%
fit_run_0.plot_final_qso_fit()
# fit_run_0.fitting_specify_class.plot_fitting_sets()
# obj_id = int(input('Which obj to measure using statmorph?\n'))
# morph = fit_run_0.cal_statmorph(obj_id=obj_id, if_plot=True)

# from statmorph.utils.image_diagnostics import make_figure
# fig = make_figure(morph)
# plt.show()
# print('xc_asymmetry =', morph.xc_asymmetry)
# print('yc_asymmetry =', morph.yc_asymmetry)
# print('ellipticity_asymmetry =', morph.ellipticity_asymmetry)
# print('elongation_asymmetry =', morph.elongation_asymmetry)
# print('orientation_asymmetry =', morph.orientation_asymmetry)
# print('C =', morph.concentration)
# print('A =', morph.asymmetry)
# print('S =', morph.smoothness)