#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.astro_tools import plt_fits
from galight.data_process import DataProcess
ID = 10004
band = 'I'
folder = 'fit_result/'
file_ = glob.glob(folder+"{0}-{1}.pkl".format(ID, band))
data_folder = '/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/'
# data_folder = './'
if file_ != []:
    file = file_[0]
    fit_run = pickle.load(open(file,'rb'))
    print(fit_run.final_result_galaxy)
    host_image = fit_run.flux_2d_out['data'] - fit_run.image_ps_list[0]
    plt_fits(host_image)
    
    ID = file.split('/')[1].split('-')[0]
    band = file.split('-')[1][0]
    fits_name = glob.glob(data_folder+'gfarm_data_download/{0}_HSC-{1}.fits'.format(ID,band)) + glob.glob(data_folder+
        'online_data_download/{0}/*cutout*-{1}-*.fits'.format(ID,band) )
    
    f = open("catalog.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    line = [lines[i] for i in range(len(lines)) if lines[i].split(' ')[0] == ID][0]
    _, Ra, Dec = line.split(' ')
    QSO_RA, QSO_DEC = float(Ra), float(Dec)
    fitsFile = pyfits.open(fits_name[0])
    file_header0 = fitsFile[0].header
    zp = 27.0
    data_process = DataProcess(fov_image = fitsFile[1].data, fov_noise_map = fitsFile[3].data ** 0.5, target_pos = [QSO_RA, QSO_DEC],
                                pos_type = 'wcs', header = fitsFile[1].header,
                                rm_bkglight = True, if_plot=False, zp = zp)
    fov_image = data_process.fov_image

    fit_run.targets_subtraction(sub_qso_list=[0], org_fov_data=fov_image, header=data_process.header, target_pos=data_process.target_pos)
    fov_image_targets_sub = fit_run.fov_image_targets_sub #Host image after sub PSF.
    plt_fits(fit_run.fov_image_targets_sub)