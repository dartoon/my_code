#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:45:49 2021

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os, copy, glob
from subprocess import call
from matplotlib.colors import LogNorm
from decomprofile.tools.measure_tools import find_loc_max, measure_FWHM, twoD_Gaussian, fit_data_twoD_Gaussian #, esti_bgkstd
from decomprofile.data_process import DataProcess
from decomprofile.fitting_specify import FittingSpeficy
from decomprofile.fitting_process import FittingProcess
from astropy.wcs import WCS

f = open("../cut_out.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

files = glob.glob('../z_over1/*')
files.sort()

run_i = 0
file = files[run_i]
image_ID = file.split('/')[-1]
line = [lines[i] for i in range(len(lines)) if image_ID in lines[i]]

_, image_RA, image_DEC = line[0].split(' ')
image_RA = float(image_RA)
image_DEC = float(image_DEC)
print("run_i = ", run_i)
print(image_ID, image_RA, image_DEC)
#%%
deep_seed = True  #Set as True to put more seed and steps to fit,
show_plot = 1
fit_data = True  #If you simply want to do the search without fitting, set False

image_folder = '../z_over1/' + image_ID + '/'
fit_folder = image_folder + 'fit_result/'

print(fit_folder)
if os.path.exists(fit_folder)==False:
    os.mkdir(fit_folder)

import shutil
shutil.copy('./fit_dual_allband.py', fit_folder)
#If only want to run I band
# band_seq = ['I'] 
# run_list = [0] 
band_seq = ['I', 'G', 'R', 'Z', 'Y']
run_list = [0, 1, 2, 3, 4]
filename_list = [image_ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(len(band_seq))]

data_process_list, zp_list = [], []
for i in range(len(band_seq)):
    # The pixel scale is all 0.168
    if len(glob.glob(image_folder+filename_list[i])) == 0:
        print(filename_list[i] + " DOES NOT EXIST!!!")
        QSO_im, err_map, PSF, _, _, qso_center, fr_c_RA_DEC = [], [], [], [], [], [], []
        run_list.remove(i)
        data_process_list.append(None)
        zp_list.append(None)
    else:
        fitsFile = pyfits.open(image_folder+filename_list[i])
        fov_image= fitsFile[1].data
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        err_data= fitsFile[3].data ** 0.5
        file_header0 = fitsFile[0].header
        zp =  27.0 
        data_process_i = DataProcess(fov_image = fov_image, fov_noise_map = err_data,
                                     target_pos = [image_RA, image_DEC],
                                     pos_type = 'wcs', header = header,
                                     rm_bkglight = True, if_plot=False, zp = zp)
        data_process_i.noise_map = err_data
        data_process_i.generate_target_materials(radius=None, radius_list = [10, 20, 30, 40],
                                                  create_mask = False, nsigma=1,
                                              exp_sz= 1.2, npixels = 9, if_plot=False)        
        PSF = pyfits.getdata(image_folder+filename_list[i].split('.fits')[0]+'_psf.fits')
        if len(PSF) != 0 and PSF.shape[0] != PSF.shape[1]:
            cut = int((PSF.shape[0] - PSF.shape[1])/2)
            if cut>0:
                PSF = PSF[cut:-cut,:]
            elif cut<0:
                PSF = PSF[:,-cut:cut]
            PSF /= PSF.sum()
            if PSF.shape[0] != PSF.shape[1]:
                raise ValueError("PSF shape is not a square.")
        data_process_i.PSF_list = [PSF]
        data_process_list.append(data_process_i)
        zp_list.append(zp) 

#%%
frame = np.max([len(data_process_list[i].target_stamp) for i in run_list])
radius = int((frame-1)/2)
for k in run_list:
    data_process_list[k].generate_target_materials(radius=20,
                                                   create_mask = False, nsigma=2.,
                                                   exp_sz= 1.2, npixels = 15, if_plot=False)
    apertures_temp = data_process_list[k].apertures
    if k == run_list[0]:
        apertures = apertures_temp
    if k != run_list[0] and len(apertures_temp)>1:
        for i in range(len(apertures_temp)):
            count = 0
            for j in range(len(apertures)):
                dis = np.sqrt(np.sum(apertures[j].positions-apertures_temp[i].positions)**2)
                if dis < 5:  #If objects is close within 5 pixels consider as a sample obj
                    count += 1
            if count == 0:
                apertures.append(apertures_temp[i])


#%%
fig, (axs) = plt.subplots(1, 5, figsize=(15, 7))
vmin = 1.e-3
vmax = data_process_list[0].target_stamp.max() * 5
color_list = ['winter', 'summer', 'afmhot', 'autumn', 'gist_heat']
plt_list = [1, 2, 0, 3, 4]
for i in range(len(plt_list)):
    p_i = plt_list[i]
    if data_process_list[p_i] != None:
        axs[i].imshow(data_process_list[p_i].target_stamp, origin='lower', cmap=color_list[i], norm=LogNorm(), vmin=vmin, vmax=vmax)
    else:
        axs[i].imshow(data_process_list[0].target_stamp * 0)
    axs[i].set_title('{0} band'.format(band_seq[p_i]))
# plt.savefig(fit_folder + 'images_5_band.pdf')
plt.show()   
from astropy.visualization import make_lupton_rgb
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval

r_i, g_i, b_i = [data_process_list[i].target_stamp for i in [0, 2, 1]]

# stretch = SqrtStretch() #+ ZScaleInterval()
# r_i = stretch(r_i)
# g_i = stretch(g_i)
# b_i = stretch(b_i)

rgb_default = make_lupton_rgb(r_i, g_i, b_i, stretch = 1)
plt.imshow(rgb_default, origin='lower')
