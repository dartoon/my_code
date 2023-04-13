#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 17:16:06 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle, glob
import warnings
warnings.filterwarnings("ignore")
from matplotlib.colors import LogNorm
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.astro_tools import plt_fits_color, plt_fits
import pickle      
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
from astropy.cosmology import FlatLambdaCDM
from astropy.visualization import make_lupton_rgb

def scale_bar(ax, d, dist=1/0.13, text='1"', text2=None, color='black', flipped=False, fontsize=20):
    p0 = d / 7.
    ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    if text2 is not None:
        ax.text(p0 + dist / 2., p0 - 0.08 * d, text2, fontsize=fontsize, color=color, ha='center')
        

# cata_list = pickle.load(open('../material/cata_list.pkl','rb'))
# check_name= 'cid_473'  #29
# check_name= 'cid_1210' #8
# check_name= 'cid_1245' #10
# check_id = [i for i in range(len(cata_list)) if cata_list[i][-1] == check_name]
# print(cata_list[check_id[0]])

f = open('fmos_alma_cosmosweb.cat','r')
string = f.read()
lines = string.split('\n')
lines = [lines[i] for i in range(len(lines)) if 'FMOS_J09' in lines[i]]


from galight.data_process import DataProcess
from scipy.ndimage import zoom
filts = ['F814W','F115W', 'F150W','F277W', 'F444W']
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
size = 90  #For LW it is 56*4 and for SW it is 56*2
sed_2d_info  = []

#%% Making cutout
image_list = [None] * len(filts)
zp_dict = {}
# idx = 0

# line = lines[idx]
for line in lines:
    for i, filt in enumerate(filts):
        name, RA, Dec, z, best_mass = line.split(' ')
        RA, Dec, z = float(RA), float(Dec), float(z)
        if filt == 'ACS':
            filename = 'mosaic_cosmos_web_2023jan_30mas_hst_acs_wfc_f814w_drz.fits'
            fitsFile = pyfits.open(filefolder+filename)
            header = fitsFile[0].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            img = fitsFile[0].data #
            #Use acstools to read HST zp:
            # from acstools import acszpt
            # q = acszpt.Query(date="2016-04-01", detector="WFC")
            # zpt_table = q.fetch()
            # print(zpt_table)
            zp = 25.937
            pixscale = read_pixel_scale(header)
            data_process = DataProcess(fov_image = img, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                                      rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = img)
            data_process.generate_target_materials(radius=size, skip = True, if_plot=False)
        
        elif filt != 'ACS':
            filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
            fitsFile = pyfits.open(filefolder+filename)
            header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            img = fitsFile[1].data #
            flux_mjsr = header['PHOTMJSR']
            pixscale = read_pixel_scale(header)
            zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
            data_process = DataProcess(fov_image = img, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                                      rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = img)
            data_process.generate_target_materials(radius=size, skip = True, if_plot=False)
        
        img_show = zoom(data_process.target_stamp, 0.5)
        img_show = img_show/np.sum(img_show) * np.sum(data_process.target_stamp)
        image_list[i] = img_show
        print(filt,'finish', img_show.shape)
        zp_dict[filt] = zp
    #plot image
    images = []
    use_filt = ''
    use_filt_id = [4, 3, 1]
    for i in use_filt_id:  #['ACS','F115W', 'F150W','F277W', 'F444W']
        images.append(image_list[i])
        use_filt += filts[i] + '+'
    use_filt = use_filt[:-1]
    deltaPix = data_process.deltaPix*2
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    kpc_per_pixel = scale_relation 
    scale = 0.5 * kpc_per_pixel
    rgb_default = make_lupton_rgb(images[0], images[1], images[2], Q=7, stretch=0.3)
    fig, ax = plt.subplots()
    plt.title(name)
    plt.imshow(rgb_default, origin='lower')
    plt.text(1,80,use_filt,fontsize=20, color = 'white')
    scale_bar(ax, len(images[0]), dist=0.5/deltaPix, text='0.5"', text2 ='~{0:.2f}kpc'.format(scale), color = 'white')
    plt.show()

    pickle.dump(image_list, open('colorimage_bin2_{0}.pkl'.format(name[5:12]), 'wb'))
    sed_image = np.zeros_like(image_list[0])
    sed_2d_info = []
    for i in range(len(sed_image[0])):
        for j in range(len(sed_image[1])):
            mag_result = {}
            for k in range(len(filts)):
                filt = filts[k]
                flux = image_list[k][i,j]
                if flux>0:
                    mag = -2.5*np.log10(flux) + zp_dict[filt]
                    mag_result[filt] = mag
            sed_2d_info.append([i, j, mag_result])
    pickle.dump(sed_2d_info, open('2d_filts_mag_bin2_{0}.pkl'.format(name[5:12]), 'wb'))


