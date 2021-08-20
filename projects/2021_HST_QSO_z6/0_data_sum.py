#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 22:48:24 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from galight.data_process import DataProcess
import galight.tools.astro_tools as astro_tools
from astropy.wcs import WCS
from galight.tools.cutout_tools import cutout
from matplotlib.colors import LogNorm

ids = glob.glob('HST_download_sh/download_sh/*')

ids = [ids[i].split('/')[-1] for i in range(len(ids))]
ids.sort()

ra_dec_file = 'HST_download_sh/ALMAz6qso_radec_for_HST_archive.list'
f = open(ra_dec_file,"r")
string = f.read()
lines = string.split('\n')   # Split in to \n
ct = 0
imgs=[]
info=[]
for ID in ids[:-1]:
# for ID in ['J0439+1634']:
    idx = [i for i in range(len(lines)) if ID in lines[i]][0]
    # print(ID, idx)
    RA, Dec = lines[idx].split(' ')[:2]
    RA, Dec = np.float(RA), np.float(Dec)
    # print(ID, RA, Dec)
    fits_files = glob.glob('HST_download_sh/download_sh/{0}/*/HST/*/*drz.fits'.format(ID))
    fits_names = [fits_files[i].split('/')[-1] for i in range(len(fits_files))]
    fits_names = list(dict.fromkeys(fits_names))
    # print(len(fits_names), fits_names)
    ct += len(fits_names)
    
    for i in range(len(fits_names)):
        fitsfile =  glob.glob('HST_download_sh/download_sh/{0}/*/HST/*/{1}'.format(ID, fits_names[i]))
        fitsFile = pyfits.open(fitsfile[0])
            # print(fitsFile[0].header['PRIMESI'])
        if fitsFile[0].header['PRIMESI'] == 'WFC3':
            fov_image = fitsFile[1].data # check the back grounp
            header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            # wht = fitsFile[2].data # The WHT map
            # exp =  astro_tools.read_fits_exp(fitsFile[0].header)  #Read the exposure time 
            # pixel_scale = astro_tools.read_pixel_scale(fitsFile[1].header)  #Read pixel scale
            # mean_wht = exp * (pixel_scale/0.135)**2
            # exp_map = exp * wht/mean_wht
            # data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
            #                       rm_bkglight = True, exptime = exp_map, if_plot=False, zp = 27.0)
            # try:
            #     data_process.generate_target_materials(radius=25, create_mask = False, if_plot=True)
            # except:
            #     print()
            wcs = WCS(header)
            target_pos = wcs.all_world2pix([[RA, Dec]], 1)[0]
            target_pos = np.int0(target_pos)
            if target_pos[0]>0 and target_pos[1]>0:
                target_stamp = cutout(image = fov_image, center = target_pos, radius=45)
                target_stamp[np.isnan(target_stamp)] = 0
                plt.imshow(target_stamp-target_stamp.min(), norm = LogNorm(), cmap='gist_heat',vmax = target_stamp.max(),
                           vmin = 1.e-4, origin='lower')
                plt.colorbar()
                plt.close()
                try:
                    filt = fitsFile[0].header['FILTER2']
                except:
                    filt = fitsFile[0].header['FILTER']
                imgs.append(target_stamp-target_stamp.min())
                info.append([ID, filt, fitsfile[0]])
del imgs[4]
del info[4]

#%%
_row = 5
fig, (axs) = plt.subplots(_row, 4, figsize=(11, 3 + 3 * (_row-1)))
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
for i in range(len(imgs)):
    _i = int(i / 4)
    _j = int(i % 4)
    axs[_i][_j].imshow(imgs[i], origin='lower', norm=LogNorm())
    frame_size = len(imgs[i])
    info_ = info[i][0]
    plttext = axs[_i][_j].text(frame_size*0.05, frame_size*0.87, "{0}".format(info_),
             fontsize=17, weight='bold', color='black')
    plttext.set_bbox(dict(facecolor='white', alpha=0.5))
    plttext = axs[_i][_j].text(frame_size*0.05, frame_size*0.13, "{0}".format(info[i][1]),
             fontsize=17, weight='bold', color='black')    
    plttext.set_bbox(dict(facecolor='white', alpha=0.5))
    axs[_i][_j].axes.xaxis.set_visible(False)
    axs[_i][_j].axes.yaxis.set_visible(False)
for i in range( 4 - len(imgs)%4 ):
    axs[-1][-(i+1)].axis('off')
plt.show()

#%%
write_file = open('fit_files_info.txt','w') 
for i in range(len(info)):
    write_file.write("{0} {1} {2}\n".format(info[i][0], info[i][2], info[i][1]) )
write_file.close()