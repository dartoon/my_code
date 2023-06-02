#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 23:19:06 2022

@author: Dartoon
"""


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

from target_info import target_info
from galight.tools.astro_tools import plt_fits, plt_many_fits
#%%
plt.rcParams["font.family"] = "sans-serif"

data_type = 'all' 
filt = 'F356W'
filt = 'F150W'
file_NO = 0

idx = 1
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']
# folder = '/Users/Dartoon/Downloads/z6JWSTNIRcam/NIRCam_{1}_stage3_{0}/bkg_removed'.format(data_type, target_id[:5])
if filt == 'F356W':
    folder = '../NIRCam_data/Nov14/bkg_removed/' 
if filt == 'F150W':
    # folder = '../NIRCam_data/Jan14/bkg_removed/' 
    folder = '../NIRCam_data/Nov14/bkg_removed/' 
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


#%%Demonstrate the local envs.
import copy, matplotlib
from matplotlib.colors import LogNorm
# cmap = ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
#             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
#             'hot', 'afmhot', 'gist_heat', 'copper']

# cmap = ['winter','summer','afmhot','spring', 'autumn', 'gist_heat','hot' ]

if filt == 'F150W' :
    cmap = 'inferno'
else:
    cmap = 'gist_heat'
# cmap = 'cubehelix'
    # cmap = 'gist_stern'
# else:
# cmap = 'gnuplot'
# cmap = 'gnuplot2'
# cmap = 'magma'
# cmap = 'CMRmap'
my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
my_cmap.set_bad('black')

def coordinate_arrows(ax, d, header, color='white', arrow_size=0.02):
    wcs = WCS(header)
    d0 = d / 12.
    deltaPix = 1
    ra0, dec0 = (d - 2*d0) / deltaPix, d0 / deltaPix
    res_ra0, res_dec0 = wcs.all_pix2world([(ra0, dec0)], 1)[0]
    xx_, yy_ = ra0, dec0
    
    xx_ra, yy_ra = wcs.all_world2pix([[res_ra0+2/3600, res_dec0]], 1)[0]
    xx_dec, yy_dec = wcs.all_world2pix([[res_ra0, res_dec0+2/3600]], 1)[0]
    xx_ra_t, yy_ra_t = wcs.all_world2pix([[res_ra0+3.2/3600, res_dec0]], 1)[0]
    xx_dec_t, yy_dec_t = wcs.all_world2pix([[res_ra0, res_dec0+3.5/3600]], 1)[0]

    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_ra - xx_) * deltaPix, (yy_ra - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc=color, ec=color, linewidth=1.2)
    ax.text(xx_ra_t * deltaPix, yy_ra_t * deltaPix, "E", color=color, fontsize=22, ha='center')
    ax.arrow(xx_ * deltaPix, yy_ * deltaPix, (xx_dec - xx_) * deltaPix, (yy_dec - yy_) * deltaPix,
             head_width=arrow_size * d, head_length=arrow_size * d, fc
             =color, ec=color, linewidth=1.2)
    ax.text(xx_dec_t * deltaPix, yy_dec_t * deltaPix, "N", color=color, fontsize=22, ha='center')
    

# def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=25):
#     if flipped:
#         p0 = d - d / 15. - dist
#         p1 = d / 15.
#         ax.plot([p0, p0 + dist], [p1, p1], linewidth=3, color=color)
#         ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
#     else:
#         p0 = d / 15.
#         ax.plot([p0, p0 + dist], [p0, p0], linewidth=3, color=color)
#         ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')

cut_kernel = None #After pos correct then, do the nearest_obj_center
# filts = ['F356W', 'F150W']
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
    expsize = 1
    exppix = 2
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], #The final cut center is dete
                               pos_type = 'wcs', header = header,rm_bkglight = False, 
                               if_plot=False, zp = zp, exptime= exp_map, 
                               fov_noise_map = None)
data_process.fov_noise_map  = data_process.fov_image
data_process.generate_target_materials(radius=400 * expsize, create_mask = False, nsigma=1.5, 
                                        cut_kernel = None, skip = True)
# plt_fits(data_process.target_stamp)
data_process.apertures = []
# data_process.plot_aperture()
fig, ax = plt.subplots(figsize=(12,12))  #!!!
vmin = 1.e-3
if filt == 'F150W':
    vmax = data_process.target_stamp.max()/3 
else:
    vmax = data_process.target_stamp.max()
plt.imshow(data_process.target_stamp, origin='lower', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax))#, vmin=vmin, vmax=vmax)
frame_size = len(data_process.target_stamp)

d = frame_size
color = 'white'
text='5"'
dist=5/data_process.deltaPix
p0 = d / 15.
ax.plot([p0 + d / 15., p0 + dist + d / 15.], [p0, p0], linewidth=3, color=color)
ax.text(p0 + dist / 2.  + d / 15., p0 + 0.02 * d , text, fontsize=25, color=color, ha='center')

angle = 0 / 180 * np.pi
coordinate_arrows(ax, frame_size, header=header, arrow_size=0.03)
# ax.set_ylim(0,800)
if filt == 'F150W':
    circle1 = plt.Circle((38, 95),32, color='white', fill=False, linewidth=2)
    circle2 = plt.Circle((400., 400),20, color='white', fill=False, linewidth=1, alpha = 1)
    # circle3 = plt.Circle((445., 365),20, color='white', fill=False, linewidth=2, alpha = 0.7)
if filt == 'F356W':
    circle1 = plt.Circle((50, 100),40, color='white', fill=False, linewidth=2)
    circle2 = plt.Circle((400., 400),35, color='white', fill=False, linewidth=1, alpha = 1)
    # circle3 = plt.Circle((445., 365),25, color='white', fill=False, linewidth=2, alpha = 0.7)
ax.text(160., 240, 'PSF-star', fontsize=24, color=color, ha='center')
ax.add_patch(circle1)
ax.set_xticks([])
ax.set_yticks([])
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

if filt == 'F356W':
    ax.arrow(148, 230, -60,-80,
             head_width=20, head_length=15, fc=color, ec=color, linewidth=1.2)
    ax.set_xlim(0,1220)
    axins = zoomed_inset_axes(ax, 2.3, loc = 'center right', bbox_to_anchor=(0,0,800,1200) )
    loc1, loc2 = 2,4
    pso1, pso2 = 350., 430
else:
    ax.arrow(148, 230, -60,-80,
             head_width=20, head_length=15, fc=color, ec=color, linewidth=1.2)
    ax.set_xlim(-420,-420+1220)
    axins = zoomed_inset_axes(ax, 2.3, loc = 'center left', bbox_to_anchor=(72,0,0,520) )
    loc1, loc2 = 2,4
    pso1, pso2 = 360., 420
    
axins.text(pso1, pso2, 'quasar', fontsize=30, color=color, ha='center')
axins.imshow(data_process.target_stamp, origin='lower', cmap=my_cmap, norm=LogNorm(vmin=vmin, vmax=vmax))#, vmin=vmin, vmax=vmax)
axins.set_xlim(300, 500)
axins.set_ylim(300, 500)

text='1"'
dist=1/data_process.deltaPix
p0 = d / 15.
axins.plot([320, 320 + dist], [320, 320], linewidth=3, color=color)
axins.text(320 + dist / 2., 320 + 0.01 * d, text, fontsize=25, color=color, ha='center')


axins.add_patch(circle2)
# axins.text(478., 380, 'obj1', fontsize=20, color=color, ha='center')
# axins.add_patch(circle3)


ax.text(d*0.2, d*0.85, filt, fontsize=35, color=color, ha='center')

if filt =='F356W':
    ax.text(d*0.475, d*1.05, target_id, fontsize=45, color='black', ha='center')
# if filt =='F150W':
#     ax.text(d*0.52, d*1.05, target_id, fontsize=45, color='black', ha='center')


axins.set_xticks([])
axins.set_yticks([])

mark_inset(ax, axins, loc1=loc1, loc2=loc2, fc="none", ec="0.6", linewidth=2.2)
ax.set(frame_on=False)  # New

plt.savefig('figures/field_overview{0}.pdf'.format(filt))
plt.show()