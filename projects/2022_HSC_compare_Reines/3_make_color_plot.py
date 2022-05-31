#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 10:35:54 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
IDs = ['J131310.12+051942.1','J233837.09-002810.4',
       'J145819.55+045451.7', 'J143450.63+033842.5', 'J141920.64+043623.3', 'J141630.82+013708.0',
       'J140018.41+050242.2', 'J134426.41+441620.0','J131305.81+012755.9','J121826.72-000750.1',
       'J104252.93+041441.1','J012159.81-010224.3','J084143.50+013149.8','J095540.47+050236.5']
# IDs = ['J084143.50+013149.8']
IDs.sort()

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
row_ = int(len(IDs)/7 -0.1) + 1
fig, (axs) = plt.subplots(row_, 7, figsize=(20, row_*3))

from astropy.visualization import make_lupton_rgb
bands = ['G', 'R', 'I', 'Z', 'Y']
# for ID in IDs:
for k,ID in enumerate(IDs):
    ra = ID[1:10]
    ra = ra[:2] + ':' + ra[2:4] + ':' +  ra[4:]
    dec = ID[10:]
    dec = dec[:3] + ':' + dec[3:5] + ':' +  dec[5:]
    images = []
    use_bands = []
    for band in bands:
        fitsname = glob.glob('./s21a/{0}/*cutout*{1}*.fits'.format(ID, band))
        if glob.glob('galight_results/'+ID+'_{0}*pkl'.format(band)) != []:
            fit_run = pickle.load(open(glob.glob('galight_results/'+ID+'_{0}*pkl'.format(band))[0],'rb')) 
            # print(ID, band,':')
            # fit_run.plot_final_qso_fit(target_ID = ID)
            obs_mag = fit_run.final_result_galaxy[0]['n_sersic']
            # print(obs_mag)
            images.append(fit_run.fitting_specify_class.kwargs_data['image_data'] - fit_run.image_ps_list[0])
            use_bands.append(band)
    r_i, g_i, b_i = images[2],images[1],images[0]
    min_s = np.min([len(r_i), len(g_i), len(b_i)])
    ct = int((len(r_i) - min_s)/2)
    if ct !=0:
        r_i = r_i[ct:-ct, ct:-ct]
    ct = int((len(g_i) - min_s)/2)
    if ct !=0:
        g_i = g_i[ct:-ct, ct:-ct]
    ct = int((len(b_i) - min_s)/2)
    if ct !=0:
        b_i = b_i[ct:-ct, ct:-ct]        
            
    rgb_default = make_lupton_rgb(r_i, g_i, b_i, stretch = 5)
    
    _i = int(k / len(axs.T))
    _j = int(k % len(axs.T))
    axs[_i][_j].imshow(rgb_default, origin='lower')
    sz = len(rgb_default)
    show_ID = ID[:4] + ID[9:14] 
    _bands = use_bands[0]+use_bands[1]+use_bands[2]
    plttext = axs[_i][_j].text(sz/30,sz/20*17.5,show_ID,color='white',fontsize=18)
    # plttext = axs[_i][_j].text(sz/30,sz/20*15,"z={0:.3f}".format(read_z(ID)),fontsize=18,color='white')
    plttext = axs[_i][_j].text(sz/30,sz/20*1, _bands,color='white',fontsize=15)
    # scale_bar(axs[_i][_j], sz, dist=1/0.168, text='1"', color='white', fontsize=18)
    # if ID == '011227.87-003151.6':
    #     plttext = axs[_i][_j].text(sz/30,sz/20*13,'not selected',color='white',fontsize=18)        
    # coordinate_arrows(axs[_i][_j], sz, arrow_size=0.03, color = 'white')
    axs[_i][_j].axes.xaxis.set_visible(False)
    axs[_i][_j].axes.yaxis.set_visible(False)
axs[1][5].axis("off")
plt.subplots_adjust(wspace=-0.02, hspace=0.04)
# plt.savefig('show_material/color_plot.pdf')
plt.show()
