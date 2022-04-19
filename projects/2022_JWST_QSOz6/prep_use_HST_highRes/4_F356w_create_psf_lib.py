#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 14:05:38 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#%%Load image and define a empty place

filt = 'f356w'
folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)
im = pyfits.open(folder+file)
data = im[1].data
header = im[1].header

# print('For flux value in unit of MJy/sr.') #https://en.wikipedia.org/wiki/AB_magnitude
value_unit = header['BUNIT']
# flux(Mjy/sr) * 2.350443 * 10**(-5) *0.03**2   #Flux to Jy  https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/18/
print("Data unit:", value_unit)
import copy

data_sb1 = copy.deepcopy(data[0:4400,0:4400])
data_sb2 = copy.deepcopy(data[:,5800:])
psfs,FWHMs = [],[]
for data_sb in [data_sb1,data_sb2]:
    zp = -2.5*np.log10(2.350443 * 10**(-5) *0.03**2/3631)
    header0 = im[0].header
    img_filter = header0['FILTER']
    img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
    exptime = header0['TEXPTIME'] #The assumed exp time.
    from galight.tools.astro_tools import plt_fits
    plt_fits(data_sb)
    # Obtain PSF stars:
    from galight.data_process import DataProcess
    # zp = 31.4 - 2.5*np.log10(header['PHOTMJSR'])  #Calculate the correspondingly zp as DN/S #This is wrong! See 3_...
    data_process_ = DataProcess(fov_image = data_sb, target_pos = [1170., 940.], pos_type = 'pixel', header = header,
                                rm_bkglight = True, exptime = np.ones_like(data_sb)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
    # data_process_.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
    #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
    data_process_.find_PSF(radius = 60, user_option = True, psf_edge=10, select_all=True)
    # data_process_.plot_overview(label = 'Example', target_label = None)
    FWHM_list = np.array(data_process_.PSF_FWHM_list)
    POS_list = np.array(data_process_.PSF_pos_list)
    psf_list = np.array(data_process_.PSF_list)
    fwhm_bools = [FWHM_list<4.5][0]
    psf_list = psf_list[fwhm_bools]
    POS_list = POS_list[fwhm_bools]
    FWHM_list = FWHM_list[fwhm_bools]
    from galight.tools.measure_tools import detect_obj
    from photutils.segmentation import SourceCatalog
    near_bools = []
    for i in range(len(psf_list)):
        _, segm_map = detect_obj(psf_list[i],if_plot=False,segm_map=True, nsigma=5)
        cat = SourceCatalog(psf_list[i], segm_map)
        tbl = cat.to_table()
        idx= segm_map.data[int(len(psf_list[i])/2), int(len(psf_list[i])/2)]
        kron_fluxes = [tbl[i]['kron_flux'] for i in range(len(tbl))]
        fluxes_ratios = np.array(kron_fluxes)/kron_fluxes[idx-1]
        if len(fluxes_ratios) == 1:
            near_bools.append(True)
        elif np.max(fluxes_ratios[fluxes_ratios!=1]) > 0.05:
            near_bools.append(False)
        else:
            near_bools.append(True)
    near_bools = np.array(near_bools)
    psf_list[near_bools]
    from galight.tools.cutout_tools import psf_clean
    psfs.append([ psf_clean(psf_list[near_bools][i])[10:-10,10:-10] for i in range(len(psf_list[near_bools]))])
    FWHMs.append( [FWHM_list[near_bools][i] for i in range(len(FWHM_list[near_bools]))])

#%%Check and remove some 'bad' ones
idx = 0
for i in range(len(psfs[idx])):
    plt_fits(psfs[idx][i])
    print(i, FWHMs[idx][i])

del psfs[1][8]
del psfs[1][1]
del psfs[0][11]

del FWHMs[1][8]
del FWHMs[1][1]
del FWHMs[0][11]

#%%
import pickle
pickle.dump([psfs,FWHMs], open(filt+'_psfs.pkl', 'wb'))   

#%%
psfs,FWHMs = pickle.load(open('f356w_psfs.pkl','rb'))