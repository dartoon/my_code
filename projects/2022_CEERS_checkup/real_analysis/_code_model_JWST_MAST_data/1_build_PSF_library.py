#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 17:43:09 2022

@author: Dartoon

To generate a PSF library for each filter. 
Filters are 
 ['f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w']
 
 For 200W and F444W, three t022 included.
 For t021, three idential fits images are included.
 
 We thus adopt only *o00* visit, each filter have 4 fits then. 
 
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob 

filters =  ['f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w']
folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'

#Check the files:
# for i in range(len(filters)):
#     filt = filters[i]
#     folder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
#     all_files= glob.glob(folder+'/*clear*/*o00*{0}_i2d.fits'.format(filt))  #For NIRCam
#     all_files.sort()
#     for j in range(len(all_files)):
#         print(all_files[j])
#     print('\n')
    
#%%

filt = filters[6]
filter_files= glob.glob(folder+'/*clear*/*o00*{0}_i2d.fits'.format(filt))  #For NIRCam
filter_files.sort()
    
# filename = filter_files[2]
for filename in filter_files:
    print("Select for", filename.split('/')[-1][8:])
    # Grab the JWST provided ERR map:
    
    fitsFile = pyfits.open(filename)
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    from galight.data_process import DataProcess
    from galight.tools.astro_tools import read_pixel_scale
        
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    wht = fitsFile[4].data # The WHT map
    # exp =  header['XPOSURE']  #Read the exposure time 
    exp = fitsFile[0].header['EFFEXPTM']
    gain_value = 2
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
        
    from galight.data_process import DataProcess
    data_process = DataProcess(fov_image = fov_image, target_pos = [-100,-100], pos_type = 'pixel', header = header,
                              rm_bkglight = True, exptime = exp_map, if_plot=False, zp = zp)
    # data_process.generate_target_materials(radius=30, create_mask = False, nsigma=2.8, if_select_obj=False,
    #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
    #PSF works.
    data_process.find_PSF(radius = 45, user_option = True, if_filter=True, nearyby_obj_filter=False, FWHM_sort=True)
    data_process.plot_overview()
    # from galight.tools.cutout_tools import psf_clean
    # PSFs = data_process.PSF_list
    # PSFs = [psf_clean(PSFs[i], print_string='clean PSF '+str(i), 
    #                   ratio_to_replace=0.01, if_plot=True) for i in range(len(PSFs))]
    # data_process.PSF_list = PSFs
    # data_process.stack_PSF()
    
    PSF_list = data_process.PSF_list
    PSF_pos_list = data_process.PSF_pos_list
    from astropy.wcs import WCS
    wcs = WCS(header)
    # sky = [wcs.pixel_to_world(PSF_pos_list[i]) for i in range(data_process.PSF_pos_list)]
    PSF_RA_DEC_list = []
    for i in range(len(PSF_pos_list)):
        RA_DEC = wcs.all_pix2world(PSF_pos_list[i][0], PSF_pos_list[i][1],0) 
        PSF_RA_DEC_list.append( [float(RA_DEC[0]), float(RA_DEC[1])] )
        
    save_name = filt + filename.split('_nircam_clear')[-2][-5:]
    import pickle
    fitsFile[1].data = data_process.fov_image
    fitsFile.writeto(folder+'/bkg_removed/'+save_name+'.fits')
    
    noselect = False
    
    if noselect == True:
        pickle.dump([[], [], []], open('material/'+save_name+'_PSF_info.pkl', 'wb'))
    else:
        pickle.dump([PSF_list, PSF_pos_list, PSF_RA_DEC_list], open('material/'+save_name+'_PSF_info.pkl', 'wb'))

# result = pickle.load(open('material/'+save_name+'_PSF_info.pkl','rb'))
    
