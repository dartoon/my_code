#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:27:48 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import glob

#%%Cutout using online cutout tool:
from galight.hsc_utils import hsc_image, hsc_psf
import os
# object_id,ra,dec='ID1', 35.2735048,  -4.6837584
# object_id,ra,dec='150200.03+444541.9', 225.5694427, 2.955502
object_id,ra,dec= '000017.88+002612.6', 0.07452999800443649, 0.4368380010128021
# object_id,ra,dec= '000017.88+002612.6_', '00:00:17.88', '00:26:12.6'

#%%
point_source_num = 0  #Number of AGN in the target.
dr='dr4'
# rerun='s21a_dud'
rerun='s21a_wide'

#%%Some settings for the fitting
fitting_level='shallow' #shallow, deep
bands = 'GRIZY'  #Band that will be download
lband = 'I' #The band fitting first and can also fit n and Re for other band.
fix_n, fix_re = True, True

#%%Mkdir folder and start downloading
if not os.path.exists(object_id):
    os.makedirs(object_id)
out_dir='./' + object_id
print('Downloading data with PSF... ... ...')
hsc_image.get_cutouts(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=120)
hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands)

#%%use galight to analyze:
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
print('Fitting using GaLight... ... ...')
if isinstance(ra, str) or isinstance(dec, str):
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    pos = SkyCoord('{0} {1}'.format(ra, dec), unit=(u.hourangle, u.deg))
    ra, dec = pos.ra.degree, pos.dec.degree

data_process_list = []
for band in bands:
    fitsFile = pyfits.open(glob.glob(object_id+'/*-cutout-HSC-{0}*.fits'.format(band))[0])
    file_header0 = fitsFile[0].header
    try:
        FLUXMAG0 = file_header0['FLUXMAG0']
        zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
    except:
        zp = 27.0
    PSF_file = glob.glob(object_id+'/*-psf*HSC-{0}*.fits'.format(band))[0]
    PSF = pyfits.getdata(PSF_file)
    data_process = DataProcess(fov_image = fitsFile[1].data, fov_noise_map = fitsFile[3].data ** 0.5, target_pos = [ra, dec],
                                pos_type = 'wcs', header = fitsFile[1].header,
                                rm_bkglight = True, if_plot=False, zp = zp)
    data_process.generate_target_materials(radius=None)
    data_process.PSF_list = [PSF]
    data_process_list.append(data_process)

#%% Determining the common settings for all bands, including cutout radius and apertures.
l_idx = [i for i in range(len(bands)) if bands[i] == lband][0]  #The first index to run
run_list = [i for i in range(len(bands))]
del(run_list[l_idx])
run_list = [l_idx] + run_list  #The list define the order to run
 
cut_radius = np.median([int(len(data_process_list[i].target_stamp)/2) for i in range(len(data_process_list))])
for i in range(len(bands)):    
    data_process_list[i].generate_target_materials(radius=cut_radius, create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
apertures = data_process_list[l_idx].apertures
for i in run_list:
    if i != l_idx:
        for j in range(len(data_process_list[i].apertures)):
            positions = [apertures[k].positions for k in range(len(apertures)) ]
            dists = [ np.sqrt((np.sum(data_process_list[i].apertures[j].positions - 
                              positions[k]))**2) for k in range(len(apertures)) ]
            if np.min(dists) > 2:
                idx = np.where(dists == np.min(dists))[0][0]
                apertures.append(data_process_list[i].apertures[idx])
fit_sepc_l, fit_run_l = [None]*5, [None]*5
for i in run_list:  
    band = bands[i]
    print("Staring fitting band-"+band+"... ... ...")
    data_process_list[i].apertures = apertures #Pass apertures to the data
    fit_sepc_l[i] = FittingSpecify(data_process_list[i])
    fix_n_list, fix_Re_list = None, None
    if i != l_idx:
        if fix_n == True:
            fix_n_list = [[0,fit_run_l[l_idx].final_result_galaxy[0]['n_sersic'] ]]
        if fix_re == True:
            fix_Re_list = [[0,fit_run_l[l_idx].final_result_galaxy[0]['R_sersic'] ]]
    fit_sepc_l[i].prepare_fitting_seq(point_source_num = point_source_num, supersampling_factor=3, 
                                      fix_n_list= fix_n_list, fix_Re_list=fix_Re_list)
    fit_sepc_l[i].plot_fitting_sets(object_id+'/fitconfig-band-{0}.png'.format(band))
    fit_sepc_l[i].build_fitting_seq()
    fit_run_l[i] = FittingProcess(fit_sepc_l[i], savename = object_id+'/result-band-{0}'.format(band), fitting_level=fitting_level)
    fit_run_l[i].run(algorithm_list = ['PSO'], setting_list=[None])
    # fit_run.plot_all(target_ID = object_id)
    if fit_run_l[i].image_ps_list != []:
        fit_run_l[i].plot_final_qso_fit(save_plot=True, target_ID= object_id +'-'+ band )
    else:
        fit_run_l[i].plot_final_galaxy_fit(save_plot=True, target_ID= object_id +'-'+ band )
        
    fit_run_l[i].cal_astrometry()
    fit_run_l[i].dump_result()
    #Script to remove the fitted target
    # fit_run_l[i].targets_subtraction(save_fitsfile=True, sub_gal_list=[1,2], sub_qso_list=[0])
