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
object_id,ra,dec='SDSS', 215.0359019, 53.00112292  # z 2.588

#%%
point_source_num = 1  #Number of AGN in the target.
ps_pix_center_list=[[0,0]] #Force the center as AGN position if there is point source to fit.
dr='dr4'
# rerun='s21a_dud'
rerun='s21a_wide'

#%%Some settings for the fitting
fitting_level='deep' #shallow, deep
bands = 'GRIZY'  #Band that will be download
lband = 'I' #The band fitting first and can also fit n and Re for other band.
fix_n, fix_re = True, True

#%%Mkdir folder and start downloading
if not os.path.exists(object_id):
    os.makedirs(object_id)
out_dir='./' + object_id
# print('Downloading data with PSF... ... ...')
# hsc_image.get_cutouts(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands,fov_arcsec=120)
# hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr=dr,rerun=rerun,filters=bands)

#For SDSS:
#hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr='dr3',rerun='s20a_wide',filters=bands)

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

#%%
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
    data_process.generate_target_materials(radius=None, #detect_tool='sep'
                                           )
    data_process.PSF_list = [PSF]
    data_process_list.append(data_process)

# % Determining the common settings for all bands, including cutout radius and apertures.
l_idx = [i for i in range(len(bands)) if bands[i] == lband][0]  #The first index to run
run_list = [i for i in range(len(bands))]
del(run_list[l_idx])
run_list = [l_idx] + run_list  #The list define the order to run

cut_radius = np.median([int(len(data_process_list[i].target_stamp)/2) for i in range(len(data_process_list))])
for i in range(len(bands)):    
    data_process_list[i].generate_target_materials(radius=cut_radius, create_mask = False, nsigma=2.8, # detect_tool='sep',
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
    data_process_list[i].checkout()
    
apertures = data_process_list[l_idx].apertures
from galight.tools.measure_tools import mask_obj   
for i in run_list:
    if i != l_idx:
        covers = mask_obj(data_process_list[i].target_stamp, apertures, if_plot=False, sum_mask = True)
        for j in range(len(data_process_list[i].apertures)):
            new_cover = mask_obj(data_process_list[i].target_stamp, [data_process_list[i].apertures[j]], if_plot=False, sum_mask = True)
            positions = [apertures[k].positions for k in range(len(apertures)) ]
            dists = [ np.sqrt((np.sum((data_process_list[i].apertures[j].positions - 
                              positions[k])**2))) for k in range(len(apertures)) ]
            cond1 = np.sum(covers - new_cover*covers) > np.sum(1-new_cover)/2 #If 1/2 of the area covered by the aperture is new)
            cond2 = np.sum( (1-new_cover) * covers ) > np.sum(1-new_cover )*4/5 #If 1/5 of the new aper is not overlap with pervious aper
            cond3 = np.min(dists) > 4  #If center offset above 4 pixel
            if cond1 and cond2 and cond3 :
                apertures.append(data_process_list[i].apertures[j])

from galight.tools.measure_tools import sort_apertures
apertures = sort_apertures(data_process_list[0].target_stamp, apertures)

#%%           
fit_sepc_l, fit_run_l = [None]*5, [None]*5
for i in run_list:  
    band = bands[i]
    if glob.glob(out_dir+'/*-{0}.pkl'.format(band)) == []:
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
                                          ps_pix_center_list=ps_pix_center_list,
                                          fix_n_list= fix_n_list, fix_Re_list=fix_Re_list)
        fit_sepc_l[i].plot_fitting_sets(out_dir+'/fitconfig-band-{0}.png'.format(band))
        fit_sepc_l[i].build_fitting_seq()
        fit_run_l[i] = FittingProcess(fit_sepc_l[i], savename = out_dir+'/result-band-{0}'.format(band), fitting_level=fitting_level)
        fit_run_l[i].run(algorithm_list = ['PSO'], setting_list=[None])
        for j in range(5):  #!!!Will be removed after Lenstrnomy debug.
            if np.sum(fit_run_l[i].image_host_list[0]) != 0:
                continue
            else:
                cut_radius = cut_radius-1
                data_process_list[i].generate_target_materials(radius=cut_radius)
                data_process_list[i].checkout()
                data_process_list[i].apertures = apertures #Pass apertures to the data
                fit_sepc_ = FittingSpecify(data_process_list[i])
                fit_sepc_.prepare_fitting_seq(point_source_num = point_source_num, 
                                                  supersampling_factor=3,
                                                  ps_pix_center_list=ps_pix_center_list)
                fit_sepc_.kwargs_params = fit_sepc_l[i].kwargs_params
                fit_sepc_.build_fitting_seq()
                fit_run_l[i] = FittingProcess(fit_sepc_, savename = out_dir+'/result-band-{0}'.format(band), 
                                              fitting_level=fitting_level)
                fit_run_l[i].run(algorithm_list = ['PSO'], setting_list=[None])
        
        if fit_run_l[i].image_ps_list != []:
            fit_run_l[i].plot_final_qso_fit(save_plot=True, target_ID= object_id +'-'+ band )
        else:
            fit_run_l[i].plot_final_galaxy_fit(save_plot=True, target_ID= object_id +'-'+ band )
            
        fit_run_l[i].cal_astrometry()
        fit_run_l[i].dump_result()
        # Script to remove the fitted target
        # fit_run_l[i].targets_subtraction(save_fitsfile=True, sub_gal_list=[1,2], sub_qso_list=[0])
