#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:15:00 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from def_functions import target_in_fits, RA_Dec_in_fit
import pickle
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.cutout_tools import psf_clean

datafolder = '/Users/Dartoon/Downloads/CEERS_JWST_data'
save_f = 'output/'
create_mask = False
# target_id = 'aegis_509'
# RA, Dec = 214.89874, 52.85251

# target_id = 'aegis_511'
# RA, Dec = 214.89561, 52.856516

filters =  ['f115w', 'f200w', 'f277w', 'f410m', 'f444w']
for filt in filters:
    filenames = glob.glob(datafolder+'/bkg_removed/*'+filt+'*.fits')
    #%% Plot the PSF library
    # PSF_lib_files = glob.glob('material/'+filt+'*_PSF_info.pkl')
    # PSF_list = []
    # for i in range(len(PSF_lib_files)):
    #     PSF_list_, _, _= pickle.load(open(PSF_lib_files[i],'rb'))
    #     PSF_list = PSF_list+PSF_list_
    # for i, psf in enumerate(PSF_list[5:8]):
    #     print('Showing PSF', i, ":")
    #     psf = psf_clean(psf, if_plot=True, if_print_fluxratio=True, 
    #                     print_string='Clean PSF{0}\n'.format(i))
    #     from galight.tools.astro_tools import plt_fits
    #     plt_fits(psf)
    
    #%%
    
    print("Fit", target_id)
    file = RA_Dec_in_fit(all_files=filenames, RA=RA, Dec=Dec)
    fitsFile = pyfits.open(file)
    
    fov_image = fitsFile[1].data # check the back grounp
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    flux_mjsr = header['PHOTMJSR']
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
    wht = fitsFile[4].data # The WHT map
    exp = fitsFile[0].header['EFFEXPTM']
    gain_value = 1.8
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    data_process = DataProcess(fov_image = fov_image, target_pos = [RA, Dec], pos_type = 'wcs', header = header,
                              rm_bkglight = False, if_plot=False, zp = zp, exptime= exp_map )
    data_process.generate_target_materials(radius=None, create_mask = create_mask, nsigma=2.8, 
                                           cut_kernel = 'center_bright', if_select_obj=False,
                                          exp_sz= 1.2, npixels = 80, contrast = 0.01, if_plot=False, )
    
    if np.sum(data_process.target_stamp ==0) >= 1:
        data_process.target_mask = data_process.target_stamp != 0
        data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
    
    PSF_lib_files = glob.glob('material/'+filt+'*_PSF_info.pkl')
    PSF_list = []
    PSF_RA_DEC_list = []
    PSF_from_list = []
    for i in range(len(PSF_lib_files)):
        PSF_list_, PSF_pos_list_, PSF_RA_DEC_list_ = pickle.load(open(PSF_lib_files[i],'rb'))
        PSF_list = PSF_list+PSF_list_
        PSF_RA_DEC_list = PSF_RA_DEC_list + PSF_RA_DEC_list_
        PSF_from_list = PSF_from_list+ [PSF_lib_files[i]] * len(PSF_list_)
    
    tbl = data_process.tbl
    ap_kron_fluxes = [float(tbl[tbl['label']==j]['kron_flux']) for j in range(len(tbl))] 
    segm_kron_dens = [tbl[np.where(tbl['label']==j)[0][0]]['segment_flux']/tbl[np.where(tbl['label']==j)[0][0]]['area'].value for j in range(len(tbl))] 
    for i in range(len(ap_kron_fluxes)-1,0, -1):
        if ap_kron_fluxes[i]/ap_kron_fluxes[0] < 0.03 or segm_kron_dens[i]/segm_kron_dens[0]<0.08:
            del data_process.apertures[i]
    fit_run_list = []
    use_PSF_RA_DEC_list = []
    psf_target_dis = np.sqrt(np.sum((np.array(PSF_RA_DEC_list) - np.array([RA, Dec]))**2, axis=1))*3600/pixscale
    PSF_list = [PSF_list[i] for i in range(len(PSF_list)) if psf_target_dis[i]>10]
    PSF_RA_DEC_list = [PSF_RA_DEC_list[i] for i in range(len(PSF_RA_DEC_list)) if psf_target_dis[i]>10]
    PSF_from_list = [PSF_from_list[i] for i in range(len(PSF_from_list)) if psf_target_dis[i]>10]
    data_process.PSF_list = []
    data_process.psf_id_for_fitting = -1
    for i in range(len(PSF_list)):
        # if i < len(PSF_list):
        use_PSF_RA_DEC_list.append(PSF_RA_DEC_list[i])
        psf = PSF_list[i]
        psf = psf_clean(psf)  
        data_process.PSF_list.append(psf)
        # if i == len(PSF_list):
        #     use_PSF_RA_DEC_list.append('stacked PSF')
        #     data_process.stack_PSF()
        fit_sepc = FittingSpecify(data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3,
                                     ps_pix_center_list = [[0,0]]) #, fix_n_list= [[0,4],[1,1]])
        fit_sepc.build_fitting_seq()
        if i == 0:
            fit_sepc.plot_fitting_sets()
        fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep'])
        fit_run.run(algorithm_list = ['PSO','PSO'])
        fit_run.plot_final_qso_fit(target_ID =target_id)
        fit_run_list.append(fit_run)
    pickle.dump([fit_run_list, use_PSF_RA_DEC_list], open(save_f+filt+'_fit_result_PsfLib'+'_'+target_id+'.pkl', 'wb'))
