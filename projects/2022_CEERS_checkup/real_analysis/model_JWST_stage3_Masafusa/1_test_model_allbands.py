#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:05:05 2022

@author: Dartoon

come form 1_build_PSF_library.py
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import sys
sys.path.insert(0,'..')
from def_functions import RA_Dec_in_fit
import pickle
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.tools.cutout_tools import psf_clean
from galight.tools.astro_tools import plt_fits
from galight.tools.measure_tools import measure_bkg
from astropy.wcs import WCS
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.cutout_tools import common_data_class_aperture
import warnings
warnings.filterwarnings("ignore")

#%%
folder = '/Volumes/Seagate_Expansion_Drive/data_backup/CEERS_data/CEERS_JWST_Masafusa'
filenames = glob.glob(folder+'/bkg_removed/'+'*.fits')
# filters = ['F356W', 'F410M', 'F444W']
# filt = filters[0]

f = open("target_info.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
result_folder = 'fit_result/'

lines = lines[1:]
#%%
# special = [31, 'aegis_471']

# fit_id = 11
create_mask = False
for idx, line in enumerate(lines[4:5]):
    target_id, RA, Dec, spec_z, photo_z = line.split(' ')
    RA, Dec, spec_z, photo_z = float(RA), float(Dec), float(spec_z), float(photo_z)
    if target_id == 'aegis_630':
        RA, Dec = 214.9421549996112, 52.946445580969026
    files = RA_Dec_in_fit(all_files=filenames, RA=float(RA), Dec=float(Dec))
    filters = [files[i].split('NIRCam')[1][2:7] for i in range(len(files))]
    files = [x for _,x in sorted(zip(filters,files))]
    filters.sort()
    # print("Fit", target_id, 'filt include', filters)
    file_list = []
    data_process_list_s = []
    data_process_list_l = []
    cut_kernel = 'nearest_obj_center'
    
    cut_RA, cut_Dec = RA, Dec
    for file in files[:]:
        print("Run next filter")
        file_list.append(file)  #recording which file is used.
        filt = file.split('NIRCam')[1][2:7]
        fitsFile = pyfits.open(file)
        fov_image = fitsFile[1].data # check the back grounp
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        flux_mjsr = header['PHOTMJSR']
        pixscale = read_pixel_scale(header)
        zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux
        wht = fitsFile[4].data # The WHT map
        fov_noise_map = fitsFile[2].data 
        
        # exp =  header['XPOSURE']  #Read the exposure time 
        exp = fitsFile[0].header['EFFEXPTM']
        if fitsFile[0].header['CHANNEL'] == 'LONG':
            gain_value = 2
            expsize = 1
        else:
            gain_value = 1.8
            expsize = 1.5
        exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
        data_process = DataProcess(fov_image = fov_image, target_pos = [cut_RA, cut_Dec], 
                                   pos_type = 'wcs', header = header,rm_bkglight = False, 
                                   if_plot=False, zp = zp, exptime= exp_map, 
                                   fov_noise_map = None)
        #estimate local bkg and remove:
        
        data_process.generate_target_materials(radius=250 * expsize, create_mask = create_mask, nsigma=2.8, 
                                                cut_kernel = cut_kernel, if_select_obj=False, skip = True,
                                                if_plot=True)
        
        bkg_std = data_process.bkg_std
        wcs = WCS(header)
        cut_RA, cut_Dec = wcs.all_pix2world([data_process.target_pos], 1)[0]
        cut_kernel = None
        
        bkglight = measure_bkg(data_process.target_stamp, if_plot=False) # Remove bkg light
        data_process.generate_target_materials(radius=70 * expsize, create_mask = create_mask, nsigma=2.8, 
                                                cut_kernel = None, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 80 * expsize, if_plot=True, bkg_std= bkg_std)
        
        data_process.noise_map[data_process.noise_map == 0] = data_process.noise_map.max()
        if np.sum(data_process.target_stamp ==0) >20:
            data_process.target_mask = data_process.target_stamp != 0
            data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
        
        ct = int((len(bkglight) - len(data_process.target_stamp ))/2)
        data_process.target_stamp = data_process.target_stamp - bkglight[ct:-ct, ct:-ct]
        # data_process.clean_aperture()
        del data_process.fov_image
        del data_process.exptime
        print('idx', idx, target_id, filt, 'apertures', len(data_process.apertures) )
        data_process.plot_aperture()
        # plot_data_apertures_point(data_process.target_stamp * data_process.target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
        #                           data_process.apertures)
        if fitsFile[0].header['CHANNEL'] == 'LONG':
            data_process_list_l.append(data_process)
        if fitsFile[0].header['CHANNEL'] == 'SHORT':
            data_process_list_s.append(data_process)
    print("Fit", target_id, 'filt include', filters)
    com_aper_s = common_data_class_aperture(data_process_list_s, l_idx=0, return_idx=0)
    com_aper_l = common_data_class_aperture(data_process_list_l, l_idx=0, return_idx=0)
    print("Common aperture:")
    plot_data_apertures_point(data_process_list_s[0].target_stamp * data_process_list_s[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                              com_aper_s)
    plot_data_apertures_point(data_process_list_l[0].target_stamp * data_process_list_l[0].target_mask, # + (self.kwargs_likelihood['image_likelihood_mask_list'][0]==0)*1.e6 , 
                              com_aper_l)
    hold = input('Hold ... OK?\n')

#%%

#     if np.sum(data_process.target_stamp ==0) > 20:
#         data_process.target_mask = data_process.target_stamp != 0
#         data_process.noise_map = np.nan_to_num(data_process.noise_map, nan=1000)
    
#     PSF_lib_files = glob.glob('material/*'+filt[:-1]+'*_PSF_Library.pkl')[0]
#     PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))
    
#     tbl = data_process.tbl
#     ap_kron_fluxes = [float(tbl[tbl['label']==j]['kron_flux']) for j in range(len(tbl))] 
#     segm_kron_dens = [tbl[np.where(tbl['label']==j)[0][0]]['segment_flux']/tbl[np.where(tbl['label']==j)[0][0]]['area'].value for j in range(len(tbl))] 
#     for i in range(len(ap_kron_fluxes)-1,0, -1):
#         if ap_kron_fluxes[i]/ap_kron_fluxes[0] < 0.05 and segm_kron_dens[i]/segm_kron_dens[0]<0.08:
#             del data_process.apertures[i]
#     fit_run_list = []
#     data_process.PSF_list = []
#     data_process.psf_id_for_fitting = -1
#     for i in range(len(PSF_list_clean)):
#         psf = PSF_list_clean[i]
#         # psf = psf_clean(psf,if_plot=True, nsigma=3, npixels=45, ratio_to_replace=0.005)
#         print("Start to model use this PSF", i)
#         plt_fits(psf)
#         data_process.PSF_list.append(psf)
#         fit_sepc = FittingSpecify(data_process)
#         fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3) #, fix_n_list= [[0,4],[1,1]])
#         fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
#         fit_sepc.build_fitting_seq()
#         if i == 0:
#             fit_sepc.plot_fitting_sets()
#         fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep'])
#         fit_run.run(algorithm_list = ['PSO','PSO'])
#         fit_run.plot_final_qso_fit(target_ID =target_id)
#         fit_run_list.append(fit_run)
#     # pickle.dump([fit_run_list, PSF_RA_DEC_list], open(result_folder+filt+'_fit_result_PsfLib'+'_'+target_id+'.pkl', 'wb'))

# #%% Test use combinined PSF by psfr:

# chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
# idx_counts = chisqs.argsort()    

# import copy
# for ct in [5, 8, 'all']: 
#     _data_process = copy.deepcopy(data_process)
#     if ct != 'all':
#         PSF_list_for_comb = [PSF_list_clean[i] for i in idx_counts[:ct]]
#     elif ct == 'all':
#         PSF_list_for_comb = [PSF_list_clean[i] for i in idx_counts]
#     _data_process.PSF_list  = copy.deepcopy(PSF_list_for_comb)
#     _data_process.stack_PSF(if_plot = True, tool = 'psfr')
#     fit_sepc = FittingSpecify(_data_process)
#     fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3) #, fix_n_list= [[0,4],[1,1]])
#     fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
#     fit_sepc.build_fitting_seq()
#     fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep'])
#     fit_run.run(algorithm_list = ['PSO','PSO'])
#     fit_run.plot_final_qso_fit(target_ID =target_id)
#     fit_run_list.append(fit_run)
#     PSF_RA_DEC_list.append('stacked {0} PSF'.format(ct))

# pickle.dump([fit_run_list, PSF_RA_DEC_list], open(result_folder+filt+'_PSFcomb_PsfLib'+'_'+target_id+'.pkl', 'wb'))

# # host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
# # AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
# # ratio = host_flux/(host_flux+AGN_flux)
# # print(filt, round(fit_run.final_result_galaxy[0]['flux_within_frame'],2),
# #       "host ratio",round(ratio,2),
# #       round(fit_run.final_result_galaxy[0]['magnitude'],2),
# #       round(fit_run.final_result_galaxy[0]['n_sersic'],2),
# #       )

# #%% Find that webbPSF are too narrow: (running using astroconda)
# # import webbpsf
# # oversample = 2
# # nc = webbpsf.NIRCam()
# # nc.detector = 'NRCA5'
# # nc.filter = filt
# # POS = [data_process.target_pos[0]/2, data_process.target_pos[1]/2]
# # if POS[0] > 2047:
# #     POS[0] = POS[0]  - 5600/2
# #     nc.detector = 'NRCB5'
# # nc.detector_position = POS
# # psf_webb_fits = nc.calc_psf(oversample=oversample)
# # psf_webb = psf_webb_fits[0].data
# # from galight.tools.measure_tools import measure_FWHM
# # print(measure_FWHM(psf_webb))
# # print(measure_FWHM(PSF_list_clean[0]))
# # # data_process.PSF_list.append(psf_webb[1:,1:])
# # # fit_sepc = FittingSpecify(data_process)
# # # fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3) #, fix_n_list= [[0,4],[1,1]])
# # # fit_sepc.build_fitting_seq()
# # # fit_run = FittingProcess(fit_sepc, savename = target_id, fitting_level=['norm','deep'])
# # # fit_run.run(algorithm_list = ['PSO','PSO'])
# # # fit_run.plot_final_qso_fit(target_ID =target_id)
