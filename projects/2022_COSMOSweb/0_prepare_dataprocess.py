#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 16:40:31 2023

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings, pickle
warnings.filterwarnings("ignore")
from astropy.wcs import WCS
from galight.tools.astro_tools import plt_many_fits
from galight.data_process import DataProcess
from galight.tools.astro_tools import read_pixel_scale
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pickle, copy, glob
from galight.tools.astro_tools import plt_fits
from galight.tools.measure_tools import mask_obj   
filefolder = '/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/'
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

#%%
filts = ['F115W', 'F150W','F277W', 'F444W']
run_i =42
for idx in range(run_i,len(cata_list)):
    pos = cata_list[idx][3:5]
    #To find the best size for 4 bands:
    _radius_l = []
    npixel = 100
    rad = '8'
    nsigma = 3.2
    if idx == 2:
        npixel = 20
    if idx == 10:
        rad = '5'
    if idx == 12:
        nsigma = 2
    if idx == 24:
        nsigma = 3.5
        npixel = 400
    if idx == 42:
        nsigma = 3.5
        npixel = 400
    if idx == 27:
        nsigma = 4.5
        npixel = 470
        
    for filt in filts:
        filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
        fitsFile = pyfits.open(filefolder+filename)
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        img = fitsFile[1].data #
        zp = -2.5*np.log10(2.350443 * 10**(-5) *read_pixel_scale(header)**2/3631)
        data_process = DataProcess(fov_image = img, target_pos = pos, pos_type = 'pixel', header = header,
                                  rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = fitsFile[2].data)
        try:
            data_process.generate_target_materials(radius=rad, #radius_list=radius_list, 
                                                   create_mask = False, nsigma=3, if_select_obj=False,
                                                   exp_sz= 1.2, npixels = 180, if_plot=False)
            _radius_l.append(data_process.radius)
        except:
            _radius_l.append(40)
    radius = np.max(_radius_l)
    data_process_list = []
    skip_i = [False] * len(filts)
    for i,filt in enumerate(filts):
        filename = 'mosaic_nircam_f{0}w_COSMOS-Web_30mas_v0_1_i2d.fits'.format(filt[1:-1])
        fitsFile = pyfits.open(filefolder+filename)
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        img = fitsFile[1].data #
        pixscale = read_pixel_scale(header)
        zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631)
        fov_noise_map = fitsFile[2].data
        data_process = DataProcess(fov_image = img, target_pos = pos, pos_type = 'pixel', header = header,
                                  rm_bkglight = False, if_plot=False, zp = zp, fov_noise_map = fov_noise_map)
        try:
            data_process.generate_target_materials(radius=radius, #radius_list=radius_list, 
                                                    create_mask = False, nsigma= nsigma, if_select_obj=False,
                                                    exp_sz= 1.2, npixels = npixel, if_plot=False,contrast=0.001)
        except:
            data_process.generate_target_materials(radius=radius, #radius_list=radius_list, 
                                                    create_mask = False, skip=True)
            skip_i[i] = True
        data_process_list.append(data_process)
    
    run_list = [3,2,1,0]
    apertures = data_process_list[run_list[0]].apertures
    for i in run_list:
            covers = mask_obj(data_process_list[i].target_stamp, apertures, if_plot=False, sum_mask = True)
            if skip_i[i] == False:
                for j in range(len(data_process_list[i].apertures)):
                    new_cover = mask_obj(data_process_list[i].target_stamp, [data_process_list[i].apertures[j]], if_plot=False, sum_mask = True)
                    if np.sum(covers - new_cover*covers) > np.sum(1-new_cover)/2 :               
                        apertures.append(data_process_list[i].apertures[j])
    rm_list = []
    for i in range(len(apertures)):
        all_cover = mask_obj(data_process_list[run_list[0]].target_stamp, apertures[:i]+apertures[i+1:], if_plot=False, sum_mask = True)
        one_cover = mask_obj(data_process_list[run_list[0]].target_stamp, [apertures[i]], if_plot=False, sum_mask = True)
        if  np.sum(all_cover) - np.sum(all_cover*one_cover) < np.sum(1-one_cover)/1.6:
            rm_list.append(i) #remove the coverred apertures
    apertures = [apertures[i] for i in range(len(apertures)) if i not in rm_list]          
    
    if idx == 18:
        del apertures[1]
    
    for i in run_list: 
        filt = filts[i]
        data_process = data_process_list[i]
        data_process.apertures = apertures
        PSF_lib_files = glob.glob('material/'+filt+'_PSF_Library_v2.pkl')[0]
        PSF_org_list, PSF_clean_list, all_PSF_pos_list, PSF_RA_DEC_list = pickle.load(open(PSF_lib_files,'rb'))
        data_process.noise_map[data_process.noise_map < 0.00025] = np.max(data_process.noise_map[data_process.noise_map!=0])
        data_process.filt = filt
        # data_process.plot_aperture()
        del data_process.fov_image
        del data_process.fov_noise_map
        # data_process.plot_aperture()
        # pickle.dump(data_process, open('material/'+'data_process_idx{0}_{1}.pkl'.format(i, filt), 'wb'))
        sz = 30
        ct = int((len(data_process.target_stamp) - 30)/2)
        test_img = data_process.target_stamp[ct:-ct,ct:-ct]
        if i == run_list[0]:
            ps_pos = np.array(np.where(test_img == test_img.max())) - len(test_img)/2
            ps_pos = ps_pos[::-1]
        target_id = cata_list[idx][7]
        for j in range(len(PSF_clean_list)):
            savename = 'fit_material/fit2_notrunyet_{2}_idx{0}_psf{1}.pkl'.format(idx,j,filt)
            psf = PSF_clean_list[j]
            ct = int((len(psf) - len(data_process.target_stamp))/2)
            if ct >10:
                psf = psf[ct:-ct, ct:-ct]
            psf[psf<0] = 0.
            data_process.PSF_list = [psf]
            fit_sepc = FittingSpecify(data_process)
            ps_pix_center_list = [copy.deepcopy(ps_pos)]
            if idx == 15:
                ps_pix_center_list = None
            fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor = 3, apertures_center_focus=True,
                                          ps_pix_center_list = ps_pix_center_list )
            # fit_sepc.kwargs_params['lens_light_model'][3][0]['R_sersic'] = 0.06
            # fit_sepc.kwargs_params['lens_light_model'][4][0]['R_sersic'] = 1.
            if j == 0:
                fit_sepc.plot_fitting_sets()
                print(idx, filt, 'apertures', len(data_process.apertures) )
                hold = input('Hold ... OK?\n')
            fit_run = FittingProcess(fit_sepc, savename = target_id)
            pickle.dump(fit_run , open(savename, 'wb'))
            # savename = savename.replace('_notrunyet_', '_run_')[:-4]+'_{0}.pkl'.format(i)
            # fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
            # fit_run.plot_final_qso_fit(target_ID =target_id)
            