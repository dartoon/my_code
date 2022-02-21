#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:36:21 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from galight.tools.astro_tools import plt_fits
from galight.tools.astro_tools import read_pixel_scale
from astropy.cosmology import FlatLambdaCDM
import sys
# sys.path.insert(0,'/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source')
from source_info import source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise
from scipy.ndimage import zoom
import copy
from scipy import signal
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.tools.astro_tools import plt_many_fits
from galight.fitting_process import FittingProcess

if_plot = False

filt_id = int(sys.argv[2])
filt = ['f150w', 'f356w'][filt_id]
# filt = 'f150w' #!!!
# filt = 'f356w'
folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)
target_info = pickle.load(open('target_info.pkl','rb'))
_psfs, _FWHMs = pickle.load(open(filt+'_psfs.pkl','rb'))
psfs, FWHMs, fluxs = [], [], []
for i in range(len(_psfs)):
    psfs = psfs + _psfs[i]
    FWHMs = FWHMs + _FWHMs[i]
    fluxs = fluxs + [np.sum(_psfs[i][j]) for j in range(len(_psfs[i]))]
    
FWHMs, fluxs = np.array(FWHMs), np.array(fluxs)

def find_close_PSF_idx(psf_id):
    sort = np.argsort(abs(FWHMs[psf_id] - FWHMs))[1:]
    for i in sort:
        # print( abs((fluxs[i]- fluxs[psf_id])/fluxs[psf_id]  ) )
        if abs( (fluxs[i]- fluxs[psf_id])/fluxs[psf_id])<0.5 and fluxs[i]>500:
            idx = i
            break
    # print(sort)
    return idx
# print(find_close_PSF_idx(0))
# 
keys = []
for key in target_info.keys():
    keys.append(key)
#%%
# ID = 0
# seed = 0
# for seed in range(0, 20):
seed = int(sys.argv[1])
for ID in range(12):
    np.random.seed(seed = seed)
    name = keys[ID] #!!! ID of target
    host_flux_ratio = np.random.uniform(0.07,0.9) #!!!
    host_Reff_kpc = np.random.uniform(1,3)   #Host effective radius, unit: Kpc #!!!
    source_id = np.random.randint(0,25) #!!!
    
    im = pyfits.open(folder+file)
    data = im[1].data
    header = im[1].header
    # print('For flux value in unit of MJy/sr.') #https://en.wikipedia.org/wiki/AB_magnitude
    # value_unit = header['BUNIT']
    # print("Data unit:", value_unit)
    # flux(Mjy/sr) * 2.350443 * 10**(-5) *0.03**2   #Flux to Jy  https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/18/
    
    psf_id = np.random.randint(0,len(psfs))
    if filt=='f356w':
        data_sbs = [data[2000:4000,2000:4000], data[100:2100,100:2100]]
        fov_cut_idx = np.random.randint(0,len(data_sbs))
        data_sb = data_sbs[fov_cut_idx]   #!!!
    elif filt=='f150w':
        pos_list = [[0,0], [0,4600], [4600,0], [4600, 4600], [250+0,11550+0], [250+0,11550+4600], [250+4600,11550+0], [250+4600, 11550+4600]]
        data_sbs = []
        for s_pos in pos_list:
            data_sbs.append(data[s_pos[0]:s_pos[0]+4400,s_pos[1]:s_pos[1]+4400])
        fov_cut_idx = np.random.randint(0,len(data_sbs))
        data_sb = data_sbs[fov_cut_idx]   #!!!
    psf_true = psfs[psf_id] #!!!
    psf_true = psf_true/np.sum(psf_true) 
    psf_id_model = find_close_PSF_idx(psf_id)
    psf_model = psfs[psf_id_model]
    
    from galight.tools.measure_tools import detect_obj
    from photutils.segmentation import SourceCatalog
    for i in range(500):
        if filt=='f150w':
            pos = [int(np.random.uniform(200, 4200)), int(np.random.uniform(200, 4200))]  #!!!
        elif filt=='f356w':
            pos = [int(np.random.uniform(100, 1900)), int(np.random.uniform(100, 1900))]  #!!!
        rad1 = 100
        cut1 = data_sb[ pos[1]-rad1:pos[1]+rad1+1, pos[0]-rad1:pos[0]+rad1+1]
        try:
            res, segm_map = detect_obj(cut1, segm_map=True)
            cat = SourceCatalog(cut1, segm_map)
            if np.max(cat.to_table()['kron_flux']) < 3:
                _ = detect_obj(cut1, if_plot=False)
                # print('max flux in empty fov:', np.max(cat.to_table()['kron_flux']))
                break
        except:
            break
        # rad2 = 40
        # cut2 = data_sb[ pos[1]-rad2:pos[1]+rad2+1, pos[0]-rad2:pos[0]+rad2+1]
        # if np.sum(cut1)<5 and np.sum(cut2)<1: #!!!Need to be updated by detect_obj()
        #     break
    
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631)
    header0 = im[0].header
    img_filter = header0['FILTER']
    img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
    exptime = header0['TEXPTIME'] #The assumed exp time.
    # plt_fits(data_sb)
    
    #%%For mock galaxy
    #Generate the QSO galaxy info
    z = target_info[name]['z']
    qso_mag = target_info[name]['mag_'+filt]
    qso_total_flux = 10**(-0.4*(qso_mag-zp))
    galaxy_flux = qso_total_flux * host_flux_ratio
    qso_flux = qso_total_flux - galaxy_flux
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    host_Reff = host_Reff_kpc/scale_relation   #In arcsec
    host_Reff_pix = host_Reff/pixscale
    
    # print('host_flux_ratio', host_flux_ratio, 'host_Reff_kpc', host_Reff_kpc, 'host_Reff_pix', host_Reff_pix)
    source_name=source_list[0][source_id]
    # print(source_name)
    
    hdu = pyfits.open('HST_real_source/{0}_fix.fits'.format(source_name))
    hd_gal_img = hdu[0].data
    hdu.close()
    
    Re = source_list[2][source_id]
    
    #The inferred are much small.
    # from galight.tools.measure_tools import flux_profile
    # r_flux, r_grids, regions = flux_profile(hd_gal_img, center = [len(hd_gal_img)/2]*2, radius=len(hd_gal_img)/2, grids=40,
    #                           if_plot=True, fits_plot=True)
    # Re = r_grids[np.sum(hd_gal_img)/2 <r_flux][0]
    
    project_gal_img = zoom(hd_gal_img, host_Reff_pix/Re)
    project_gal_img = project_gal_img/np.sum(project_gal_img)*galaxy_flux
    galaxy_mag = -2.5*np.log10(galaxy_flux) + zp
    
    #convolve:
    conv_project_gal_img = signal.fftconvolve(project_gal_img, psf_true, mode='full')
    
    conv_project_qso_img = copy.deepcopy(conv_project_gal_img)
    psf_true = psf_true/np.sum(psf_true) #Make sure PSF is normalized
    cut = (len(conv_project_gal_img) - len(psf_true))/2
    if cut == int(cut):
        cut = int(cut)
        conv_project_qso_img[cut:-cut,cut:-cut] = conv_project_qso_img[cut:-cut,cut:-cut] + qso_flux*psf_true
    else:
        cut = int(cut)
        conv_project_qso_img[cut:-cut-1,cut:-cut-1] = conv_project_qso_img[cut:-cut-1,cut:-cut-1] + qso_flux*psf_true
    
    #Add Poisson Noise:
    noise_conv_project_qso_img=np.random.poisson(lam=abs(conv_project_qso_img)*exptime)/(exptime)
    #Add to fov image
    rad = len(noise_conv_project_qso_img)/2
    data_mock = copy.deepcopy(data_sb)
    if rad != int(rad):
        rad = int(rad)
        data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += noise_conv_project_qso_img
    elif rad == int(rad):
        rad = int(rad)
        data_mock[ pos[1]-rad:pos[1]+rad, pos[0]-rad:pos[0]+rad] += noise_conv_project_qso_img
    
    # plt_fits(hd_gal_img)    
    # plt_fits(project_gal_img)
    # plt_fits(conv_project_qso_img)
    # plt_fits(noise_conv_project_qso_img)
    folder_save = 'sim_results/'
    filename = folder_save+ 'qsoID'+str(ID)+'_filt_'+filt+'_seed'+str(seed)
    #Plot and save sim details
    plot_sim_name = filename + '_sim.pdf'
    labels = ['org_gal_img', 'project_gal_img', 'conv_gal_img', 'add_ps_noise']
    plt_many_fits([hd_gal_img, project_gal_img, conv_project_gal_img, noise_conv_project_qso_img], labels = labels,
                  savename=plot_sim_name, if_plot=False)
    #%%Obtain PSF stars:
    
    data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header = header,
                                rm_bkglight = True, exptime = np.ones_like(data_sb)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
    data_process.generate_target_materials(radius=rad*0.8, create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
    data_process.plot_overview(label = filt+'_'+str(fov_cut_idx)+'_FOV', target_label=name[:7],
                               ifsave=True, filename = filename+'_FOV', if_plot=if_plot)
    # data_process.find_PSF(radius = 50, user_option = True, psf_edge=10)
    data_process.apertures = [data_process.apertures[0]]
    info = {}
    for psf, name in [[psf_true,'same_psf'], [psf_model,'diff_psf']]:
    # for psf, name in [[psf_model,'diff_psf']]:    
        plot_fit_name = filename + name+ '_fit'
        data_process.PSF_list = [psf]
        #Start to produce the class and params for lens fitting.
        fit_sepc = FittingSpecify(data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
        fit_sepc.build_fitting_seq()
        #Plot the initial settings for fittings. 
        # fit_sepc.plot_fitting_sets()
        #
        fit_run = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='norm')
        fit_run.run(algorithm_list = ['PSO', 'PSO'])
        if fit_run.reduced_Chisq > 1.3:
            fit_run = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='deep')
            fit_run.run(algorithm_list = ['PSO', 'PSO'])
        info['inferred_host_flux_'+name] = fit_run.final_result_galaxy[0]['flux_within_frame']
        info['inferred_magnitude_'+name] = fit_run.final_result_galaxy[0]['magnitude']
        info['inferred_R_sersic_'+name] = fit_run.final_result_galaxy[0]['R_sersic']
        info['inferred_n_sersic_'+name] = fit_run.final_result_galaxy[0]['n_sersic']
        info['plot_fit_name_'+name] = plot_fit_name +'_qso_final_plot.pdf'
        print('inferred galaxy flux, mag, Re (arcsec)):\n\t', round(fit_run.final_result_galaxy[0]['flux_within_frame'],2), 
              round(fit_run.final_result_galaxy[0]['magnitude'],2), 
              round(fit_run.final_result_galaxy[0]['R_sersic'],2))
        fit_run.plot_final_qso_fit(target_ID=name[:5]+'_HostRto_'+str(round(host_flux_ratio,3)), save_plot=True, 
                                   show_plot=if_plot)
    #%%
    print('True galaxy flux, mag, Re (arcsec)):\n\t', round(galaxy_flux,2), round(galaxy_mag,2),
          round(host_Reff,2))
    info['PSF_id_true'] = psf_id
    info['PSF_id_model'] = psf_id_model
    info['ID'] = ID
    info['filter'] = filt
    info['zp'] = zp
    info['target_name'] = name
    info['source_id'] = source_id
    info['insert_fov'] = cut1
    info['fov_cut_idx'] = fov_cut_idx
    info['plot_sim_name'] = plot_sim_name
    info['true_host_flux_ratio'] = host_flux_ratio
    info['true_host_flux'] = galaxy_flux
    info['true_host_mag'] = galaxy_mag
    info['assumed_host_Re_kpc'] = host_Reff_kpc
    info['galfit_Re'] = host_Reff
    pickle.dump(info, open(filename+'_result.pkl', 'wb'))   
    
    # infos = ID, filt, zp, target_name, host_flux_ratio, host_Reff_kpc, host_Reff, source_id, fov_cut_idx, psf_id, add_pos
    # save_result = true_flux, true_mag, true_galfit_Re, inferred_flux, inferred_mag, inferred_Re, inferred_n
    # save_plot_name
