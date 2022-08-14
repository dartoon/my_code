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
from source_info import source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise
from scipy.ndimage import zoom
import copy
from scipy import signal
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.tools.astro_tools import plt_many_fits
from galight.fitting_process import FittingProcess
from galight.tools.measure_tools import detect_obj
import sys

if_plot = False
filt_id = int(sys.argv[2])
filt = ['f150w', 'f356w'][filt_id]
folder = '../prep_use_HST_highRes/JWST_CEERS/'
ceers_file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)
target_info = pickle.load(open('../prep_use_HST_highRes/target_info.pkl','rb'))

psfs, FWHMs = pickle.load(open(filt+'_psfs_star.pkl','rb'))
fluxs= []
for i in range(len(psfs)):
    fluxs.append(np.sum(psfs[i]))
fluxs = np.array(fluxs)
    
def find_close_PSF_idx(psf_id):
    sort = np.argsort(abs(FWHMs[psf_id] - FWHMs))[1:]
    idx = sort[0] #Setting a initial value
    for i in sort:
        # print( abs((fluxs[i]- fluxs[psf_id])/fluxs[psf_id]  ) )
        if abs( (fluxs[i]- fluxs[psf_id])/fluxs[psf_id])<0.5 and fluxs[i]>500:
            idx = i
            break
    return idx

keys = []
for key in target_info.keys():
    keys.append(key)
#%%
import glob
# TNG_files = glob.glob('TNG50_img/*photo.fits')
# TNG_files.sort()
TNG_ids = ['101491', '101499', '119454', '140982', '15', '219845', '272230', '321717', '561512', '579945']
TNG_files = [ 'TNG50_img/shalo_091-{0}_v0_photo.fits'.format(i) for i in TNG_ids  ]

#!!!
#U band 364 nm ~ 368 * 7.2 /1.1 = 2408 nm for F150w
#G band 464 nm ~ 480 * 7.2 /1.1 = 3141 nm for F356w
#R band 622 nm ~ 622 * 7.2 /1.1 = 4071 nm
TNG_band = ['CFHT_MegaCam.u', 'SUBARU_HSC.G'][filt_id] #
from galight.tools.measure_tools import flux_profile

seed = int(sys.argv[1])
for ID in range(len(keys)):
    np.random.seed(seed = seed)
    name = keys[ID] # ID of target
    host_flux_ratio = np.random.uniform(0.07,0.95) #
    # host_Reff_kpc = np.random.uniform(1,3)   #Host effective radius, unit: Kpc #
    source_id = np.random.randint(0,len(TNG_ids)) #
    im = pyfits.open(folder+ceers_file)
    data = im[1].data
    header = im[1].header
    flux_mjsr = header['PHOTMJSR']
    data = data/flux_mjsr # To change MJ/sr to flux
    
    psf_id = np.random.randint(0,len(psfs))
    if filt=='f356w':
        data_sbs = [data[2000:4000,2000:4000], data[100:2100,100:2100]]
        fov_cut_idx = np.random.randint(0,len(data_sbs))
        data_sb = data_sbs[fov_cut_idx]   #
    elif filt=='f150w':
        pos_list = [[0,0], [0,4600], [4600,0], [4600, 4600], [250+0,11550+0], [250+0,11550+4600], [250+4600,11550+0], [250+4600, 11550+4600]]
        data_sbs = []
        for s_pos in pos_list:
            data_sbs.append(data[s_pos[0]:s_pos[0]+4400,s_pos[1]:s_pos[1]+4400])
        fov_cut_idx = np.random.randint(0,len(data_sbs))
        data_sb = data_sbs[fov_cut_idx]   #
    psf_true = psfs[psf_id] #
    psf_true = psf_true/np.sum(psf_true) 
    psf_id_model = find_close_PSF_idx(psf_id)
    psf_model = psfs[psf_id_model]
    
    for i in range(200):
        if filt=='f150w':
            pos = [int(np.random.uniform(400, 4000)), int(np.random.uniform(400, 4000))]  #
        elif filt=='f356w':
            pos = [int(np.random.uniform(300, 1700)), int(np.random.uniform(300, 1700))]  #
        rad1 = 100
        cut1 = data_sb[ pos[1]-rad1:pos[1]+rad1+1, pos[0]-rad1:pos[0]+rad1+1]
        _, _, _, tbl = detect_obj(cut1, nsigma=1, npixels=6, use_moments=False)
        tbl_kron_flux = np.nan_to_num(tbl['kron_flux'])
        if np.max(tbl_kron_flux) < 4:
            break
    
    pixscale = read_pixel_scale(header)
    zp = -2.5*np.log10(2.350443 * 10**(-5) *pixscale**2/3631) - 2.5*np.log10(flux_mjsr)  #zp for flux
    header0 = im[0].header
    img_filter = header0['FILTER']
    img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
    exptime = header0['TEXPTIME'] #The assumed exp time.
    
    #%%For mock galaxy
    #Generate the QSO galaxy info
    z = target_info[name]['z']
    qso_mag = target_info[name]['mag_'+filt]
    qso_total_flux = 10**(-0.4*(qso_mag-zp))
    galaxy_flux = qso_total_flux * host_flux_ratio
    qso_flux = qso_total_flux - galaxy_flux
    # cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    # scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    # host_Reff = host_Reff_kpc/scale_relation   #In arcsec
    # host_Reff_pix = host_Reff/pixscale
    # print('host_flux_ratio', host_flux_ratio, 'host_Reff_kpc', host_Reff_kpc, 'host_Reff_pix', host_Reff_pix)
    source_name=source_list[0][source_id]
    # print(source_name)
    
    # hdu = pyfits.open('HST_real_source/{0}_fix.fits'.format(source_name))
    # hd_gal_img = hdu[0].data
    # hdu.close()
    # Re = source_list[2][source_id]
    # project_gal_img = zoom(hd_gal_img, host_Reff_pix/Re)
    # project_gal_img = project_gal_img/np.sum(project_gal_img)*galaxy_flux
    galaxy_mag = -2.5*np.log10(galaxy_flux) + zp
    
    # TNG_file = 'TNG50_img/shalo_091-540258_v0_photo.fits'
    # file = 'TNG50_img/shalo_091-542669_v0_photo.fits'
    # file = 'TNG50_img/shalo_091-572599_v0_photo.fits'
    TNG_file = TNG_files[source_id]
    im = pyfits.open(TNG_file)
    fits = im[TNG_band].data  
    flux_hd = 10**(-0.4*(fits-zp))
    z_s = 6.2
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/"
    scale = 0.1/scale_relation/ pixscale #Project as 0.1 #!!!
    flux_zoom = zoom(flux_hd, scale)
    flux_zoom = flux_zoom/np.sum(flux_zoom) * galaxy_flux
    
    _fluxs, rad, _ = flux_profile(flux_zoom, center=[len(flux_zoom)/2]*2 , radius=len(flux_zoom)/2, if_plot=False, #if_annuli=(True), 
                      fits_plot=(False), grids=50, x_gridspace=None)
    Reff_rad = rad[_fluxs<_fluxs[-1]/2][-1]
    Reff_kpc = Reff_rad / scale * 0.1
    Reff_arcsec = Reff_rad * pixscale
    
    #convolve:
    conv_project_gal_img = signal.fftconvolve(flux_zoom, psf_true, mode='full')
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
    
    #Get the Noise QSO - AGN image.
    cut = (len(conv_project_gal_img) - len(psf_true))/2
    _noise_conv_project_qsosub_img = copy.deepcopy(noise_conv_project_qso_img)
    if cut == int(cut):
        cut = int(cut)
        _noise_conv_project_qsosub_img[cut:-cut,cut:-cut] = _noise_conv_project_qsosub_img[cut:-cut,cut:-cut] - qso_flux*psf_true
    else:
        cut = int(cut)
        _noise_conv_project_qsosub_img[cut:-cut-1,cut:-cut-1] = _noise_conv_project_qsosub_img[cut:-cut-1,cut:-cut-1] - qso_flux*psf_true
    # plt_fits(hd_gal_img)    
    # plt_fits(project_gal_img)
    # plt_fits(conv_project_qso_img)
    # plt_fits(noise_conv_project_qso_img)
    folder_save = 'sim_results/'
    filename = folder_save+ 'qsoID'+str(ID)+'_filt_'+filt+'_seed'+str(seed) #+'_smallcutout'
    #Plot and save sim details
    plot_sim_name = filename + '_sim.pdf'
    labels = ['org_gal_img', 'project_gal_img', 'conv_gal_img', 'add_ps+Poss_noi.', 'Host image (ps sub.)']
    plt_many_fits([flux_hd, flux_zoom, conv_project_gal_img, noise_conv_project_qso_img, _noise_conv_project_qsosub_img], labels = labels,
                  savename=plot_sim_name, if_plot=if_plot)
    
    #%%Obtain PSF stars:
    data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header = header,
                                rm_bkglight = True, exptime = np.ones_like(data_sb)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
    data_process.generate_target_materials(radius=np.max([Reff_rad*3.5, 80]), 
                                           create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 15, if_plot=if_plot)
    data_process.plot_overview(label = filt+'_'+str(fov_cut_idx)+'_FOV', target_label=name[:7],
                               ifsave=True, filename = filename+'_FOV', if_plot=if_plot)
    # data_process.find_PSF(radius = 50, user_option = True, psf_edge=10)
    data_process.apertures = [data_process.apertures[0]]
    data_process.deltaPix = pixscale
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
        info['inferred_R_sersic_'+name] = fit_run.final_result_galaxy[0]['R_sersic'] * fit_run.final_result_galaxy[0]['q']
        info['inferred_n_sersic_'+name] = fit_run.final_result_galaxy[0]['n_sersic']
        info['plot_fit_name_'+name] = plot_fit_name +'_qso_final_plot.pdf'
        print('inferred galaxy flux, mag, Re (arcsec)):\n\t', round(fit_run.final_result_galaxy[0]['flux_within_frame'],2), 
              round(fit_run.final_result_galaxy[0]['magnitude'],2), 
              round(fit_run.final_result_galaxy[0]['R_sersic'],2))
        fit_run.plot_final_qso_fit(target_ID=name[:5]+'_HostRto_'+str(round(host_flux_ratio,3)), save_plot=True, 
                                   show_plot=if_plot)
        # fit_run.dump_result()
    #%%
    # print('True galaxy flux, mag, Re (arcsec)):\n\t', round(galaxy_flux,2), round(galaxy_mag,2),
    #       round(host_Reff,2))
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
    info['true_host_Reff_arcsec'] = Reff_arcsec
    info['true_host_Reff_kpc'] = Reff_kpc
    info['TNG_ID'] = TNG_file
    # info['assumed_host_Re_kpc'] = host_Reff_kpc
    # info['galfit_Re'] = host_Reff
    # pickle.dump(info, open(filename+'_result.pkl', 'wb'))   
    # infos = ID, filt, zp, target_name, host_flux_ratio, host_Reff_kpc, host_Reff, source_id, fov_cut_idx, psf_id, add_pos
    # save_result = true_flux, true_mag, true_galfit_Re, inferred_flux, inferred_mag, inferred_Re, inferred_n
    # save_plot_name
