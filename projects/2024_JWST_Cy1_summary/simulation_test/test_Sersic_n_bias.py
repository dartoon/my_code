#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:06:55 2024

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.plot_tools import plot_data_apertures_point
from galight.tools.astro_tools import plt_fits
from photutils.aperture import EllipticalAperture
import sys
import copy, matplotlib
import random

import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '../../2022_JWST_QSOz6/model_z6_data_id0/')
from target_info import target_info
run_folder = '../material/fit_result/'
idx = 6
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

filt = 'F150W'
#%%Grab PSF
add_cond = '_fixn1'
fit_run_list_sp1 = []
psf_sp = 1
load_files_sp1 = glob.glob(run_folder+'fit_run_{0}_idx{1}_psfidx*_psfsf{2}{3}.pkl'.format(filt, idx, psf_sp, add_cond))
load_files_sp1.sort()
chisqs_idx = []
for file in load_files_sp1:
    fit_run_list_sp1.append(pickle.load(open(file,'rb')))
chisqs = np.array([fit_run_list_sp1[i].reduced_Chisq for i in range(len(fit_run_list_sp1))])
sort_Chisq_sp1 = chisqs.argsort()  

diff_n = []
diff_Reff = []
for seed in range(10):
    PSF1, PSF2 = [fit_run_list_sp1[sort_Chisq_sp1[i]].fitting_specify_class.data_process_class.PSF_list[0] for i in random.sample(range(9), 2) ] 
    fit_run = fit_run_list_sp1[sort_Chisq_sp1[0]]
    
    #%% Make simulation taking 
    folder = '../../2022_JWST_QSOz6/NIRCam_data/*/bkg_removed/'   #!!!
    jwst_all_filenames = glob.glob(folder+'*{0}*{1}*{2}*_rmbkg.fits'.format(target_id[:5],target_id[-4:], filt))  #For NIRCam
    jwst_all_filenames.sort()
    file = jwst_all_filenames[0]
    im = pyfits.open(file)
    data = im[1].data
    header = im[1].header
    zp = fit_run.zp
    
    #%% Use grab a empty sky
    rad = 80                #The radius to cutout, final simulation in 2*rad+1
    for i in range(1000):
        pos = [int(np.random.uniform(400, len(data)-400)), int(np.random.uniform(400, len(data)-400))]  #!!!
        cut1 = data[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1]
        if np.sum(abs(cut1)>0.02) <10 and filt == 'F356W' or np.sum(abs(cut1)>0.1) <10 and filt == 'F150W':
            break
    from galight.tools.astro_tools import plt_fits
    print("The region will be added with our our simulation using 'pos' values")
    plt_fits(cut1)  
    
    #%% Start to make the mocks:
    import warnings
    warnings.filterwarnings("ignore")
    from lenstronomy.ImSim.image_model import ImageModel
    import lenstronomy.Util.simulation_util as sim_util
    from lenstronomy.Data.imaging_data import ImageData
    from lenstronomy.Data.psf import PSF
    from lenstronomy.PointSource.point_source import PointSource
    from lenstronomy.LightModel.light_model import LightModel
    import lenstronomy.Util.param_util as param_util
    from galight.data_process import DataProcess
    from galight.fitting_specify import FittingSpecify
    from galight.fitting_process import FittingProcess
        
    seed = 0
    psf = PSF1 #!!! Use PSF 
    psf[psf<0] = 0.  #!!! Only allow positive values in our PSF
    
    #Use Lenstronomy to make the simulation
    pix_s = fit_run.fitting_specify_class.deltaPix
    
    kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': psf, 'pixel_size': pix_s}
    psf_class = PSF(**kwargs_psf)
    kwargs_data = sim_util.data_configure_simple(2*rad+1, pix_s, inverse=True)
    data_class = ImageData(**kwargs_data)
    point_source_list = ['UNLENSED']
    
    kwargs_numerics = {'supersampling_factor': 5, 'supersampling_convolution': False} # 'point_source_supersampling_factor': 2}
    
    kwargs_sersic = {'R_sersic': 0.3, #fit_run.final_result_galaxy[0]['R_sersic'],   #The Reff size of the Sersic to use
                     'n_sersic': 1.3}
    
    light_model_list = ['SERSIC_ELLIPSE']
    lightModel = LightModel(light_model_list=light_model_list,sersic_major_axis=True) #sersic_major_axis = True in Galight and Galfit
    
    host_mag = 23 #fit_run.final_result_galaxy[0]['magnitude']  #The host magnitude #!!!
    _host_flux = 10**(0.4*(zp-host_mag))  
    q = fit_run.final_result_galaxy[0]['q']
    phi = np.random.uniform(0.,2*np.pi)   
    e1, e2 = param_util.phi_q2_ellipticity(phi=phi, q=q)
    kwargs_sersic['e1'] = e1
    kwargs_sersic['e2'] = e2
    kwargs_sersic['amp'] = 1
    center_x, center_y = np.random.uniform(-1.5, 1.5) * pix_s, np.random.uniform(-1.5,1.5)* pix_s
    kwargs_sersic['center_x'] = center_x
    kwargs_sersic['center_y'] = center_y
    kwargs_host = [kwargs_sersic]
    imageModel_host = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                    point_source_class=None, kwargs_numerics=kwargs_numerics)
    medi_host_flux = np.sum(imageModel_host.image(kwargs_lens_light=kwargs_host, unconvolved=True))
    amp = 1. / medi_host_flux * _host_flux        
    kwargs_sersic['amp'] = amp
    kwargs_host = [kwargs_sersic]
    host_highres = imageModel_host.image(kwargs_lens_light=kwargs_host, unconvolved=False)
    if im[0].header['CHANNEL'] == 'LONG':
        gain_value = 2
    else:
        gain_value = 1.8
    exp = header['XPOSURE']
    flux_mjsr = header['PHOTMJSR']
    #Adding Possion noise to the host mock
    host_highres_noise = np.random.poisson(lam=host_highres/header['PHOTMJSR']*exp*gain_value) / (1/header['PHOTMJSR']*exp*gain_value)  #Non-drizzled imaged
    
    pointSource = PointSource(point_source_type_list=point_source_list)
    imageModel_ps = ImageModel(data_class, psf_class, lens_light_model_class=None,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
    ps_mag = 21.5   #The PS magnitude #!!!
    _ps_amp = 10**(0.4*(zp-ps_mag))  
    PS_pos = np.random.uniform(-1.5, 1.5) * pix_s, np.random.uniform(-1.5,1.5)* pix_s
    kwargs_pointsource= {'ra_image': np.array([PS_pos[0]]),
                         'dec_image': np.array([PS_pos[1]]),
                         'point_amp': np.array([_ps_amp])}
    kwargs_ps = [kwargs_pointsource]
    ps_highres = imageModel_ps.image(kwargs_ps=kwargs_ps, unconvolved=False)
    
    # print("The simulation for Sersic host, without adding noise ")
    # plt_fits(host_highres_noise)  
    # print("The simulation for PS, without adding noise ")
    # plt_fits(ps_highres)  
    # print("The simulation for Sersic+PS, without adding noise ")
    # plt_fits(host_highres_noise+ps_highres)  
    
    #%% Adding simulated Host and PS to the mock data
    import copy 
    data_mock = copy.deepcopy(data)
    data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += host_highres_noise
    # data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += ps_highres
    
    #%%Fitting the mock data
    wht = im[4].data # The WHT map
    exp = im[0].header['EFFEXPTM']
    if im[0].header['CHANNEL'] == 'LONG':
        gain_value = 2
        expsize = 1
    else:
        gain_value = 1.8
        expsize = 1 
    flux_mjsr = header['PHOTMJSR']
    # exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
    data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header=header,
                                rm_bkglight = False, if_plot=False, zp = zp, exptime=exp_map)  #Gain value assuming as 1
    data_process.generate_target_materials(radius='7', 
                                           create_mask = False, nsigma=2.8, if_select_obj=False,
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
    data_process.apertures = [data_process.apertures[0]]
    use_psf = PSF2
    data_process.PSF_list = [use_psf]
    fit_sepc = FittingSpecify(data_process)
    fit_sepc.prepare_fitting_seq(point_source_num = 0) #, fix_n_list= [[0,4],[1,1]])
    fit_sepc.build_fitting_seq()
    # fit_sepc.kwargs_params['lens_light_model'][0] = [fit_run_.kwargs_result['kwargs_lens_light'][0]]
    plot_fit_name = 'sim_result_{1}_seed{0}'.format(seed,filt)
    fit_run_sim = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='norm')
    fit_run_sim.run(algorithm_list = ['PSO','PSO','PSO'], fitting_level=['norm','deep', 'deep'])
    fit_run_sim.plot_final_galaxy_fit()
    diff_n.append(fit_run_sim.final_result_galaxy[0]['n_sersic'] - 1.3)
    diff_Reff.append(fit_run_sim.final_result_galaxy[0]['R_sersic'] - 0.3)
    
