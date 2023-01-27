#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:12:51 2022

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import sys
import copy
sys.path.insert(0,'../model_z6_data_id0/')
from target_info import target_info

#%%Input number
filt = 'F356W'
idx = 0
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx)

#%% Load data
# folder = '../NIRCam_data/Jan14/bkg_removed'.format(idx)
folder = '../NIRCam_data/Nov14/bkg_removed'.format(idx)
info = target_info[str(idx)]
files = glob.glob(folder+'/*{1}*{2}*{0}*.fits'.format(filt, info['target_id'][:5], info['target_id'][-4:]))
if len(files)==1:
    file = files[0]
else:
    raise ValueError('More than one files found')
im = pyfits.open(file)
data = im[1].data
header = im[1].header

#Load PSF information:
PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
PSF_list, PSF_list_clean, PSF_RA_DEC_list, PSF_from_file_list = pickle.load(open(PSF_lib_files,'rb'))

#Load top PSF result from the overall fittings
if idx == 1:
    fit_files = glob.glob(run_folder+'*fit_material*/fit_run*_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
elif idx == 0:
    fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
fit_files.sort()
fit_run_list = []
for i in range(len(fit_files)):
    fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
sort_Chisq = chisqs.argsort() 
fit_run_ = fit_run_list[sort_Chisq[0]]  # use top PSF result to run the simulation.

#%%
from galight.tools.measure_tools import detect_obj
rad = 50
for i in range(1000):
    pos = [int(np.random.uniform(400, len(data)-400)), int(np.random.uniform(400, len(data)-400))]  #!!!
    cut1 = data[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1]
    if np.sum(abs(cut1)>0.02) <10 and filt == 'F356W' or np.sum(abs(cut1)>0.1) <10 and filt == 'F150W':
        # apertures, segm_deblend, mask_apertures, tbl = detect_obj(cut1, nsigma=1, npixels = 8)
        # if np.max(np.nan_to_num(tbl['kron_flux']))<0.3:
        break
from galight.tools.astro_tools import plt_fits
# print(i, np.max(np.nan_to_num(tbl['kron_flux'])))
plt_fits(cut1)

cut1 = cut1 / np.sqrt(580/fit_run_.fitting_specify_class.data_process_class.header['XPOSURE'])

#%% Make the simulation/
#One random part is the background. DONE!!!
#!!! One random part the parameter

import warnings
warnings.filterwarnings("ignore")
from lenstronomy.ImSim.image_model import ImageModel
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LightModel.light_model import LightModel
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess

chisqs_PSF = np.array([fit_run_list[i].reduced_Chisq for i in range(3, len(fit_run_list))])
sort_PSF_Chisq = chisqs_PSF.argsort()

Based_PSF_list= [PSF_list_clean[i] for i in sort_PSF_Chisq]

# if idx == 1:
#     del Based_PSF_list[1]
import lenstronomy.Util.param_util as param_util
zp = fit_run_.zp
# for seed in range(30,100):
for seed in range(100,101):
    sim_psf = np.random.randint(5)
    for fit_psf in range(4):
        psf = Based_PSF_list[sim_psf] #!!! Use PSF 
        psf[psf<0] = 0.
        use_psf_list = Based_PSF_list[:sim_psf]+Based_PSF_list[sim_psf+1:]
        pix_s = fit_run_.fitting_specify_class.deltaPix
        kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': psf, 'pixel_size': pix_s}
        psf_class = PSF(**kwargs_psf)
        kwargs_data = sim_util.data_configure_simple(2*rad+1, pix_s, inverse=True)
        data_class = ImageData(**kwargs_data)
        point_source_list = ['UNLENSED']
        pointSource = PointSource(point_source_type_list=point_source_list)
        # kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
        light_model_list = ['SERSIC_ELLIPSE']
        lightModel = LightModel(light_model_list=light_model_list,sersic_major_axis=fit_run_.sersic_major_axis)
        kwargs_numerics = {'supersampling_factor': 5, 'supersampling_convolution': False}
        imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                        point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
        kwargs_sersic = fit_run_.kwargs_result['kwargs_lens_light'][0]
        q = np.random.uniform(0.5,0.9)
        phi = np.random.uniform(0.,2*np.pi)   
        # host_flux = 15
        host_flux = 60  #!!!
        e1, e2 = param_util.phi_q2_ellipticity(phi=phi, q=q)
        kwargs_sersic['e1'] = e1
        kwargs_sersic['e2'] = e2
        kwargs_sersic['amp'] = 1
        center_x, center_y = np.random.uniform(-1.5, 1.5) * pix_s, np.random.uniform(-1.5,1.5)* pix_s
        kwargs_sersic['center_x'] = center_x
        kwargs_sersic['center_y'] = center_y
        kwargs_host = [kwargs_sersic]
        medi_host_flux = np.sum(imageModel.image(kwargs_lens_light=kwargs_host, unconvolved=True))
        amp = 1. / medi_host_flux * host_flux        
        kwargs_sersic['amp'] = amp
        kwargs_host = [kwargs_sersic]
        
        kwargs_pointsource= fit_run_.kwargs_result['kwargs_ps'][0]
        kwargs_ps = [kwargs_pointsource]
        host_highres = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=kwargs_ps, 
                                        unconvolved=False)
        wht = im[4].data # The WHT map
        exp = im[0].header['EFFEXPTM']
        if im[0].header['CHANNEL'] == 'LONG':
            gain_value = 2
            expsize = 1
        else:
            gain_value = 1.8
            expsize = 1 
            
        host_highres= np.random.poisson(lam=host_highres*580*gain_value)/(580*gain_value)  #Non-drizzled imaged
        
        data_mock = copy.deepcopy(data)
        data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += host_highres
            
        flux_mjsr = header['PHOTMJSR']
        exp_map = 580 * wht/wht.max() / flux_mjsr * gain_value
        data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header=fit_run_.fitting_specify_class.data_process_class.header,
                                    rm_bkglight = False, if_plot=False, zp = zp, exptime=exp_map)  #Gain value assuming as 1
        data_process.generate_target_materials(radius=int(len(fit_run_.image_host_list[0])/2), 
                                                create_mask = False, nsigma=2.8, if_select_obj=False,
                                              exp_sz= 1.2, npixels = 15, if_plot=False)
        data_process.apertures = [data_process.apertures[0]]
        use_psf = use_psf_list[fit_psf]
        use_psf[use_psf<0] = 0.
        data_process.PSF_list = [use_psf]
        fit_sepc = FittingSpecify(data_process)
        if idx !=1:
            fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
        else:
            fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[0,1]])
        fit_sepc.build_fitting_seq()
        fit_sepc.kwargs_params['lens_light_model'][0] = [fit_run_.kwargs_result['kwargs_lens_light'][0]]
        plot_fit_name = 'sim_result_less_exp/sim_idx{3}_{4}_seed{0}BasedPSF{1}_FitBasedPSF{2}'.format(seed,sim_psf,fit_psf,idx, filt)
        fit_run = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='norm')
        fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
        fit_run.plot_final_qso_fit()
        print(fit_run.final_result_galaxy)
        fit_run.dump_result()
    
    #%% recombin use top-2 PSF:
    #Load top PSF result from the overall fittings
    fit_files = glob.glob('sim_result_less_exp/sim_idx{1}_{2}_seed{0}B*.pkl'.format(seed,idx,filt))#+\
    fit_files.sort()
    fit_run_list = []
    for i in range(len(fit_files)):
        fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
    chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
    idx_counts = chisqs.argsort() 
    ct = 2
    _data_process_list = [fit_run_list[i].fitting_specify_class.data_process_class for i in idx_counts[:ct]] 
    PSF_list_for_comb = [_data_process_list[i].PSF_list[0] for i in range(len(_data_process_list))]
    data_process = copy.deepcopy(_data_process_list[0])
    data_process.PSF_list  = copy.deepcopy(PSF_list_for_comb)
    data_process.stack_PSF(if_plot = False, tool = 'psfr')
    fit_sepc = FittingSpecify(data_process)
    if idx !=1:
        fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
    else:
        fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[0,1]])
    fit_sepc.build_fitting_seq()
    fit_sepc.kwargs_params['lens_light_model'][0] = [fit_run_.kwargs_result['kwargs_lens_light'][0]]
    plot_fit_name = 'sim_result_less_exp/sim_idx{2}_{3}_seed{0}BasedPSF{1}_FitBasedCombPSF'.format(seed,sim_psf,idx,filt)
    fit_run = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='norm')
    fit_run.run(algorithm_list = ['PSO','PSO', 'PSO'], fitting_level=['norm','deep', 'deep'])
    fit_run.plot_final_qso_fit()
    fit_run.dump_result()

# #%% Read Best fitting result:
# fit_files = glob.glob('sim_result/sim_idx0_seed{0}B*.pkl'.format(seed))#+\
# fit_files.sort()
# fit_run_list = []
# for i in range(len(fit_files)):
#     fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
# chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
# sort_Chisq = chisqs.argsort() 
# best_sim_run = fit_run_list[sort_Chisq[0]]

