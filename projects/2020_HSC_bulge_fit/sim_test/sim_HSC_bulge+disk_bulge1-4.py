#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 23:19:01 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from matplotlib.colors import LogNorm
from lenstronomy.ImSim.image_model import ImageModel
import lenstronomy.Util.param_util as param_util
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.psf import PSF
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LightModel.light_model import LightModel
import copy
from decomprofile.tools.plot_tools import profile_plots
import pickle
    
zp = 27.0 #
deltaPix = 0.167 #arcsec

def condition_bulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
    logL = 0
    phi0, q0 = param_util.ellipticity2phi_q(kwargs_source[0]['e1'], kwargs_source[0]['e2'])
    phi1, q1 = param_util.ellipticity2phi_q(kwargs_source[1]['e1'], kwargs_source[1]['e2'])
    cond_0 = (kwargs_source[0]['R_sersic'] > kwargs_source[1]['R_sersic'] * 0.9)
    cond_1 = (kwargs_source[0]['R_sersic'] < kwargs_source[1]['R_sersic']*0.15)
    cond_2 = (q0 < q1)
    if cond_0 or cond_1 or cond_2:
        logL -= 10**15
    return logL

array_l_means = ['id', 'z', 'ra', 'dec', 'fix_sersic_n', 'sersic_n_fitted', 'sersic_re_fitted', 'sersic_n_corrected',
         'sersic_re_corrected', 'host_mag_g', 'host_mag_r', 'host_mag_i', 'host_mag_z', 'host_mag_y',
         'ps_mag_g', 'ps_mag_r', 'ps_mag_i', 'ps_mag_z', 'ps_mag_y', 'decomposition_chisq', 'stellar_mass', 
         'sed_chisq', 'logMBH', 'logMBH_err']
infers  = np.loadtxt('../../2020_HSC_to_numerical_simulation/HSC_fitting/sdss_quasar_decomposition_v1.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_ = infers[:,1].astype(np.float)
HSC_host_mag_ = infers[:,11].astype(np.float)
HSC_AGN_mag_ = infers[:,16].astype(np.float)
HSC_re_ = infers[:,6].astype(np.float)

flags_  = np.loadtxt('../../2020_HSC_to_numerical_simulation/HSC_fitting/sdss_quasar_decomposition_v1_catalog_flag.txt', dtype=str)
flags = flags_[:,0]
IDs, HSC_z, HSC_AGN_mag, HSC_host_mag, HSC_re=[], [], [], [], []
for i in range(len(IDs_)):
    idx = np.where(IDs_[i] == flags)[0][0]
    if flags_[idx][1] == 'y' and HSC_AGN_mag_[i]>0 and HSC_z_[i]<0.3:
        IDs.append(IDs_[i])
        HSC_z.append(HSC_z_[i])
        HSC_AGN_mag.append(HSC_AGN_mag_[i])
        HSC_host_mag.append(HSC_host_mag_[i])
        HSC_re.append(HSC_re_[i])

HSC_host_mag = np.array(HSC_host_mag)
HSC_AGN_mag = np.array(HSC_AGN_mag)
HSC_re = np.array(HSC_re)

#%%
bulge_n = 4
for j in range(78*4):
    i = j%78
    host_mag = HSC_host_mag[i] 
    AGN_mag = HSC_AGN_mag[i]
    disk_reff =  HSC_re[i]*1.2 #1 # 
    B2T = np.random.uniform(0.3,0.7) #Bulge to total ratio
    bulge_reff = np.max([ disk_reff * np.random.uniform(0.2,0.5) , 0.2 ])
    bulge_reff = np.min([ bulge_reff , disk_reff])
    host_flux = 10**(0.4*(zp - host_mag))
    bulge_flux = host_flux*B2T
    disk_flux = host_flux* (1-B2T)
    AGN_flux = 10**(0.4*(zp - AGN_mag))
    psf_take_id = 0
    stdd = 0.04  #!!!
    exptim= 2.1 * 60 * 60  #units of seconds, assuming deep.  https://hsc.mtk.nao.ac.jp/ssp/survey/
    files = glob.glob("../SDSS_0.2-0.3/*_psf.fits")
    psf_data = pyfits.getdata(files[psf_take_id])
    if len(psf_data) != 0 and psf_data.shape[0] != psf_data.shape[1]:
        cut = int((psf_data.shape[0] - psf_data.shape[1])/2)
        if cut>0:
            psf_data = psf_data[cut:-cut,:]
        elif cut<0:
            psf_data = psf_data[:,-cut:cut]
        psf_data /= psf_data.sum()
    #%%
    numPix = 89
    light_model_list = ['SERSIC_ELLIPSE']
    kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': deltaPix}
    kwargs_data_high_res = sim_util.data_configure_simple(numPix, deltaPix)
    data_class = ImageData(**kwargs_data_high_res)
    psf_class = PSF(**kwargs_psf_high_res)
    center_x, center_y = np.random.uniform(-0.2, 0.2) * deltaPix, np.random.uniform(-0.2,0.2)* deltaPix
    point_amp = AGN_flux
    
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
    light_model_list = ['SERSIC_ELLIPSE']
    lightModel = LightModel(light_model_list=light_model_list)
    
    q0 = np.random.uniform(0.5,0.9)
    phi0 = np.random.uniform(0.,2*np.pi)            
    q1 = np.random.uniform(0.5,0.9)
    phi1 = np.random.uniform(0.,2*np.pi)            
    e1_0, e2_0 = param_util.phi_q2_ellipticity(phi=phi0, q=q0)
    e1_1, e2_1 = param_util.phi_q2_ellipticity(phi=phi1, q=q1)
    
    kwargs_numerics = {'supersampling_factor': 3, 'supersampling_convolution': False}
    kwargs_bulge = {'amp': 1. , 'n_sersic': bulge_n, 'R_sersic': bulge_reff, 'e1': e1_0, 'e2': e2_0,
                      'center_x': center_x + np.random.uniform(-0.1, 0.1)*deltaPix,
                      'center_y': center_y + np.random.uniform(-0.1, 0.1)*deltaPix} #!!!
    
    kwargs_disk = {'amp': 1. , 'n_sersic': 1, 'R_sersic': disk_reff, 'e1': e1_1, 'e2': e2_1,
                      'center_x': center_x + np.random.uniform(-0.1, 0.1)*deltaPix,
                      'center_y': center_y + np.random.uniform(-0.1, 0.1)*deltaPix} #!!!
    
    # kwargs_host_medi = [kwargs_bulge, kwargs_disk]
    imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
    
    medi_bluge_flux = np.sum(imageModel.image(kwargs_lens_light=[kwargs_bulge], unconvolved=True))
    medi_disk_flux = np.sum(imageModel.image(kwargs_lens_light=[kwargs_disk], unconvolved=True))
    
    kwargs_bulge['amp'] = 1. / medi_bluge_flux * bulge_flux
    kwargs_disk['amp'] = 1. / medi_disk_flux * disk_flux
    
    # ## simulate image with the parameters we have defined above #
    bulge_image = imageModel.image(kwargs_lens_light=[kwargs_bulge], 
                                   kwargs_ps= [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [0]}], 
                                   unconvolved=False)
    disk_image = imageModel.image(kwargs_lens_light=[kwargs_disk], 
                                   kwargs_ps= [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [0]}], 
                                   unconvolved=False)
    disk_PS_image = imageModel.image(kwargs_lens_light=[kwargs_disk], kwargs_ps=kwargs_ps, unconvolved=False)
    #%%
    AGN_image = disk_PS_image - disk_image
    sim_image = bulge_image + disk_PS_image
    kwargs_bulge['flux_within_frame'] =  np.sum(bulge_image)
    kwargs_disk['flux_within_frame'] =  np.sum(disk_image)
    
    # sim_image = disk_image
    # plt.imshow(bulge_image, origin='lower',cmap='gist_heat', norm=LogNorm())
    # plt.colorbar()
    # plt.show()
    rms=(sim_image/(exptim)+stdd**2)**0.5 #RMS not saved, thus its value not used here
    bkg_noise= stdd
    noiz=np.random.normal(0, bkg_noise, size=rms.shape)
    sim_image[sim_image<0] = 0
    sim_image_noise= noiz + np.random.poisson(lam=sim_image*exptim)/(exptim)  #Non-drizzled imaged
    plt.imshow(sim_image_noise, origin='lower',cmap='gist_heat', norm=LogNorm())
    plt.colorbar()
    plt.show()
    
    flux_list_2d = [bulge_image, disk_image, AGN_image]
    label_list_2d = ['Bulge', 'Disk', 'nuclei']
    flux_list_1d = [bulge_image, disk_image, AGN_image]
    label_list_1d = ['Bulge', 'Disk', 'nuclei']
    
    profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
                  deltaPix = deltaPix,
                  target_ID = 'Simulation', if_annuli=True)
    #%%Model:
    from decomprofile.data_process import DataProcess
    from decomprofile.fitting_specify import FittingSpeficy
    from decomprofile.fitting_process import FittingProcess
    
    save_name = 'sim_result_bulge_n{1}/round0_ID{0}_'.format(j,bulge_n)
    data_process_0 = DataProcess(fov_image = sim_image_noise, fov_noise_map = rms, 
                                  target_pos = [len(sim_image_noise)/2, len(sim_image_noise)/2],
                                  pos_type = 'pixel', header = None,
                                  rm_bkglight = False, if_plot=False, zp = zp)
    data_process_0.deltaPix = deltaPix
    data_process_0.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 15, if_plot=False)
    
    data_process_0.PSF_list = [psf_data]
    data_process_0.checkout() #Check if all the materials is known.
    #%%Start to produce the class and params for lens fitting.
    # #Manually input another component:
    # apertures_0 = copy.deepcopy(data_process_0.apertures)
    # add_aperture1 = copy.deepcopy(apertures_0[0])
    # add_pos = [60, 60]   #The position of the component.
    # if isinstance(add_aperture1.positions[0],float): 
    #     add_aperture1.positions = np.array(add_pos)
    # elif isinstance(add_aperture1.positions[0],np.ndarray):
    #     add_aperture1.positions = np.array([add_pos])
    # add_aperture1.a, add_aperture1.b = 2, 2  #define the a, b value of this component, i.e., Reff = sqrt(a^2 +b^2)
    # apertures_0 = apertures_0 + [add_aperture1]  #attach the added aperture into the group.    
    # data_process_0.apertures = apertures_0 #Pass apertures to the data
    fit_sepc_0 = FittingSpeficy(data_process_0)
    fit_sepc_0.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    # fit_sepc_0.plot_fitting_sets()
    fit_sepc_0.build_fitting_seq()
    #Setting the fitting method and run.
    fit_run_0 = FittingProcess(fit_sepc_0, savename = save_name+'single_Sersic')
    fit_run_0.run(algorithm_list = ['PSO'], setting_list = [None])
    fit_run_0.plot_final_qso_fit()
    bic_0 = fit_run_0.fitting_seq.bic
    fit_run_0.dump_result()
    
    #%%Fitting as disk + bulge:
    data_process_1 = copy.deepcopy(data_process_0)
    apertures = copy.deepcopy(data_process_1.apertures)
    comp_id = 0 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    add_aperture0 = copy.deepcopy(apertures[comp_id])
    add_aperture0.a, add_aperture0.b = add_aperture0.a/3, add_aperture0.b/3
    apertures = apertures[:comp_id] + [add_aperture0] + apertures[comp_id:]
    data_process_1.apertures = apertures #Pass apertures to the data
    fit_sepc_1 = FittingSpeficy(data_process_1)
    fit_sepc_1.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[1,1]],  #First component fix n = 4 (bluge), second one fix to 1 (disk).
                                  fix_center_list = [[0,0]], condition = condition_bulgedisk)
    # fit_sepc_1.kwargs_params['source_model']
    use_true = 0
    if use_true == 0:
        fit_sepc_1.plot_fitting_sets()
        fit_sepc_1.build_fitting_seq()
        fit_sepc_1.kwargs_params['source_model'][0][1] = fit_run_0.final_result_galaxy[0]
        fit_sepc_1.kwargs_params['source_model'][4][0]['R_sersic'] = fit_run_0.final_result_galaxy[0]['R_sersic']
        fit_sepc_1.kwargs_params['source_model'][4][1]['R_sersic'] = fit_run_0.final_result_galaxy[0]['R_sersic']*2
    elif use_true == 1:
        fit_sepc_1.source_params = [kwargs_bulge, kwargs_disk] #!!!
        
    fit_run_1 = FittingProcess(fit_sepc_1, savename = save_name+'bulge+disk')
    fit_run_1.run(algorithm_list = ['PSO'],
                  setting_list = [ None ])
    # fit_run_1.run(algorithm_list = ['PSO']*5, 
    #               setting_list= [ {'sigma_scale': 0.8, 'n_particles': 300, 
    #                                 'n_iterations': 300}] *5)
    fit_run_1.plot_final_qso_fit()
    bic_1 = fit_run_1.fitting_seq.bic
    fit_run_1.dump_result()
    
    bulge = fit_run_1.image_host_list[0]
    disk = fit_run_1.image_host_list[1]
    B2T_inf = np.sum(bulge)/np.sum(bulge+disk)
    AGN = fit_run_1.image_ps_list[0]
    
    bulge_Re_inf = fit_run_1.final_result_galaxy[0]['R_sersic']
    disk_Re_inf = fit_run_1.final_result_galaxy[1]['R_sersic']
    
    flux_list_2d = [bulge, disk, AGN]
    label_list_2d = ['Bulge', 'Disk', 'nuclei']
    flux_list_1d = [bulge, disk, AGN] 
    label_list_1d = ['Bulge', 'Disk', 'nuclei']
    
    profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
                  deltaPix = fit_run_1.fitting_specify_class.deltaPix,
                  target_ID =  "Fitting inference", if_annuli=True)
    print("BIC compare:", round(bic_0,1), round(bic_1,1))
    print("One VS two Sersics Chisq:", round(fit_run_0.reduced_Chisq,1), round(fit_run_1.reduced_Chisq,1))
    print("True VS inf B/T:", round(B2T,2), round(B2T_inf,2))
    print('True Reff: bulge, disk: ', round(bulge_reff,2), round(disk_reff,2) )
    print('Inf Reff: bulge, disk: ', round(bulge_Re_inf,2), round(disk_Re_inf,2) )

    picklename = save_name + 'QuickResult.pkl'
    bic_result = ['bic_0 VS bic_1', bic_0, bic_1]
    chisq_result = ['Chisq0 VS Chisq1', fit_run_0.reduced_Chisq, fit_run_1.reduced_Chisq]
    B2T_result = ['B/T True VS inf', B2T, B2T_inf ]
    True_param = ['kwargs_bulge, kwargs_disk, kwargs_ps', kwargs_bulge, kwargs_disk, kwargs_ps]
    infer_param_single = ['inferred singel sersic: kwargs_sersic, kwargs_ps', fit_run_0.final_result_galaxy[0], fit_run_0.final_result_ps] 
    infer_param_DB = ['inferred bulge disk: kwargs_bulge, kwargs_disk, kwargs_ps', fit_run_1.final_result_galaxy[0], 
                      fit_run_1.final_result_galaxy[1], fit_run_1.final_result_ps] 
    savepkl = [bic_result, chisq_result, B2T_result, True_param, infer_param_single, infer_param_DB]
    pickle.dump(savepkl, open(picklename, 'wb'))
    
    # if_file = glob.glob(filename)
    # if if_file == []:
    #     write_file =  open(filename,'w') 
    # else:
    #     write_file =  open(filename,'r+') 
    #     write_file.read()
    # write_file.write("#################################\n")
    # write_file.write("BIC compare:{0} {1}\n".format(round(bic_0,1), round(bic_1,1)) )
    # write_file.write("One VS two Sersics Chisq:{0} {1}\n".format(round(fit_run_0.reduced_Chisq,1), round(fit_run_1.reduced_Chisq,1)))
    # write_file.write("True VS inf B/T:{0} {1}\n".format(round(B2T,2), round(B2T_inf,2)))
    # write_file.write("True Reff: bulge, disk:{0} {1}\n".format(round(bulge_reff,2), round(disk_reff,2)))
    # write_file.write("Inf Reff: bulge, disk:{0} {1}\n".format(round(bulge_Re_inf,2), round(disk_Re_inf,2)))    
    # write_file.close()
    
