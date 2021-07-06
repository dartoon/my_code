#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 14:54:38 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import copy
import glob
import lenstronomy.Util.param_util as param_util
import pickle

fitsFile_ = glob.glob('SDSS_0.2-0.3/*_HSC-I.fits')

def condition_bulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
    logL = 0
    phi0, q0 = param_util.ellipticity2phi_q(kwargs_lens_light[0]['e1'], kwargs_lens_light[0]['e2'])
    phi1, q1 = param_util.ellipticity2phi_q(kwargs_lens_light[1]['e1'], kwargs_lens_light[1]['e2'])
    cond_0 = (kwargs_lens_light[0]['R_sersic'] > kwargs_lens_light[1]['R_sersic'] * 0.9)
    cond_1 = (kwargs_lens_light[0]['R_sersic'] < kwargs_lens_light[1]['R_sersic']*0.15)
    cond_2 = (q0 < q1)
    if cond_0 or cond_1 or cond_2:
        logL -= 10**15
    return logL

def condition_barbulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
    logL = 0
    phi0, q0 = param_util.ellipticity2phi_q(kwargs_lens_light[0]['e1'], kwargs_lens_light[0]['e2'])
    phi1, q1 = param_util.ellipticity2phi_q(kwargs_lens_light[1]['e1'], kwargs_lens_light[1]['e2'])
    phi2, q2 = param_util.ellipticity2phi_q(kwargs_lens_light[2]['e1'], kwargs_lens_light[2]['e2'])
    cond_0 = (kwargs_lens_light[1]['R_sersic'] > kwargs_lens_light[2]['R_sersic'] * 0.9)
    cond_1 = (kwargs_lens_light[1]['R_sersic'] < kwargs_lens_light[2]['R_sersic']*0.15)
    cond_2 = (q1 < q2)
    cond_3 = (q2 < q0)
    if cond_0 or cond_1 or cond_2 or cond_3:
        logL -= 10**15
    return logL


# NO = 1
# for NO in range(21):
# for NO in range(21,42):
# for NO in range(0,63):
# for NO in [8]: 
for NO in [13]:     
# for NO in range(63, len(fitsFile_)):    
# for NO in [66]:    
    fitsFile  = fitsFile_[NO]
    ID = fitsFile.split('/')[1].split('_')[0]
    PSF_filename = fitsFile.split('.fits')[0]+ '_psf.fits'
    save_name = 'fit_result/' + fitsFile.split('.fits')[0].split('/')[1]
    f = open("SDSS_0.2-0.3/SDSS_0.2-0.3.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    line = [i for i in range(len(lines)) if lines[i].split(' ')[0] == ID][0]
    QSO_RA, QSO_DEC = lines[line].split(' ')[1:]
    
    
    fits = pyfits.open(fitsFile)
    fov_image= fits[1].data
    header = fits[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    err_data= fits[3].data ** 0.5
    
    file_header0 = fits[0].header
    zp = 27.0# This is something Xuheng can't make sure.
    PSF = pyfits.getdata(PSF_filename)
    
    #Start to use galight
    data_process_0 = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [float(QSO_RA), float(QSO_DEC)],
                               pos_type = 'wcs', header = header,
                              rm_bkglight = True, if_plot=False, zp = zp)
    
    data_process_0.noise_map = err_data
    
    data_process_0.generate_target_materials(radius=65, create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 15, if_plot=True)
    
    data_process_0.apertures = [data_process_0.apertures[0]]
    
    data_process_0.PSF_list = [PSF]
    
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
    fit_sepc_0 = FittingSpecify(data_process_0)
    fit_sepc_0.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    fit_sepc_0.plot_fitting_sets()
    fit_sepc_0.build_fitting_seq()
    #Setting the fitting method and run.
    fit_run_0 = FittingProcess(fit_sepc_0, savename = save_name+'single_Sersic')
    fit_run_0.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}]*2)
    fit_run_0.plot_final_qso_fit()
    bic_0 = fit_run_0.fitting_seq.bic
    fit_run_0.dump_result()
    
    
    #%%Fitting as disk + bulge:
    data_process_1 = copy.deepcopy(data_process_0)
    apertures = copy.deepcopy(data_process_1.apertures)
    comp_id = 0 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    add_aperture0 = copy.deepcopy(apertures[comp_id])
    add_aperture0.a, add_aperture0.b = add_aperture0.a/2, add_aperture0.a/2
    apertures = apertures[:comp_id] + [add_aperture0] + apertures[comp_id:]
    data_process_1.apertures = apertures #Pass apertures to the data
    fit_sepc_1 = FittingSpecify(data_process_1)
    fit_sepc_1.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[1,1]],  #First component fix n = 4 (bluge), second one fix to 1 (disk).
                                  fix_center_list = [[0,0]], condition = condition_bulgedisk)
    
    fit_sepc_1.plot_fitting_sets()
    
    picklename = save_name+'single_Sersic.pkl'
    fit_run_0 = pickle.load(open(picklename,'rb'))
    
    fit_sepc_1.build_fitting_seq()
    fit_sepc_1.kwargs_params['lens_light_model'][0][1] = fit_run_0.final_result_galaxy[0]
    fit_sepc_1.kwargs_params['lens_light_model'][4][0]['R_sersic'] = fit_run_0.final_result_galaxy[0]['R_sersic']
    fit_sepc_1.kwargs_params['lens_light_model'][4][1]['R_sersic'] = fit_run_0.final_result_galaxy[0]['R_sersic']*2
    
    
    fit_run_1 = FittingProcess(fit_sepc_1, savename = save_name+'bulge+disk')
    
    fit_run_1.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}]*2)
    fit_run_1.plot_final_qso_fit()
    bic_1 = fit_run_1.fitting_seq.bic
    fit_run_1.dump_result()

    #%%Fitting as disk + bulge + bar:
    data_process_2 = copy.deepcopy(data_process_1)
    apertures = copy.deepcopy(data_process_2.apertures)
    comp_id = 0 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    add_aperture0 = copy.deepcopy(apertures[comp_id])
    add_aperture0.a, add_aperture0.b = add_aperture0.a, add_aperture0.b/3
    apertures = apertures[:comp_id] + [add_aperture0] + apertures[comp_id:]
    data_process_2.apertures = apertures #Pass apertures to the data
    fit_sepc_2 = FittingSpecify(data_process_2)
    fit_sepc_2.prepare_fitting_seq(point_source_num = 1, fix_n_list= [[2,1], [0,0.5]],  #First component fix n = 4 (bluge), second one fix to 1 (disk).
                                  fix_center_list = [[0,1]], condition = condition_barbulgedisk)
    
    fit_sepc_2.plot_fitting_sets()
    
    picklename = save_name+'bulge+disk.pkl'
    fit_run_2 = pickle.load(open(picklename,'rb'))
    
    fit_sepc_2.build_fitting_seq()
    fit_sepc_2.kwargs_params['lens_light_model'][0][1] = fit_run_1.final_result_galaxy[0]
    fit_sepc_2.kwargs_params['lens_light_model'][0][2] = fit_run_1.final_result_galaxy[1]
    
    fit_run_2 = FittingProcess(fit_sepc_2, savename = save_name+'bar+bulge+disk')
    
    fit_run_2.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}]*2)
    fit_run_2.plot_final_qso_fit()
    bic_2 = fit_run_2.fitting_seq.bic
    fit_run_2.dump_result()

    #%% 
    print("BIC compare:", round(fit_run_0.fitting_seq.bic,1), round(bic_1,1), round(bic_2,1), )
    print("Chisq:", round(fit_run_0.reduced_Chisq,1), round(fit_run_1.reduced_Chisq,1), round(fit_run_2.reduced_Chisq,1))
    filename = 'fit_result/' + 'BIC_result_2nd.txt'
    if_file = glob.glob(filename)
    if if_file == []:
        write_file =  open(filename,'w') 
    else:
        write_file =  open(filename,'r+') 
        write_file.read()
    if fit_run_0.fitting_seq.bic < fit_run_1.fitting_seq.bic:
        print(ID+" is a singel Sersic model; "+"glob Number: " + str(NO))
        write_file.write(ID+" is a singel Sersic model;" +"glob Number: " + str(NO))
    else:
        print(ID+" is a disk+bulge!!! "+"glob Number: " + str(NO))
        write_file.write(ID+"is a disk+bulge!!! "+"glob Number: " + str(NO))
    write_file.write("\n")
    write_file.close()