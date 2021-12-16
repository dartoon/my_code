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
from galight.tools.plot_tools import profile_plots, total_compare

#%%
"""
# NOTES:


#a Bulge (two comp) feature:
#b Bar feature:    
#c Spiral feature:
  
0 000557.23+002837.7 a  #Mask will be needed
4 010224.54-001501.5 a 
5 010335.35-005527.3 a
8 011304.19+001918.9 a
10 012059.12+003024.2 c
13 020301.63-033140.8 a
16 020427.68-043739.2 a
18 020559.66-063736.6 a
20 020852.34-063403.4 a
21 021122.59-021129.7 a
23 021235.98-053211.6 a  #Extended feature
25 021342.89-060758.4 a
26
27
28
"""  



#%%

fitsFile_ = glob.glob('SDSS_0.2-0.3/*_HSC-I.fits')

def condition_diskbulge(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
    logL = 0
    phi0, q0 = param_util.ellipticity2phi_q(kwargs_lens_light[0]['e1'], kwargs_lens_light[0]['e2'])
    phi1, q1 = param_util.ellipticity2phi_q(kwargs_lens_light[1]['e1'], kwargs_lens_light[1]['e2'])
    if (kwargs_lens_light[1]['R_sersic'] > kwargs_lens_light[0]['R_sersic'] * 0.9):
        logL -= 10**15
    if (kwargs_lens_light[1]['R_sersic'] < kwargs_lens_light[0]['R_sersic']*0.15):
        logL -= 10**15
    if (q1 < 0.9):
        logL -= 10**15
    return logL

def condition_diskbulgebar(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
    logL = 0
    phi0, q0 = param_util.ellipticity2phi_q(kwargs_lens_light[0]['e1'], kwargs_lens_light[0]['e2'])
    phi1, q1 = param_util.ellipticity2phi_q(kwargs_lens_light[1]['e1'], kwargs_lens_light[1]['e2'])
    phi2, q2 = param_util.ellipticity2phi_q(kwargs_lens_light[2]['e1'], kwargs_lens_light[2]['e2'])
    if (kwargs_lens_light[0]['R_sersic'] < kwargs_lens_light[1]['R_sersic']):
        logL -= 10**15
    if (kwargs_lens_light[0]['R_sersic']*0.15 > kwargs_lens_light[1]['R_sersic']):
        logL -= 10**15     
    if (kwargs_lens_light[0]['R_sersic'] > kwargs_lens_light[2]['R_sersic']):
        logL -= 10**15
    if (kwargs_lens_light[0]['R_sersic']*0.15 < kwargs_lens_light[2]['R_sersic']):
        logL -= 10**15 
    if (q0 < q2) or (q1 < q0):
        logL -= 10**15
    if abs(phi0 - (-apertures[0].theta)) > 5/180*np.pi:
        logL -= 10**15
    if abs(phi2 - (-apertures[2].theta)) > 5/180*np.pi:
        logL -= 10**15
    return logL

#%%
fitsFile_.sort()
# NO = 1
# for NO in [21]:
for NO in range(len(fitsFile_)):    
    fitsFile  = fitsFile_[NO]
    ID = fitsFile.split('/')[1].split('_')[0]
    print(NO, ID)
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
    
    # #Start to use galight
    # data_process_0 = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [float(QSO_RA), float(QSO_DEC)],
    #                            pos_type = 'wcs', header = header,
    #                           rm_bkglight = True, if_plot=False, zp = zp)
    
    # data_process_0.noise_map = err_data
    
    # data_process_0.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
    #                                       exp_sz= 1.2, npixels = 15, if_plot=True)
    # data_process_0.PSF_list = [PSF]
    # data_process_0.checkout() #Check if all the materials is known.
    # fit_sepc_0 = FittingSpecify(data_process_0)
    # fit_sepc_0.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    # fit_sepc_0.plot_fitting_sets()
    # fit_sepc_0.build_fitting_seq()
    # #Setting the fitting method and run.
    # fit_run_0 = FittingProcess(fit_sepc_0, savename = save_name+'single_Sersic')
    # fit_run_0.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 150}]*2)
    # fit_run_0.plot_final_qso_fit()
    # bic_0 = fit_run_0.fitting_seq.bic
    # fit_run_0.dump_result()
    
    # Read 1 fitting and check residual
    picklename = save_name+'single_Sersic.pkl'
    fit_run_0 = pickle.load(open(picklename,'rb'))
    # fit_run_0.fitting_specify_class.plot_fitting_sets()
    cutout_data = fit_run_0.fitting_specify_class.kwargs_data['image_data']
    noise = fit_run_0.fitting_specify_class.kwargs_data['noise_map']
    model = np.zeros_like(cutout_data)
    for i in range(len(fit_run_0.image_host_list)):
        model = model+fit_run_0.image_host_list[i]    
    for i in range(len(fit_run_0.image_ps_list)):
        model = model+fit_run_0.image_ps_list[i]    
    residual = (cutout_data-model)/noise
    removed = copy.deepcopy(cutout_data)
    removed = cutout_data-model + fit_run_0.image_host_list[0]  
    print(fit_run_0.reduced_Chisq)
    flux_list_2d = [cutout_data, model, removed, residual]
    label_list_2d = ['data', 'model', 'data-PS-nearby', 'normalized residual']
    flux_list_1d = [cutout_data, model, fit_run_0.image_ps_list[0], model - fit_run_0.image_ps_list[0]]
    label_list_1d = ['data', 'model', 'Point Source', '{0} galaxy(s)'.format(len(fit_run_0.image_host_list))] 
    fig = total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, deltaPix = fit_run_0.fitting_specify_class.deltaPix,
                  zp=fit_run_0.zp, target_ID=ID)
    plt.show()
    
    # #%%Fitting as disk + bulge:
    # data_process_1 = copy.deepcopy(data_process_0)
    # data_process_1.target_stamp -= fit_run_0.image_ps_list[0]
    # apertures = copy.deepcopy(data_process_1.apertures)
    # comp_id = 0 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    # add_aperture0 = copy.deepcopy(apertures[comp_id])
    # add_aperture0.a, add_aperture0.b = add_aperture0.a/3, add_aperture0.a/3
    # apertures = apertures[:comp_id+1] + [add_aperture0] + apertures[comp_id+1:]
    # data_process_1.apertures = apertures #Pass apertures to the data
    # fit_sepc_1 = FittingSpecify(data_process_1)
    # fit_sepc_1.prepare_fitting_seq(point_source_num = 0, fix_n_list= [[0,1]],  #First component fix n = 4 (bluge), second one fix to 1 (disk).
    #                               condition = condition_diskbulge)
    
    # fit_sepc_1.plot_fitting_sets()
    
    # picklename = save_name+'single_Sersic.pkl'
    # fit_run_0 = pickle.load(open(picklename,'rb'))
    
    # fit_sepc_1.build_fitting_seq()
    # fit_sepc_1.kwargs_params['lens_light_model'][0][0] = fit_run_0.final_result_galaxy[0] #Input to bulge
    
    # fit_run_1 = FittingProcess(fit_sepc_1, savename = save_name+'bulge+disk')
    
    # fit_run_1.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 150}]*2)
    # bic_1 = fit_run_1.fitting_seq.bic
    # fit_run_1.dump_result()
    # print(fit_run_1.final_result_galaxy[0]['q'], fit_run_1.final_result_galaxy[1]['q'])
    # fit_run_1.plot_final_galaxy_fit(target_ID =  ID)
    # bulge1 = fit_run_1.image_host_list[1]
    # disk1 = fit_run_1.image_host_list[0]
    # B2T = np.sum(bulge1)/np.sum(bulge1+disk1)
    # AGN1 = fit_run_0.image_ps_list[0]
    # disk_Re1 = fit_run_1.final_result_galaxy[0]['R_sersic']
    # bulge_Re1 = fit_run_1.final_result_galaxy[1]['R_sersic']
    # flux_list_2d = [bulge1, disk1, AGN1]
    # label_list_2d = ['Bulge', 'Disk', 'nuclei']
    # flux_list_1d = [bulge1, disk1, AGN1] 
    # label_list_1d = ['Bulge', 'Disk', 'nuclei']
    # profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
    #               deltaPix = fit_run_1.fitting_specify_class.deltaPix,
    #               target_ID =  ID, if_annuli=True)   
    # print("B/T:", round(B2T,2))

    # #%%Fitting as disk + bulge:
    # data_process_1b = copy.deepcopy(data_process_0)
    # data_process_1b.target_stamp -= fit_run_0.image_ps_list[0]
    # apertures = copy.deepcopy(data_process_1b.apertures)
    # comp_id = 0 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    # add_aperture0 = copy.deepcopy(apertures[comp_id])
    # add_aperture0.a, add_aperture0.b = add_aperture0.a/3, add_aperture0.a/9
    # add_aperture0.theta = 0
    # apertures = apertures[:comp_id+1] + [add_aperture0] + apertures[comp_id+1:]
    # data_process_1b.apertures = apertures #Pass apertures to the data
    # fit_sepc_1b = FittingSpecify(data_process_1b)
    # fit_sepc_1b.prepare_fitting_seq(point_source_num = 0, fix_n_list= [[0,1], [1,0.5]])  #First component fix n = 4 (bluge), second one fix to 1 (disk).
    #                               # condition = condition_diskbulge)
    # fit_sepc_1b.plot_fitting_sets()
    
    # picklename = save_name+'single_Sersic.pkl'
    # fit_run_0 = pickle.load(open(picklename,'rb'))
    
    # fit_sepc_1b.build_fitting_seq()
    # fit_sepc_1b.kwargs_params['lens_light_model'][0][0] = fit_run_0.final_result_galaxy[0] #Input to bulge
    
    # fit_run_1b = FittingProcess(fit_sepc_1b, savename = save_name+'bar+disk')
    
    # fit_run_1b.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 150}]*2)
    # bic_1b = fit_run_1b.fitting_seq.bic
    # fit_run_1b.dump_result()
    # print(fit_run_1b.final_result_galaxy[0]['q'], fit_run_1b.final_result_galaxy[1]['q'])
    # fit_run_1b.plot_final_galaxy_fit(target_ID =  ID)
    # bar1b = fit_run_1b.image_host_list[1]
    # disk1b = fit_run_1b.image_host_list[0]
    # B2T = np.sum(bulge1)/np.sum(bulge1+disk1)
    # AGN1b = fit_run_0.image_ps_list[0]
    # disk_Re1b = fit_run_1b.final_result_galaxy[0]['R_sersic']
    # bulge_Re1b = fit_run_1b.final_result_galaxy[1]['R_sersic']
    # flux_list_2d = [bar1b, disk1b, AGN1b]
    # label_list_2d = ['Bar', 'Disk', 'nuclei']
    # flux_list_1d = [bar1b, disk1b, AGN1b] 
    # label_list_1d = ['Bar', 'Disk', 'nuclei']
    # profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
    #               deltaPix = fit_run_1b.fitting_specify_class.deltaPix,
    #               target_ID =  ID, if_annuli=True)   
    # print("B/T:", round(B2T,2))

    # #%%Fitting as disk + bulge + bar:
    # data_process_2 = copy.deepcopy(data_process_0)
    # # fit_run_0 = pickle.load(open('fit_result/092256.42+042733.3_HSC-Isingle_Sersic.pkl','rb'))
    # data_process_2.target_stamp -= fit_run_0.image_ps_list[0]
    # apertures = copy.deepcopy(data_process_1.apertures)
    # comp_id = 1 #Change the component (galaxy) id = 0 into to components (i.e., bulge + disk)
    # add_aperture0 = copy.deepcopy(apertures[comp_id])
    # add_aperture0.a, add_aperture0.b = add_aperture0.a, add_aperture0.b/3
    # add_aperture0.theta = 0
    # apertures = apertures[:comp_id+1] + [add_aperture0] + apertures[comp_id+1:]
    
    # data_process_2.apertures = apertures #Pass apertures to the data
    # fit_sepc_2 = FittingSpecify(data_process_2)
    # fit_sepc_2.prepare_fitting_seq(point_source_num = 0, fix_n_list= [[0,1],[2,0.5]],  #First component fix n = 4 (bluge), second one fix to 1 (disk).
    #                               condition = condition_diskbulgebar)
    #                               # ,
    #                               # apertures_center_focus = True)
    # fit_sepc_2.build_fitting_seq()
    # # fit_sepc_2.kwargs_params['lens_light_model'][4][1]['R_sersic'] = 1.5
    # # fit_sepc_2.kwargs_params['lens_light_model'][0][1] = fit_run_1.final_result_galaxy[0]
    # # fit_sepc_2.kwargs_params['lens_light_model'][0][2] = fit_run_1.final_result_galaxy[1]
    # fit_sepc_2.plot_fitting_sets()
    # # fit_sepc_2.kwargs_params['lens_light_model'][0][0] = fit_run_1.final_result_galaxy[0] #Input to bulge
    # # fit_sepc_2.kwargs_params['lens_light_model'][0][1] = fit_run_1.final_result_galaxy[1] #Input to bulge
    # fit_run_2 = FittingProcess(fit_sepc_2, savename = save_name+'bar+bulge+disk')
    # fit_run_2.run(algorithm_list = ['PSO', 'PSO'], setting_list = [{'sigma_scale': 1., 'n_particles': 150, 'n_iterations': 150}]*2)
    # bic_2 = fit_run_2.fitting_seq.bic
    # fit_run_2.dump_result()
    # print("Fitting as Bar + Bulge + Disk")  
    # fit_run_2.plot_final_galaxy_fit(target_ID =  ID)
    # disk2 = fit_run_2.image_host_list[0]
    # bulge2 = fit_run_2.image_host_list[1]
    # bar2 = fit_run_2.image_host_list[2]
    # B2T = np.sum(bulge2)/np.sum(bulge2+disk2+bar2)
    # AGN2 = fit_run_0.image_ps_list[0]
    # disk_Re2 = fit_run_2.final_result_galaxy[0]['R_sersic']
    # bulge_Re2 = fit_run_2.final_result_galaxy[1]['R_sersic']
    # bar_Re2 = fit_run_2.final_result_galaxy[2]['R_sersic']
    # # bar2_resi = fit_run_2.fitting_specify_class.kwargs_data['image_data'] -AGN2
    # flux_list_2d = [disk2, bulge2, bar2, AGN2]
    # label_list_2d = ['Disk', 'Bulge', 'Bar', 'nuclei']
    # flux_list_1d = [disk2, bulge2, bar2, AGN2] 
    # label_list_1d = ['Disk', 'Bulge', 'Bar', 'nuclei']      
    # profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
    #               deltaPix = fit_run_0.fitting_specify_class.deltaPix,
    #               target_ID =  ID, if_annuli=True)
    
    # print("B/T:", round(B2T,2))
    # print('Reff:disk, bulge, bar: ', round(disk_Re2,2), round(bulge_Re2,2), round(bar_Re2) )
    # hold = input("hold:")

    # print("BIC compare:", round(bic_0,1), round(bic_1,1), round(bic_2,1) )
    # print("Chisq:", round(fit_run_0.reduced_Chisq,1), round(fit_run_1.reduced_Chisq,1),
    #         round(fit_run_2.reduced_Chisq,1))


    #%% 
    # print("BIC compare:", round(fit_run_0.fitting_seq.bic,1), round(bic_1,1), round(bic_2,1), )
    # print("Chisq:", round(fit_run_0.reduced_Chisq,1), round(fit_run_1.reduced_Chisq,1), round(fit_run_2.reduced_Chisq,1))
    # filename = 'fit_result_202110/' + 'BIC_result_2nd.txt'
    # if_file = glob.glob(filename)
    # if if_file == []:
    #     write_file =  open(filename,'w') 
    # else:
    #     write_file =  open(filename,'r+') 
    #     write_file.read()
    # if fit_run_0.fitting_seq.bic < fit_run_1.fitting_seq.bic:
    #     print(ID+" is a singel Sersic model; "+"glob Number: " + str(NO))
    #     write_file.write(ID+" is a singel Sersic model;" +"glob Number: " + str(NO))
    # else:
    #     print(ID+" is a disk+bulge!!! "+"glob Number: " + str(NO))
    #     write_file.write(ID+"is a disk+bulge!!! "+"glob Number: " + str(NO))
    # write_file.write("\n")
    # write_file.close()