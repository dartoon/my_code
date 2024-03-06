#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:49:21 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os, glob
import shutil 
import pickle

# name_list = {1:'SDSS1420+5300A', 2: 'SDSS1420+5300B', 0: 'SDSS1419+5254', 
#               51: 'AEGIS 477', 35: 'AEGIS 482'}
name_list = {1:'SDSS1420A', 2: 'SDSS1420B', 0: 'SDSS1419', 
              51: 'AEGIS 477', 35: 'AEGIS 482'}

def load_info(idx):
    f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
    
    f1 = open("../model_JWST_stage3_Masafusa/material/target_info.txt","r")
    string1 = f1.read()
    lines1 = string1.split('\n')   # Split in to \n
    spec_z = [lines1[i].split(' ')[3] for i in range(1, len(lines1)) if lines1[i].split(' ')[0] == target_id][0]
    return target_id, spec_z

def load_prop(idx, root_folder = './', if_plot=False, prop_name = 'magnitude', path = None):
    result = {}
    files = glob.glob(root_folder+'*fit_material/data_process_idx{0}_*_psf*.pkl'.format(idx))
    files.sort()
    _collect_info = []
    for i in range(len(files)):
        _file = files[i]
        idx_info = _file.split('idx')[1].split('_')[0]
        filt_info = _file.split('_psf')[0].split('_')[-1]
        this_info = [idx_info, filt_info]
        if this_info not in _collect_info:
            _collect_info.append(this_info)
    
    filters = [_collect_info[i][1] for i in range(len(_collect_info))]
    if 'F814W' in filters: #Move F606W and F814W to top
        filters = ['F814W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F814W' ]
    if 'F606W' in filters:
        filters = ['F606W'] + [filters[i] for i in range(len(filters)) if filters[i]!='F606W' ]
        
    if path is None:
        f = open("../model_JWST_stage3_Masafusa/target_idx_info.txt","r")
    else:
        f = open(path,"r")
    string = f.read()
    lines = string.split('\n')   # Split in to \n
    target_id = [lines[i].split(' ')[1] for i in range(len(lines)) if lines[i].split(' ')[0] == str(idx)][0]
    for count in range(len(filters)):
        fit_run_list = []
        idx = idx_info
        filt = filters[count]
        # idx, filt= item
        fit_files = glob.glob(root_folder+'*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))+\
                    glob.glob('../model_HST_material/fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))
        fit_files.sort()
        warn_strs = ['F115W_psf6', 'F150W_psf7', 'F277W_psf2']
        for warn_str in warn_strs:
            fit_files = [fit_files[i] for i in range(len(fit_files)) if warn_str not in fit_files[i]]
        
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, filt, "Total PSF NO.", len(sort_Chisq))
        if 'HST_material' in fit_files[0]:
            weight = np.zeros(len(chisqs))
            weight[sort_Chisq[0]] = 1
        elif 'JWST_stage3' in fit_files[0]:
            count_n = 5
            Chisq_best = chisqs[sort_Chisq[0]]
            Chisq_last= chisqs[sort_Chisq[count_n-1]]
            inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
            weight = np.zeros(len(chisqs))
            for i in sort_Chisq[:count_n]:
                weight[i] = np.exp(-1/2. * (chisqs[i]-Chisq_best)/(Chisq_best* inf_alp))
        fit_run = fit_run_list[sort_Chisq[0]]
        if idx == '35' and filt == 'F444W':  #!!! Manually use one.
            fit_run = fit_run_list[sort_Chisq[3]]
        if idx == '0' and filt == 'F444W':  #!!!
            fit_run = fit_run_list[sort_Chisq[1]]
        if if_plot == True:
            fit_run.plot_final_qso_fit(target_ID = target_id+'$-$'+filt)
        # # prop_name = 'R_sersic'
        # all_values = [fit_run_list[i].final_result_galaxy[0][prop_name] for i in range(len(fit_run_list))]
        # weighted_value = np.sum(np.array(all_values)*weight) / np.sum(weight)
        # rms_value = np.sqrt(np.sum((np.array(all_values)-weighted_value)**2*weight) / np.sum(weight))
        # # result.append([filt, fit_run.fitting_specify_class.zp, weighted_value, rms_value])
        # host_flux = fit_run.final_result_galaxy[0]['flux_within_frame']
        # AGN_flux = fit_run.final_result_ps[0]['flux_within_frame']
        # ratio = host_flux/(host_flux+AGN_flux)
        if prop_name == 'magnitude':
            mag_correct = 0
            if filt == 'F444W':
                if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000: #module A
                    correct = 0.44157708 / 0.343
                    mag_correct = +2.5*np.log10(correct)
                if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000: #module B
                    correct = 0.3899884 / 0.335
                    mag_correct = +2.5*np.log10(correct)
            elif filt == 'F410M':
                if fit_run.fitting_specify_class.data_process_class.target_pos[0] < 5000:
                    correct = 0.9355298 / 0.832
                    mag_correct = +2.5*np.log10(correct)
                if fit_run.fitting_specify_class.data_process_class.target_pos[0] > 5000:
                    correct = 0.9272488 / 0.811
                    mag_correct = +2.5*np.log10(correct)
            result[filt] = fit_run.final_result_galaxy[0]['magnitude'] + mag_correct
        elif prop_name == 'chisq':
            result[filt] = fit_run.reduced_Chisq
        elif prop_name == 'fit_run':
            result[filt] = fit_run
        else:
            result[filt] = fit_run.final_result_galaxy[0][prop_name]
    return result

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362', 'F090W': '351'}

def esti_smass(ID, mags_dict, z, folder = 'esti_smass/', flag = 0, if_run_gsf=True, band_as_upper = [],
               mag_err = [], just_run = False, metallicity = 0.0):
    from gsf import gsf
    if just_run == False:
        ID = ID
        z = z
        #%reate a cat file
        folder_path = folder + ID + '/'
        if glob.glob(folder_path) != []:   
            shutil.rmtree(folder_path)
        os.mkdir(path = folder_path)
        # text_temp = "# id F352 E352 F353 E353 F354 E354 F355 E355 F356 E356 F362 E362 F357 E357\n"
        text_temp = "# id "
        mags = []
        filterIDs=''
        filt_id  = jwst_filt_id | hst_filt_id
        # print(filt_id)
        if_hst = []
        for key in mags_dict.keys():
            text_temp = text_temp + ' F{0} E{0}'.format(filt_id[key])
            filterIDs = filterIDs + ','+filt_id[key]
            mags.append(mags_dict[key])
            if key in hst_filt_id.keys() or key in band_as_upper:
                if_hst.append(True)
            else:
                if_hst.append(False)
        filterIDs = filterIDs[1:]
        text_temp = text_temp + "\n"
        if mag_err == []:
            mag_err = [0.2] * len(mags)
        fnu = [10 ** ((mags[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_up = [10 ** ((mags[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_dw = [10 ** ((mags[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mags))]
            
        # for i in range(7):
        #     if mags[i] <0:
        #         fnu[i] = 100
        #         fnu_err[i] = 1000000
        write_file = open(folder_path+'sample.cat','w') 
        write_file.write(text_temp)
        # _string = str(int(ID)) + " {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}".format(round(fnu[0],3), round(fnu_err[0],3), round(fnu[1],3), round(fnu_err[1],3), round(fnu[2],3), round(fnu_err[2],3),
        #       round(fnu[3],3), round(fnu_err[3],3), round(fnu[4],3), round(fnu_err[4],3), round(fnu[5],3), round(fnu_err[5],3), round(fnu[6],3), round(fnu_err[6],3)) 
        _string = str(int(ID))
        for i in range(len(fnu)):
            _string = _string + " {0:.8f} {1:.8f}".format(fnu[i], fnu_err[i])
        # for i in range(len(fnu)):
        #     if if_hst[i] == False:
        #         _string = _string + " {0:.8f} {1:.8f}".format(fnu[i], fnu_err[i])
        #     else:
                # _string = _string + " {0:.8f} {1:.8f}".format(0, fnu[i]/2)
        write_file.write(_string)
        write_file.close()
    
        #Create a input file
        f = open("../../../../template/gsf_JWSTNIRCam_temp/sample_template.input","r")
        string = f.read()
        string = string.replace("idname", str(int(ID)))
        string = string.replace("zinfo", str(z))
        string = string.replace("folder/", folder_path)
        string = string.replace("filterIDs", filterIDs)
        string = string.replace("metaltemp", str(metallicity))
        write_file = open(folder_path+'sample.input','w') 
        write_file.write(string)
        write_file.close()
        
    if if_run_gsf == True:
        gsf.run_gsf_all(folder+'{0}/sample.input'.format(ID), flag, idman=None)
        gsf.run_gsf_all(folder+'{0}/sample.input'.format(ID), 6, idman=None)
        #Move things in position.
        mv_files = glob.glob('*_{0}_*'.format(ID)) + glob.glob('*_{0}.*'.format(ID))
        for mv_file in mv_files:
            shutil.move(mv_file, folder+'{0}/'.format(ID)+mv_file)
            
            
            
            