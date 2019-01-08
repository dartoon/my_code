#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:14 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from gen_fit_ingredient import gen_fit_ingredient
from photutils import make_source_mask
import os

from fit_qso import fit_qso
from mask_objects import detect_obj
from matplotlib.colors import LogNorm
import copy
from subprocess import call

deep_seed = True  #Set as True to put more seed and steps to fit,

for seq in range(1,59):
    ID_seq = seq
    file_name_seq = [ i + (ID_seq-1)*5 + 2 for i in range(5)]  #The 'NO.' in the filename.
    band_seq = ['G', 'R', 'I', 'Z', 'Y']
    QSO_im_list, err_map_list, PSF_list=[], [], []
    zp_list, filename_list = [], []
    qso_center_list, frm_c_RA_DEC_list = [], []
    
    for i in range(len(band_seq)):
        # The pixel scale is all 0.168
        file_seq = file_name_seq[i]-1   # if count from 2, i.e. 2-*.fits is file_seq 1
        QSO_im, err_map, PSF, pix_scale, zp, qso_ID, qso_fr_center, fr_c_RA_DEC, filename = gen_fit_ingredient(file_seq, cut_frame=120)
        if '-{0}-'.format(band_seq[i]) not in filename:
            raise ValueError("The band information is not matched.")
        QSO_im_list.append(QSO_im)
        err_map_list.append(err_map)
        qso_center_list.append(qso_fr_center)
        frm_c_RA_DEC_list.append(fr_c_RA_DEC)
        if PSF.shape[0] != PSF.shape[1]:
            cut = ((PSF.shape[0] - PSF.shape[1])/2)
            if cut>0:
                PSF = PSF[cut:-cut,:]
            elif cut<0:
                PSF = PSF[:,-cut:cut]
            PSF /= PSF.sum()
            if PSF.shape[0] != PSF.shape[1]:
                raise ValueError("PSF shape is not a square.")
        PSF_list.append(PSF)
        zp_list.append(zp)
        if i ==0:
            qso_ID_last = qso_ID
        if i>0 and qso_ID_last != qso_ID:
            raise ValueError("The QSO ID is not consistent.")
        qso_ID_last = qso_ID
        filename_list.append(filename.split('/')[-1])
    
    #==============================================================================
    # quickly evaluate the background rms
    #==============================================================================
    background_rms_list = []
    for i in range(len(band_seq)):
        mask = make_source_mask(QSO_im_list[i], snr=3, npixels=5, dilate_size=11)
        background_rms_list.append(np.std(QSO_im_list[i]* (1-mask*1)))
    #    plt.imshow(QSO_im_list[i]*(1-mask*1), origin='low')
    #    plt.show()
    
    fit_frame_size = 81
    ct = (len(QSO_im)-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]
    
    #==============================================================================
    # Start set up for fitting:
    #==============================================================================
    if os.path.exists('fit_result_detect')==False:
        os.mkdir('fit_result_detect')
    
    psf_l, QSO_img_l, QSO_std_l = [], [], []
    for k in range(len(band_seq)):
        psf_l.append(PSF_list[k])
        QSO_img_l.append(QSO_im_list[k][ct:-ct,ct:-ct])
        QSO_std_l.append(err_map_list[k][ct:-ct,ct:-ct])
    
    for k in range(len(band_seq)):
        objs, Q_index = detect_obj(QSO_img_l[k], pltshow=0)
        qso_info = objs[Q_index]
        obj_temp = [objs[i] for i in range(len(objs)) if i != Q_index]
        if k == 0:
            obj = obj_temp
        if k>0 and len(obj_temp)>0:
            for i in range(len(obj_temp)):
                count = 0
                for j in range(len(obj)):
                    dis = np.sqrt(np.sum((np.asarray(obj[j][0])-np.asarray(obj_temp[i][0]))**2))
                    if i ==1:
                        print dis, (obj[j][1]+obj_temp[i][1])/2
                    if dis < (obj[j][1]+obj_temp[i][1])/2:
                        count += 1
                if count == 0:
                    obj.append(obj_temp[i])
            print "the number of nearby objs:", len(obj)
    
    from mask_objects import find_loc_max
    for k in range(0,5):  #['G', 'R', 'I', 'Z', 'Y']
        print "fiting the filename {0}, {1} band, QSO ID: ".format(file_name_seq[k],band_seq[k]) + qso_ID
        psf, QSO_img, QSO_std = psf_l[k], QSO_img_l[k], QSO_std_l[k]
            
        #for a quick glimsp of the number of image's local maximum
        x, y = find_loc_max(QSO_img)
        arr_x, arr_y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        center = len(QSO_img)/2
        bool_x, bool_y = (arr_x>(center-7))*(arr_x<(center+7)), (arr_y>(center-7))*(arr_y<(center+7))
        arr_x = arr_x[bool_x*bool_y]
        arr_y = arr_y[bool_x*bool_y]
        if len(arr_x)>=2:
            print "This image is likely to be a BHBH system!!!"
            print "Comparing the fitting Chisq:"
            if os.path.exists('fit_result_detect/{0}/'.format(file_name_seq[k]))==False:
                os.mkdir('fit_result_detect/{0}/'.format(file_name_seq[k]))
            reduced_Chisq_list_0, reduced_Chisq_list_1 = [], []
            pixels=len(QSO_img)**2
            for ft in range(1):     #The fitting rounds for each band
                fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
                #==============================================================================
                # input the objects components and parameteres
                #==============================================================================
                # input the parameter for the extended source
                fixed_source = []
                kwargs_source_init = []
                kwargs_source_sigma = []
                kwargs_lower_source = []
                kwargs_upper_source = []      
                fixed_source.append({})  
                kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
                kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
                kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
                
                if len(obj) >= 1:
                    for i in range(len(obj)):
                        fixed_source.append({})  
                        kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-10, 'center_y': obj[i][0][1]*pix_scale-10})
                        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+10, 'center_y': obj[i][0][1]*pix_scale+10})
                source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
                # input the parameter for the point source
                fixed_ps = []
                kwargs_ps_init = []
                kwargs_ps_sigma = []
                kwargs_lower_ps = []
                kwargs_upper_ps = []
                #==============================================================================
                # to fit and save the inference
                #==============================================================================
                fixcenter = False
                tag = 'fit_result_detect/{0}/fit_image_{0}_fittime-{1}'.format(file_name_seq[k],fit_time+1)
                source_result_0, ps_result_0, image_ps_0, image_host_0, error_map_0=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                                  background_rms=background_rms_list[k],
                                                                                  source_params=source_params, QSO_msk = None, fixcenter=fixcenter,
                                                                                  pix_sz = pix_scale, no_MCMC =True,
                                                                                  QSO_std =QSO_std, tag=tag, deep_seed= deep_seed)
                if len(image_host_0)>1:
                    image_host_0 = np.sum(image_host_0,axis=0)
                chiq_map_0 = ((QSO_img-image_ps_0-image_host_0)/(error_map_0))**2
                red_Chisq = chiq_map_0.sum()/pixels
                tag_name = tag + "_fitted_image"
                print call("mv {0} {1}".format(tag_name+'.pdf', tag_name+"_chisq_"+repr(round(red_Chisq,1)))+'.pdf', shell=True)
    #            reduced_Chisq_list_0.append(chiq_map_0.sum()/pixels)
            #==============================================================================
            # fitting the QSO as a BHBH        
            #==============================================================================
            for ft in range(1):     #The fitting rounds for each band            
                fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
                tag = 'fit_result_detect/{0}/fit_image_BHBH_{0}_fittime-{1}'.format(file_name_seq[k],fit_time+1)
                ps_x = ((arr_x - 40) * pix_scale) * -1
                ps_y = (arr_y - 40) * pix_scale
                fixed_ps = []
                kwargs_ps_init = []
                kwargs_ps_sigma = []
                kwargs_lower_ps = []
                kwargs_upper_ps = []
                point_amp = QSO_im.sum()/(len(x)+1.)
                for i in range(len(arr_x)):
                    fixed_ps.append({})
                    kwargs_ps_init.append({'ra_image': [ps_x[i]], 'dec_image': [ps_y[i]], 'point_amp': [point_amp]})
                    kwargs_ps_sigma.append({'ra_image': [0.05], 'dec_image': [0.05]})
                    kwargs_lower_ps.append({'ra_image': [-0.6], 'dec_image': [-0.6]})
                    kwargs_upper_ps.append({'ra_image': [0.6], 'dec_image': [0.6]})
                ps_param_2 = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]            
                source_result_1, ps_result_1, image_ps_1, image_host_1, error_map_1=fit_qso(QSO_img, psf_ave=psf, psf_std = None, ps_param = ps_param_2,
                                                                                  background_rms=background_rms_list[k],
                                                                                  source_params=source_params, QSO_msk = None, fixcenter=fixcenter,
                                                                                  pix_sz = pix_scale, no_MCMC =True,
                                                                                  QSO_std =QSO_std, tag=tag, deep_seed= deep_seed)                       
                if len(image_host_1)>1:
                    image_host_1 = np.sum(image_host_1,axis=0)
                chiq_map_1 = ((QSO_img-image_ps_1-image_host_1)/(error_map_1))**2
                red_Chisq = chiq_map_1.sum()/pixels
                tag_name = tag + "_fitted_image"                                      
                print call("mv {0} {1}".format(tag_name+'.pdf', tag_name+"_chisq_"+repr(round(red_Chisq,1)))+'.pdf', shell=True)
    #            reduced_Chisq_list_1.append(chiq_map_1.sum()/pixels) 
                                      
