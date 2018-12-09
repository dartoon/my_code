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

from fit_qso import fit_qso_multiband
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare
from matplotlib.colors import LogNorm
import copy

deep_seed = True  #Set as True to put more seed and steps to fit,
import time
t1 = time.time()
for QSO_i in range(1,2):  # run the rest sample
    ID_seq = QSO_i
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
    
    psf_l, QSO_img_l, QSO_std_l = [], [], []
    for k in range(len(band_seq)):
        psf_l.append(PSF_list[k])
        QSO_img_l.append(QSO_im_list[k][ct:-ct,ct:-ct])
        QSO_std_l.append(err_map_list[k][ct:-ct,ct:-ct])
    
    for k in range(len(band_seq)):
        objs, Q_index = detect_obj(QSO_img_l[k])
        qso_info = objs[Q_index]
        obj_temp = [objs[i] for i in range(len(objs)) if i != Q_index]
        if k == 0:
            obj = obj_temp
        if k>0 and len(obj_temp)>0:
            for i in range(len(obj_temp)):
                count = 0
                for j in range(len(obj)):
                    dis = np.sqrt(np.sum((np.asarray(obj[j][0])-np.asarray(obj_temp[i][0]))**2))
                    if dis < (obj[j][1]+obj_temp[i][1])/2:
                        count += 1
                if count == 0:
                    obj.append(obj_temp[i])
            print "the number of nearby objs:", len(obj)
    
    print "fiting the filename {0}-{1}*.fits, {2} band, ID: ".format(file_name_seq[0],file_name_seq[-1],band_seq) + qso_ID
    psfs, QSO_imgs, QSO_stds = psf_l, QSO_img_l, QSO_std_l
#    #==============================================================================
#    #Compare the PSF and QSO profiles     
#    #==============================================================================
#    profiles_compare([psf, QSO_img], [1,1], ['PSF', 'QSO'], if_annuli=True)
#    plt.show()
    #==============================================================================
    # input the objects components and parameteres
    #==============================================================================
    fixed_source = []
    kwargs_source_init = []
    kwargs_source_sigma = []
    kwargs_lower_source = []
    kwargs_upper_source = []      
    fixed_source.append({})  
    kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -qso_info[0][0]*pix_scale, 'center_y': qso_info[0][1]*pix_scale})
    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.5, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
    if len(obj) >= 1:
        for i in range(len(obj)):
            fixed_source.append({})  
            kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-10, 'center_y': obj[i][0][1]*pix_scale-10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+10, 'center_y': obj[i][0][1]*pix_scale+10})
    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    
    
    #==============================================================================
    # Start set up for fitting file:
    #==============================================================================
    if os.path.exists('fit_result_joint_band')==False:
        os.mkdir('fit_result_joint_band')
    
    filename = 'fit_result_joint_band/fit_result_joint_{0}-{1}_ID{2}.txt'.format(file_name_seq[0],file_name_seq[-1],qso_ID)
    if_file = glob.glob(filename)   
    if if_file == []:
        fit_result =  open(filename,'w') 
    elif if_file is not []:
        fit_result = open(filename,"r+")
        fit_result.read()
    
    #==============================================================================
    # to fit and save the inference
    #==============================================================================
    fixcenter = False
    tag = 'fit_result_joint_band/fit_image_{0}-{1}'.format(file_name_seq[0],file_name_seq[-1])
    source_result_list, ps_result_list, image_ps_list, image_host_list, error_map_list, shift_RADEC_list, fitting_seq = fit_qso_multiband(QSO_imgs, psf_ave_list=psfs, psf_std_list = None,
                                                              background_rms_list=background_rms_list,
                                                              source_params=source_params, QSO_msk_list = None, fixcenter=fixcenter,
                                                              pix_sz = pix_scale, no_MCMC =True,
                                                              QSO_std_list =QSO_stds, tag=tag, deep_seed= deep_seed)
    for k in range(len(band_seq)):
        result_k = transfer_to_result(data=QSO_imgs[k], pix_sz = pix_scale,
                                    source_result=source_result_list[k], ps_result=ps_result_list[k], image_ps=image_ps_list[k], image_host=image_host_list[k], error_map=error_map_list[k],
                                    zp=zp_list[k], fixcenter=fixcenter,ID=qso_ID+band_seq[k], QSO_msk = None, tag=tag)
        del result_k['amp']
        fitsFile = pyfits.open('../IMAGES/'+filename_list[k])
        file_header = copy.deepcopy(fitsFile[1].header)
        file_header['CRPIX1'] = file_header['CRPIX1']-qso_center_list[k][0]+len(QSO_imgs[k])/2
        file_header['CRPIX2'] = file_header['CRPIX2']-qso_center_list[k][1]+len(QSO_imgs[k])/2
        fit_result.write("#QSO_img intensity in the band {1}: {0} \n".format(round(np.sum(QSO_imgs[k]),2), band_seq[k]))
        fit_result.write("#fit result of the QSO by in band {0}: \n".format(band_seq[k]))
        fit_result.write(repr(result_k) + "\n")
        if len(source_result_list[k])>1:
            for i in range(1,len(source_result_list[k])):
                result_i = transfer_obj_to_result(source_result_list[k][i],image_host_list[k][i], zp)
                del result_i['amp']
                fit_result.write("#obj {0}: \n".format(i))
                fit_result.write(repr(result_i) + "\n")
                plt.imshow(QSO_imgs[k], origin='low', norm=LogNorm())
                obj_x, obj_y = len(QSO_imgs[k])/2 - source_result_list[k][i]['center_x']/pix_scale, len(QSO_imgs[k])/2+source_result_list[k][i]['center_y']/pix_scale
                plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
            plt.savefig('fit_result_joint_band/fit_image_{0}_objID'.format(file_name_seq[k]))
            plt.show()
            objs_img = np.zeros_like(image_host_list[k][1])
            for i in range(1,len(image_host_list[k])):
                objs_img += image_host_list[k][i]
            pyfits.PrimaryHDU(QSO_imgs[k]-image_ps_list[k]-image_host_list[k][0]-objs_img,header=file_header).writeto(tag+'_data-qso-host-objs(residual)_band{0}.fits'.format(k),overwrite=True)
        pyfits.PrimaryHDU(QSO_imgs[k]-image_ps_list[k],header=file_header).writeto(tag+'_data-qso_band{0}.fits'.format(k),overwrite=True)
        pyfits.PrimaryHDU(QSO_imgs[k]-image_host_list[k][0],header=file_header).writeto(tag+'_data-host_band{0}.fits'.format(k),overwrite=True)
        pyfits.PrimaryHDU(QSO_imgs[k]-image_ps_list[k]-image_host_list[k][0],header=file_header).writeto(tag+'_data-qso-host_band{0}.fits'.format(k),overwrite=True)
    fit_result.close()
    
    t2 = time.time()
    time_sp = t2-t1
    time_ave = time_sp/((QSO_i-1) + 1)
    time_total = time_ave * 59
    t_left= time_total - time_sp
    print "Finish percent %:",time_sp/time_total ,"total time needed :", time_total/60, "mins", "time_left", t_left/60, 'mins'
#import os
#os.system('say "your program has finished"')
