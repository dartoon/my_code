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
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile
from matplotlib.colors import LogNorm
import copy

deep_seed = False  #Set as True to put more seed and steps to fit,
import time
t1 = time.time()
fix_on_I_band = input("Do you want, based on I band's inference, fix the other band's Reff and n???\n input 1 means yes, input others means no. \n")
pltshow = input("Do you want to plot while fitting?input 0 means no, input others means yes. \n")

filename_ascii = 'fit_result/ascii_result.txt'
if_file = glob.glob(filename_ascii)   
if if_file == []:
    ascii =  open(filename_ascii,'w') 
    ascii.write("#ID, RA, DEC, zp, pixel_scale, \
host_flux_ratio(%), host_x, host_y, host_flux, host_mag, \
Reff(arcsec), n_sersic, host_q, Kron_radius(arcsec), host_Kron_flux, host_Kron_mag, \
qso_x, qso_y, qso_flux, qso_mag, host_qso_center_mismatch(arcsec)\
\n")
elif if_file is not []:
    ascii = open(filename_ascii,"r+")
    ascii.read()

#ID, RA, DEC, zp, pixel_scale, host_flux_ratio(%),host_x, host_y, host_flux, host_mag, Kron_radius(arcsec), host_Kron_flux, host_Kron_mag, Reff(arcsec), n_sersic, host_q,qso_x, qso_y, qso_flux, qso_mag, host_qso_center_mismatch(arcsec)                               

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
    
    fit_frame_size = 81
    ct = (len(QSO_im)-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]
    
    #==============================================================================
    # Start set up for fitting:
    #==============================================================================
    if os.path.exists('fit_result')==False:
        os.mkdir('fit_result')
    
    filename = 'fit_result/fit_result_from_{0}-{1}_ID{2}.txt'.format(file_name_seq[0],file_name_seq[-1],qso_ID)
    if_file = glob.glob(filename)   
    if if_file == []:
        fit_result =  open(filename,'w') 
    elif if_file is not []:
        fit_result = open(filename,"r+")
        fit_result.read()


    psf_l, QSO_img_l, QSO_std_l = [], [], []
    for k in range(len(band_seq)):
        psf_l.append(PSF_list[k])
        QSO_img_l.append(QSO_im_list[k][ct:-ct,ct:-ct])
        QSO_std_l.append(err_map_list[k][ct:-ct,ct:-ct])
    
    for k in range(len(band_seq)):
        objs, Q_index = detect_obj(QSO_img_l[k],pltshow = pltshow)
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
    
    data_host_list = []
    for k in [2,0,1,3,4]:#len(band_seq)):
        print "fiting the filename {0}-{1}*.fits, {2} band, ID: ".format(file_name_seq[0],file_name_seq[-1],band_seq[k]) + qso_ID
        psf, QSO_img, QSO_std = psf_l[k], QSO_img_l[k], QSO_std_l[k]
        #==============================================================================
        #Compare the PSF and QSO profiles     
        #==============================================================================
        profiles_compare([psf, QSO_img], [1,1], ['PSF', 'QSO'], if_annuli=True)
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
        #==============================================================================
        # input the objects components and parameteres
        #==============================================================================
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        if band_seq[k] == 'I' or fix_on_I_band != 1: 
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
        elif band_seq[k] != 'I' and fix_on_I_band == 1: 
            fixed_source.append({'R_sersic': Iband_inf[0]['R_sersic'],'n_sersic': Iband_inf[0]['n_sersic']})  
            kwargs_source_init.append({'R_sersic': Iband_inf[0]['R_sersic'], 'n_sersic': Iband_inf[0]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -qso_info[0][0]*pix_scale, 'center_y': qso_info[0][1]*pix_scale})
            kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -0.5, 'center_y': -0.5})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': 0.5, 'center_y': 0.5})
            if len(obj) >= 1:
                for i in range(len(obj)):
                    fixed_source.append({'R_sersic': Iband_inf[i+1]['R_sersic'],'n_sersic': Iband_inf[i+1]['n_sersic']})  
                    kwargs_source_init.append({'R_sersic': Iband_inf[i+1]['R_sersic'], 'n_sersic': Iband_inf[i+1]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                    kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -obj[i][0][0]*pix_scale-10, 'center_y': obj[i][0][1]*pix_scale-10})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': -obj[i][0][0]*pix_scale+10, 'center_y': obj[i][0][1]*pix_scale+10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
        #==============================================================================
        # to fit and save the inference
        #==============================================================================
        fixcenter = False
        tag = 'fit_result/fit_image_{0}'.format(file_name_seq[k])
        source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                          background_rms=background_rms_list[k],
                                                                          source_params=source_params, QSO_msk = None, fixcenter=fixcenter,
                                                                          pix_sz = pix_scale, no_MCMC =True,
                                                                          QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow)
        if band_seq[k] == 'I':
            Iband_inf = source_result # in order to save the I band source_reslt to fix Reff and n for other band.

        if pltshow == 0:
            plot_compare=False
            fits_plot =False
        else:
            plot_compare=True
            fits_plot =True
        result = transfer_to_result(data=QSO_img, pix_sz = pix_scale,
                                    source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=error_map,
                                    zp=zp_list[k], fixcenter=fixcenter,ID=qso_ID+band_seq[k], QSO_msk = None, tag=tag, plot_compare = plot_compare)
        
        fitsFile = pyfits.open('../IMAGES/'+filename_list[k])
        file_header = copy.deepcopy(fitsFile[1].header)
        file_header['CRPIX1'] = file_header['CRPIX1']-qso_center_list[k][0]+len(QSO_img)/2
        file_header['CRPIX2'] = file_header['CRPIX2']-qso_center_list[k][1]+len(QSO_img)/2
        fit_result.write("#QSO_img intensity in the band {1}: {0} \n".format(round(np.sum(QSO_img),2), band_seq[k]))
        fit_result.write("#fit result of the QSO by in band {0}: \n".format(band_seq[k]))
        objs_img = np.zeros_like(image_host[0])
        if len(source_result)>1:
            for i in range(1,len(image_host)):
                objs_img += image_host[i]
        data_host = QSO_img - image_ps - objs_img
        data_host_list.append(data_host)
        center = np.asarray(data_host.shape) /2
        r_flux, r_grids, regions=flux_profile(data_host, center,
                                              radius=center.min(),
                                              grids=50, ifplot=False,
                                              fits_plot= fits_plot, mask_list=None)
        plt.plot(r_grids*pix_scale, r_flux)
        plt.title('Kron curve')
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('Kron flux')
        plt.savefig('fit_result/fit_image_host_Kron_curve_objID-{0}.pdf'.format(file_name_seq[k]))
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
        Kron_radius = 4.
        print "Kron_radius take {0} arcsec".format(Kron_radius)
        Kron_flux = r_flux[r_grids>4./pix_scale][0]
        Kron_mag = -2.5*np.log10(Kron_flux) + zp_list[k]
        result['Kron_radius'] = Kron_radius
        result['Kron_flux'] = round(Kron_flux,3)
        result['Kron_mag'] = round(Kron_mag,3)
        fit_result.write(repr(result) + "\n")
        if len(source_result)>1:
            for i in range(1,len(source_result)):
                result_i = transfer_obj_to_result(source_result[i],image_host[i], zp)
                fit_result.write("#obj {0}: \n".format(i))
                fit_result.write(repr(result_i) + "\n")
                plt.imshow(QSO_img, origin='low', norm=LogNorm())
                obj_x, obj_y = len(QSO_img)/2 - source_result[i]['center_x']/pix_scale, len(QSO_img)/2+source_result[i]['center_y']/pix_scale
                plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
            plt.savefig('fit_result/fit_image_{0}_objID'.format(file_name_seq[k]))
            if pltshow == 0:
                plt.close()
            else:
                plt.show()
            pyfits.PrimaryHDU(QSO_img-image_ps-image_host[0]-objs_img,header=file_header).writeto(tag+'_data-qso-host-objs(residual).fits',overwrite=True)
        pyfits.PrimaryHDU(QSO_img-image_ps,header=file_header).writeto(tag+'_data-qso.fits',overwrite=True)
        pyfits.PrimaryHDU(QSO_img-image_host[0],header=file_header).writeto(tag+'_data-host.fits',overwrite=True)
        pyfits.PrimaryHDU(QSO_img-image_ps-image_host[0],header=file_header).writeto(tag+'_data-qso-host.fits',overwrite=True)
        pyfits.PrimaryHDU(data_host,header=file_header).writeto(tag+'_data-qso-objs(host).fits',overwrite=True)
        
        ascii.write("{0} {1} {2} {3} {4} ".format(str(file_name_seq[k])+'-'+band_seq[k]+'-'+qso_ID, frm_c_RA_DEC_list[k][0][0][0], frm_c_RA_DEC_list[k][1][0][0], zp_list[k],0.168))
        ascii.write("{0}% {1} {2} {3} {4} ".format(result['host_flux_ratio_percent'], result['center_x'], result['center_y'],result['host_amp'],result['host_mag']))
        ascii.write("{0} {1} {2} {3} {4} {5} ".format(result['R_sersic'],result['n_sersic'],result['q'],result['Kron_radius'], result['Kron_flux'], result['Kron_mag']))
        qso_mag = -2.5*np.log10(result['QSO_amp']) + zp_list[k]
        if fixcenter == False:
            c_miss =np.sqrt((result['center_x']-result['qso_x'])**2+(result['center_y']-result['qso_y'])**2)
            ascii.write("{0} {1} {2} {3} {4}\n".format(result['qso_x'],result['qso_y'],result['QSO_amp'], qso_mag, c_miss))
        else:
            ascii.write("{0} {1} {2} {3} {4}\n".format(result['center_x'],result['center_y'],result['QSO_amp'], qso_mag, 0))
        t2 = time.time()
        time_sp = t2-t1
        time_ave = time_sp/((QSO_i-1)*5 + k + 1)
        time_total = time_ave * 59 * 5
        t_left= time_total - time_sp
        print "Finish percent %:",time_sp/time_total ,"total time needed :", time_total/60, "mins", "time_left", t_left/60, 'mins'
    fit_result.close()
    
#    ['G', 'R', 'I', 'Z', 'Y']
    if len(data_host_list) == 5:
        tag_c = 'fit_result/fit_image_{0}-{1}'.format(file_name_seq[0],file_name_seq[4])
        file_name = file_name_seq[0]
        host_g_i = data_host_list[0] -data_host_list[2]
        pyfits.PrimaryHDU(host_g_i,header=file_header).writeto(tag_c+'_host_g-i.fits',overwrite=True)
        host_r_z = data_host_list[1] -data_host_list[3]
        pyfits.PrimaryHDU(host_r_z,header=file_header).writeto(tag_c+'_host_r-z.fits',overwrite=True)
        host_i_y = data_host_list[2] -data_host_list[4]
        pyfits.PrimaryHDU(host_i_y,header=file_header).writeto(tag_c+'_host_i-y.fits',overwrite=True)
ascii.close()
    
import os
os.system('say "your program has finished"')
