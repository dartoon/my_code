#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:14 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from gen_fit_id_binary import gen_fit_id_binary
from photutils import make_source_mask
import os

from fit_qso import fit_qso, fit_galaxy
from mask_objects import detect_obj
from matplotlib.colors import LogNorm
import copy
from subprocess import call
import sys

#image_ID ='084710.40-001302.6' 
#image_RA = 131.7934293
#image_DEC = -0.2174615916

#image_ID ='092158.92+034235.7' 
#image_RA = 140.495501
#image_DEC = 3.709936749

#image_ID ='094943.34+000536.2' 
#image_RA = 147.4305489
#image_DEC = 0.09334982316

#image_ID ='121405.12+010205.1' 
#image_RA = 183.5213596
#image_DEC = 1.034757599

image_ID = sys.argv[1] #'141637.44+003352.2' 
image_RA = float(sys.argv[2]) #214.15602111816406
image_DEC = float(sys.argv[3]) #0.5645210146903992

#image_ID ='100043.13+020637.2' 
#image_RA = 150.1797789
#image_DEC = 2.110369603

#image_ID = "000017.88+002612.6"
#image_RA = 0.07452999800443649
#image_DEC = 0.4368380010128021


print image_ID, image_RA, image_DEC

deep_seed = True  #Set as True to put more seed and steps to fit,
pltshow = 1

image_folder = './images_directory/'
    
if os.path.exists('fit_result_detect')==False:
    os.mkdir('fit_result_detect')

filename_ascii = 'RESULTS/' + image_ID + '_result.txt'

#band_run_list = [2,0,1,3,4]  #run I band first
band_run_list=[0]
band_seq = ['I']
filename_list = [image_ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(len(band_run_list))]

run_list = copy.deepcopy(band_run_list)
QSO_im_list, err_map_list, PSF_list=[], [], []
zp_list = []
for i in range(len(band_seq)):
    # The pixel scale is all 0.168
    if len(glob.glob(image_folder+filename_list[i])) == 0:
        print filename_list[i] + " DOES NOT EXIST!!!"
        QSO_im, err_map, PSF, pix_scale, zp, qso_center, fr_c_RA_DEC = [], [], [], [], [], [], [], []
        run_list.remove(i)
    else:
        QSO_im, err_map, PSF, pix_scale, zp, qso_center, fr_c_RA_DEC = gen_fit_id_binary(image_folder, image_RA, image_DEC, filename_list[i], cut_frame=80)
    QSO_im_list.append(QSO_im)
    err_map_list.append(err_map)
    if len(PSF) != 0 and PSF.shape[0] != PSF.shape[1]:
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

#==============================================================================
# quickly evaluate the background rms
#==============================================================================
background_rms_list = []
for i in range(len(band_seq)):
    if i in run_list:
        mask = make_source_mask(QSO_im_list[i], snr=3, npixels=5, dilate_size=11)
        background_rms_list.append(np.std(QSO_im_list[i]* (1-mask*1)))
    else:
        background_rms_list.append([])

fit_frame_size = 81
ct = (len(QSO_im)-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

#==============================================================================
# Start set up for fitting:
#==============================================================================
psf_l, QSO_img_l, QSO_std_l = [], [], []
for k in range(len(band_seq)):
    if k in run_list:
        psf_l.append(PSF_list[k])
        QSO_img_l.append(QSO_im_list[k][ct:-ct,ct:-ct])
        QSO_std_l.append(err_map_list[k][ct:-ct,ct:-ct])
    else:
        psf_l.append([])
        QSO_img_l.append([])
        QSO_std_l.append([])            

for k in run_list:
    objs, Q_index = detect_obj(QSO_img_l[k], pltshow=pltshow)
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

from mask_objects import find_loc_max, measure_FWHM

for k in run_list:  #['G', 'R', 'I', 'Z', 'Y']
    print "Fiting the: "+ filename_list[k]
    if_dual = False
    psf, QSO_img, QSO_std = psf_l[k], QSO_img_l[k], QSO_std_l[k]
        
    #for a quick glimsp of the number of image's local maximum
    x, y = find_loc_max(QSO_img, neighborhood_size = 3, threshold = 1)
    arr_x, arr_y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    center = len(QSO_img)/2
    bool_x, bool_y = (arr_x>(center-18))*(arr_x<(center+18)), (arr_y>(center-18))*(arr_y<(center+18))
    arr_x = arr_x[bool_x*bool_y]
    arr_y = arr_y[bool_x*bool_y]
    qsoid = filename_list[k].split('.fits')[0]
    
    if len(arr_x)>=2:
        if_dual = True
        claim = "This {0} is likely to be a {1} system!!!".format(filename_list[k], 'BH'*len(arr_x))
    elif len(arr_x)==1:
        FWHM_ver, FWHM_hor, FWHM_xy, FWHM_xy_ = measure_FWHM(QSO_img)    #Measure FWHM from ver, hor, up_right, up_left
        FWHM_array = np.array([FWHM_ver, FWHM_hor, FWHM_xy, FWHM_xy_])
        if (np.max(FWHM_array) - np.median(FWHM_array))/np.median(FWHM_array) > 0.2 :
            if_dual = True
            direction = np.where(FWHM_array == FWHM_array.max())[0][0]
            arr_x = np.concatenate([arr_x,arr_x])
            arr_y = np.concatenate([arr_y,arr_y])
            claim = "This {0} is likely to be very closed dual AGN system!!!".format(filename_list[k])
    if if_dual == True:
        print claim
        print "Comparing the fitting Chisq:"
        if os.path.exists('fit_result_detect/{0}/'.format(qsoid))==False:
            os.mkdir('fit_result_detect/{0}/'.format(qsoid))
        plt.imshow(QSO_img, origin='low', norm=LogNorm())
        for i in range(len(arr_x)):
            plt.text(arr_x[i], arr_y[i],'BH{0}'.format(i))
            plt.plot(arr_x[i], arr_y[i],'ro')
        plt.savefig('fit_result_detect/{0}/highlight_BHBH.pdf'.format(qsoid))
        if pltshow == 1:
            plt.show()
        else:
            plt.close()
        write_result =  open('fit_result_detect/{0}/fit_result.txt'.format(qsoid),'w') 
        write_result.write("#The fitting information:\n")
            
        pixels=len(QSO_img)**2
        #==============================================================================
        # fitting the QSO as a BH + Sersic       
        #==============================================================================
        for ft in range(1):     #The fitting rounds for each band
            fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
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
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
            source_params_0 = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
            #==============================================================================
            # to fit and save the inference
            #==============================================================================
            fixcenter = False
            tag = 'fit_result_detect/{0}/fit_image0_PS+Sersic_fittime-{1}'.format(qsoid,fit_time+1)
            print "fitting the QSO as one BH + Sersic "
            source_result_0, ps_result_0, image_ps_0, image_host_0, error_map_0, reduced_Chisq_0=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_0, QSO_msk = None, fixcenter=fixcenter,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow, return_Chisq=True)
            host_mag = -2.5*np.log10(image_host_0[0].sum()) + zp_list[k]
            c_miss = np.sqrt((source_result_0[0]['center_x']-ps_result_0[0]['ra_image'])**2+(source_result_0[0]['center_y']-ps_result_0[0]['dec_image'])**2)
            AGN_mag = -2.5*np.log10(ps_result_0[0]['point_amp'].sum()) + zp_list[k]
            write_result.write("1. Fitting as a regular QSO,i.e. one PS + Sersic:\n")
            write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_0,3)))
            write_result.write("\nHost mag: "+repr(round(host_mag,3)))
            write_result.write("\nAGN mag: "+repr(round(AGN_mag,3)))
            write_result.write("\nPS Sersic center offset (arcsec): "+repr(round(c_miss,3)) + "; ")
            write_result.write("\n=======================================================\n")
            tag_name = tag + "_fitted_image"
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_0,1)))+'.pdf', shell=True)
        #==============================================================================
        # fitting the QSO as a BHBH        
        #==============================================================================
        for ft in range(1):     #The fitting rounds for each band            
            fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
            tag = 'fit_result_detect/{0}/fit_image1_PSPS_fittime-{1}'.format(qsoid,fit_time+1)
            ps_x = ((arr_x - fit_frame_size/2) * pix_scale) * -1
            ps_y = (arr_y - fit_frame_size/2) * pix_scale
            fixed_source = []
            kwargs_source_init = []
            kwargs_source_sigma = []
            kwargs_lower_source = []
            kwargs_upper_source = []     
            if len(obj) >= 1:
                for i in range(len(obj)):
                    fixed_source.append({})  
                    kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
                source_params_1 = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
            else:
                source_params_1 = []  #This means no Sersic source
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
                kwargs_lower_ps.append({'ra_image': [ps_x[i]-0.6], 'dec_image': [ps_y[i]-0.6]})
                kwargs_upper_ps.append({'ra_image': [ps_x[i]+0.6], 'dec_image': [ps_y[i]+0.6]})
            ps_param_1 = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]  
            print "fitting the QSO as {0} point sources".format(len(arr_x))
            source_result_1, ps_result_1, image_ps_1, image_host_1, error_map_1, reduced_Chisq_1=fit_qso(QSO_img, psf_ave=psf, psf_std = None, ps_param = ps_param_1,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_1, QSO_msk = None,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow,return_Chisq=True)                  
            if image_host_1 == []:
                image_host_1 = np.zeros_like(image_ps_1)
            if len(ps_result_1) == 2:
                c_miss = np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)
            elif len(ps_result_1) > 2:
                c_miss = [np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)]
                c_miss.append(np.sqrt((ps_result_1[1]['ra_image']-ps_result_1[2]['ra_image'])**2+(ps_result_1[1]['dec_image']-ps_result_1[2]['dec_image'])**2))
                c_miss.append(np.sqrt((ps_result_1[2]['ra_image']-ps_result_1[0]['ra_image'])**2+(ps_result_1[2]['dec_image']-ps_result_1[0]['dec_image'])**2))
                c_miss = np.average(c_miss)
            AGN_mags = [-2.5*np.log10(ps_result_1[i]['point_amp'].sum()) + zp_list[k] for i in range(len(ps_result_1))]
            write_result.write("2. Fitting as {0}PS:\n".format(len(ps_result_1)))
            write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_1,3)))
            write_result.write("\nAGN mag: ")
            for i in range(len(ps_result_1)):
                write_result.write(repr(round(AGN_mags[i],3))+' ')
            write_result.write("\n")
            for i in range(len(ps_result_1)):
                write_result.write("AGN{0} position: ".format(i))
                write_result.write("x: "+repr(round(ps_result_1[i]['ra_image'],3))+' y: '+repr(round(ps_result_1[i]['dec_image'],3))+ "; ")            
            write_result.write("\nPS PS center offset (arcsec): "+repr(round(c_miss,3)))
            write_result.write("\n=======================================================\n")
            tag_name = tag + "_fitted_image"  
            tag_name = tag + "_fitted_image"                                      
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_1,1)))+'.pdf', shell=True)
        #==============================================================================
        # fitting the QSO as a BHBH + Sersic       
        #==============================================================================
        for ft in range(1):     #The fitting rounds for each band            
            fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
            tag = 'fit_result_detect/{0}/fit_image2_PSPS+Sersic_fittime-{1}'.format(qsoid,fit_time+1)
            ps_x = ((arr_x - fit_frame_size/2) * pix_scale) * -1
            ps_y = (arr_y - fit_frame_size/2) * pix_scale
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
                kwargs_lower_ps.append({'ra_image': [ps_x[i]-0.6], 'dec_image': [ps_y[i]-0.6]})
                kwargs_upper_ps.append({'ra_image': [ps_x[i]+0.6], 'dec_image': [ps_y[i]+0.6]})
            ps_param_2 = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps] 
            print "fitting the QSO as {0} point sources + Sersic".format(len(arr_x))
            source_result_2, ps_result_2, image_ps_2, image_host_2, error_map_2, reduced_Chisq_2=fit_qso(QSO_img, psf_ave=psf, psf_std = None, ps_param = ps_param_2,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_0, QSO_msk = None, fixcenter=fixcenter,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow, return_Chisq=True)                       
            host_mag = -2.5*np.log10(image_host_2[0].sum()) + zp_list[k]
            if len(ps_result_2) == 2:
                c_miss = np.sqrt((ps_result_2[0]['ra_image']-ps_result_2[1]['ra_image'])**2+(ps_result_2[0]['dec_image']-ps_result_2[1]['dec_image'])**2)
            elif len(ps_result_2) > 2:
                c_miss = [np.sqrt((ps_result_2[0]['ra_image']-ps_result_2[1]['ra_image'])**2+(ps_result_2[0]['dec_image']-ps_result_2[1]['dec_image'])**2)]
                c_miss.append(np.sqrt((ps_result_2[1]['ra_image']-ps_result_2[2]['ra_image'])**2+(ps_result_2[1]['dec_image']-ps_result_2[2]['dec_image'])**2))
                c_miss.append(np.sqrt((ps_result_2[2]['ra_image']-ps_result_2[0]['ra_image'])**2+(ps_result_2[2]['dec_image']-ps_result_2[0]['dec_image'])**2))
                c_miss = np.average(c_miss)
            AGN_mags = [-2.5*np.log10(ps_result_2[i]['point_amp'].sum()) + zp_list[k] for i in range(len(ps_result_2))]
            write_result.write("3. Fitting as {0}PS + Sersic:\n".format(len(ps_result_2)))
            write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_2,3)))
            write_result.write("\nHost mag: "+repr(round(host_mag,3)))
            write_result.write("\nAGN mag: ")
            for i in range(len(ps_result_2)):
                write_result.write(repr(round(AGN_mags[i],3))+' ')
            write_result.write("\n")
            for i in range(len(ps_result_2)):
                write_result.write("AGN{0} position: ".format(i))
                write_result.write("x: "+repr(round(ps_result_2[i]['ra_image'],3))+' y: '+repr(round(ps_result_2[i]['dec_image'],3))+ "; ")
            write_result.write("\nPS PS center offset (arcsec): "+repr(round(c_miss,3)))
            write_result.write("\n=======================================================\n")
            tag_name = tag + "_fitted_image"  
            objs_img = np.zeros_like(image_host_2[0])
            if len(image_host_2)>1:
                for i in range(1,len(image_host_2)):
                    objs_img += image_host_2[i]
            fitsFile = pyfits.open(image_folder+filename_list[k])
            file_header = copy.deepcopy(fitsFile[1].header)
            file_header['CRPIX1'] = file_header['CRPIX1']-qso_center[0]+len(QSO_img)/2
            file_header['CRPIX2'] = file_header['CRPIX2']-qso_center[1]+len(QSO_img)/2
            pyfits.PrimaryHDU(QSO_img-image_ps_2-objs_img,header=file_header).writeto('fit_result_detect/{0}/data-BHBH(host).fits'.format(qsoid),overwrite=True)
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_2,1)))+'.pdf', shell=True)   
        #==============================================================================
        # fitting the QSO as a sSersic + lSersic        
        #==============================================================================
        for ft in range(1):     #The fitting rounds for each band            
            fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
            tag = 'fit_result_detect/{0}/fit_image3_sSersic+Sersic_fittime-{1}'.format(qsoid,fit_time+1)
            ps_x = ((arr_x - fit_frame_size/2) * pix_scale) * -1
            ps_y = (arr_y - fit_frame_size/2) * pix_scale
            fixed_source = []
            kwargs_source_init = []
            kwargs_source_sigma = []
            kwargs_lower_source = []
            kwargs_upper_source = []   
            for i in range(len(arr_x)):
                fixed_source.append({'e1': 0., 'e2': 0.})  
                kwargs_source_init.append({'R_sersic': 0.05, 'n_sersic': 1., 'e1': 0., 'e2': 0., 'center_x':ps_x[i], 'center_y': ps_y[i]})
                kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.01, 'center_x': 0.1, 'center_y': 0.1})
                kwargs_lower_source.append({'R_sersic': 0.01, 'n_sersic': 0.5, 'center_x': ps_x[i]-0.6, 'center_y': ps_y[i]-0.6})
                kwargs_upper_source.append({'R_sersic': 0.3, 'n_sersic': 4., 'center_x':ps_y[i]+0.6, 'center_y': ps_y[i]+0.6})            
            
            fixed_source.append({})  
            kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.2, 'center_y': 0.2})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
            if len(obj) >= 1:
                for i in range(len(obj)):
                    fixed_source.append({})  
                    kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
            source_params_3 = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
            print "fitting the QSO as {0} small Sersic + Sersic".format(len(arr_x))
            source_result_3, image_host_3, error_map_3, reduced_Chisq_3=fit_galaxy(QSO_img, psf_ave=psf, psf_std = None,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_3, galaxy_msk = None,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              galaxy_std =QSO_std, tag=tag, deep_seed= 'very_deep', pltshow=pltshow, return_Chisq=True)
#            host_mag = -2.5*np.log10(image_host_3[len(arr_x)+1].sum()) + zp_list[k]
            host_mag = -2.5*np.log10(image_host_3[len(arr_x)].sum()) + zp_list[k]
            AGN_mags = [-2.5*np.log10(image_host_3[i].sum()) + zp_list[k] for i in range(len(arr_x))]
            if len(arr_x) == 2:
                c_miss = np.sqrt((source_result_3[0]['center_x']-source_result_3[1]['center_x'])**2+(source_result_3[0]['center_y']-source_result_3[1]['center_y'])**2)
            elif len(arr_x) > 2:
                c_miss = [np.sqrt((source_result_3[0]['center_x']-source_result_3[1]['center_x'])**2+(source_result_3[0]['center_y']-source_result_3[1]['center_y'])**2)]
                c_miss.append(np.sqrt((source_result_3[1]['center_x']-source_result_3[2]['center_x'])**2+(source_result_3[1]['center_y']-source_result_3[2]['center_y'])**2))
                c_miss.append(np.sqrt((source_result_3[2]['center_x']-source_result_3[0]['center_x'])**2+(source_result_3[2]['center_y']-source_result_3[0]['center_y'])**2))
                c_miss = np.average(c_miss)
            write_result.write("4. Fitting as {0}_small_Sersic + larger_Sersic:\n".format(len(arr_x)))
            write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_3,3)))
            write_result.write("\nHost mag: "+repr(round(host_mag,3)))
            write_result.write("\nAGN mag: ")
            for i in range(len(arr_x)):
                write_result.write(repr(round(AGN_mags[i],3))+' ')
            write_result.write("\n")
            for i in range(len(arr_x)):
                write_result.write("AGN{0} position: ".format(i))
                write_result.write("x: "+repr(round(source_result_3[i]['center_x'],3))+' y: '+repr(round(source_result_3[i]['center_y'],3))+ "; ")
            write_result.write("\nPS PS center offset (arcsec): "+repr(round(c_miss,3)))
            write_result.write("\n====================================================end\n")
            tag_name = tag + "_fitted_image"                                      
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_3,1)))+'.pdf', shell=True)       
        write_result.close()                             
#os.system('say "your program has finished"')
print("Program has finished")
#
