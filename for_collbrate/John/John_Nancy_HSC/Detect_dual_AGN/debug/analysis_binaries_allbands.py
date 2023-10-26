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
from astropy.wcs import WCS

import sys
sys.path.insert(0, '../../../../py_tools/')

from fit_qso import fit_qso, fit_galaxy
from mask_objects import detect_obj
from matplotlib.colors import LogNorm
import copy
from subprocess import call

#image_ID, image_RA, image_DEC = '011935.29-002033.5' ,19.897051, -0.342664
#image_ID, image_RA, image_DEC = '090654.53+021315.2' ,136.727234,  2.220898
#image_ID, image_RA, image_DEC = '100935.02+004536.5' ,152.395935,  0.760158
#image_ID, image_RA, image_DEC = '125216.06+003141.1' ,193.066925,  0.528084
#image_ID, image_RA, image_DEC = '220501.19+003122.8' ,331.254974,  0.523022
#image_ID, image_RA, image_DEC = '221101.45+001449.0' ,332.756073,  0.246971
#image_ID, image_RA, image_DEC = '220228.65+004901.9' ,330.619415,  0.817198
#image_ID, image_RA, image_DEC = '233718.07+002550.6' ,354.325317,  0.430744
#image_ID, image_RA, image_DEC = '010748.49-000740.4' ,16.952066, -0.127890
#image_ID, image_RA, image_DEC = '020318.87-062321.3' ,30.828645, -6.389267
#image_ID, image_RA, image_DEC = '022105.64-044101.5' ,35.273505, -4.683758
#image_ID, image_RA, image_DEC = '094132.90+000731.1' ,145.387085,  0.125333
#image_ID, image_RA, image_DEC = '100719.31+030554.5' ,151.830475,  3.098492
#image_ID, image_RA, image_DEC = '120417.10+003653.7' ,181.071285,  0.614940
#image_ID, image_RA, image_DEC = '123821.66+010518.6' ,189.590281,  1.088519
#image_ID, image_RA, image_DEC = '124604.03-010954.6' ,191.516818, -1.165177
#image_ID, image_RA, image_DEC = '220132.66+030733.9' ,330.386099,  3.126103
#image_ID, image_RA, image_DEC = '220910.38-001601.5' ,332.293274, -0.267085
#image_ID, image_RA, image_DEC = '221101.45+001449.0' ,332.756073,  0.246971
#image_ID, image_RA, image_DEC = '222929.45+010438.4' ,337.372742,  1.077335
#image_ID, image_RA, image_DEC = '232853.40+011221.8' ,352.222507,  1.206083
#image_ID, image_RA, image_DEC = '235619.42+020907.7' ,359.080947,  2.152151

image_ID, image_RA, image_DEC = '001401.62-002600.5' , 3.506767211,  -0.4334980001

#image_ID = sys.argv[1] #'141637.44+003352.2' 
#image_RA = float(sys.argv[2]) #214.15602111816406
#image_DEC = float(sys.argv[3]) #0.5645210146903992

print image_ID, image_RA, image_DEC

deep_seed = False  #Set as True to put more seed and steps to fit,
pltshow = 0
fix_on_I_band = 1 #Set as 1 to fix the Sersic in as I band.

image_folder = '../images_directory/'
#image_folder = '../dual_AGNs/images/'
    
if os.path.exists('fit_result_detect')==False:
    os.mkdir('fit_result_detect')

filename_ascii = 'RESULTS/' + image_ID + '_result.txt'

band_run_list = [2,0,1]  #run I band first
#band_run_list=[0]


band_seq = ['G', 'R', 'I']
filename_list = [image_ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(len(band_seq))]

run_list = copy.deepcopy(band_run_list)
QSO_im_list, err_map_list, PSF_list=[], [], []
zp_list = []
for i in range(len(band_seq)):
    # The pixel scale is all 0.168
    if len(glob.glob(image_folder+filename_list[i])) == 0:
        print filename_list[i] + " DOES NOT EXIST!!!"
        QSO_im, err_map, PSF, _, zp, qso_center, fr_c_RA_DEC = [], [], [], [], [], [], []
        if i in band_run_list:
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
    if k == run_list[0]:
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
for k in run_list:  #['G', 'R', 'I', 'Z', 'Y']
    print "Testing the: "+ filename_list[k]
    psf, QSO_img, QSO_std = psf_l[k], QSO_img_l[k], QSO_std_l[k]
     
    #for a quick glimsp of the number of image's local maximum in I band
    if k == run_list[0]:
#        x, y = find_loc_max(QSO_img_l[2], neighborhood_size=4, threshold = 3)
#        threshold = background_rms_list[k]/0.1*3
        x, y = find_loc_max(QSO_img_l[k], neighborhood_size=3, threshold = 1)
        arr_x, arr_y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        center = len(QSO_img_l[2])/2
        bool_x, bool_y = (arr_x>(center-18))*(arr_x<(center+18)), (arr_y>(center-18))*(arr_y<(center+18))
        arr_x = arr_x[bool_x*bool_y]
        arr_y = arr_y[bool_x*bool_y]
    if len(arr_x)>=2:
        print "This {0} is likely to be a {1} system!!!".format(filename_list[k], 'BH'*len(arr_x))
        print "Comparing the fitting Chisq:"
        if os.path.exists('fit_result_detect/{0}/'.format(image_ID))==False:
            os.mkdir('fit_result_detect/{0}/'.format(image_ID))
        plt.imshow(QSO_img, origin='low', norm=LogNorm())
        for i in range(len(arr_x)):
            plt.text(arr_x[i], arr_y[i],'BH{0}'.format(i))
            plt.plot(arr_x[i], arr_y[i],'ro')
        plt.savefig('fit_result_detect/{0}/{1}band_highlight_BHBH.pdf'.format(image_ID,band_seq[k]))
        if pltshow == 1:
            plt.show()
        else:
            plt.close()

#Delete the obj[i] that is covered by the BH:
        obj_copy = copy.deepcopy(obj)
        obj = []
        for i in range(len(obj_copy)):
            pos_i = np.array(obj_copy[i][0])+center
            dis = [np.sqrt((pos_i[0] - arr_x[ii])**2 + (pos_i[1] - arr_y[ii])**2) for ii in range(len(arr_x))]
            if np.min(dis)>4:
                obj.append(obj_copy[i])
#            print dis 
#            print [np.sqrt((pos_i[0] - arr_x[ii])**2 + (pos_i[1] - arr_y[ii])**2) for ii in range(len(arr_x))]
        
        filename = 'fit_result_detect/{0}/fit_result.txt'.format(image_ID)
        if_file = glob.glob(filename)   
        if if_file == []:
            write_result =  open(filename,'w') 
            write_result.write("#The fitting information:\n")
        elif if_file is not []:
            write_result = open(filename,"r+")
            write_result.read()
        
        write_result.write("\n====================================================\n")
        write_result.write("#Fitting {0}-band :\n".format(band_seq[k]))
        pixels=len(QSO_img)**2
        fitsFile = pyfits.open(image_folder+filename_list[k])
        file_header = copy.deepcopy(fitsFile[1].header)
        file_header['CRPIX1'] = file_header['CRPIX1']-qso_center[0]+len(QSO_img)/2
        file_header['CRPIX2'] = file_header['CRPIX2']-qso_center[1]+len(QSO_img)/2
        wcs = WCS(file_header)   
        qsoid = filename_list[k].split('.fits')[0]
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
            if band_seq[k] == 'I' or fix_on_I_band != 1: 
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
            else:
                fixed_source.append({'R_sersic': Iband0_inf[0]['R_sersic'],'n_sersic': Iband0_inf[0]['n_sersic']})  
                kwargs_source_init.append({'R_sersic': Iband0_inf[0]['R_sersic'], 'n_sersic': Iband0_inf[0]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': 0, 'center_y': 0})
                kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -0.5, 'center_y': -0.5})
                kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': 0.5, 'center_y': 0.5})
                if len(obj) >= 1:
                    for i in range(len(obj)):
                        fixed_source.append({'R_sersic': Iband0_inf[i+1]['R_sersic'],'n_sersic': Iband0_inf[i+1]['n_sersic']})  
                        kwargs_source_init.append({'R_sersic': Iband0_inf[i+1]['R_sersic'], 'n_sersic': Iband0_inf[i+1]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                        kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
            source_params_0 = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
            #==============================================================================
            # to fit and save the inference
            #==============================================================================
            fixcenter = False
            tag = 'fit_result_detect/{0}/{2}band_fit_image0_PS+Sersic_fittime-{1}'.format(image_ID,fit_time+1,band_seq[k])
            print "fitting the QSO as one BH + Sersic "
            source_result_0, ps_result_0, image_ps_0, image_host_0, error_map_0, reduced_Chisq_0=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_0, QSO_msk = None, fixcenter=fixcenter,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow, return_Chisq=True)
            if band_seq[k] == 'I':
                Iband0_inf = source_result_0 # in order to save the I band source_reslt to fix Reff and n for other band.
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
            tag = 'fit_result_detect/{0}/{2}band_fit_image1_PSPS_fittime-{1}'.format(image_ID,fit_time+1,band_seq[k])
            ps_x = ((arr_x - fit_frame_size/2) * pix_scale) * -1
            ps_y = (arr_y - fit_frame_size/2) * pix_scale
            fixed_source = []
            kwargs_source_init = []
            kwargs_source_sigma = []
            kwargs_lower_source = []
            kwargs_upper_source = []
            
            
            if len(obj) >= 1:
                if band_seq[k] == 'I' or fix_on_I_band != 1: 
                    for i in range(len(obj)):
                        fixed_source.append({})  
                        kwargs_source_init.append({'R_sersic': obj[i][1] * pix_scale, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': obj[i][1] * pix_scale/5, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
                else:
                    for i in range(len(obj)):
                        fixed_source.append({'R_sersic': Iband1_inf[i]['R_sersic'],'n_sersic': Iband1_inf[i]['n_sersic']})  
                        kwargs_source_init.append({'R_sersic': Iband1_inf[i]['R_sersic'], 'n_sersic': Iband1_inf[i]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                        kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
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
            if band_seq[k] == 'I':
                Iband1_inf = source_result_1 # in order to save the I band source_reslt to fix Reff and n for other band.
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
            for i in range(len(ps_result_1)):
                write_result.write("\nAGN{0} position: ".format(i))
                write_result.write("x: "+repr(round(ps_result_1[i]['ra_image'],3))+' y: '+repr(round(ps_result_1[i]['dec_image'],3))+ "; ")            
            for i in range(len(ps_result_1)):  
                x_pixel = -ps_result_1[i]['ra_image']/pix_scale +center+1
                y_pixel = ps_result_1[i]['dec_image']/pix_scale +center+1
                Ra, Dec = wcs.all_pix2world([x_pixel],[y_pixel],1)
                write_result.write("\nAGN{0} WCS position: ".format(i))
                write_result.write("RA: "+repr(float(Ra))+' Dec: '+repr(float(Dec))+ "; ")            
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
            tag = 'fit_result_detect/{0}/{2}band_fit_image2_PSPS+Sersic_fittime-{1}'.format(image_ID,fit_time+1,band_seq[k])
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
            if band_seq[k] == 'I' or fix_on_I_band != 1: 
                source_params_2 = source_params_0
            else:
                fixed_source = []
                kwargs_source_init = []
                kwargs_source_sigma = []
                kwargs_lower_source = []
                kwargs_upper_source = []
                fixed_source.append({'R_sersic': Iband2_inf[0]['R_sersic'],'n_sersic': Iband2_inf[0]['n_sersic']})  
                kwargs_source_init.append({'R_sersic': Iband2_inf[0]['R_sersic'], 'n_sersic': Iband2_inf[0]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': 0, 'center_y': 0})
                kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -0.5, 'center_y': -0.5})
                kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': 0.5, 'center_y': 0.5})
                if len(obj) >= 1:
                    for i in range(len(obj)):
                        fixed_source.append({'R_sersic': Iband2_inf[i+1]['R_sersic'],'n_sersic': Iband2_inf[i+1]['n_sersic']})  
                        kwargs_source_init.append({'R_sersic': Iband2_inf[i+1]['R_sersic'], 'n_sersic': Iband2_inf[i+1]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                        kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -obj[i][0][0]*pix_scale-0.5, 'center_y': obj[i][0][1]*pix_scale-0.5})
                        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': -obj[i][0][0]*pix_scale+0.5, 'center_y': obj[i][0][1]*pix_scale+0.5})
                source_params_2 = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
            print "fitting the QSO as {0} point sources + Sersic".format(len(arr_x))
            source_result_2, ps_result_2, image_ps_2, image_host_2, error_map_2, reduced_Chisq_2=fit_qso(QSO_img, psf_ave=psf, psf_std = None, ps_param = ps_param_2,
                                                                              background_rms=background_rms_list[k],
                                                                              source_params=source_params_2, QSO_msk = None, fixcenter=fixcenter,
                                                                              pix_sz = pix_scale, no_MCMC =True,
                                                                              QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow, return_Chisq=True)                       
            if band_seq[k] == 'I':
                Iband2_inf = source_result_2                
            
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
            for i in range(len(ps_result_2)):  
                x_pixel = -ps_result_2[i]['ra_image']/pix_scale +center+1
                y_pixel = ps_result_2[i]['dec_image']/pix_scale +center+1
                Ra, Dec = wcs.all_pix2world([x_pixel],[y_pixel],1)
                write_result.write("\nAGN{0} WCS position: ".format(i))
                write_result.write("RA: "+repr(float(Ra))+' Dec: '+repr(float(Dec))+ "; ")  
            write_result.write("\nPS PS center offset (arcsec): "+repr(round(c_miss,3)))
            write_result.write("\n=======================================================\n")
            tag_name = tag + "_fitted_image"  
            objs_img = np.zeros_like(image_host_2[0])
            if len(image_host_2)>1:
                for i in range(1,len(image_host_2)):
                    objs_img += image_host_2[i]
            pyfits.PrimaryHDU(QSO_img-image_ps_2-objs_img,header=file_header).writeto('fit_result_detect/{0}/{1}band_data-BHBH(host).fits'.format(image_ID, band_seq[k]),overwrite=True)
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_2,1)))+'.pdf', shell=True)   
        #==============================================================================
        # fitting the QSO as a sSersic + lSersic        
        #==============================================================================
        for ft in range(1):     #The fitting rounds for each band            
            fit_time = ft #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
            tag = 'fit_result_detect/{0}/{2}band_fit_image3_sSersic+Sersic_fittime-{1}'.format(image_ID,fit_time+1,band_seq[k])
            ps_x = ((arr_x - fit_frame_size/2) * pix_scale) * -1
            ps_y = (arr_y - fit_frame_size/2) * pix_scale
            fixed_source = []
            kwargs_source_init = []
            kwargs_source_sigma = []
            kwargs_lower_source = []
            kwargs_upper_source = [] 
            if band_seq[k] == 'I' or fix_on_I_band != 1:             
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
            else:
                for i in range(len(arr_x)):
                    fixed_source.append({'R_sersic': Iband3_inf[i]['R_sersic'],'n_sersic': Iband3_inf[i]['n_sersic'], 'e1': 0., 'e2': 0.})  
                    kwargs_source_init.append({'R_sersic': Iband3_inf[i]['R_sersic'], 'n_sersic': Iband3_inf[i]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x':ps_x[i], 'center_y': ps_y[i]})
                    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.01, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'R_sersic': 0.01, 'n_sersic': 0.5, 'center_x': ps_x[i]-0.6, 'center_y': ps_y[i]-0.6})
                    kwargs_upper_source.append({'R_sersic': 0.3, 'n_sersic': 4., 'center_x':ps_y[i]+0.6, 'center_y': ps_y[i]+0.6})            
                
                fixed_source.append({'R_sersic': Iband3_inf[len(arr_x)]['R_sersic'],'n_sersic': Iband3_inf[len(arr_x)]['n_sersic']})  
                kwargs_source_init.append({'R_sersic': Iband3_inf[len(arr_x)]['R_sersic'], 'n_sersic': Iband3_inf[len(arr_x)]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
                kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.2, 'center_y': 0.2})
                kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
                kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
                if len(obj) >= 1:
                    for i in range(len(obj)):
                        fixed_source.append({'R_sersic': Iband3_inf[len(arr_x)+1+i]['R_sersic'],'n_sersic': Iband3_inf[len(arr_x)+1+i]['n_sersic']})  
                        kwargs_source_init.append({'R_sersic': Iband3_inf[len(arr_x)+1+i]['R_sersic'], 'n_sersic': Iband3_inf[len(arr_x)+1+i]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
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
            if band_seq[k] == 'I':
                Iband3_inf = source_result_3
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
            for i in range(len(arr_x)):  
                x_pixel = -source_result_3[i]['center_x']/pix_scale +center+1
                y_pixel = source_result_3[i]['center_y']/pix_scale +center+1
                Ra, Dec = wcs.all_pix2world([x_pixel],[y_pixel],1)
                write_result.write("\nAGN{0} WCS position: ".format(i))
                write_result.write("RA: "+repr(float(Ra))+' Dec: '+repr(float(Dec))+ "; ")             
            write_result.write("\nPS PS center offset (arcsec): "+repr(round(c_miss,3)))
            write_result.write("\n====================================================end\n")
            tag_name = tag + "_fitted_image"                                      
            print call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_3,1)))+'.pdf', shell=True)       
        write_result.close()                             
#os.system('say "your program has finished"')
#print("Program has finished")
#
