#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:14 2018

@author: Dartoon

NOTE:
Put this script in folder “analysis_directory/“.
To fit the QSO, you can command in the terminal as follows:
    
python analysis_id_follow.py 084702.80+013001.4

This will allowing you to run the QSO ID 084702.80+013001.4 
Also, you can run two QSO one by one by put two QSO IDs:
    
python analysis_id_follow.py 084702.80+013001.4 084940.34+014638.4
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from gen_fit_id import gen_fit_id
from photutils import make_source_mask
import os
#plt.ion()

import sys
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result, transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile
from matplotlib.colors import LogNorm
import copy
import time

#Setting the fitting condition:
deep_seed = False  #Set as True to put more seed and steps to fit.
fix_on_I_band = 1 # input("Do you want, based on I band's inference, fix the other band's Reff and n???\n input 1 means yes, input others means no. \n")
pltshow = 0 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
sub_bkg = True
fit_frame_size = 71
image_folder = '../images_directory/'
run_MCMC = True

QSO_id = sys.argv[1]
QSO_RA = float(sys.argv[2])
QSO_DEC = float(sys.argv[3])
#print isinstance(QSO_id,str), isinstance(QSO_RA,float), isinstance(QSO_DEC,float)

#QSO_id = "000017.88+002612.6"
#QSO_RA = 0.07452999800443649
#QSO_DEC = 0.4368380010128021

#QSO_id = '084702.80+013001.4'
#QSO_RA = 131.7616729736328
#QSO_DEC = 1.5004160404205322

#QSO_id = '084940.34+014638.4' 
#QSO_RA = 132.4180908203125
#QSO_DEC = 1.7773410081863403

#QSO_id = '085056.94+014711.0' 
#QSO_RA = 132.73727416992188
#QSO_DEC = 1.7863889932632446

t1 = time.time()
if os.path.exists('fit_result')==False:
    os.mkdir('fit_result')
filename_ascii = 'fit_result/ascii_result.txt'
if_file = glob.glob(filename_ascii)   
if if_file == []:
    ascii =  open(filename_ascii,'w') 
    ascii.write("#ID, RA, DEC, zp, pixel_scale, \
host_flux_ratio(%), host_x, host_y, host_flux, host_mag, \
Reff(arcsec), n_sersic, host_q, Kron_radius(arcsec), host_Kron_flux, host_Kron_mag, \
qso_x, qso_y, qso_flux, qso_mag, host_qso_center_mismatch(arcsec), qso_local_bkg_light(pixel/flux)\
\n")
    
filename_image_not_exit = 'fit_result/image_not_exit.txt'
if_file = glob.glob(filename_image_not_exit) 
if if_file == []:
    image_not_exit = open(filename_image_not_exit,'w') 
    image_not_exit.write("#The following images are not exist:\n")
    
#ID, RA, DEC, zp, pixel_scale, host_flux_ratio(%),host_x, host_y, host_flux, host_mag, Kron_radius(arcsec), host_Kron_flux, host_Kron_mag, Reff(arcsec), n_sersic, host_q,qso_x, qso_y, qso_flux, qso_mag, host_qso_center_mismatch(arcsec)                               
band_run_list = [2,0,1,3,4]  #run I band first
run_list = copy.deepcopy(band_run_list)
band_seq = ['G', 'R', 'I', 'Z', 'Y']
filename_list = [QSO_id+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(5)]
QSO_im_list, err_map_list, QSO_bkg_list, PSF_list=[], [], [], []
zp_list = []
qso_center_list, frm_c_RA_DEC_list = [], []

for i in range(len(band_seq)):
    # The pixel scale is all 0.168
    if len(glob.glob(image_folder+filename_list[i])) == 0:
        image_not_exit = open(filename_image_not_exit,"r+")
        image_not_exit.read()
        print filename_list[i] + " DOES NOT EXIST!!!"
        image_not_exit.write(filename_list[i] + ' DOES NOT EXIST!!!'+ '\n')
        QSO_im, err_map, QSO_bkg, PSF, pix_scale, zp, qso_fr_center, fr_c_RA_DEC = [], [], [], [], [], [], [], []
        run_list.remove(i)
    else:
        QSO_im, err_map, QSO_bkg, PSF, pix_scale, zp, qso_fr_center, fr_c_RA_DEC = gen_fit_id(image_folder, QSO_RA, QSO_DEC, filename_list[i],cut_frame=120, subbkl=sub_bkg)
    QSO_im_list.append(QSO_im)
    err_map_list.append(err_map)
    QSO_bkg_list.append(QSO_bkg)
    qso_center_list.append(qso_fr_center)
    frm_c_RA_DEC_list.append(fr_c_RA_DEC)
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

ct = (len(QSO_im_list[2])-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

#==============================================================================
# Start set up for fitting:
#==============================================================================

filename = 'fit_result/fit_result_for_{0}.txt'.format(QSO_id)
if_file = glob.glob(filename)   
if if_file == []:
    fit_result =  open(filename,'w') 
elif if_file is not []:
    fit_result = open(filename,"r+")
    fit_result.read()

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
    objs, Q_index = detect_obj(QSO_img_l[k],pltshow = pltshow)
    qso_info = objs[Q_index]
    obj_temp = [objs[i] for i in range(len(objs)) if i != Q_index]
    if k == run_list[0]:
        obj = obj_temp
    if k != run_list[0] and len(obj_temp)>0:
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
Iband_inf = []
count = 0
for k in run_list:#len(band_seq)):
    print "Fiting the: "+ filename_list[k]
    psf, QSO_img, QSO_std = psf_l[k], QSO_img_l[k], QSO_std_l[k]
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
        kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0, 'center_y': 0})
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
    tag = 'fit_result/fit_image_{0}'.format(filename_list[k].split('.fits')[0])
    source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                      background_rms=background_rms_list[k],
                                                                      source_params=source_params, QSO_msk = None, fixcenter=fixcenter,
                                                                      pix_sz = pix_scale, no_MCMC = (run_MCMC==False),
                                                                      QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow,
                                                                      corner_plot=False, flux_ratio_plot=True, dump_result=run_MCMC)
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
                                zp=zp_list[k], fixcenter=fixcenter,ID=QSO_id+band_seq[k], QSO_msk = None, tag=tag, plot_compare = plot_compare)
    
    fitsFile = pyfits.open(image_folder+filename_list[k])
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
    plt.savefig('fit_result/fit_image_{0}_host_Kron_curve_objID.pdf'.format(filename_list[k].split('.fits')[0]))
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
        plt.savefig('fit_result/fit_image_{0}_objID.pdf'.format(filename_list[k].split('.fits')[0]))
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
        pyfits.PrimaryHDU(QSO_img-image_ps-image_host[0]-objs_img,header=file_header).writeto(tag+'_data-qso-host-objs(residual).fits',overwrite=True)
    pyfits.PrimaryHDU(QSO_img-image_ps,header=file_header).writeto(tag+'_data-qso.fits',overwrite=True)
    pyfits.PrimaryHDU(QSO_img-image_host[0],header=file_header).writeto(tag+'_data-host.fits',overwrite=True)
    pyfits.PrimaryHDU(QSO_img-image_ps-image_host[0],header=file_header).writeto(tag+'_data-qso-host.fits',overwrite=True)
    pyfits.PrimaryHDU(data_host,header=file_header).writeto(tag+'_data-qso-objs(host).fits',overwrite=True)
    
    ascii = open(filename_ascii,"r+")
    ascii.read()
    ascii.write("{0} {1} {2} {3} {4} ".format(QSO_id+'-'+band_seq[k], frm_c_RA_DEC_list[k][0][0][0], frm_c_RA_DEC_list[k][1][0][0], zp_list[k],0.168))
    ascii.write("{0}% {1} {2} {3} {4} ".format(result['host_flux_ratio_percent'], result['center_x'], result['center_y'],result['host_amp'],result['host_mag']))
    ascii.write("{0} {1} {2} {3} {4} {5} ".format(result['R_sersic'],result['n_sersic'],result['q'],result['Kron_radius'], result['Kron_flux'], result['Kron_mag']))
    qso_mag = -2.5*np.log10(result['QSO_amp']) + zp_list[k]
    if fixcenter == False:
        c_miss =np.sqrt((result['center_x']-result['qso_x'])**2+(result['center_y']-result['qso_y'])**2)
        ascii.write("{0} {1} {2} {3} {4} ".format(result['qso_x'],result['qso_y'],result['QSO_amp'], round(qso_mag,3), round(c_miss,3)))
    else:
        ascii.write("{0} {1} {2} {3} {4} ".format(result['center_x'],result['center_y'],result['QSO_amp'], round(qso_mag,3), 0))
    if isinstance(QSO_bkg_list[k],str): 
        ascii.write("{0}\n".format('NotMeasure'))
    else:
        ascii.write("{0}\n".format(round(np.mean(QSO_bkg_list[2][ct:-ct,ct:-ct]),5)))
    ascii.close()
    
    t2 = time.time()
    time_sp = t2-t1
    time_ave = time_sp/(count + 1)
    time_total = time_ave  * len(run_list)
    t_left= time_total - time_sp
    count = count+1
    print "Finish percent:",time_sp/time_total*100,"%" ,"total time needed :", time_total/60, "mins", "time_left", t_left/60, 'mins'
fit_result.close()

#    ['G', 'R', 'I', 'Z', 'Y']
if len(data_host_list) == 5:
    tag_c = 'fit_result/fit_image_{0}'.format(QSO_id)
    host_g_i = data_host_list[0] -data_host_list[2]
    pyfits.PrimaryHDU(host_g_i,header=file_header).writeto(tag_c+'_host_g-i.fits',overwrite=True)
    host_r_z = data_host_list[1] -data_host_list[3]
    pyfits.PrimaryHDU(host_r_z,header=file_header).writeto(tag_c+'_host_r-z.fits',overwrite=True)
    host_i_y = data_host_list[2] -data_host_list[4]
    pyfits.PrimaryHDU(host_i_y,header=file_header).writeto(tag_c+'_host_i-y.fits',overwrite=True)

