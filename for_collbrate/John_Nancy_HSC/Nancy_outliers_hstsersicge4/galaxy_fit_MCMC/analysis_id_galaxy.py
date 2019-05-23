#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:14 2018

@author: Dartoon

NOTE:
Put this script in folder “analysis_directory/“.
To fit the QSO, you can command in the terminal as follows:
    
python analysis_id.py 084702.80+013001.4

This will allowing you to run the QSO ID 084702.80+013001.4 
Also, you can run two QSO one by one by put two QSO IDs:
    
python analysis_id.py 084702.80+013001.4
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import os
from photutils import make_source_mask
from matplotlib.colors import LogNorm
import copy
import time
import pickle

from gen_fit_id import gen_fit_id
from fit_qso import fit_galaxy
from transfer_to_result import transfer_obj_to_result
from mask_objects import detect_obj
from flux_profile import flux_profile, concentration_profile


t1 = time.time()

deep_seed = 'very_deep'  #Set as True to put more seed and steps to fit,
fix_on_I_band = 1 # input("Do you want, based on I band's inference, fix the other band's Reff and n???\n input 1 means yes, input others means no. \n")
pltshow = 1 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
fitting_strategy =  None # you can set as 'add_galaxy_center_mask' or 'boost_galaxy_center_rms'
sub_bkg = True
image_folder = '../data/'
run_MCMC = True
fit_frame_size =41

galaxy_ID = '18184' 
galaxy_RA = 150.121811
galaxy_DEC = 2.394572

#galaxy_ID = '6310'
#galaxy_RA = 150.154602
#galaxy_DEC = 2.251147

#%%
if os.path.exists('RESULTS')==False:
    os.mkdir('RESULTS')

filename_ascii = 'RESULTS/' + galaxy_ID + '_result.txt'

if_file = glob.glob(filename_ascii)   
if if_file == []:
    ascii =  open(filename_ascii,'w') 
    ascii.write("#ID, RA, DEC, zp, pixel_scale, \
galaxy_x, galaxy_y, galaxy_flux, galaxy_mag, \
Reff(arcsec), n_sersic, galaxy_q, Kron_radius(arcsec), galaxy_Kron_flux, galaxy_Kron_mag, concentration, reduced_Chisq , local_bkg_light(pixel/flux)\
\n")
    ascii.close()
    
band_run_list = [0]  #run I band first
run_list = copy.deepcopy(band_run_list)
band_seq = ['I']
filename_list = [galaxy_ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(1)]
galaxy_im_list, err_map_list, PSF_list, galaxy_bkg_list=[], [], [], []
zp_list = []
galaxy_center_list, frm_c_RA_DEC_list = [], []
    
for i in range(len(band_seq)):
        # The pixel scale is all 0.168
        galaxy_im, err_map, galaxy_bkg, PSF, pix_scale, zp, galaxy_fr_center, fr_c_RA_DEC = gen_fit_id(image_folder, galaxy_RA, galaxy_DEC, filename_list[i],cut_frame=120, subbkl=sub_bkg)
        galaxy_im_list.append(galaxy_im)
        err_map_list.append(err_map)
        galaxy_center_list.append(galaxy_fr_center)
        frm_c_RA_DEC_list.append(fr_c_RA_DEC)
        galaxy_bkg_list.append(galaxy_bkg)        
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
            mask = make_source_mask(galaxy_im_list[i], snr=3, npixels=5, dilate_size=11)
            background_rms_list.append(np.std(galaxy_im_list[i]* (1-mask*1)))
        else:
            background_rms_list.append([])
    
ct = (len(galaxy_im_list[0])-fit_frame_size)/2     # If want to cut to 61, galaxy_im[ct:-ct,ct:-ct]
    
    #==============================================================================
    # Start set up for fitting:
    #==============================================================================
    
filename = 'RESULTS/fit_result_for_{0}.txt'.format(galaxy_ID)
if_file = glob.glob(filename)   
if if_file == []:
        fit_result =  open(filename,'w') 
        fit_result.close()

psf_l, galaxy_img_l, galaxy_std_l = [], [], []
for k in range(len(band_seq)):
    if k in run_list:
        psf_l.append(PSF_list[k])
        galaxy_img_l.append(galaxy_im_list[k][ct:-ct,ct:-ct])
        galaxy_std_l.append(err_map_list[k][ct:-ct,ct:-ct])
    else:
        psf_l.append([])
        galaxy_img_l.append([])
        galaxy_std_l.append([])            

#%%    
from astropy.wcs import WCS
import re
reg_file = glob.glob(image_folder+'*id{0}.dat'.format(galaxy_ID))
if reg_file != []:
    fitsFile = pyfits.open(image_folder+filename_list[run_list[0]])
    wcs = WCS(fitsFile[1].header)  
    reg_name = reg_file[0]
    with open(reg_name, 'r') as input_file:
        reg_string=input_file.read()
    regs_list = reg_string.split('\n')
    regs_list = [regs_list[i] for i in range(len(regs_list)) if '.' in regs_list[i]] 
    regs_list = [re.findall("\d+\.\d+", regs_list[i]) for i in range(len(regs_list))] 
    regs_x_list = [float(regs_list[i][0]) for i in range(len(regs_list))]
    regs_y_list = [float(regs_list[i][1]) for i in range(len(regs_list))]
    re_list = np.asarray([float(regs_list[i][2]) for i in range(len(regs_list))])
    regs_pixpos = np.asarray(wcs.all_world2pix(regs_x_list,regs_y_list,1)).T
    regs_pix = regs_pixpos - galaxy_fr_center
    regs_pix_target = regs_pix[(abs(regs_pix[:,0])<fit_frame_size/2) & (abs(regs_pix[:,1])<fit_frame_size/2)]  # derive all the reg in the frame
    re_target = re_list[(abs(regs_pix[:,0])<fit_frame_size/2) & (abs(regs_pix[:,1])<fit_frame_size/2)]  # derive all the reg in the frame
    selet_galaxy = ((abs(regs_pix_target[:,0])<2) & (abs(regs_pix_target[:,1])<2))
    galaxy_info = [regs_pix_target[selet_galaxy][0], re_target[selet_galaxy][0]] # The galaxy pixel position
    selet_obj = ((abs(regs_pix_target[:,0])>2) & (abs(regs_pix_target[:,1])>2))
    obj_pos = regs_pix_target[selet_obj] # The objects pixel position
    obj_re = re_target[selet_obj]
    obj = [[obj_pos[i], obj_re[i]] for i in range(len(obj_pos))]  #Change the objects position list to right format
    print "the number of nearby objs:", len(obj)
else:
    for k in run_list:
            objs, Q_index = detect_obj(galaxy_img_l[k],pltshow = 1, snr=1.2)
            galaxy_info = objs[Q_index]
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
            print "the number of nearby objs:", len(obj)
plt.imshow(galaxy_img_l[0], origin='low', norm=LogNorm())
rad = len(galaxy_img_l[0])/2
for i in range(len(obj)):
    plt.text(obj[i][0][0]+rad, obj[i][0][1]+rad, "obj{0}".format(i), fontsize=15, color='k')    
plt.show()
#%%    
data_host_list = []
Iband_inf = []
count = 0
for k in run_list:#len(band_seq)):
        print "Fiting the: "+ filename_list[k]
        psf, galaxy_img, galaxy_std = psf_l[k], galaxy_img_l[k], galaxy_std_l[k]
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
            kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -galaxy_info[0][0]*pix_scale, 'center_y': galaxy_info[0][1]*pix_scale})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
            if len(obj) >= 1:
                for i in range(len(obj)):
                    fixed_source.append({})  
                    kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -obj[i][0][0]*pix_scale-5*pix_scale, 'center_y': obj[i][0][1]*pix_scale-5*pix_scale})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': -obj[i][0][0]*pix_scale+5*pix_scale, 'center_y': obj[i][0][1]*pix_scale+5*pix_scale})
        elif band_seq[k] != 'I' and fix_on_I_band == 1: 
            fixed_source.append({'R_sersic': Iband_inf[0]['R_sersic'],'n_sersic': Iband_inf[0]['n_sersic']})  
            kwargs_source_init.append({'R_sersic': Iband_inf[0]['R_sersic'], 'n_sersic': Iband_inf[0]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -galaxy_info[0][0]*pix_scale, 'center_y': galaxy_info[0][1]*pix_scale})
            kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -0.5, 'center_y': -0.5})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': 0.5, 'center_y': 0.5})
            if len(obj) >= 1:
                for i in range(len(obj)):
                    fixed_source.append({'R_sersic': Iband_inf[i+1]['R_sersic'],'n_sersic': Iband_inf[i+1]['n_sersic']})  
                    kwargs_source_init.append({'R_sersic': Iband_inf[i+1]['R_sersic'], 'n_sersic': Iband_inf[i+1]['n_sersic'], 'e1': 0., 'e2': 0., 'center_x': -obj[i][0][0]*pix_scale, 'center_y': obj[i][0][1]*pix_scale})
                    kwargs_source_sigma.append({'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'center_x': -obj[i][0][0]*pix_scale-5*pix_scale, 'center_y': obj[i][0][1]*pix_scale-5*pix_scale})
                    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'center_x': -obj[i][0][0]*pix_scale+5*pix_scale, 'center_y': obj[i][0][1]*pix_scale+5*pix_scale})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
        #==============================================================================
        # to fit and save the inference
        #==============================================================================
        img_c = fit_frame_size/2
        galaxy_msk = None
        if fitting_strategy == 'boost_galaxy_center_rms':
            galaxy_std[galaxy_img == galaxy_img[img_c-3:img_c+4,img_c-3:img_c+4].max()] = 10**6
        elif fitting_strategy == 'add_galaxy_center_mask':
            galaxy_msk = np.ones_like(galaxy_img)
            galaxy_msk[galaxy_img == galaxy_img[img_c-3:img_c+4,img_c-3:img_c+4].max()] = 0
        tag = 'RESULTS/fit_image_{0}'.format(filename_list[k].split('.fits')[0])
        source_result, image_host, error_map, reduced_Chisq=fit_galaxy(galaxy_img, psf_ave=psf, psf_std = None,
                                                      background_rms=background_rms_list[k],
                                                      source_params=source_params, galaxy_msk = galaxy_msk,
                                                      pix_sz = pix_scale, no_MCMC = (run_MCMC==False),
                                                      galaxy_std =galaxy_std, tag=tag, deep_seed= deep_seed,
                                                      pltshow=pltshow, return_Chisq=True,  corner_plot=False, 
                                                      dump_result=run_MCMC, flux_corner_plot = True)
        if band_seq[k] == 'I':
            Iband_inf = source_result # in order to save the I band source_reslt to fix Reff and n for other band.

        if pltshow == 0:
            plot_compare=False
            fits_plot =False
        else:
            plot_compare=True
            fits_plot =True
        result = transfer_obj_to_result(source_result[0],image_host[0], zp)
        result['reduced_Chisq'] = reduced_Chisq
        
        fitsFile = pyfits.open(image_folder+filename_list[k])
        file_header = copy.deepcopy(fitsFile[1].header)
        file_header['CRPIX1'] = file_header['CRPIX1']-galaxy_center_list[k][0]+len(galaxy_img)/2
        file_header['CRPIX2'] = file_header['CRPIX2']-galaxy_center_list[k][1]+len(galaxy_img)/2
        
        fit_result = open(filename,"r+")
        fit_result.read()
        fit_result.write("#galaxy_img intensity in the band {1}: {0} \n".format(round(np.sum(galaxy_img),2), band_seq[k]))
        fit_result.write("#fit result of the galaxy by in band {0}: \n".format(band_seq[k]))
        objs_img = np.zeros_like(image_host[0])
        if len(source_result)>1:
            for i in range(1,len(image_host)):
                objs_img += image_host[i]
        from flux_profile import galaxy_total_compare
        flux_list = [galaxy_img, image_host[0]*0, image_host[0]+objs_img, error_map]
        fig = galaxy_total_compare(label_list = ['data', 'model', 'normalized residual'], flux_list = flux_list, target_ID = galaxy_ID, pix_sz=pix_scale,
                  zp = zp, plot_compare=plot_compare, msk_image = np.ones_like(image_host[0]))
        
        data_host = galaxy_img - objs_img
        data_host_list.append(data_host)
        center = np.asarray(data_host.shape) /2
        r_flux, r_grids, regions=flux_profile(data_host, center,
                                              radius=center.min(),
                                              grids=100, ifplot=False,
                                              fits_plot= False, mask_list=None)
        plt.plot(r_grids*pix_scale, r_flux)
        plt.title('Kron curve')
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('Kron flux')
        plt.savefig('RESULTS/fit_image_{0}_host_Kron_curve_objID.pdf'.format(filename_list[k].split('.fits')[0]))
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
        Kron_radius = 3.
        print "Kron_radius take {0} arcsec".format(Kron_radius)
        if Kron_radius/pix_scale < fit_frame_size/2:
            Kron_flux = r_flux[r_grids>Kron_radius/pix_scale][0]
        else:
            print "Frame size smaller than the Kron_radius!!!".format(Kron_radius)
            Kron_flux = r_flux[-1]
        Kron_mag = -2.5*np.log10(Kron_flux) + zp_list[k]
        result['Kron_radius'] = Kron_radius
        result['Kron_flux'] = round(Kron_flux,3)
        result['Kron_mag'] = round(Kron_mag,3)
        r20, r80 = concentration_profile(data_host, center, total_flux = image_host[0].sum(),
                                    radius=center.min(), grids=100, ifplot=pltshow,fits_plot= pltshow)
        concentration = round(5 * np.log10(r80/r20),3)
        result['concentration'] = concentration
        fit_result.write(repr(result) + "\n")
        if len(source_result)>1:
            for i in range(1,len(source_result)):
                result_i = transfer_obj_to_result(source_result[i],image_host[i], zp)
                fit_result.write("#obj {0}: \n".format(i))
                fit_result.write(repr(result_i) + "\n")
                plt.imshow(galaxy_img, origin='low', norm=LogNorm())
                obj_x, obj_y = len(galaxy_img)/2 - source_result[i]['center_x']/pix_scale, len(galaxy_img)/2+source_result[i]['center_y']/pix_scale
                plt.text(obj_x, obj_y, "obj{0}".format(i), fontsize=15, color='k')
            plt.savefig('RESULTS/fit_image_{0}_objID.pdf'.format(filename_list[k].split('.fits')[0]))
            if pltshow == 0:
                plt.close()
            else:
                plt.show()
            pyfits.PrimaryHDU(galaxy_img-image_host[0]-objs_img,header=file_header).writeto(tag+'_data-host-objs(residual).fits',overwrite=True)
        pyfits.PrimaryHDU(galaxy_img-image_host[0],header=file_header).writeto(tag+'_data-host.fits',overwrite=True)
        pyfits.PrimaryHDU(data_host,header=file_header).writeto(tag+'_data-objs(host).fits',overwrite=True)
        fit_result.close()
        ascii = open(filename_ascii,"r+")
        ascii.read()
        ascii.write("{0} {1} {2} {3} {4} ".format(str(filename_list[k].split('.fits')[0]),frm_c_RA_DEC_list[k][0][0][0], frm_c_RA_DEC_list[k][1][0][0], zp_list[k],0.168))
        ascii.write("{0} {1} {2} {3} ".format(result['center_x'], result['center_y'],result['host_amp'],result['host_mag']))
        ascii.write("{0} {1} {2} {3} {4} {5} {6} {7} ".format(result['R_sersic'],result['n_sersic'],result['q'],result['Kron_radius'], result['Kron_flux'],
                    result['Kron_mag'], result['concentration'], result['reduced_Chisq'],))
        if isinstance(galaxy_bkg_list[k],str): 
            ascii.write("{0}\n".format('NotMeasure'))
        else:
            ascii.write("{0}\n".format(round(np.mean(galaxy_bkg_list[k][ct:-ct,ct:-ct]),5)))
        ascii.close()
        t2 = time.time()
        time_sp = t2-t1
        time_ave = time_sp/(k + 1)
        time_total = time_ave * len(run_list)
        t_left= time_total - time_sp
        count = count+1
        print "Finish percent:",time_sp/time_total*100,"%" ,"total time needed :", time_total/60, "mins", "time_left", t_left/60, 'mins'

#%%    
#    ['G', 'R', 'I', 'Z', 'Y']
if len(data_host_list) == 5:
        tag_c = 'RESULTS/fit_image_{0}'.format(galaxy_ID)
        host_g_i = data_host_list[0] -data_host_list[2]
        pyfits.PrimaryHDU(host_g_i,header=file_header).writeto(tag_c+'_host_g-i.fits',overwrite=True)
        host_r_z = data_host_list[1] -data_host_list[3]
        pyfits.PrimaryHDU(host_r_z,header=file_header).writeto(tag_c+'_host_r-z.fits',overwrite=True)
        host_i_y = data_host_list[2] -data_host_list[4]
        pyfits.PrimaryHDU(host_i_y,header=file_header).writeto(tag_c+'_host_i-y.fits',overwrite=True)

#%%
        #Save the mag together with error level into ascii file:
filename_ascii = 'RESULTS/ascii_mag_error.txt'
if_file = glob.glob(filename_ascii)   
if if_file == []:
    ascii_err =  open(filename_ascii,'w') 
    ascii_err.write("#ID, mag, mag_err\n")
    ascii_err.close()
ascii_err = open(filename_ascii,"r+")
ascii_err.read()
for band in band_seq:
    if_file = glob.glob("RESULTS/*{0}*{1}.pkl".format(galaxy_ID,band))
    m_m, m_err = -99, -99
    if if_file == []:
        print "MCMC for {0}-{1} not performed!".format(galaxy_ID,band)
    else:
        picklename = 'RESULTS/'+'fit_image_'+galaxy_ID+ '_HSC-{0}.pkl'.format(band)
        result = pickle.load(open(picklename,'rb'))
        [source_result, image_host, samples_mcmc, param_mcmc, paras, chain_list, param_list] = result
        [source_params_2, mcmc_new_list, labels_new] = paras
        mcmc_new_list = np.asarray(mcmc_new_list)
        idx = 0     #The translated flux for the host
        v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
        v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
        v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
        #print labels_new[idx], ":", v_l, v_m, v_h
        m_l = -2.5 * np.log10(v_h) + zp
        m_m = -2.5 * np.log10(v_m) + zp
        m_h = -2.5 * np.log10(v_l) + zp
        m_err = ((m_h-m_m)**2 + (m_m-m_l)**2)**0.5
    ascii_err.write("{0} {1} {2}\n".format(galaxy_ID+'-'+band, round(m_m,3), round(m_err,3)))
ascii_err.close()