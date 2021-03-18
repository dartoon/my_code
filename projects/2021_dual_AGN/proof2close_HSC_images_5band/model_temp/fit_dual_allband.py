#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:14 2018

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os, copy, glob
from subprocess import call
from matplotlib.colors import LogNorm
from decomprofile.tools.measure_tools import find_loc_max, measure_FWHM, twoD_Gaussian, fit_data_twoD_Gaussian #, esti_bgkstd
from decomprofile.data_process import DataProcess
from decomprofile.fitting_specify import FittingSpeficy
from decomprofile.fitting_process import FittingProcess
from astropy.wcs import WCS

# f = open("../cut_out.txt","r")
# string = f.read()
# lines = string.split('\n')   # Split in to \n

# files = glob.glob('../z_over1/*')
# files.sort()

# run_i = 23
# file = files[run_i]
# image_ID = file.split('/')[-1]
# line = [lines[i] for i in range(len(lines)) if image_ID in lines[i]]

# _, image_RA, image_DEC = line[0].split(' ')
# image_RA = float(image_RA)
# image_DEC = float(image_DEC)
# print("run_i = ", run_i)
# print(image_ID, image_RA, image_DEC)

image_ID = '144034.78+441520.5'
f = open("../../_pdfs_2close/DR144.4_short.asc","r")
string = f.read()
lines_ = string.split('\n')   # Split in to \n
line= [lines_[i] for i in range(len(lines_)) if image_ID in lines_[i]]
image_RA, image_DEC = float(line[0].split(' ')[1]), float(line[0].split(' ')[2]), 
# 
#%%
deep_seed = True  #Set as True to put more seed and steps to fit,
show_plot = 1
fit_data = True  #If you simply want to do the search without fitting, set False

image_folder = '../z_below1/' + image_ID + '/'
fit_folder = image_folder + 'fit_result/'

print(fit_folder)
if os.path.exists(fit_folder)==False:
    os.mkdir(fit_folder)

import shutil
shutil.copy('./fit_dual_allband.py', fit_folder)
#If only want to run I band
# band_seq = ['I'] 
# run_list = [0] 
band_seq = ['I', 'G', 'R', 'Z', 'Y']
run_list = [0, 1, 2, 3, 4]
filename_list = [image_ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(len(band_seq))]

data_process_list, zp_list = [], []
for i in range(len(band_seq)):
    # The pixel scale is all 0.168
    if len(glob.glob(image_folder+filename_list[i])) == 0:
        print(filename_list[i] + " DOES NOT EXIST!!!")
        QSO_im, err_map, PSF, _, _, qso_center, fr_c_RA_DEC = [], [], [], [], [], [], []
        run_list.remove(i)
        data_process_list.append(None)
        zp_list.append(None)
    else:
        fitsFile = pyfits.open(image_folder+filename_list[i])
        fov_image= fitsFile[1].data
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        err_data= fitsFile[3].data ** 0.5
        file_header0 = fitsFile[0].header
        zp =  27.0 
        data_process_i = DataProcess(fov_image = fov_image, fov_noise_map = err_data,
                                     target_pos = [image_RA, image_DEC],
                                     pos_type = 'wcs', header = header,
                                     rm_bkglight = True, if_plot=False, zp = zp)
        data_process_i.noise_map = err_data
        data_process_i.generate_target_materials(radius=None, radius_list = [10, 20, 30, 40],
                                                  create_mask = False, nsigma=1,
                                              exp_sz= 1.2, npixels = 9, if_plot=False)        
        PSF = pyfits.getdata(image_folder+filename_list[i].split('.fits')[0]+'_psf.fits')
        if len(PSF) != 0 and PSF.shape[0] != PSF.shape[1]:
            cut = int((PSF.shape[0] - PSF.shape[1])/2)
            if cut>0:
                PSF = PSF[cut:-cut,:]
            elif cut<0:
                PSF = PSF[:,-cut:cut]
            PSF /= PSF.sum()
            if PSF.shape[0] != PSF.shape[1]:
                raise ValueError("PSF shape is not a square.")
        data_process_i.PSF_list = [PSF]
        data_process_list.append(data_process_i)
        zp_list.append(zp) 

#%%
frame = np.max([len(data_process_list[i].target_stamp) for i in run_list])
radius = int((frame-1)/2)
for k in run_list:
    data_process_list[k].generate_target_materials(radius=radius,
                                                   create_mask = False, nsigma=2.,
                                                   exp_sz= 1.2, npixels = 15, if_plot=False)
    apertures_temp = data_process_list[k].apertures
    if k == run_list[0]:
        apertures = apertures_temp
    if k != run_list[0] and len(apertures_temp)>1:
        for i in range(len(apertures_temp)):
            count = 0
            for j in range(len(apertures)):
                dis = np.sqrt(np.sum(apertures[j].positions-apertures_temp[i].positions)**2)
                if dis < 5:  #If objects is close within 5 pixels consider as a sample obj
                    count += 1
            if count == 0:
                apertures.append(apertures_temp[i])
# del apertures[-1]
for k in run_list:    
    data_process_list[k].apertures = apertures #Pass apertures to the data
fix_n_list_0, fix_Re_list_0, fix_n_list_1, fix_Re_list_1, fix_n_list_2, fix_Re_list_2 = [None] * 6
ps_param0, ps_param1, ps_param2 = None, None, None

# from decomprofile.fitting_specify import ps_params_generator
# ps_param1 = ps_params_generator([[0,0], [-1,-2]], [30,30], deltaPix = data_process_list[0].deltaPix)
# ps_param2 = ps_param1

#%%
fig, (axs) = plt.subplots(1, 5, figsize=(15, 7))
vmin = 1.e-3
vmax = data_process_list[0].target_stamp.max() * 5
color_list = ['winter', 'summer', 'afmhot', 'autumn', 'gist_heat']
plt_list = [1, 2, 0, 3, 4]
for i in range(len(plt_list)):
    p_i = plt_list[i]
    if data_process_list[p_i] != None:
        axs[i].imshow(data_process_list[p_i].target_stamp, origin='lower', cmap=color_list[i], norm=LogNorm(), vmin=vmin, vmax=vmax)
    else:
        axs[i].imshow(data_process_list[0].target_stamp * 0)
    axs[i].set_title('{0} band'.format(band_seq[p_i]))
plt.savefig(fit_folder + 'images_5_band.pdf')
plt.show()   

#%%
if deep_seed ==True:
    pso_setting = None
else:
    pso_setting = {'sigma_scale': 0.8, 'n_particles': 60, 'n_iterations': 60}

num_BHBH = 2
neighborhood_size, threshold = 4, 4

for k in run_list:  #Search for Dual Based on I band.
    print("Fiting the: "+ filename_list[k])
    # print("Comparing the fitting Chisq:")
    qso_center = data_process_list[k].target_pos
    file_header = copy.deepcopy(fitsFile[1].header)
    file_header['CRPIX1'] = file_header['CRPIX1']-qso_center[0]+int(len(data_process_list[k].target_stamp)/2)
    file_header['CRPIX2'] = file_header['CRPIX2']-qso_center[1]+int(len(data_process_list[k].target_stamp)/2)        
    wcs = WCS(file_header)
    write_result =  open(fit_folder+'fit_result_{0}-band.txt'.format(band_seq[k]),'w') 
    write_result.write("#The fitting information for {0}:\n".format(filename_list[k]))

    #==============================================================================
    # fitting the QSO as a BH (PS) + Sersic       
    #==============================================================================
    print("fitting the QSO as one BH + Sersic ")
    tag = fit_folder + 'fit_{0}-band_fit0_PS+Sersic'.format(band_seq[k])
    _fit_sepc_0 = FittingSpeficy(data_process_list[k])
    _fit_sepc_0.prepare_fitting_seq(point_source_num = 1, fix_n_list= fix_n_list_0, fix_Re_list= fix_Re_list_0, ps_params = ps_param0)
    _fit_sepc_0.build_fitting_seq()
    if k == run_list[0]:
        _fit_sepc_0.plot_fitting_sets(savename = fit_folder + 'fitting0_used_aper.pdf'.format(image_ID),
                                    show_plot=show_plot)         
    
    _fit_run_0 = FittingProcess(_fit_sepc_0, savename = tag)
    _fit_run_0.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
    _fit_run_0.translate_result()
    _fit_run_0.plot_final_qso_fit(target_ID = image_ID, save_plot = True, show_plot = show_plot)
    source_result_0, ps_result_0 = _fit_run_0.final_result_galaxy, _fit_run_0.final_result_ps
    host_mag, AGN_mag = source_result_0[0]['magnitude'], ps_result_0[0]['magnitude']
    c_miss = np.sqrt((source_result_0[0]['center_x']-ps_result_0[0]['ra_image'])**2+(source_result_0[0]['center_y']-ps_result_0[0]['dec_image'])**2)
    reduced_Chisq_0 = _fit_run_0.reduced_Chisq
    write_result.write("1. Fitting as a regular QSO,i.e. one PS + Sersic:\n")
    write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_0,3)))
    write_result.write("\nHost mag: "+repr(round(host_mag,3)))
    write_result.write("\nAGN mag: "+repr(round(AGN_mag,3)))
    write_result.write("\nPS Sersic center offset (arcsec): "+repr(round(float(c_miss),3)) + "; ")
    write_result.write("\nmodel_Sersic_result: "+repr(source_result_0))
    write_result.write("\nmodel_PS_result: "+repr(ps_result_0))
    write_result.write("\n=======================================================\n")
    tag_name = tag + "_qso_final_plot"
    print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_0,1)))+'.pdf', shell=True))
    if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.        
        fix_Re_list_0 = [[i, source_result_0[i]['R_sersic']] for i in range(len(source_result_0))] 
        fix_n_list_0 = [[i, source_result_0[i]['n_sersic']] for i in range(len(source_result_0))]
        ps_param0 =  _fit_sepc_0.kwargs_params['point_source_model']    
    #==============================================================================
    # fitting the QSO as a BHBH (PS+PS)        
    #==============================================================================
    print("fitting the QSO as {0} point sources".format(num_BHBH))
    tag = fit_folder +'fit_{0}-band_fit1_PSPS'.format(band_seq[k])
    _fit_sepc_1 = FittingSpeficy(data_process_list[k])
    del _fit_sepc_1.apertures[0]
    _fit_sepc_1.prepare_fitting_seq(point_source_num = num_BHBH, fix_n_list= fix_n_list_1, 
                                    fix_Re_list= fix_Re_list_1, ps_params = ps_param1,
                                    neighborhood_size = neighborhood_size, threshold = threshold)
    _fit_sepc_1.build_fitting_seq()
    if k == run_list[0]:
        _fit_sepc_1.plot_fitting_sets(savename = fit_folder + 'fitting1_used_aper.pdf'.format(image_ID),
                                    show_plot=show_plot)        
    
    _fit_run_1 = FittingProcess(_fit_sepc_1, savename = tag)
    _fit_run_1.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
    _fit_run_1.translate_result()
    _fit_run_1.plot_final_qso_fit(target_ID = image_ID, save_plot = True, show_plot = show_plot)            
    source_result_1, ps_result_1 = _fit_run_1.final_result_galaxy, _fit_run_1.final_result_ps 
    AGN_mags = [ps_result_1[i]['magnitude'] for i in range(len(ps_result_1))]
    if len(ps_result_1) == 2:
        c_miss = np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)
    elif len(ps_result_1) > 2:
        c_miss = [np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)]
        c_miss.append(np.sqrt((ps_result_1[1]['ra_image']-ps_result_1[2]['ra_image'])**2+(ps_result_1[1]['dec_image']-ps_result_1[2]['dec_image'])**2))
        c_miss.append(np.sqrt((ps_result_1[2]['ra_image']-ps_result_1[0]['ra_image'])**2+(ps_result_1[2]['dec_image']-ps_result_1[0]['dec_image'])**2))
        c_miss = np.average(c_miss)
    reduced_Chisq_1 = _fit_run_1.reduced_Chisq
    write_result.write("2. Fitting as {0}PS:\n".format(len(ps_result_1)))
    write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_1,3)))
    write_result.write("\nAGN mag: ")
    for i in range(len(ps_result_1)):
        write_result.write(repr(round(AGN_mags[i],3))+' ')
    write_result.write("\n")
    for i in range(len(ps_result_1)):
        write_result.write("AGN{0} position: ".format(i))
        ps_x = -ps_result_1[i]['ra_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        ps_y = ps_result_1[i]['dec_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        ra, dec = wcs.all_pix2world([[ps_x, ps_y]], 1)[0]
        write_result.write("x: "+repr(round(ps_x,3))+' y: '+repr(round(ps_y,3))+ " RA: "+repr(ra)+' DEC: '+repr(dec)+ "; ")         
    write_result.write("\nPS PS center offset (arcsec): "+repr(round(float(c_miss),3)))
    write_result.write("\nmodel_Sersic_result: "+repr(source_result_1))
    write_result.write("\nmodel_PS_result: "+repr(ps_result_1))        
    write_result.write("\n=======================================================\n")
    tag_name = tag + "_qso_final_plot"  
    print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_1,1)))+'.pdf', shell=True))
    if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.     
        fix_Re_list_1 = [[i, source_result_1[i]['R_sersic']] for i in range(len(source_result_1))] 
        fix_n_list_1 = [[i, source_result_1[i]['n_sersic']] for i in range(len(source_result_1))]
        ps_param1 =  _fit_sepc_1.kwargs_params['point_source_model']    
    #==============================================================================
    # fitting the QSO as a BHBH (PS+PS) + Sersic       
    #==============================================================================
    print("fitting the QSO as {0} point sources + Sersic".format(num_BHBH))
    tag = fit_folder +'fit_{0}-band_fit2_PSPS+Sersic'.format(band_seq[k])
    _fit_sepc_2 = FittingSpeficy(data_process_list[k])
    _fit_sepc_2.prepare_fitting_seq(point_source_num = num_BHBH, fix_n_list= fix_n_list_2,
                                    fix_Re_list= fix_Re_list_2, ps_params = ps_param2,
                                    neighborhood_size = neighborhood_size, threshold = threshold)
    _fit_sepc_2.build_fitting_seq()
    if k == run_list[0]:
        _fit_sepc_2.plot_fitting_sets(savename = fit_folder + 'fitting2_used_aper.pdf'.format(image_ID),
                                    show_plot=show_plot)        
    _fit_run_2 = FittingProcess(_fit_sepc_2, savename = tag)
    _fit_run_2.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
    _fit_run_2.translate_result()
    _fit_run_2.plot_final_qso_fit(target_ID = image_ID, save_plot = True, show_plot = show_plot)     
    source_result_2, ps_result_2 = _fit_run_2.final_result_galaxy, _fit_run_2.final_result_ps 
    host_mag = source_result_2[0]['magnitude']
    AGN_mags = [ps_result_2[i]['magnitude'] for i in range(len(ps_result_2))]
    if len(ps_result_2) == 2:
        c_miss = np.sqrt((ps_result_2[0]['ra_image']-ps_result_2[1]['ra_image'])**2+(ps_result_2[0]['dec_image']-ps_result_2[1]['dec_image'])**2)
    elif len(ps_result_2) > 2:
        c_miss = [np.sqrt((ps_result_2[0]['ra_image']-ps_result_2[1]['ra_image'])**2+(ps_result_2[0]['dec_image']-ps_result_2[1]['dec_image'])**2)]
        c_miss.append(np.sqrt((ps_result_2[1]['ra_image']-ps_result_2[2]['ra_image'])**2+(ps_result_2[1]['dec_image']-ps_result_2[2]['dec_image'])**2))
        c_miss.append(np.sqrt((ps_result_2[2]['ra_image']-ps_result_2[0]['ra_image'])**2+(ps_result_2[2]['dec_image']-ps_result_2[0]['dec_image'])**2))
        c_miss = np.average(c_miss)
    reduced_Chisq_2 = _fit_run_2.reduced_Chisq
    write_result.write("3. Fitting as {0}PS + Sersic:\n".format(len(ps_result_2)))
    write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_2,3)))
    write_result.write("\nHost mag: "+repr(round(host_mag,3)))
    write_result.write("\nAGN mag: ")
    for i in range(len(ps_result_2)):
        write_result.write(repr(round(AGN_mags[i],3))+' ')
    write_result.write("\n")
    for i in range(len(ps_result_2)):
        write_result.write("AGN{0} position: ".format(i))
        ps_x = -ps_result_2[i]['ra_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        ps_y = ps_result_2[i]['dec_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        ra, dec = wcs.all_pix2world([[ps_x, ps_y]], 1)[0]
        write_result.write("x: "+repr(round(ps_x,3))+' y: '+repr(round(ps_y,3))+ " RA: "+repr(ra)+' DEC: '+repr(dec)+ "; ")
    write_result.write("\nPS PS center offset (arcsec): "+repr(round(float(c_miss),3)))
    write_result.write("\nmodel_Sersic_result: "+repr(source_result_2))
    write_result.write("\nmodel_PS_result: "+repr(ps_result_2))            
    write_result.write("\n=======================================================\n")
    tag_name = tag + "_qso_final_plot"  
    image_host_2, image_ps_2 = _fit_run_2.image_host_list, _fit_run_2.image_ps_list
    _image_ps_2 = np.zeros_like(image_ps_2[0])
    for i in range(len(image_ps_2)):
        _image_ps_2 += image_ps_2[i]
    objs_img = np.zeros_like(image_host_2[0])
    if len(image_host_2)>1:
        for i in range(1,len(image_host_2)):
            objs_img += image_host_2[i]
    fitsFile = pyfits.open(image_folder+filename_list[k])
    pyfits.PrimaryHDU(data_process_list[k].target_stamp-_image_ps_2-objs_img,header=file_header).writeto(fit_folder+'data-BHBH(host image)_{1}-band.fits'.format(image_ID, band_seq[k]),overwrite=True)
    print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_2,1)))+'.pdf', shell=True))
    if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.
        fix_Re_list_2 = [[i, source_result_2[i]['R_sersic']] for i in range(len(source_result_2))] 
        fix_n_list_2 = [[i, source_result_2[i]['n_sersic']] for i in range(len(source_result_2))]
        ps_param2 = _fit_sepc_2.kwargs_params['point_source_model']  
    write_result.close() 
#os.system('say "your program has finished"')
print("Program has finished")
#
