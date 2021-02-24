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

image_ID = sys.argv[1] #'141637.44+003352.2' 
image_RA = float(sys.argv[2]) #214.15602111816406
image_DEC = float(sys.argv[3]) #0.5645210146903992

# image_ID ='100043.13+020637.2' 
# image_RA = 150.1797789
# image_DEC = 2.110369603

# image_ID ='010048.81+021604.0'
# image_RA = 15.2033809 
# image_DEC = 2.2677925

# image_ID ='141637.44+003352.2' 
# image_RA = 214.156021
# image_DEC = 0.564521

# # run_line = [5069, 5070, 5086, 5093, 5096, 5110][1]
# run_line = 0
# targets_info = open("QSO_id_file.txt", "r")
# targets_info = targets_info.read().split('\n')
# info = copy.deepcopy(targets_info[run_line])
# info = info.split(' ')
# info = [info[i] for i in range(len(info)) if info[i] != '']
# image_ID, image_RA, image_DEC = info[:3]
# image_RA, image_DEC = float(image_RA), float(image_DEC)

#%%
print(image_ID, image_RA, image_DEC)

deep_seed = False  #Set as True to put more seed and steps to fit,
show_plot = 0
fit_data = True  #If you simply want to do the search without fitting, set False

image_folder = '../images_directory/'
    
if os.path.exists('fit_result_detect')==False:
    os.mkdir('fit_result_detect')

filename_ascii = 'RESULTS/' + image_ID + '_result.txt'

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
        # FLUXMAG0 = file_header0['FLUXMAG0']
        zp =  27.0 #2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
        
        data_process_i = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [image_RA, image_DEC],
                                   pos_type = 'wcs', header = header,
                                  rm_bkglight = True, if_plot=False, zp = zp)
        data_process_i.noise_map = err_data
        data_process_i.generate_target_materials(radius=None, create_mask = False, nsigma=2.8,
                                              exp_sz= 1.2, npixels = 15, if_plot=False)        
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
k = 0  # k = 0 means use I band to select Dual
QSO_img = data_process_list[k].target_stamp

neighborhood_size, threshold = 4, 4
x, y = find_loc_max(QSO_img, neighborhood_size = neighborhood_size, threshold = threshold)
arr_x, arr_y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
center = len(QSO_img)/2
bool_x, bool_y = (arr_x>(center-18))*(arr_x<(center+18)), (arr_y>(center-18))*(arr_y<(center+18))
arr_x = arr_x[bool_x*bool_y]
arr_y = arr_y[bool_x*bool_y]
#Remove arr_x's elememt if the corresponding pixel is too faint:
arr_bool = np.array([True] * len(arr_x), dtype=bool)
for i in range(len(arr_x)):
    if QSO_img[int(arr_y[i]), int(arr_x[i])] < 3.0 :
        arr_bool[i] = False
arr_x = arr_x[arr_bool]

qsoid = filename_list[k].split('_')[0]
claim = ''
if len(arr_x)>=2:
    if os.path.exists('fit_result_detect/{0}/'.format(qsoid))==False:
        os.mkdir('fit_result_detect/{0}/'.format(qsoid))
    claim = claim+"\nThis {0} is likely to be a {1} system (based on multi-peaks)!!!".format(filename_list[k], 'BH'*len(arr_x))
    plt.imshow(QSO_img, origin='lower', norm=LogNorm())
    for i in range(len(arr_x)):
        plt.text(arr_x[i], arr_y[i],'BH{0}'.format(i))
        plt.plot(arr_x[i], arr_y[i],'ro')
    plt.savefig('fit_result_detect/{0}/proof-BHBH.pdf'.format(qsoid))
    if show_plot == 1:
        plt.show()
    else:
        plt.close()

# #Fit central as twoD Gaussian Anyway.            
# twoD_Gau_p_PSF =  fit_data_twoD_Gaussian(data_process_list[k].PSF_list[0])
# frz = int(center/2)
# twoD_Gau_p_data = fit_data_twoD_Gaussian(QSO_img[frz:-frz,frz:-frz])
# q_PSF = twoD_Gau_p_PSF[3]/twoD_Gau_p_PSF[4]
# q_PSF = min(q_PSF, 1/q_PSF)
# q_data = twoD_Gau_p_data[3]/twoD_Gau_p_data[4]
# q_data = min(q_data, 1/q_data)
# if abs((q_data-q_PSF)/q_PSF) > 0.15 :   #!!! Set the level as 15% mismatch to PSF
#     if os.path.exists('fit_result_detect/{0}/'.format(qsoid))==False:
#         os.mkdir('fit_result_detect/{0}/'.format(qsoid))
#     claim = claim+"\nThis {0} is also likely to have closed dual AGN pair (based on FWHM)!!!".format(filename_list[k])
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 10))
#     ax1.imshow(QSO_img[frz:-frz,frz:-frz], origin='lower', cmap='gist_heat', norm=LogNorm(), vmin=1.e-4 , vmax = np.max(QSO_img) )
#     ax1.set_title('Data')
#     QSO_2D_fitted = twoD_Gaussian(len(QSO_img[frz:-frz,frz:-frz]), *twoD_Gau_p_data)
#     ax2.imshow(QSO_2D_fitted.reshape(len(QSO_img[frz:-frz,frz:-frz]), len(QSO_img[frz:-frz,frz:-frz])), origin='lower', cmap='gist_heat', norm=LogNorm(), vmin=1.e-4 , vmax = np.max(QSO_img) )
#     ax2.set_title('fitted Gaussian Image')
#     plt.savefig('fit_result_detect/{0}/proof-2close.pdf'.format(qsoid))
#     if show_plot == 1:
#         plt.show()
#     else:
#         plt.close()

if os.path.exists('fit_result_detect/{0}/'.format(qsoid)) == True and fit_data == True:  #Start to fit the data:
    #Define the right frame to fit:
    print(claim)
    
    #Re-define the radius to make the target the have the same frame size.
    frame = np.max([len(data_process_list[i].target_stamp) for i in run_list])
    radius = int((frame-1)/2)
    for k in run_list:
        data_process_list[k].generate_target_materials(radius=radius, create_mask = False, nsigma=2.8,
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
        data_process_list[k].apertures = apertures #Pass apertures to the data
    fix_n_list_0, fix_Re_list_0, fix_n_list_1, fix_Re_list_1, fix_n_list_2, fix_Re_list_2 = [None] * 6
    ps_param0, ps_param1, ps_param2 = None, None, None
    
    if deep_seed ==True:
        pso_setting = None
    else:
        pso_setting = {'sigma_scale': 0.8, 'n_particles': 60, 'n_iterations': 60}
    
    for k in run_list:  #Search for Dual Based on I band.
        print("Fiting the: "+ filename_list[k])
        # print("Comparing the fitting Chisq:")
        qso_center = data_process_list[k].target_pos
        file_header = copy.deepcopy(fitsFile[1].header)
        file_header['CRPIX1'] = file_header['CRPIX1']-qso_center[0]+int(len(data_process_list[k].target_stamp)/2)
        file_header['CRPIX2'] = file_header['CRPIX2']-qso_center[1]+int(len(data_process_list[k].target_stamp)/2)        
        wcs = WCS(file_header)
        
        write_result =  open('fit_result_detect/{0}/fit_result_{1}-band.txt'.format(qsoid,band_seq[k]),'w') 
        write_result.write("#The fitting information for {0}:\n".format(filename_list[k]))
        # #==============================================================================
        # # fitting the QSO as a BH (PS) + Sersic       
        # #==============================================================================
        # print("fitting the QSO as one BH + Sersic ")
        # fit_time = 1
        # tag = 'fit_result_detect/{0}/fit_{2}-band_image0_PS+Sersic_fittime-{1}'.format(qsoid,fit_time+1, band_seq[k])
        # _fit_sepc_0 = FittingSpeficy(data_process_list[k])
        # _fit_sepc_0.prepare_fitting_seq(point_source_num = 1, fix_n_list= fix_n_list_0, fix_Re_list= fix_Re_list_0, ps_params = ps_param0)
        # _fit_sepc_0.build_fitting_seq()
        # _fit_run_0 = FittingProcess(_fit_sepc_0, savename = tag)
        # _fit_run_0.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
        # _fit_run_0.translate_result()
        # _fit_run_0.plot_final_qso_fit(target_ID = qsoid, save_plot = True, show_plot = show_plot)
        # source_result_0, ps_result_0 = _fit_run_0.final_result_galaxy, _fit_run_0.final_result_ps
        # host_mag, AGN_mag = source_result_0[0]['magnitude'], ps_result_0[0]['magnitude']
        # c_miss = np.sqrt((source_result_0[0]['center_x']-ps_result_0[0]['ra_image'])**2+(source_result_0[0]['center_y']-ps_result_0[0]['dec_image'])**2)
        # reduced_Chisq_0 = _fit_run_0.reduced_Chisq
        # write_result.write("1. Fitting as a regular QSO,i.e. one PS + Sersic:\n")
        # write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_0,3)))
        # write_result.write("\nHost mag: "+repr(round(host_mag,3)))
        # write_result.write("\nAGN mag: "+repr(round(AGN_mag,3)))
        # write_result.write("\nPS Sersic center offset (arcsec): "+repr(round(float(c_miss),3)) + "; ")
        # write_result.write("\nmodel_Sersic_result: "+repr(source_result_0))
        # write_result.write("\nmodel_PS_result: "+repr(ps_result_0))
        # write_result.write("\n=======================================================\n")
        # tag_name = tag + "_qso_final_plot"
        # print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_0,1)))+'.pdf', shell=True))
        # if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.        
        #     fix_Re_list_0 = [[i, source_result_0[i]['R_sersic']] for i in range(len(source_result_0))] 
        #     fix_n_list_0 = [[i, source_result_0[i]['n_sersic']] for i in range(len(source_result_0))]
        #     ps_param0 =  _fit_sepc_0.ps_params        
        # #==============================================================================
        # # fitting the QSO as a BHBH (PS+PS)        
        # #==============================================================================
        # print("fitting the QSO as {0} point sources".format(len(arr_x)))
        # num_BHBH = max(len(arr_x), 2)
        # fit_time = 1 #len(glob.glob("fit_result_detect/{0}/fit_image_*_SB_profile_annuli*.pdf".format(file_name_seq[k])))
        # tag = 'fit_result_detect/{0}/fit_{2}-band_image1_PSPS_fittime-{1}'.format(qsoid,fit_time+1, band_seq[k])
        # _fit_sepc_1 = FittingSpeficy(data_process_list[k])
        # del _fit_sepc_1.apertures[0]
        # _fit_sepc_1.prepare_fitting_seq(point_source_num = num_BHBH, fix_n_list= fix_n_list_1, fix_Re_list= fix_Re_list_1, ps_params = ps_param1)
        # _fit_sepc_1.build_fitting_seq()
        # _fit_run_1 = FittingProcess(_fit_sepc_1, savename = tag)
        # _fit_run_1.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
        # _fit_run_1.translate_result()
        # _fit_run_1.plot_final_qso_fit(target_ID = qsoid, save_plot = True, show_plot = show_plot)            
        # source_result_1, ps_result_1 = _fit_run_1.final_result_galaxy, _fit_run_1.final_result_ps 
        # AGN_mags = [ps_result_1[i]['magnitude'] for i in range(len(ps_result_1))]
        # if len(ps_result_1) == 2:
        #     c_miss = np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)
        # elif len(ps_result_1) > 2:
        #     c_miss = [np.sqrt((ps_result_1[0]['ra_image']-ps_result_1[1]['ra_image'])**2+(ps_result_1[0]['dec_image']-ps_result_1[1]['dec_image'])**2)]
        #     c_miss.append(np.sqrt((ps_result_1[1]['ra_image']-ps_result_1[2]['ra_image'])**2+(ps_result_1[1]['dec_image']-ps_result_1[2]['dec_image'])**2))
        #     c_miss.append(np.sqrt((ps_result_1[2]['ra_image']-ps_result_1[0]['ra_image'])**2+(ps_result_1[2]['dec_image']-ps_result_1[0]['dec_image'])**2))
        #     c_miss = np.average(c_miss)
        # reduced_Chisq_1 = _fit_run_1.reduced_Chisq
        # write_result.write("2. Fitting as {0}PS:\n".format(len(ps_result_1)))
        # write_result.write("Reduced Chisq: "+repr(round(reduced_Chisq_1,3)))
        # write_result.write("\nAGN mag: ")
        # for i in range(len(ps_result_1)):
        #     write_result.write(repr(round(AGN_mags[i],3))+' ')
        # write_result.write("\n")
        # for i in range(len(ps_result_1)):
        #     write_result.write("AGN{0} position: ".format(i))
        #     ps_x = -ps_result_1[i]['ra_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        #     ps_y = ps_result_1[i]['dec_image'][0]/data_process_list[k].deltaPix + int(len(data_process_list[k].target_stamp)/2)
        #     ra, dec = wcs.all_pix2world([[ps_x, ps_y]], 1)[0]
        #     write_result.write("x: "+repr(round(ps_x,3))+' y: '+repr(round(ps_y,3))+ " RA: "+repr(ra)+' DEC: '+repr(dec)+ "; ")         
        # write_result.write("\nPS PS center offset (arcsec): "+repr(round(float(c_miss),3)))
        # write_result.write("\nmodel_Sersic_result: "+repr(source_result_1))
        # write_result.write("\nmodel_PS_result: "+repr(ps_result_1))        
        # write_result.write("\n=======================================================\n")
        # tag_name = tag + "_qso_final_plot"  
        # print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_1,1)))+'.pdf', shell=True))
        # if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.     
        #     fix_Re_list_1 = [[i, source_result_1[i]['R_sersic']] for i in range(len(source_result_1))] 
        #     fix_n_list_1 = [[i, source_result_1[i]['n_sersic']] for i in range(len(source_result_1))]
        #     ps_param1 =  _fit_sepc_1.ps_params        
        #==============================================================================
        # fitting the QSO as a BHBH (PS+PS) + Sersic       
        #==============================================================================
        print("fitting the QSO as {0} point sources + Sersic".format(len(arr_x)))
        num_BHBH = max(len(arr_x), 2)
        fit_time = 1
        tag = 'fit_result_detect/{0}/fit_{2}-band_image2_PSPS+Sersic_fittime-{1}'.format(qsoid,fit_time+1, band_seq[k])
        _fit_sepc_2 = FittingSpeficy(data_process_list[k])
        _fit_sepc_2.prepare_fitting_seq(point_source_num = num_BHBH, fix_n_list= fix_n_list_2, fix_Re_list= fix_Re_list_2, ps_params = ps_param2, neighborhood_size = neighborhood_size, threshold = threshold)
        _fit_sepc_2.build_fitting_seq()
        if k == run_list[0]:
            _fit_sepc_2.plot_fitting_sets(savename = 'fit_result_detect/{0}/fitting_used_aper.pdf'.format(qsoid),
                                        show_plot=show_plot)        
        _fit_run_2 = FittingProcess(_fit_sepc_2, savename = tag)
        _fit_run_2.run(algorithm_list = ['PSO'], setting_list= [pso_setting]) 
        _fit_run_2.translate_result()
        _fit_run_2.plot_final_qso_fit(target_ID = qsoid, save_plot = True, show_plot = show_plot)     
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
        pyfits.PrimaryHDU(data_process_list[k].target_stamp-_image_ps_2-objs_img,header=file_header).writeto('fit_result_detect/{0}/data-BHBH(host image)_{1}-band.fits'.format(qsoid, band_seq[k]),overwrite=True)
        print(call("mv {0} {1}".format(tag_name+'.pdf', tag+"_chisq_"+repr(round(reduced_Chisq_2,1)))+'.pdf', shell=True))
        if band_seq[k] == 'I':   #Get the parameters of I band to fix to the other band.
            fix_Re_list_2 = [[i, source_result_2[i]['R_sersic']] for i in range(len(source_result_2))] 
            fix_n_list_2 = [[i, source_result_2[i]['n_sersic']] for i in range(len(source_result_2))]
            ps_param2 = _fit_sepc_2.ps_params
        write_result.close() 
#os.system('say "your program has finished"')
print("Program has finished")
#
