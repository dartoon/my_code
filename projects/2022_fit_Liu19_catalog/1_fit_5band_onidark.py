#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:47:31 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#Select which data you plan to download
dr='dr4'
rerun='s21a_wide'  #Wide 
bands = 'GRIZY'  #Band that will be download
f = open("catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.tools.measure_tools import mask_obj   
import sys

point_source_num = 1  #Number for Point sources (AGNs, or stars) in the target. 0 means no PS.
lband = 'I' #The band fitting first and can also fit n and Re for other band.
fix_n, fix_re = False, False  #Fix sersic n and Re based on I band fitting.
idx = int(sys.argv[1])

#%%
line = lines[idx-1]
ID, Ra, Dec = line.split(' ')
QSO_RA, QSO_DEC = float(Ra), float(Dec)
data_process_list = []
for i, band in enumerate(bands):
    img_globname = glob.glob('gfarm_data_download/{0}_HSC-{1}.fits'.format(ID,band)) + glob.glob('online_data_download/{0}/*cutout*-{1}-*.fits'.format(ID,band) )
    psf_globname = glob.glob('gfarm_data_download/{0}_HSC-{1}_psf.fits'.format(ID,band)) + glob.glob('online_data_download/{0}/*psf*-{1}-*.fits'.format(ID,band) )
    # print(img_globname, psf_globname)
    if len(img_globname) == 1 and len(psf_globname) == 1:
        fitsFile = pyfits.open(img_globname[0])
        file_header0 = fitsFile[0].header
        try:
            FLUXMAG0 = file_header0['FLUXMAG0']
            zp =  2.5 * np.log10(FLUXMAG0)   # This is something Xuheng can't make sure.
        except:
            zp = 27.0
        PSF = pyfits.getdata(psf_globname[0])
        data_process = DataProcess(fov_image = fitsFile[1].data, fov_noise_map = fitsFile[3].data ** 0.5, target_pos = [QSO_RA, QSO_DEC],
                                    pos_type = 'wcs', header = fitsFile[1].header,
                                    rm_bkglight = True, if_plot=False, zp = zp)
        data_process.generate_target_materials(radius=None, radius_list=[30, 35, 40, 50, 60])
        data_process.PSF_list = [PSF]
        data_process_list.append(data_process)
    else:
        data_process_list.append(None)
        
#%% Determining the common settings for all bands, including cutout radius and apertures.
run_list = [i for i in range(len(data_process_list)) if data_process_list[i] is not None]
l_idx = run_list[0]
if data_process_list[2] is not None:
    l_idx = [i for i in range(len(bands)) if bands[i] == lband][0]  #The first index to run
    run_list = [run_list[i] for i in range(len(run_list)) if run_list[i]!=l_idx]
    run_list = [l_idx] + run_list  #The list define the order to run     
cut_radius = np.median([int(len(data_process_list[i].target_stamp)/2) for i in run_list])
for i in run_list:    
    data_process_list[i].generate_target_materials(radius=cut_radius, create_mask = False, nsigma=2.8,
                                          exp_sz= 1.2, npixels = 25, if_plot=False)
    data_process_list[i].clean_aperture()
    data_process_list[i].checkout()
apertures = data_process_list[l_idx].apertures

for i in run_list[1:]:
    covers = mask_obj(data_process_list[i].target_stamp, apertures, if_plot=False, sum_mask = True)
    for j in range(len(data_process_list[i].apertures)):
        new_cover = mask_obj(data_process_list[i].target_stamp, [data_process_list[i].apertures[j]], if_plot=False, sum_mask = True)
        if np.sum(covers - new_cover*covers) > np.sum(1-new_cover)/2 :               
            apertures.append(data_process_list[i].apertures[j])
rm_list = []
for i in range(1,len(apertures)):
    other_apertures = [apertures[j] for j in range(len(apertures)) if i!=j]
    all_cover = mask_obj(data_process_list[l_idx].target_stamp, other_apertures, if_plot=False, sum_mask = True)
    one_cover = mask_obj(data_process_list[l_idx].target_stamp, [apertures[i]], if_plot=False, sum_mask = True)
    if  np.sum(all_cover) - np.sum(all_cover*one_cover) < np.sum(1-one_cover)/1.6:
        rm_list.append(i)
apertures = [apertures[i] for i in range(len(apertures)) if i not in rm_list]
                
fit_sepc_l, fit_run_l = [None]*5, [None]*5
for i in run_list:  
    band = bands[i]
    print("Staring fitting band-"+band+"... ... ...")
    data_process_list[i].apertures = apertures #Pass apertures to the data
    fit_sepc_l[i] = FittingSpecify(data_process_list[i])
    fix_n_list, fix_Re_list = None, None
    if i != l_idx:
        if fix_n == True:
            fix_n_list = [[0,fit_run_l[l_idx].final_result_galaxy[0]['n_sersic'] ]]
        if fix_re == True:
            fix_Re_list = [[0,fit_run_l[l_idx].final_result_galaxy[0]['R_sersic'] ]]
    fit_sepc_l[i].prepare_fitting_seq(point_source_num = point_source_num, supersampling_factor=3, 
                                      fix_n_list= fix_n_list, fix_Re_list=fix_Re_list,
                                      ps_pix_center_list = [[0,0]])
    fit_sepc_l[i].plot_fitting_sets(savename='fit_result/'+ID+'-{0}_set.png'.format(band), show_plot=False)
    fit_sepc_l[i].build_fitting_seq()
    fit_run_l[i] = FittingProcess(fit_sepc_l[i], savename = 'fit_result/'+ID+'-{0}'.format(band), fitting_level=['norm', 'deep'])
    fit_run_l[i].run(algorithm_list = ['PSO', 'PSO'])
    fit_run_l[i].plot_final_qso_fit(save_plot=True, target_ID= 'fit_result/'+ID +'-'+ band, show_plot=False)
    fit_run_l[i].dump_result()
    


        