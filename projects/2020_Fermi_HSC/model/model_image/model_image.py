#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 15:17:23 2021

@author: Xuheng Ding
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

import glob
from decomprofile.data_process import DataProcess
from decomprofile.fitting_specify import FittingSpeficy
from decomprofile.fitting_process import FittingProcess

Bla_type = 'BL_Lac'
# Bla_type = 'FSRQ'
for band in ['G', 'R', 'I' ,'Z' , 'Y']:
    ids_file = glob.glob('../../data/{1}/*-{0}.fits'.format(band, Bla_type))
    ids = [ids_file[i].split('/')[-1].split('_')[0] for i in range(len(ids_file))]
    for i in range(len(ids_file)):    
        fitsFile = pyfits.open(ids_file[i])
        fov_image= fitsFile[1].data
        header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        err_data= fitsFile[3].data ** 0.5
        PSF = pyfits.getdata(ids_file[1].split('.fits')[0] + '_psf.fits')
        zp =  27.0
        
        if 'BL_Lac' in ids_file[0]:
            pos_file = '../../data/Fermi_bll.txt'
        elif 'FSRQ' in ids_file[0]:
            pos_file = '../../data/Fermi_fsrq.txt'
        f =  open(pos_file,"r")
        pos_info = f.read() 
        pos_info = pos_info.split('\n')   # Split in to \n
        line = [j for j in range(len(pos_info)) if ids[i] in pos_info[j]][0]
        test_id, QSO_RA, QSO_DEC = pos_info[line].split(' ')
        QSO_RA = float(QSO_RA)
        QSO_DEC = float(QSO_DEC)
        if test_id!=ids[i]:
            raise ValueError("When load pos for the target, line for ID not right")
        data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                                   pos_type = 'wcs', header = header,
                                   rm_bkglight = True, if_plot=False, zp = zp)
        
        
        data_process.noise_map = err_data
        try:
            data_process.generate_target_materials(radius=None, create_mask = False, nsigma=1.5,
                                                  exp_sz= 1.2, npixels = 15, if_plot=False)
        except:
            continue
        data_process.PSF_list = [PSF]
        data_process.checkout() #Check if all the materials is known.
        
        #Start to produce the class and params for lens fitting.
        fit_sepc = FittingSpeficy(data_process)
        fit_sepc.prepare_fitting_seq(point_source_num = 1)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
        fit_sepc.plot_fitting_sets()
        fit_sepc.build_fitting_seq()
        
        #Setting the fitting method and run.
        fit_run = FittingProcess(fit_sepc, savename = 'model_result/' + ids_file[i].split('/')[-1][:-5])
        fit_run.run(algorithm_list = ['PSO'], setting_list = [None])
        fit_run.plot_final_qso_fit(save_plot=True, target_ID = ids_file[i].split('/')[-1][:-5])
        fit_run.dump_result()
        # print(fit_run.final_result_galaxy[0])
        filename = 'model_result/{0}_host_result.txt'.format(Bla_type)
        if_file = glob.glob(filename)
        if if_file == []:
            write_file =  open(filename,'w') 
        else:
            write_file =  open(filename,'r+')
            write_file.read()
        write_file.write(ids_file[i].split('/')[-1][:-5]+' '+str(fit_run.final_result_galaxy[0])+'\n')
        write_file.close()        