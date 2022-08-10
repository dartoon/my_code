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
bands = 'I'  #Band that will be download
f = open("catalog.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
import glob
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
from galight.tools.measure_tools import mask_obj   
# import sys

point_source_num = 1  #Number for Point sources (AGNs, or stars) in the target. 0 means no PS.
lband = 'I' #The band fitting first and can also fit n and Re for other band.

#%%
for _, line in enumerate(lines):
    ID, Ra, Dec = line.split(' ')
    # if ID in ['10547', '13120', '13338', '14503', '1956', '2831', '5328', '7007', '8457', '8917', '9009', '9209']:
    if ID == '68':        
        print(_, ID)
        QSO_RA, QSO_DEC = float(Ra), float(Dec)
        data_process_list = []
        for i, band in enumerate(bands):
            img_globname = glob.glob('/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/gfarm_data_download/{0}_HSC-{1}.fits'.format(ID,band)) + glob.glob(
                '/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/online_data_download/{0}/*cutout*-{1}-*.fits'.format(ID,band) )
            psf_globname = glob.glob('/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/gfarm_data_download/{0}_HSC-{1}_psf.fits'.format(ID,band)) + glob.glob(
                '/Volumes/Seagate_Expansion_Drive/data_backup/Liu19_catalog/online_data_download/{0}/*psf*-{1}-*.fits'.format(ID,band) )
            print(img_globname, psf_globname)
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
                data_process.generate_target_materials(radius=35, radius_list=[30, 35, 40, 50, 60], if_plot=True,
                                                       nsigma=1.0, npixels = 4 , contrast=0.0001)
                PSF = PSF[1:-1,:]
                data_process.PSF_list = [PSF]
                
            # aperture = data_process.apertures[1]
            # import copy
            # aperture_1 = copy.deepcopy(aperture)
            # aperture_1.positions = np.array([41,47])
            # aperture_1.a, aperture_1.b = aperture_1.a/3, aperture_1.b/3
            # aperture_2 = copy.deepcopy(aperture)
            # aperture_2.positions = np.array([44,43])
            # aperture_2.a, aperture_2.b = aperture_2.a/3, aperture_2.b/3    
            # del data_process.apertures[1]
            data_process.apertures = aper_bk
            fit_sepc = FittingSpecify(data_process)
            fit_sepc.prepare_fitting_seq(point_source_num = point_source_num, supersampling_factor=3, 
                                              ps_pix_center_list = [[0,0]])
            fit_sepc.plot_fitting_sets('fit_result/'+ID+'-{0}_set.png'.format(band))
            fit_sepc.build_fitting_seq()
            fit_run = FittingProcess(fit_sepc, savename = 'fit_result/'+ID+'-{0}'.format(band), fitting_level=['norm', 'deep'])
            fit_run.run(algorithm_list = ['PSO', 'PSO'])
            fit_run.plot_final_qso_fit(save_plot=False, target_ID= 'fit_result/'+ID +'-'+ band )
            # fit_run_l[i].dump_result()
            
        
        
                