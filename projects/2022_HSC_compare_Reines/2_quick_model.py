#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:03:33 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.tools.plot_tools import plt_fits
Reines_t1 = np.loadtxt('2021_previous/Reines_2015_table_1.txt', dtype=str)

f = open("2021_previous/Reines_ID_RaDec.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n
# IDs = glob.glob('s21a/*')
# IDs.sort()
IDs = ['J131310.12+051942.1','J233837.09-002810.4',
       'J145819.55+045451.7', 'J143450.63+033842.5', 'J141920.64+043623.3', 'J141630.82+013708.0',
       'J140018.41+050242.2', 'J134426.41+441620.0','J131305.81+012755.9','J121826.72-000750.1',
       'J104252.93+041441.1','J012159.81-010224.3','J084143.50+013149.8','J095540.47+050236.5']
IDs.sort()
# names = glob.glob('*.pkl')
# IDs =  [name.split('.pkl')[0] for name in names]

bands = ['G', 'R', 'I', 'Z', 'Y']
from astropy.cosmology import FlatLambdaCDM
# for ID in IDs[:5]:
for ID in IDs[5:10]:
# for ID in IDs[10:14]:    
    for band in bands:
        idx = np.where(Reines_t1[:,2] == ID)[0][0]
        Reines_iMag = float(Reines_t1[idx, 5])
        z = float(Reines_t1[idx, 4])
        ra = ID[1:10]
        ra = ra[:2] + ':' + ra[2:4] + ':' +  ra[4:]
        dec = ID[10:]
        dec = dec[:3] + ':' + dec[3:5] + ':' +  dec[5:]
        
        fitsname = glob.glob('./s21a/{0}/*cutout*{1}*.fits'.format(ID, band))
        if fitsname != []:
            fitsFile = pyfits.open(fitsname[0])
            fov_image= fitsFile[1].data
            header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            err_data= fitsFile[3].data ** 0.5
            
            file_header0 = fitsFile[0].header
            FLUXMAG0 = file_header0['FLUXMAG0']
            zp =  2.5 * np.log10(FLUXMAG0)
            psfname = glob.glob('./s21a/{0}/*psf*{1}*.fits'.format(ID, band))
            PSF = pyfits.getdata(psfname[0])
            radius = None
            if glob.glob('galight_results/'+ID+'_{0}*pkl'.format(band)) != []:
                fit_run = pickle.load(open(glob.glob('galight_results/'+ID+'_{0}*pkl'.format(band))[0],'rb')) 
            else:
                print(ID, band)
                data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [ra, dec],
                                           pos_type = 'wcs', header = header,
                                           rm_bkglight = True, if_plot=False, zp = zp)
                
                data_process.generate_target_materials(radius=radius, create_mask = False, nsigma=2.8,
                                                       radius_list = [80, 100, 120, 140, 160, 180], 
                                                      exp_sz= 1.2, npixels = 15, if_plot=False)
                radius = int(len(data_process.target_stamp)/2)
                data_process.PSF_list = [PSF]
                data_process.checkout() #Check if all the materials is known.
            
                data_process.apertures = [data_process.apertures[0]]
            
                fit_sepc = FittingSpecify(data_process)
                fit_sepc.prepare_fitting_seq(point_source_num = 1, supersampling_factor=3)#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
                #Plot the initial settings for fittings. 
                fit_sepc.plot_fitting_sets()
                
                fit_sepc.build_fitting_seq()
                
                from galight.fitting_process import FittingProcess
                fit_run = FittingProcess(fit_sepc, savename = 'galight_results/'+ID+'_'+band, fitting_level=['shallow','deep'])
                fit_run.run(algorithm_list = ['PSO', 'PSO']) 
                
                fit_run.dump_result(savedata=True)
            print(ID, band)
            fit_run.plot_final_qso_fit(target_ID = ID)
            print(fit_run.final_result_galaxy[0]['magnitude'])
            obs_mag = fit_run.final_result_galaxy[0]['magnitude']
            # obs_mag = -2.5*np.log10(np.sum(fit_run.flux_2d_out['data'])) + zp
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            dl = cosmo.luminosity_distance(z).value  # Mpc
            abs_mag = obs_mag - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
            Reines_iMag = float(Reines_t1[idx, 5])
            print('abs_mag, Reines_iMag: ', abs_mag, Reines_iMag)
        
