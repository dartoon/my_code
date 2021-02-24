#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:01:48 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

Reines_t1 = np.loadtxt('Reines_2015_table_1.txt', dtype=str)
# ID  Running identification number
# NSAID NASA-Sloan Atlas identification number
# SDSS SDSS identification
# P-MJD-F Plate-MJD-Fiber of SDSS spectrum
# z  Redshift (1)
# mag  iMag Absolute i-band magnitude (2)
# mag  g-i  The (g-i) color (2)
# [solMass] logM* Log stellar mass; corrected for AGN contribution (3)
# [solMass] logMBH Log black hole mass (3)
#%%
# J012159.81-010224.3 20.49921333 -1.040074815
# J084143.50+013149.8 130.4312678 1.530471828
# J095540.47+050236.5 148.9187161 5.043512791
# J100207.03+030327.6 150.5293397 3.057698011
# J104252.93+041441.1 160.7205967 4.244771374
# J110032.29+020656.0 165.1355266 2.116022222
# J120257.81+045045.0 180.7409093 4.845859658
# J121826.72-000750.1 184.6113252 -0.1305986132 ?
# J131305.81+012755.9 198.2742239 1.465547128
# J134426.41+441620.0 206.1100925 44.27223612
# J134952.60+020440.0 207.4701858 2.079184649
# J141630.82+013708.0 214.1283952 1.618930223
# J141920.64+043623.3 214.8360752 4.606470097
# J143450.63+033842.5 218.7109935 3.645158632
# J145819.55+045451.7 224.5815106 4.914388788 only-1-band
# J160558.11+440319.5 241.4921964 44.05540825
# J224824.63+000920.6 342.1026171 0.1557374447
# J231815.67+001540.2 349.5653755 0.2611389517
# J233837.09-002810.4 354.6545807 -0.4695773356

import glob
info = np.loadtxt('./Reines/Reines_ID_RaDec.txt', dtype=str)
ID_ = info[:,0]
RA_ = info[:,1].astype(np.float)
Dec_ = info[:,2].astype(np.float)

files = glob.glob('./Reines/*-I.fits')
IDs = [files[i].split('/')[-1].split('_HSC')[0] for i in range(len(files))]
IDs.sort()

mag_mis = []

for i in range(len(IDs)):
    ID = IDs[i]
    RA = RA_[(np.where(ID == ID_)[0][0])]
    Dec =  Dec_[(np.where(ID == ID_)[0][0])]
    #%%
    fitsFile = pyfits.open('./Reines/{0}_HSC-I.fits'.format(ID))
    fov_image= fitsFile[1].data
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    err_data= fitsFile[3].data ** 0.5
    
    file_header0 = fitsFile[0].header
    zp =  27.0
    PSF = pyfits.getdata('./Reines/{0}_HSC-I_psf.fits'.format(ID))
    if len(PSF) != 0 and PSF.shape[0] != PSF.shape[1]:
        cut = int((PSF.shape[0] - PSF.shape[1])/2)
        if cut>0:
            PSF = PSF[cut:-cut,:]
        elif cut<0:
            PSF = PSF[:,-cut:cut]
        PSF /= PSF.sum()
    
    #%%Start to use decomprofile
    from decomprofile.data_process import DataProcess
    QSO_RA = RA
    QSO_DEC = Dec
    data_process = DataProcess(fov_image = fov_image, fov_noise_map = err_data, target_pos = [QSO_RA, QSO_DEC],
                                pos_type = 'wcs', header = header,
                              rm_bkglight = False, if_plot=True, zp = zp)
    
    data_process.noise_map = err_data
    
    data_process.generate_target_materials(radius=None, create_mask = False, nsigma=3,
                                           radius_list = [120, 140, 160,180],
                                          exp_sz= 1.2, npixels = 55, if_plot=True)
    
    data_process.PSF_list = [PSF]
    
    # data_process.checkout() #Check if all the materials is known.
    
    # #Start to produce the class and params for lens fitting.
    # from decomprofile.fitting_specify import FittingSpeficy
    # fit_sepc = FittingSpeficy(data_process)
    # fit_sepc.apertures = [fit_sepc.apertures[0]]
    
    # fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_center_list = [[0,0]])#, fix_n_list= [[0,4]], fix_center_list = [[0,0]])
    # fit_sepc.plot_fitting_sets()
    # fit_sepc.build_fitting_seq()
    
    # #%%Setting the fitting method and run.
    # from decomprofile.fitting_process import FittingProcess
    # fit_run = FittingProcess(fit_sepc, savename = '{0}_I'.format(ID))
    # fit_run.run(['PSO'], [None])
    # fit_run.plot_all()
    # fit_run.dump_result()
    # fit_run.translate_result()
    # # # print(fit_run.final_result_galaxy[0])
    
    #%%
    # # Test load pkl
    # import pickle
    # picklename = '{0}_I'.format(ID) + '.pkl'
    # fit_run = pickle.load(open(picklename,'rb'))
    # fit_run.plot_all()
    # # fitting_run_class.run_diag()
    
    #%% Calculate abs magnitude
    idx = np.where(Reines_t1[:,2] == ID)[0][0]
    Reines_iMag = float(Reines_t1[idx, 5])
    z = float(Reines_t1[idx, 4])
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dl = cosmo.luminosity_distance(z).value  # Mpc
    # obs_mag = fit_run.final_result_galaxy[0]['magnitude']
    
    obs_mag =  -2.5*np.log10(np.sum(data_process.target_stamp)) + 27
    abs_mag = obs_mag - 5*(np.log10(dl * 10**6)-1)   #dl is the luminosity distance which is a function of redshift:
    print("mag miss:", abs_mag - Reines_iMag)
    mag_mis.append(abs_mag - Reines_iMag)