#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:52:07 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from decomprofile.fitting_process import FittingProcess
import pickle, glob
import copy, matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
from matplotlib.colors import LogNorm
from decomprofile.tools.plot_tools import profile_plots
import lenstronomy.Util.param_util as param_util

fitsFile_ = glob.glob('SDSS_0.2-0.3/*_HSC-I.fits')
# NO = 
# for NO in range(1):
# for NO in range(21,42):
# for NO in range(42,63):
# for NO in range(55):    
for NO in range(len(fitsFile_)):  
# for NO in [66]:
    fitsFile  = fitsFile_[NO]
    ID = fitsFile.split('/')[1].split('_')[0]
    PSF_filename = fitsFile.split('.fits')[0]+ '_psf.fits'
    save_name = 'fit_result/' + fitsFile.split('.fits')[0].split('/')[1]
    ## Test load pkl
    
    picklename = save_name+'single_Sersic.pkl'
    fitting_run_class_0 = pickle.load(open(picklename,'rb'))

    
    picklename = save_name+'bulge+disk_2nd.pkl'
    fitting_run_class_1 = pickle.load(open(picklename,'rb'))

    
    bic_0 = fitting_run_class_0.fitting_seq.bic
    bic_1 = fitting_run_class_1.fitting_seq.bic
    if bic_0 < bic_1:
        print(ID+" is a singel Sersic model; "+"glob Number: " + str(NO))
    else:
        print(ID+" is a disk+bulge!!! "+"glob Number: " + str(NO))

        print(ID+"Single Sersic fit:")
        fitting_run_class_0.plot_final_qso_fit(target_ID =  ID)
        print(ID+"Disk + Bulge fit:")    
        fitting_run_class_1.plot_final_qso_fit(target_ID =  ID)
        
        print("BIC compare:", round(bic_0,1), round(bic_1,1))
        print("Chisq:", round(fitting_run_class_0.reduced_Chisq,1), round(fitting_run_class_1.reduced_Chisq,1))

        bulge = fitting_run_class_1.image_host_list[0]
        disk = fitting_run_class_1.image_host_list[1]
        B2T = np.sum(bulge)/np.sum(bulge+disk)
        
        AGN = fitting_run_class_1.image_ps_list[0]
        
        bulge_Re = fitting_run_class_1.final_result_galaxy[0]['R_sersic']
        disk_Re = fitting_run_class_1.final_result_galaxy[1]['R_sersic']
        
        flux_list_2d = [bulge, disk, AGN]
        label_list_2d = ['Bulge', 'Disk', 'nuclei']
        flux_list_1d = [bulge, disk, AGN] 
        label_list_1d = ['Bulge', 'Disk', 'nuclei']
        
        profile_plots(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d,
                      deltaPix = fitting_run_class_1.fitting_specify_class.deltaPix,
                      target_ID =  ID, if_annuli=True)
        
        print("B/T:", round(B2T,2))
        print('Reff: bulge, disk: ', round(bulge_Re,2), round(disk_Re,2) )
        hold = input("hold:")


# =============================================================================
# Sec with better prior
# =============================================================================
# 092256.42+042733.3 is a disk+bulge!!! glob Number: 13  # Good dis bulge 1 D shape.
# 140325.35+011338.6 is a disk+bulge!!! glob Number: 14
# 021122.59-021129.7 is a disk+bulge!!! glob Number: 20
# 101326.05-000136.4 is a disk+bulge!!! glob Number: 21
# 143348.21+012937.6 is a disk+bulge!!! glob Number: 25 # Clear host feature
# 151109.93+423434.5 is a disk+bulge!!! glob Number: 42 # Merger feature
# 134510.33+001852.2 is a disk+bulge!!! glob Number: 43 # Merger feature

# =============================================================================
# First 
# =============================================================================
# [21,25, 29,40,55]    
# 144747.02+005518.0 is a disk+bulge!!! glob Number: 11 # intersting case, can be fitted more carefully
# 144257.79+444301.5 is a singel Sersic model; glob Number: 18 # larger frame is needed 
# 101326.05-000136.4 is a disk+bulge!!! glob Number: 21 # Bulge can be seen? #!!!
# 143348.21+012937.6 is a singel Sersic model; glob Number: 25 # larger frame is needed, bulge featuer is clear #!!!
#!!! 022021.75-014809.3 is a disk+bulge!!! glob Number: 29 #Bulge can be seen from the redisual map! Inetersting!!!
# 120312.13+015321.3 is a singel Sersic model; glob Number: 36 #Bar image is clear
# 012344.00-003554.7 is a disk+bulge!!! glob Number: 38 # Central image residual imporved.
# 020559.66-063736.6 is a disk+bulge!!! glob Number: 40 #Bulge Clearly #!!!
# 235441.54-000448.6 is a disk+bulge!!! glob Number: 41 #Residual imporved significant 
# 235441.54-000448.6 is a singel Sersic model; glob Number: 41 #Why BIC not make sense? so much?
# 151109.93+423434.5 is a singel Sersic model; glob Number: 42 #Merging system? Can be fitting better.
# 134510.33+001852.2 is a disk+bulge!!! glob Number: 43 #Merging system?
# 023344.40-015042.8 is a disk+bulge!!! glob Number: 48 #Merging system?
# 134113.83-005626.3 is a singel Sersic model; glob Number: 51 # Interesting system. Jet?
# 090054.40+012605.3 is a singel Sersic model; glob Number: 55 #larger frame is needed, bulge featuer is clear #!!!
# 134826.86-005943.9 is a singel Sersic model; glob Number: 56 #Host is  very extended
# 000557.23+002837.7 is a singel Sersic model; glob Number: 66 #Merging system!!!
# 085431.28-003650.6 is a disk+bulge!!! glob Number: 70 #Host clear 
# 145857.84-013419.3 is a disk+bulge!!! glob Number: 73 #Host clear 
# 010335.35-005527.3 is a disk+bulge!!! glob Number: 76 #Host clear 
# 144359.75+012956.2 is a disk+bulge!!! glob Number: 82 #Clear image, beautiful!
# 101804.46+002059.9 is a disk+bulge!!! glob Number: 83 #Clear image and bulge
