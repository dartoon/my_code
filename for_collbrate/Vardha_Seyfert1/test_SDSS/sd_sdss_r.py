#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:34:48 2019

@author: Dartoon

Fit the Seyfert as Point source + two Sersic (Bluge + Disk).
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../fitting_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result, transfer_obj_to_result
import glob
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from flux_profile import total_compare
import lenstronomy.Util.param_util as param_util
from fit_qso import condition_bulgedisk

import numpy as np
import corner
from lenstronomy.Util import kernel_util
from lenstronomy.Util import image_util

from lenstronomy.Plots import model_plot
import math
pi = math.pi

namelist = ['103'] 
#namelist = ['138'] 
# namelist = ['2','9','15','19','20','22','24','26','31','35','39','44','45','51','52','53','58','61','62','63','70','71','73','78','91','99','100','102','108','109','143','155','156','180','187','202']
#disks (36)
for i in range(len(namelist)):
    ID = namelist[i]
        
    rf = open("./sdss_input.dat",'r')
    for line in rf:
        ob = line.split()[0]
        if ob == ID:
            bulgereff = float(line.split()[1])
            bulgen = float(line.split()[2])
            bulgee1 = float(line.split()[3])
            bulgee2 = float(line.split()[4])
            diskreff = float(line.split()[5])
            diske1 = float(line.split()[6])
            diske2 = float(line.split()[7])
            barreff = float(line.split()[8])
            bare1 = float(line.split()[9])
            bare2 = float(line.split()[10])
            phi0, q0 = param_util.ellipticity2phi_q(bulgee1, bulgee2)
            phi1, q1 = param_util.ellipticity2phi_q(diske1, diske2)
            phi2, q2 = param_util.ellipticity2phi_q(bare1, bare2)
            phi0 = -1*np.degrees(phi0)+90
            phi1 = -1*np.degrees(phi1)+90
            phi2 = -1*np.degrees(phi2)+90
            
    rf = open("./orientation_hst.dat",'r')
    for line in rf:
        ob = line.split()[0]
        if ob == ID:
            pahst = float(line.split()[1])
            if pahst < 0:
                pahst = pahst+360
            
    rf = open("./orientation_sdss.dat",'r')
    for line in rf:
        ob = line.split()[0]
        if ob == ID:
            pasdss = float(line.split()[1])
            if pasdss < 0:
                pasdss = pasdss+360
    
    #orientdiff = pasdss-pahst
    orientdiff = pahst-pasdss
    if orientdiff < 0:
        orientdiff = orientdiff+180
    print (pahst, pasdss, orientdiff)
    phi0new = (phi0-orientdiff)
    if phi0new > 360:
        phi0new = phi0new-360
    phi1new = (phi1-orientdiff)
    if phi1new > 360:
        phi1new = phi1new-360
    phi2new = (phi2-orientdiff)
    if phi2new > 360:
        phi1new = phi2new-360
    phi0final = np.radians(-1*(phi0new-90))
    phi1final = np.radians(-1*(phi1new-90))
    phi2final = np.radians(-1*(phi2new-90))
    print (phi1, phi1new, phi1final, orientdiff)
    e1bulge, e2bulge = param_util.phi_q2_ellipticity(phi0final, q0)
    e1disk, e2disk = param_util.phi_q2_ellipticity(phi1final, q1)
    e1bar, e2bar = param_util.phi_q2_ellipticity(phi2final, q2)
    
    obj = '{0}'.format(ID)
    ID = 'L{0}'.format(ID)             # Object ID
    deep_seed = False                  # False = trial run
    no_MCMC = True                     # True = without MCMC 
    mpi_para = False
    pltshow = 0  #Change to 1 to show the plot while fitting.
    pix_sz = 0.396  # The pixel scale of the SDSS images
    exp_time = 53.91 #Exposure time in seconds of the SDSS images
    lowerpos = -3*pix_sz #lower limit offset in x
    upperpos = 3*pix_sz #lower limit offset in y
    spreadx = 2*pix_sz #spread in x
    spready = 2*pix_sz #spread in y
    #kernel_size = 125                   # psf kernel size (odd number)

    #mask_image = pyfits.getdata('./fits/{0}_mask.fits'.format(ID))
    #QSO_msk = 1-mask_image


    rf = open("./start.dat",'r')
    for line in rf:
        ob = line.split()[0]
        if ob==obj:
            centerxg = float(line.split()[4])
            centeryg = float(line.split()[5])
            centerxr = float(line.split()[6])
            centeryr = float(line.split()[7])
            centerxi = float(line.split()[8])
            centeryi = float(line.split()[9])
            centerxz = float(line.split()[10])
            centeryz = float(line.split()[11])

    name_save = '{0}_sersicdisk_r'.format(ID)
    band = 'r'
    agn_image = pyfits.getdata('./{0}{1}_final.fits'.format(ID,band))
    frsize = agn_image.shape[1]
    QSO_msk = np.ones_like(agn_image)

    center_x = round((frsize/2-centerxr)*pix_sz,2) #offset in arcsec; object on the left side should have center_x as positive, since the east are on the left direction.
    center_y = round(-(frsize/2-centeryr)*pix_sz,2) #offset in arcsec
    #print (centerxr, frsize, center_x, centeryr, frsize, center_y)

    
    agn_stdd = pyfits.getdata('./{0}{1}_noi.fits'.format(ID,band))
    psf = pyfits.getdata('./{0}_psf2_{1}.fits'.format(ID,band))
    #kernel_util.cut_psf(psf, psf_size=kernel_size)
        
    fixed_source = []
    kwargs_source_init = []
    kwargs_source_sigma = []
    kwargs_lower_source = []
    kwargs_upper_source = []      
    fixed_source.append({'R_sersic': bulgereff, 'n_sersic': bulgen, 'e1': e1bulge, 'e2': e2bulge})  
    kwargs_source_init.append({'R_sersic': bulgereff, 'n_sersic': bulgen, 'e1': e1bulge, 'e2': e2bulge, 'center_x': center_x, 'center_y': center_y})
    kwargs_source_sigma.append({'center_x': spreadx, 'center_y': spready})
    kwargs_lower_source.append({'center_x': lowerpos, 'center_y': lowerpos})
    kwargs_upper_source.append({'center_x': upperpos, 'center_y': upperpos})
    fixed_source.append({'R_sersic': diskreff, 'n_sersic': 1, 'e1': e1disk, 'e2': e2disk})  
    kwargs_source_init.append({'R_sersic': diskreff, 'n_sersic': 1., 'e1': e1disk, 'e2': e2disk, 'center_x': center_x, 'center_y': center_y})
    kwargs_source_sigma.append({'center_x': spreadx, 'center_y': spready})
    kwargs_lower_source.append({'center_x': lowerpos, 'center_y': lowerpos})
    kwargs_upper_source.append({'center_x': upperpos, 'center_y': upperpos})

    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    

    fixed_ps = []
    kwargs_ps_init = []
    kwargs_ps_sigma = []
    kwargs_lower_ps = []
    kwargs_upper_ps = []
    fixed_ps.append({})
    kwargs_ps_init.append({'ra_image': [center_x], 'dec_image': [center_y]})
    kwargs_ps_sigma.append({'ra_image': [spreadx], 'dec_image': [spready]})
    kwargs_lower_ps.append({'ra_image': [lowerpos], 'dec_image': [lowerpos]})
    kwargs_upper_ps.append({'ra_image': [upperpos], 'dec_image': [upperpos]})
    ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]         

    source_result, ps_result, image_ps, image_host, error_map=fit_qso(agn_image, psf_ave=psf, ps_param = ps_param, 
                                                                  pix_sz = pix_sz, exp_time = exp_time, QSO_std = agn_stdd,  QSO_msk = QSO_msk,
                                                                  source_params=source_params, fixcenter=False, no_MCMC = no_MCMC,
                                                                  tag=name_save, deep_seed= deep_seed, pltshow = pltshow, flux_ratio_plot=True,
                                                                      dump_result = True, condition=condition_bulgedisk)

  
