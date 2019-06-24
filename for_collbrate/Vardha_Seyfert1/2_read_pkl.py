#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 21:41:53 2019

@author: Dartoon

Read the pickle result
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy

ID = 'l106'
picklename = ID+'_fitting.pkl'
result = pickle.load(open(picklename,'rb'))
[best_fit,pso_fit,mcmc_fit, trans_paras] = result

source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, param_list, _ = pso_fit
samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
_, _, mcmc_new_list, labels_new, _ = trans_paras
mcmc_new_list = np.asarray(mcmc_new_list)

pix_sz = 0.04
zp = 25.0985
#%% Readout the fitted parameters for each component:
Reff_idxs = [i for i in range(len(param_mcmc)) if param_mcmc[i] == 'R_sersic_source_light']
#Read out the Point source magnitude:
ps_id = 0
flux_l, flux_m, flux_h = [np.percentile(mcmc_new_list[:,ps_id],percent,axis=0) for percent in [16, 50, 84]]
mag_l, mag_m , mag_h = [ (-2.5*np.log10(flux) + zp) for flux in [flux_h, flux_m, flux_l]]
print "PointSource mag :", round(mag_m,3), '-',round(mag_m-mag_l,3), '+',round(mag_h-mag_m,3)

#Read the fitting result for bulge, which is the first fitting component:
Bulge_id = 0
Reff_idx = Reff_idxs[Bulge_id] #The index for Reff
n_idx = Bulge_id+1 #The index for Sersic n 
flux_idx = Bulge_id+1 #The first is QSO
Reff_l, Reff_m, Reff_h=[np.percentile(samples_mcmc[:,Reff_idx],percent,axis=0) for percent in [16, 50, 84]]
print "Bulge", param_mcmc[Reff_idx], ":",round(Reff_m,3), '-',round(Reff_m-Reff_l,3), '+',round(Reff_h-Reff_m,3)
n_l, n_m, n_h= [np.percentile(samples_mcmc[:,n_idx],percent,axis=0) for percent in [16, 50, 84]]
print "Bulge", param_mcmc[n_idx], "(arcsec):", round(n_m,3), '-',round(n_m-n_l,3), '+',round(n_h-n_m,3)
flux_l, flux_m, flux_h = [np.percentile(mcmc_new_list[:,flux_idx],percent,axis=0) for percent in [16, 50, 84]]
mag_l, mag_m , mag_h = [ (-2.5*np.log10(flux) + zp) for flux in [flux_h, flux_m, flux_l]]
print "Bulge mag :", round(mag_m,3), '-',round(mag_m-mag_l,3), '+',round(mag_h-mag_m,3)

#Read the fitting result for disk, which is the second fitting component:
Disk_id = 1
Reff_idx = Reff_idxs[Disk_id] #The index for Reff
if source_result[Disk_id]['n_sersic'] !=1:
    print "Warning: the Disk's Sersic index is not n=1 !!!"
flux_idx = Disk_id+1 #The first is QSO
Reff_l, Reff_m, Reff_h=[np.percentile(samples_mcmc[:,Reff_idx],percent,axis=0) for percent in [16, 50, 84]]
print "Disk", param_mcmc[Reff_idx], "(arcsec):",round(Reff_m,3), '-',round(Reff_m-Reff_l,3), '+',round(Reff_h-Reff_m,3)
flux_l, flux_m, flux_h = [np.percentile(mcmc_new_list[:,flux_idx],percent,axis=0) for percent in [16, 50, 84]]
mag_l, mag_m , mag_h = [ (-2.5*np.log10(flux) + zp) for flux in [flux_h, flux_m, flux_l]]
print "Disk mag:", round(mag_m,3), '-',round(mag_m-mag_l,3), '+',round(mag_h-mag_m,3)

#Read the fitting result for disk, which is the second fitting component:
Bar_id = 2
Reff_idx = Reff_idxs[Bar_id] #The index for Reff
if source_result[Bar_id]['n_sersic'] !=0.5:
    print "Warning: the Bar's Sersic index is not n=0.5 !!!"
flux_idx = Bar_id+1 #The first is QSO
Reff_l, Reff_m, Reff_h=[np.percentile(samples_mcmc[:,Reff_idx],percent,axis=0) for percent in [16, 50, 84]]
print "Bar", param_mcmc[Reff_idx], "(arcsec):",round(Reff_m,3), '-',round(Reff_m-Reff_l,3), '+',round(Reff_h-Reff_m,3)
flux_l, flux_m, flux_h = [np.percentile(mcmc_new_list[:,flux_idx],percent,axis=0) for percent in [16, 50, 84]]
mag_l, mag_m , mag_h = [ (-2.5*np.log10(flux) + zp) for flux in [flux_h, flux_m, flux_l]]
print "Bar mag:", round(mag_m,3), '-',round(mag_m-mag_l,3), '+',round(mag_h-mag_m,3)

#%% To save the result fits file:
#[0] original image
#[1] total model
#[2] PSF model
#[3] bulge model
#[4] disk model (if present)
#[5] bar model (if present)
#[6] residual (original - total model)
center_QSO = np.array([2053, 2475])   #!!! The is the position array for the QSO that introduced in 0_cutout.py
fitsFile = pyfits.open(ID+"_sci.fits")
file_header = copy.deepcopy(fitsFile[0].header)

flux_hdulist = pyfits.open(ID+'_flux_list.fits')  #includes [agn_image,image_ps, extended_source, error_map, QSO_msk]
flux_list = [flux_hdulist[i].data for i in range(5)]
file_header['CRPIX1'] = file_header['CRPIX1']-center_QSO[0]+len(flux_list[0])/2
file_header['CRPIX2'] = file_header['CRPIX2']-center_QSO[1]+len(flux_list[0])/2

data_image_hdu = pyfits.PrimaryHDU(flux_list[0])
data_image_hdu.header['targname'] = "original image"
data_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
data_image_hdu.header['CRPIX2'] = file_header['CRPIX2']
model_image_hdu = pyfits.ImageHDU((flux_list[1]+flux_list[2]))
model_image_hdu.header['targname'] = 'total model'
model_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
model_image_hdu.header['CRPIX2'] = file_header['CRPIX2']
residual_image_hdu = pyfits.ImageHDU((flux_list[0]-(flux_list[1]+flux_list[2])))
residual_image_hdu.header['targname'] = 'residual (original - total model)'
residual_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
residual_image_hdu.header['CRPIX2'] = file_header['CRPIX2']
PS_hdu = pyfits.ImageHDU(image_ps)
PS_hdu.header['targname'] = 'Fitted Point Source image'
PS_hdu.header['CRPIX1'] = file_header['CRPIX1']
PS_hdu.header['CRPIX2'] = file_header['CRPIX2']
bulge_image_hdu = pyfits.ImageHDU(image_host[0])
bulge_image_hdu.header['targname'] = 'Fitted Bulge image'
bulge_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
bulge_image_hdu.header['CRPIX2'] = file_header['CRPIX2']
if len(source_result)>1:
    if source_result[1]['n_sersic'] ==1: #Disk Sersic n must == 1
        disk_image_hdu = pyfits.ImageHDU(image_host[1])
        disk_image_hdu.header['targname'] = 'Fitted Disk image'
        disk_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
        disk_image_hdu.header['CRPIX2'] = file_header['CRPIX2']        
else:
    disk_image_hdu = pyfits.ImageHDU(image_host[0] * 0)
if len(source_result)>2:
    if source_result[2]['n_sersic'] ==0.5: #Disk Sersic n must == 1
        bar_image_hdu = pyfits.ImageHDU(image_host[2])
        bar_image_hdu.header['targname'] = 'Fitted Disk image'
        bar_image_hdu.header['CRPIX1'] = file_header['CRPIX1']
        bar_image_hdu.header['CRPIX2'] = file_header['CRPIX2']
else:
    bar_image_hdu = pyfits.ImageHDU(image_host[0] * 0)
thdu_fluxlist = pyfits.HDUList([data_image_hdu, model_image_hdu,PS_hdu,
                                bulge_image_hdu, disk_image_hdu, bar_image_hdu,residual_image_hdu])
thdu_fluxlist.writeto(ID+'_result.fits', overwrite=True)
#%% Plot the updated SB profile:
import sys
sys.path.insert(0,'./fitting_tools/')
from flux_profile import total_compare
flux_hdulist = pyfits.open(ID+'_flux_list.fits')  # includes [agn_image,image_ps, extended_source, error_map, QSO_msk, psf]
flux_list = [flux_hdulist[i].data for i in range(4)]
QSO_msk = flux_hdulist[4].data
label = ['data', 'QSO', 'extended sources', 'model', 'normalized residual']
fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz, zp = zp,
                    plot_compare = False, msk_image = QSO_msk, host_comp_name=['bulge', 'disk', 'bar'], host_comp=[image_host[0], image_host[1], image_host[2]])
plt.show()
