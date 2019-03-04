#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:11:32 2019

@author: Dartoon

1. Remove l106 background light.
2. Cutout l106 image.
3. Select potential PSFs and compare the CID1174 ones.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from cut_image import cut_center_bright, cut_image, save_loc_png
from mask_objects import mask_obj
from est_bkg import sub_bkg
from flux_profile import flux_profile, profiles_compare
ID = 'l106'

#==============================================================================
# Cut out the lens image
#==============================================================================
fitsFile = pyfits.open('l106_sci.fits')
img = fitsFile[0].data # check the back grounp
#img[img==0.] = np.nan
#img_sub, bkg_light = sub_bkg(img)
#pyfits.PrimaryHDU(img_sub).writeto('l106_sci_sublight.fits',overwrite=True)
img_sub = pyfits.getdata('l106_sci_sublight.fits')

#plt.figure(figsize=(10, 10))
#plt.imshow(img, norm=LogNorm(), origin='lower', vmax = 1)
#plt.colorbar()
#plt.show()
#plt.figure(figsize=(10, 10))
#plt.imshow(img_sub, norm=LogNorm(), origin='lower', vmax = 1)
#plt.colorbar()
#plt.show()

img_sub = np.nan_to_num(img_sub)
center_QSO = np.array([2053, 2475]) #np.where(img_sub==img_sub.max()) #c_psf_list[QSO_loc]

fr = 400  # The radius of the frame

agn_image = cut_image(img_sub, center_QSO, fr)
plt.imshow(agn_image, norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()

#lens_bkg_light = cut_image(bkg_light, center_QSO, 100)
#plt.imshow(lens_bkg_light, norm=LogNorm(),origin='low')
#plt.colorbar()
#plt.show()
##==============================================================================
## #Estimate the exposure time and error map:
##==============================================================================
wht = pyfits.getdata('l106_wht.fits')
exp = 800
mean_wht = exp * (0.04/0.04)**2
exp_map = exp * wht/mean_wht
from photutils import make_source_mask
mask = make_source_mask(agn_image, snr=1, npixels=5, dilate_size=11)
plt.imshow(agn_image * (1-mask*1), norm=LogNorm(), origin='low')
plt.show()
stdd = np.std(agn_image[mask ==False])  # TOO SMALL
print "stdd:",stdd
#stdd = 0.0064   # Read out from the empty place

agn_exp_map =  cut_image(exp_map, center_QSO, fr)
agn_exp_map[agn_exp_map==0] = agn_exp_map[agn_exp_map!=0].min()
agn_rms =(abs(agn_image/agn_exp_map)+stdd**2)**0.5
plt.imshow(agn_rms, norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()

#==============================================================================
# Cutout and compare the PSF
#==============================================================================
psfs_pos = np.array([[3599,1513], [3456,798]])
dist_psfs = (psfs_pos-center_QSO)[:,0]**2+(psfs_pos-center_QSO)[:,1]**2
psfs_pos = psfs_pos[dist_psfs.argsort()]
count = 0
PSF_gauss_centers, PSF_bright_centers=[],[]
PSFs_ = []
for i in range(len(psfs_pos)):
    print 'PSF',count
    PSF, PSF_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=75,  return_center=True, plot=False, center_float=True)
    PSF_gauss_centers.append(PSF_center)
    _, PSF_br_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=75, kernel = 'center_bright', return_center=True, plot=False)
    PSF_bright_centers.append(PSF_br_center)
    #Create and save PSFs mask
    _, PSF_mask = mask_obj(img=PSF, exp_sz=1.4)
    if len(PSF_mask) > 1:
        PSF_mask = np.sum(np.asarray(PSF_mask),axis=0)
    elif len(PSF_mask) == 1:
        PSF_mask = PSF_mask[0]
    PSF_mask = (1 - (PSF_mask != 0)*1.)
    PSFs_.append([PSF,PSF_center,PSF_mask])
    count += 1
center_match = (np.sum(abs(np.asarray(PSF_gauss_centers)-np.asarray(PSF_bright_centers)<0.7),axis = 1) == 2)
psf_pos_list = [PSF_gauss_centers[i] for i in range(len(PSFs_))] # if center_match[i]==True]
PSFs = [PSFs_[i] for i in range(len(PSFs_))] # if center_match[i]==True]
save_loc_png(img_sub,center_QSO,psf_pos_list, ID=ID, label='PSF-' ,reg_ty = None)
PSF_list = [PSFs[i][0] for i in range(len(PSFs))]
fig = profiles_compare(PSF_list, scal_list=np.ones(len(PSFs)),
                       prf_name_list=['PSF'+str(i) for i in range(len(PSFs))],
                       gridspace = 'log', if_annuli=True)
# =============================================================================
# Test the bkg level for lensed QSO and PSF
# =============================================================================
for i in range(len(PSF_list)):
    _, _, _ = flux_profile(PSF_list[i], [len(PSF_list[i])/2]*2,
                            radius=len(PSF_list[i])/2, grids=50, ifplot=True, fits_plot= True)

#==============================================================================
# Save the lens fitting ingredients including:
#    lensing image and noise map
#    PSFs
#==============================================================================
pyfits.PrimaryHDU(agn_image).writeto('{0}_cutout.fits'.format(ID),overwrite=True) 
pyfits.PrimaryHDU(agn_rms).writeto('{0}_stdd.fits'.format(ID),overwrite=True) 
pyfits.PrimaryHDU(agn_exp_map).writeto('{0}_exp_map.fits'.format(ID),overwrite=True) 
for i in range(len(PSFs)):
    pyfits.PrimaryHDU(PSF_list[i]).writeto('PSF{0}.fits'.format(i),overwrite=True) 

