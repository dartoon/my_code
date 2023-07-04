#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:31:18 2023

@author: Xuheng Ding
"""
#%%
import numpy as np
import astropy.io.fits as pyfits
import glob, copy
from astropy.wcs import WCS
import warnings
warnings.filterwarnings("ignore")
from astropy.wcs.utils import proj_plane_pixel_scales

#%% Settings for the cutout
folder = '/lustre/work/JWST/'

cat_file = 'CEERS_cutout.cat'
fclts = []
filts = []

#%%
def RA_Dec_in_fit(all_files, RA,Dec):
    filename = []
    for i in range(len(all_files)):
        file = all_files[i]
        fitsFile = pyfits.open(file)
        try:
            fov_image = fitsFile['SCI'].data # check the back grounp
            header = fitsFile['SCI'].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        except:
            fov_image = fitsFile['PRIMARY'].data # check the back grounp
            header = fitsFile['PRIMARY'].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
        wcs_i = WCS(header)
        pos = wcs_i.all_world2pix([[RA, Dec]], 1)[0]
        try:
            if pos[0]>0 and pos[1]>0 and fov_image[int(pos[1])-1, int(pos[0])-1 ] > 0 :
                filename.append(file)
        except:
            None
    return filename

folder = '/lustre/work/JWST/grizli_v6/'

_fitsfiles = glob.glob(folder+'*sci.fits')

_fitsfiles = [_fitsfiles[i] for i in range(len(_fitsfiles)) if 'wht' not in _fitsfiles[i]]
fclt_fitsfiles = []
for file in _fitsfiles:
    for fclt in fclts:
        if fclt in file:
            fclt_fitsfiles.append(file)
if fclts == []:
    fclt_fitsfiles = _fitsfiles
            
fitsfiles = []
for filt in filts:
    fitsfiles = fitsfiles + [fclt_fitsfiles[i] for i in range(len(fclt_fitsfiles)) if filt in fclt_fitsfiles[i]]
if filts == []:
    fitsfiles = fclt_fitsfiles

dtype = np.dtype([("name", "U30"), ("RA", float), ("Dec", float),("rad", float)])
targets = np.loadtxt(cat_file, dtype=dtype)

if targets.shape == ():
    targets = [[str(targets['name']), float(targets['RA']), float(targets['Dec']),float(targets['rad'])]]


#%%
for target in targets:
    name, RA, Dec, cut_rad = target
    # print(len(fitsfiles))
    files = RA_Dec_in_fit(fitsfiles, RA, Dec)
    # print(len(files))
    i = 0
    for file in files:
        cut_SCI_image, cut_Err_image, cut_WHT_image = None, None, None
        fitsFile = pyfits.open(file)
        # layer_names = [fitsFile[i].name for i in range(len(fitsFile))]
        # print(layer_names)
        fitsFile.info()
        SCI_image = fitsFile['PRIMARY'].data
        SCI_header = fitsFile['PRIMARY'].header
        
        whtfitsFile = pyfits.open(file.replace('sci', 'wht') )
        WHT_image = whtfitsFile['PRIMARY'].data
        WHT_header = whtfitsFile['PRIMARY'].header
        try:
            filt = fitsFile[0].header['FILTER']
        except KeyError:
            filt = fitsFile[0].header['FILTER1'] 
        fac = fitsFile[0].header['INSTRUME'] 
        fov = file.split('/')[-1].split('-')[0]
            
        wcs = WCS(SCI_header)
        scales = proj_plane_pixel_scales(wcs) * 3600  #From degree to arcsec
        if abs(scales[0]-scales[1])/scales[0]>1.e-5:
            print('Warning: Pixel scale is not the same along x and y!!!')
        pix_scale = scales[0] 
        
        ct = int(cut_rad/pix_scale)
        pos = WCS(SCI_header).all_world2pix([[RA, Dec]], 1)[0]
        
        file_header = copy.deepcopy(SCI_header)
        if ct<pos[0]:
            file_header['CRPIX1'] = file_header['CRPIX1']-pos[0]+ct
        if ct<pos[1]:
            file_header['CRPIX2'] = file_header['CRPIX2']-pos[1]+ct
        
        savename = copy.deepcopy(name)
        if 'nircam_nircam' in file:
            savename = name +'_'+ file.split('_nircam_')[1].split('_')[0]
        
        hdul = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU(None,header=fitsFile['PRIMARY'].header)
        hdul.insert(0, hdu)    
        
        cuts_reg = np.array([int(pos[1])-ct,int(pos[1])+ct, int(pos[0])-ct,int(pos[0])+ct])
        cuts_reg[cuts_reg<0] = 0
        x1, x2, y1, y2 = cuts_reg
        print(np.shape(SCI_image))
        cut_SCI_image = SCI_image[x1:x2, y1:y2]
        if 'EXTNAME' not in file_header:
            file_header['EXTNAME'] = 'SCI'
        hdu = pyfits.ImageHDU(cut_SCI_image,header=file_header)
        hdul.insert(1, hdu)
        cut_WHT_image = WHT_image[x1:x2, y1:y2]
        if 'EXTNAME' not in WHT_header:
            WHT_header['EXTNAME'] = 'WHT'
        hdu = pyfits.ImageHDU(cut_WHT_image,header=WHT_header) #.writeto(savename+'_' +fac + '_' + filt + '_' +'WHT' + '.fits',overwrite=True)
        hdul.insert(2, hdu)
        
        print(hdul.info())
        hdul.writeto(savename+ '_' + fov + '_' +fac + '_' + filt+'_'+'cutout_' + str(i) + '.fits',overwrite=True)
        i = i+1


