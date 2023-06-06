#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:31:18 2023

@author: Xuheng Ding
"""

import numpy as np
import astropy.io.fits as pyfits
import glob, copy
from astropy.wcs import WCS
import warnings
warnings.filterwarnings("ignore")
from astropy.wcs.utils import proj_plane_pixel_scales

#%% Settings for the cutout
folder = '/lustre/work/JWST/'

field = 'CEERS' #all, CEERS, COSMOSweb, PRIMER_COSMOS
cat_file = 'CEERS_cutout.cat'

# field = 'COSMOSweb' #all, CEERS, COSMOSweb, PRIMER_COSMOS
# cat_file = 'COSMOS_coutout.cat'

# name, RA, Dec, cut_rad 
#fclts = ['acs', 'wfc3', 'nircam', 'miri']  #
#filts = ['f814w']
fclts = []
filts = []

#%%
def RA_Dec_in_fit(all_files, RA,Dec):
    filename = []
    for i in range(len(all_files)):
        file = all_files[i]
        # print(file)
        fitsFile = pyfits.open(file)
        # fitsFile.info()
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

#%% To check which file cover the position 
if field == 'all':
    _fitsfiles = glob.glob(folder+'/*/'+'*fits') + glob.glob(folder+'/*/*/'+'*fits')
elif field != 'CEERS':
    _fitsfiles = glob.glob(folder+field+'/'+'*fits')
else:
    _fitsfiles = glob.glob(folder+field+'/*/'+'*fits')

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

#%% To get cutout for each targets
dtype = np.dtype([("name", "U12"), ("RA", float), ("Dec", float),("rad", float)])
targets = np.loadtxt(cat_file, dtype=dtype)

if targets.shape == ():
    targets = [[str(targets['name']), float(targets['RA']), float(targets['Dec']),float(targets['rad'])]]

for target in targets:
    name, RA, Dec, cut_rad = target
    # print(len(fitsfiles))
    files = RA_Dec_in_fit(fitsfiles, RA, Dec)
    # print(len(files))
    for file in files:
        cut_SCI_image, cut_Err_image, cut_WHT_image = None, None, None
        fitsFile = pyfits.open(file)
        layer_names = [fitsFile[i].name for i in range(len(fitsFile))]
        SCI_BKSUB_image = None
        if 'CEERS' in file:
                SCI_BKSUB_image = fitsFile['SCI_BKSUB'].data
                SCI_BKSUB_header = fitsFile['SCI_BKSUB'].header
                SCI_image = fitsFile['SCI'].data
                SCI_header = fitsFile['SCI'].header
                
                try:
                    Err_image = fitsFile['RMS'].data
                    Err_header = fitsFile['RMS'].header
                except:
                    Err_image = fitsFile['ERR'].data
                    Err_header = fitsFile['ERR'].header
                    WHT_image = fitsFile['WHT'].data
                    WHT_header = fitsFile['WHT'].header
                    SCI_header['EFFEXPTM'] = fitsFile[0].header['EFFEXPTM'] 
                try:
                    filt = fitsFile[0].header['FILTER'] 
                except:
                    filt = fitsFile[0].header['FILTER1'] 
                    if 'CLEAR' in filt:
                        filt = fitsFile[0].header['FILTER2'] 
                fac = fitsFile[0].header['INSTRUME'] 
                        
        else:
            if 'SCI' in layer_names:
                SCI_image = fitsFile['SCI'].data
                SCI_header = fitsFile['SCI'].header
                Err_image = fitsFile['ERR'].data
                Err_header = fitsFile['ERR'].header
                WHT_image = fitsFile['WHT'].data
                WHT_header = fitsFile['WHT'].header
            else:
                SCI_image = fitsFile['PRIMARY'].data
                SCI_header = fitsFile['PRIMARY'].header
                whtfitsFile = pyfits.open(file.replace('drz', 'wht') )
                WHT_image = whtfitsFile['PRIMARY'].data
                WHT_header = whtfitsFile['PRIMARY'].header
            try:
                filt = SCI_header['FILTER2'] 
                fac = SCI_header['INSTRUME'] 
            except:
                filt = fitsFile[0].header['FILTER'] 
                fac = fitsFile[0].header['INSTRUME'] 
                SCI_header['EFFEXPTM'] = fitsFile[0].header['EFFEXPTM'] 
        
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
        
        cuts_reg = np.array([int(pos[1])-1-ct,int(pos[1])-1+ct, int(pos[0])-1-ct,int(pos[0])-1+ct])
        cuts_reg[cuts_reg<0] = 0
        x1, x2, y1, y2 = cuts_reg
        
        cut_SCI_image = SCI_image[x1:x2, y1:y2]
        if 'EXTNAME' not in file_header:
            file_header['EXTNAME'] = 'SCI'
        hdu = pyfits.ImageHDU(cut_SCI_image,header=file_header)
        hdul.insert(1, hdu)
        
        if 'WHT' in layer_names:
            cut_WHT_image = WHT_image[x1:x2, y1:y2]
            if 'EXTNAME' not in WHT_header:
                WHT_header['EXTNAME'] = 'WHT'
            hdu = pyfits.ImageHDU(cut_WHT_image,header=WHT_header) #.writeto(savename+'_' +fac + '_' + filt + '_' +'WHT' + '.fits',overwrite=True)
            hdul.insert(2, hdu)
        
        if 'ERR' in layer_names or 'RMS' in layer_names:
            cut_Err_image = Err_image[x1:x2, y1:y2]
            hdu = pyfits.ImageHDU(cut_Err_image,header=Err_header) #.writeto(savename+ '_' +fac + '_' + filt+'_'+'ERR' + '.fits',overwrite=True)
            hdul.insert(3, hdu)
            
        if SCI_BKSUB_image is not None:
            cut_SCI_BKSUB_image = SCI_BKSUB_image[x1:x2, y1:y2]
            file_header = copy.deepcopy(SCI_BKSUB_header)
            if ct<pos[0]:
                file_header['CRPIX1'] = file_header['CRPIX1']-pos[0]+ct
            if ct<pos[1]:
                file_header['CRPIX2'] = file_header['CRPIX2']-pos[1]+ct
            hdu = pyfits.ImageHDU(cut_SCI_image,header=file_header)
            hdul.insert(4, hdu)
        
        print("Making cutout for:", savename, filt)
        hdul.writeto(savename+ '_' +fac + '_' + filt+'_'+'cutout' + '.fits',overwrite=True)

