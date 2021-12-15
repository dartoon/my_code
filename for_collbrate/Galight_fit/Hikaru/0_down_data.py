#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:27:48 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import glob

#%%Cutout using online cutout tool:
from galight.hsc_utils import hsc_image, hsc_psf
import os
import pandas as pd

part = 0
if part == 0:
    qso_type = 'Radio_DOGs'
elif part == 1:
    qso_type = 'all_DOGs'
sample = pd.read_csv('S19A_{0}_information.csv'.format(qso_type))
ID_list = sample['object_id_1']
RA_list = sample['ra_1']
Dec_list = sample['dec_1']

#%%
point_source_num = 1  #Number of AGN in the target.
ps_pix_center_list=[[0,0]] #Force the center as AGN position
dr='dr3'
rerun='s19a'
show_plot = False
if rerun!='s19a':
    psf_rerun = rerun
else:
    psf_rerun = 's20a'
    
#%%Some settings for the fitting
bands = 'GRIZY'  #Band that will be download
for i_ in range(0, len(ID_list)): 
    object_id,ra,dec= '{0}_'.format(i_)+str(ID_list[i_]), RA_list[i_], Dec_list[i_]
    print(object_id)
    #%%Mkdir folder and start downloading
    out_dir='./' + qso_type +'/' + object_id
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if glob.glob(out_dir+'/*psf*fits') == []:
        print('Downloading data with PSF... ... ...')
        try:
            hsc_image.get_cutouts(object_id,ra,dec,out_dir,dr=dr,rerun=rerun+'_dud',filters=bands,fov_arcsec=120)
            hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr=dr,rerun=psf_rerun+'_dud',filters=bands)
        except:
            hsc_image.get_cutouts(object_id,ra,dec,out_dir,dr=dr,rerun=rerun+'_wide',filters=bands,fov_arcsec=120)
            hsc_psf.get_psfs(object_id,ra,dec,out_dir,dr=dr,rerun=psf_rerun+'_wide',filters=bands)
    else:
        continue
