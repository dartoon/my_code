#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:02:18 2022

@author: Dartoon
"""

#To make sure the used PSFs are stars.

import numpy as np
import pickle
# from photutils.segmentation import SourceCatalog
# from source_info import source_list

filt_id = 0 #int(sys.argv[2])
filt = ['f150w', 'f356w'][filt_id]
folder = '../prep_use_HST_highRes/JWST_CEERS/'

from psf_if_star import if_star
# #Before select:
# _psfs, _FWHMs = pickle.load(open('../prep_use_HST_highRes/'+filt+'_psfs.pkl','rb'))
# psfs, FWHMs, fluxs = [], [], []
# for i in range(len(_psfs)):
#     psfs = psfs + _psfs[i]
#     FWHMs = FWHMs + _FWHMs[i]
#     fluxs = fluxs + [np.sum(_psfs[i][j]) for j in range(len(_psfs[i]))]
# FWHMs, fluxs = np.array(FWHMs), np.array(fluxs)


# select_ids = []
# for i in range(len(psfs)):
#     if if_star(psfs[i], filt) == True:
#         select_ids.append(i)
        
# select_psfs = [psfs[i] for i in select_ids]
# select_FWHMs = [FWHMs[i] for i in select_ids]

# import pickle
# pickle.dump([select_psfs,select_FWHMs], open(filt+'_psfs_star.pkl', 'wb'))   


# Check if OK:
psfs, FWHMs = pickle.load(open(filt+'_psfs_star.pkl','rb'))
select_ids = []
for i in range(len(psfs)):
    if if_star(psfs[i], filt) == True:
        select_ids.append(i)