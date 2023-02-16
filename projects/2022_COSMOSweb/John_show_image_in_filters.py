#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:50:02 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
from matplotlib.colors import LogNorm
# filt_i = 0
# filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
# cata_list = pickle.load(open('material/cata_list.pkl','rb'))
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

check_name= 'cid_473'  #29
check_name= 'cid_1210' #8
check_name= 'cid_1245' #10
check_id = [i for i in range(len(cata_list)) if cata_list[i][-1] == check_name]

print(cata_list[check_id[0]])

#%%
filts = ['F115W', 'F150W','F277W', 'F444W']
image_list = []
fit_file_folder ='fit_material'

ifACS = True

for idx in check_id:
    if ifACS == True:
        filename = glob.glob('otherfiles/{0}_*sci.fits'.format(check_name))[0]
        image = pyfits.getdata(filename)
        image_list.append(image)
    for filt in filts:
        fit_run_list = []
        fit_files = glob.glob(fit_file_folder+'/fit_notrunyet_{0}_idx{1}_psf0.pkl'.format(filt,idx))
        fit_files.sort()
        z = cata_list[idx][6]
        if z >0:
            zinfo = 'Zspec'+str(z)
        elif z <0:
            zinfo = 'Zphot'+str(cata_list[idx][5])
        
        fit_run = pickle.load(open(fit_files[0],'rb'))
        image_list.append(fit_run.fitting_specify_class.data_process_class.target_stamp)
if ifACS==True:
    filts  = ['ACS']+filts
        
import copy, matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')

fig, axs = plt.subplots(len(image_list), figsize=(5,18))
for i in range(len(image_list)):
    norm = LogNorm(vmin = 0.001, vmax =image_list[i].max() )
    axs[i].imshow(image_list[i], norm=norm, origin='lower',cmap = my_cmap) 
    axs[i].set_ylabel(filts[i],fontsize=20)
    axs[i].tick_params(labelsize=15)
fig.suptitle('{0}'.format(check_name),fontsize=35)
fig.tight_layout()
plt.show()
