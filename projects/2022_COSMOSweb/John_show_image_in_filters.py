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

def scale_bar(ax, d, dist=1/0.13, text='1"', text2=None, color='black', flipped=False, fontsize=20):
    p0 = d / 7.
    ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
    ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    if text2 is not None:
        ax.text(p0 + dist / 2., p0 - 0.08 * d, text2, fontsize=fontsize, color=color, ha='center')
        


# filt_i = 0
# filt = ['F115W', 'F150W','F277W', 'F444W'][filt_i]
# cata_list = pickle.load(open('material/cata_list.pkl','rb'))
cata_list = pickle.load(open('material/cata_list.pkl','rb'))

check_name= 'cid_473'  #29
# check_name= 'cid_1210' #8
# check_name= 'cid_1245' #10
check_id = [i for i in range(len(cata_list)) if cata_list[i][-1] == check_name]

print(cata_list[check_id[0]])

if check_name == 'cid_473':
    WCS = 149.97944107802815, 2.3090658535647433
    size = 120
if check_name == 'cid_1210':
    WCS = 149.95357779220066, 2.383045149701304
    size = 80
if check_name == 'cid_1245':
    WCS = 149.99728884735094, 2.44924997883808
    size = 140


#%%
filts = ['F115W', 'F150W','F277W', 'F444W']
image_list = []
fit_file_folder ='fit_result'
ifACS = True

for idx in check_id:
    if ifACS == True:
        filename = glob.glob('otherfiles/{0}_cutout.fits'.format(check_name))[0]
        image = pyfits.getdata(filename)
        # fitsfile = pyfits.open('/Volumes/Seagate_Expansion_Drive/data_backup/JWST_COSMOS/mosaic_cosmos_web_2023jan_30mas_hst_acs_wfc_f814w_drz.fits')
        # header = fitsfile[0].header
        # fov_image = fitsfile[0].data
        # from galight.data_process import DataProcess
        # data_process = DataProcess(fov_image = fov_image, target_pos = WCS, pos_type = 'wcs', header = header,
        #                   rm_bkglight = False, if_plot=False, zp = 27.0)  #zp use 27.0 for convinence.
        # data_process.generate_target_materials(radius=120, skip=True)
        # image = data_process.target_stamp
        image_list.append(image)
    for filt in filts:
        fit_run_list = []
        fit_files = glob.glob(fit_file_folder+'/fit2_run_{0}_idx{1}_psf*.pkl'.format(filt,idx))
        fit_files.sort()
        z = cata_list[idx][6]
        if z >0:
            zinfo = 'Zspec'+str(z)
        elif z <0:
            zinfo = 'Zphot'+str(cata_list[idx][5])
        fit_run = pickle.load(open(fit_files[0],'rb'))
        # fit_run.plot_final_qso_fit(target_ID = check_name+'_'+filt)
        image_list.append(fit_run.fitting_specify_class.data_process_class.target_stamp)

image_list_ct = []
for image in image_list:
    ct = int((len(image) - size)/2)
    image = image[ct:-ct, ct:-ct]
    image_list_ct.append(image)

if ifACS==True:
    filts  = ['ACS']+filts
#%%        
import copy, matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')

deltaPix = fit_run.fitting_specify_class.deltaPix
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
if z <0:
    z = cata_list[idx][5]
scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
kpc_per_pixel = scale_relation 
scale = 0.5 * kpc_per_pixel

fig, axs = plt.subplots(len(image_list_ct), figsize=(5,18))
for i in range(len(image_list_ct)):
    norm = LogNorm(vmin = np.std(image_list_ct[i][:,:2])*0.6, vmax =image_list_ct[i].max()/1.5 )
    axs[i].imshow(image_list_ct[i], norm=norm, origin='lower',cmap = my_cmap) 
    axs[i].set_ylabel(filts[i],fontsize=20)
    axs[i].tick_params(labelsize=15)
    if i == 0:
        scale_bar(axs[i], len(image_list_ct[i]), dist=0.5/deltaPix, text='0.5"', text2 ='{0:.2f}kpc'.format(scale), color = 'white')
show_name = check_name.replace('cid_', 'CID ')
fig.suptitle('{0}'.format(show_name),fontsize=35)
fig.tight_layout()
fig.savefig('/Users/Dartoon/Downloads/{0}_filters_image.pdf'.format(check_name))
plt.show()

print(z)