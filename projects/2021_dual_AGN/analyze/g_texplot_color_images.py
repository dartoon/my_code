#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:48:40 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob 

from decomprofile.data_process import DataProcess
from decomprofile.fitting_specify import FittingSpeficy
from decomprofile.fitting_process import FittingProcess
from photutils import EllipticalAperture

from ast import literal_eval


from decomprofile.tools.plot_tools import scale_bar, coordinate_arrows
def read_string_list(string):
    """
    Translate a string-list-dict to a list-dict. 
    Not allow array inside the dict...
    """
    string = ''.join(string.split('array(['))
    string = ''.join(string.split('])'))    
    string = string.replace('nan', '-99')
    string = string.replace('inf', '-99')
    string_list = string.split('{')[1:]
    string_list = [literal_eval("{"+string_list[i].split('}')[0]+"}") for i in range(len(string_list))]
    return string_list

f = open("../_pdfs_2close/DR144.4_short.asc","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
def read_z(ID):
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    if line != []:
        z = float(line[0].split(' ')[-1])
    else:
        z = -99
    return z

from ID_list import ID_list

# ID_list.append('011227.87-003151.6')
import pandas as pd
sample = pd.read_csv('material/Gaia_catalog.csv')

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
# row_ = int(len(ID_list)/6 -0.1) + 1
fig, (axs) = plt.subplots(1, 3, figsize=(9, 4))

for k in range(len(ID_list)):
    string = ID_list[k]
    string = string.replace(':', '')
    string = string.replace(' ', '')
    ID = [string][0]
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    pos = SkyCoord('{0} {1}'.format(ID[:2]+':'+ID[2:4]+':'+ID[4:9], ID[9:12]+':'+ID[12:14]+':'+ID[14:]), unit=(u.hourangle, u.deg))
    RA, Dec = pos.ra.degree, pos.dec.degree
    image_RA = float(RA)
    image_DEC = float(Dec)
    # print(image_RA, image_DEC)
    
    files_1 = glob.glob('../proof2close_HSC_images_5band/*/' + ID + '/fit_result/')
    files_2 = glob.glob('../extra/*/fit_result*/' + ID + '/')
    
    files = files_1 + files_2
    
    image_folder_0 = glob.glob('../proof2close_HSC_images_5band/*/' + ID +'/')
    image_folder_1 = glob.glob('../extra/*/' + ID + '*')
    image_folder = image_folder_0 + image_folder_1
    if 'extra' in image_folder[0]:
        image_folder = image_folder_1[0].split(ID)[0]
    else:
        image_folder = image_folder[0]
    
    #%%
    band_seq = ['I', 'G', 'R', 'Z', 'Y']
    run_list = [0, 1, 2, 3, 4]
    filename_list = [ID+'_HSC-{0}.fits'.format(band_seq[i]) for i in range(len(band_seq))]
    data_process_list, zp_list = [], []
    for i in range(len(band_seq)):
        # The pixel scale is all 0.168
        if len(glob.glob(image_folder+filename_list[i])) == 0:
            print(filename_list[i] + " DOES NOT EXIST!!!")
            QSO_im, err_map, PSF, _, _, qso_center, fr_c_RA_DEC = [], [], [], [], [], [], []
            run_list.remove(i)
            data_process_list.append(None)
            zp_list.append(None)
        else:
            fitsFile = pyfits.open(image_folder+filename_list[i])
            fov_image= fitsFile[1].data
            header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
            err_data= fitsFile[3].data ** 0.5
            file_header0 = fitsFile[0].header
            zp =  27.0 
            data_process_i = DataProcess(fov_image = fov_image, fov_noise_map = err_data,
                                         target_pos = [image_RA, image_DEC],
                                         pos_type = 'wcs', header = header,
                                         rm_bkglight = False, if_plot=False, zp = zp)
            data_process_i.noise_map = err_data
            data_process_list.append(data_process_i)
    
    try:
        fit_size = len(pyfits.getdata(glob.glob(files[0]+'data-BHBH(host image)_I-band.fits')[0]))
    except:
        fit_size = 61 
    radius = int((fit_size-1)/2)
    for j in run_list:
        data_process_list[j].generate_target_materials(radius=radius,
                                                       create_mask = False, nsigma=1,
                                                       exp_sz= 1.2, npixels = 5, if_plot=False)
    vmin = 1.e-3
    vmax = data_process_list[0].target_stamp.max() * 5
    color_list = ['winter', 'summer', 'afmhot', 'autumn', 'gist_heat']
    plt_list = [1, 2, 0, 3, 4]
    from astropy.visualization import make_lupton_rgb
    Band = ['G','R','I','Z','Y',]
    data_process_list[0], data_process_list[1], data_process_list[2] = data_process_list[1], data_process_list[2], data_process_list[0]
    Band = [Band[i] for i in range(len(Band)) if data_process_list[i] !=None ]
    data_process_list = [data_process_list[i] for i in range(len(data_process_list)) if data_process_list[i] !=None ]
    try:
        r_i, g_i, b_i = [data_process_list[i].target_stamp for i in [2,1,0]]
    except:
        continue
    rgb_default = make_lupton_rgb(r_i, g_i, b_i, stretch = 1)
    for i in [0, 1, 2, 3, 4]:
        if len(rgb_default) > 50:
            show_size = 47
            cut =  int(  (len(rgb_default) - 47) / 2 )
            rgb_default = rgb_default[cut:-cut,cut:-cut,:]
    
    # sz = len(rgb_default)
    # plt.text(sz/20,sz/20*18,ID,color='white',fontsize=15)
    # plt.text(sz/20,sz/20*16.5,"I-band offset: {0:.2f} arcsec".format(offset_list[2]),fontsize=15,color='white')
    # plt.text(sz/20,sz/20*16.5,"z={0:.3f}".format(read_z(ID)),fontsize=15,color='white')
    # plt.text(sz/20,sz/20*1, Band[0]+Band[1]+Band[2]+'-band Used',color='white',fontsize=15)
    # print(ID, "z=",read_z(ID))
    # plt.imshow(rgb_default, origin='lower')
    # _i = int(k / len(axs.T))
    _j = k
    axs[_j].imshow(rgb_default, origin='lower')
    sz = len(rgb_default)
    show_ID = ID[:4] + ID[9:14] 
    bands = Band[0]+Band[1]+Band[2]
    
    file = glob.glob(files[0]+'fit_result_{0}-band.txt'.format('I'))
    if file != []:
        f = open(file[0],"r")    
    string = f.read()    
    lines = string.split('\n')   # Split in to \n
    trust = 2    
    l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
    AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
    AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
    offset = np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5 
    
    plttext = axs[_j].text(sz*5/20,sz/20*3,"offset: {0:.2f}'""'  ".format(offset),fontsize=15,color='white')
    plttext = axs[_j].text(sz/30,sz/20*17.5,show_ID,color='white',fontsize=18)
    plttext = axs[_j].text(sz/30,sz/20*15,"z={0:.3f}".format(read_z(ID)),fontsize=18,color='white')
    plttext = axs[_j].text(sz*15/20,sz/20*17, bands,color='white',fontsize=15)
    
    cent = rgb_default.shape[0]/2
    
    if _j in [10]:
        offs = 1.5
    elif _j in [1]:
        offs = 1
    else:
        offs = 0.5
    plttext = axs[_j].scatter(AGN_pos[0][0]/0.168 + cent -offs, AGN_pos[0][1]/0.168+ cent -offs, 
                              marker="x",c='red',  s = 100, 
                              zorder=100 ,alpha=1)
    plttext = axs[_j].scatter(AGN_pos[1][0]/0.168 + cent -offs, AGN_pos[1][1]/0.168+ cent -offs,
                              marker="x",c='red',  s = 100, 
                              zorder=100 ,alpha=1)
    
    scale_bar(axs[_j], sz, dist=1/0.168, text='1"', color='white', fontsize=18)
    # if ID == '011227.87-003151.6':
    #     plttext = axs[_i][_j].text(sz/30,sz/20*13,'not selected',color='white',fontsize=18)        
    coordinate_arrows(axs[_j], sz, arrow_size=0.03, color = 'white')
    axs[_j].axes.xaxis.set_visible(False)
    axs[_j].axes.yaxis.set_visible(False)

# axs[1][4].axes.xaxis.set_visible(False)
# axs[1][4].axes.yaxis.set_visible(False)    
# axs[4].axis("off")
plt.tight_layout()    
# plt.subplots_adjust(wspace=-0.005)
plt.savefig('show_material/color_plot.pdf')
plt.show()

