#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:15:53 2021

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
# ID_list =  ['00:00:50.56 -01:30:55.2', '00:14:59.72 +00:23:19.2', #'00:29:38.47 +00:20:39.7',   #00:29:38.47 +00:20:39.7 Removed because too faint
#             '00:36:59.26 -00:19:22.7', '00:36:59.43 -00:18:50.2', '01:12:27.87 -00:31:51.6', 
#             '01:21:10.93 +01:07:03.3', '01:37:36.57 +00:57:42.3', '01:38:34.18 -00:05:09.3', 
#             '01:39:30.80 -00:21:31.6',  #01:39:30.80 -00:21:31.6' remove because too faint
#             '02:24:04.85 +01:49:41.9', '02:29:06.04 -05:14:28.9', 
#             '02:36:00.28 -01:04:32.3', '02:38:29.90 -01:12:24.2', '09:06:54.53 +02:13:15.2', 
#             '09:25:32.13 -02:08:06.1', '09:52:18.04 -00:04:59.1', '10:46:44.31 +00:03:29.7', 
#             '10:54:58.01 +04:33:10.6', '12:46:18.51 -00:17:50.2', '13:15:12.46 +01:50:21.6', 
#             '13:24:41.58 -01:54:01.8', '13:42:57.16 -01:39:12.9', '15:02:16.66 +02:57:19.8', 
#             '15:30:08.91 +42:56:34.8', '16:25:01.98 +43:09:31.6', '22:06:42.82 +00:30:16.2', 
#             '22:10:11.62 -00:16:54.9', '22:11:01.45 +00:14:49.0', '23:04:02.77 -00:38:55.4', 
#             '23:37:18.07 +00:25:50.6', '12:04:17.10 +00:36:53.7', '12:46:04.03 -01:09:54.6', 
#             '12:52:16.06 +00:31:41.1', '14:43:08.16 -00:49:13.4', '14:53:47.46 +00:39:27.0', 
#             '22:04:22.46 +07:31:38.3', '22:09:10.38 -00:16:01.5', '15:21:12.96 +44:14:52.5']

from ID_list import ID_list

import pandas as pd
sample = pd.read_csv('material/Gaia_catalog.csv')

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
    
    fit_size = len(pyfits.getdata(glob.glob(files[-1]+'data-BHBH(host image)_I-band.fits')[0]))
    radius = int((fit_size-1)/2)
    for j in run_list:
        data_process_list[j].generate_target_materials(radius=radius,
                                                       create_mask = False, nsigma=2.,
                                                       exp_sz= 1.2, npixels = 15, if_plot=False)
    
    #%%
    from matplotlib.colors import LogNorm
    fig, (axs) = plt.subplots(1, 5, figsize=(15, 7))
    vmin = 1.e-3
    vmax = data_process_list[0].target_stamp.max() * 5
    color_list = ['winter', 'summer', 'afmhot', 'autumn', 'gist_heat']
    plt_list = [1, 2, 0, 3, 4]
    AGN_mag_list = [[-99,-99]] * 5
    AGN_pos_list = []
    host_mag_list = [[-99]] * 5
    AGN_RA_DEC_list = [[-99,-99]] * 5
    for i in range(len(plt_list)):
        p_i = plt_list[i]
        if data_process_list[p_i] != None:
            axs[i].imshow(data_process_list[p_i].target_stamp, origin='lower', cmap=color_list[i], norm=LogNorm(), vmin=vmin, vmax=vmax)
            file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format(band_seq[p_i]))
            # print('fit_result_{0}-band.txt'.format(band_seq[p_i]))
            if file != []:
                f = open(file[0],"r")    
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            trust = 2
            if ID == '120417.10+003653.7':
                trust = 1
            l0 = [i for i in range(len(lines)) if 'AGN0 position:' in lines[i]]
            AGN_RA_DEC_list[i] = [[ float(lines[l0[trust-1]].split('RA: ')[1][:11]), float(lines[l0[trust-1]].split('DEC: ')[1][:11])  ],
                          [ float(lines[l0[trust-1]].split('RA: ')[2][:11]), float(lines[l0[trust-1]].split('DEC: ')[2][:11])  ]]
            # AGN_RA_DEC_list.append(AGN_RA_DEC)
            l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
            AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])
            AGN_mag = [AGN_dic[i]['magnitude'] for i in range(2)]
            AGN_mag_list[i] = AGN_mag #In order of G, R, I, Z, Y
            # AGN_0_pos = np.array([AGN_dic[0]['ra_image'], AGN_dic[0]['dec_image']])
            # AGN_1_pos = np.array([AGN_dic[1]['ra_image'], AGN_dic[1]['dec_image']])
            AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
            AGN_pos_list.append(AGN_pos)
            AGN_pos = AGN_pos/0.168 + radius
            l2 = [i for i in range(len(lines)) if 'model_Sersic_result:' in lines[i]]
            galaxy_dic = read_string_list(string = lines[l2[trust]].split('model_Sersic_result: ')[1])
            galaxy_mag = [galaxy_dic[i]['magnitude'] for i in range(len(galaxy_dic))]
            host_mag_list[i] = galaxy_mag
            galaxy_pos = np.array([[-1*galaxy_dic[i]['center_x'], galaxy_dic[i]['center_y']] for i in range(len(galaxy_dic))])    
            galaxy_pos = galaxy_pos/0.168 + radius
            apertures = []
            for j in range(len(galaxy_dic)):
                a = galaxy_dic[j]['R_sersic']/0.168
                b = galaxy_dic[j]['R_sersic']/0.168 * galaxy_dic[j]['q']
                theta = galaxy_dic[j]['phi_G'] 
                apertures.append(EllipticalAperture(galaxy_pos[j], a, b, theta=theta))            
            np.random.seed(seed = 4)
            for j in range(len(AGN_pos)):
                label = None
                if i == 0:
                    label = 'PS {0}'.format(j)            
                axs[i].scatter(AGN_pos[j][0], AGN_pos[j][1], 
                            color=(np.random.uniform(0, 1), np.random.uniform(0, 1), np.random.uniform(0, 1)),
                            s=180, marker=".",label = label)
            np.random.seed(seed = 3)
            for j in range(len(apertures)):
                label = None
                if i == 0:
                    label = 'comp {0}'.format(j)
                aperture = apertures[j]
                if host_mag_list[i][j] < 24.5:
                    aperture.plot(color= (np.random.uniform(0, 1), np.random.uniform(0, 1), np.random.uniform(0, 1)),
                                  axes=axs[i], lw=3.5, label = label)     
            axs[0].legend()
        else:
            axs[i].imshow(data_process_list[0].target_stamp * 0)
            AGN_pos_list.append([-99,-99])
        axs[i].set_title('{0} band'.format(band_seq[p_i]))
    # plt.savefig('savefig/{0}_fiveband.png'.format(ID))
    # plt.close()
    from astropy.visualization import make_lupton_rgb
    Band = ['G','R','I','Z','Y',]
    data_process_list[0], data_process_list[1], data_process_list[2] = data_process_list[1], data_process_list[2], data_process_list[0]
    Band = [Band[i] for i in range(len(Band)) if data_process_list[i] !=None ]
    data_process_list = [data_process_list[i] for i in range(len(data_process_list)) if data_process_list[i] !=None ]
    r_i, g_i, b_i = [data_process_list[i].target_stamp for i in [2,1,0]]
    rgb_default = make_lupton_rgb(r_i, g_i, b_i, stretch = 1)
    
    # #%%
    # blue_line = np.loadtxt('Shenli_materials/stellar_blueline.txt')
    # plt.figure(figsize=(11, 11))
    # plt.plot(blue_line[:,0], blue_line[:,1],linewidth = 5, zorder = -1)
    
    # g_r_0 = AGN_mag_list[1][0] - AGN_mag_list[2][0]
    # r_i_0 = AGN_mag_list[2][0] - AGN_mag_list[0][0]
    # g_r_1 = AGN_mag_list[1][1] - AGN_mag_list[2][1]
    # r_i_1 = AGN_mag_list[2][1] - AGN_mag_list[0][1]
    # plt.scatter([g_r_0,g_r_1], [r_i_0, r_i_1],color = '#C942C7')
    # plt.plot([g_r_0,g_r_1], [r_i_0, r_i_1],color = '#C942C7')
    # # plt.text([g_r_0,g_r_1], [r_i_0, r_i_1], ['0', '1'],color = '#C942C7')
    
    # plt.text(g_r_0,r_i_0, '0',fontsize=25)
    # plt.text(g_r_1,r_i_1, '1',fontsize=25)
    # plt.xlabel('HSC g-r',fontsize=27)
    # plt.ylabel('HSC r-i',fontsize=27)
    # if np.min([g_r_0,g_r_1]) > -2:
    #     plt.xlim([-0.5, 2.5])
    # if np.min([r_i_0,r_i_1]) > -2:    
    #     plt.ylim([-0.5, 3.0])
    # plt.tick_params(labelsize=20)
    # plt.show()
    
    #%%print information
    print("\n\nID order: ",  ID, 'z=', read_z(ID), RA, Dec )  
    offset_list = [-99] * 5
    for i in [0, 1, 2, 3, 4]:
        offset = np.sum( (AGN_pos_list[i][0] -  AGN_pos_list[i][1])**2)**0.5 
        # print('POS {1} {2}  Offset {0:.2f} '.format(offset, AGN_RA_DEC_list[i][0], AGN_RA_DEC_list[i][1]) )
        offset_list[i] = offset
        print(['G','R','I','Z','Y'][i] + ' band PS Mag: {0:.2f}, {1:.2f} \t Host0 Mag: {2:.2f} Offset {3:.2f}'.format( AGN_mag_list[i][0], 
                                                                                                                   AGN_mag_list[i][1], 
                                                                                                                   host_mag_list[i][0], 
                                                                                                                   offset))
        if len(rgb_default) > 50:
            show_size = 47
            cut =  int(  (len(rgb_default) - 47) / 2 )
            rgb_default = rgb_default[cut:-cut,cut:-cut,:]
    info = ['None'] * 2    
    for i in range(2):
        _ra, _dec  = AGN_RA_DEC_list[0][i]
        _dis = [ np.sqrt((_ra - sample['ra'][k])**2 + (_dec - sample['dec'][k])**2) for k in range(len(sample))]
        if np.min(_dis)< 1.e-04:
            k = np.where(_dis == np.min(_dis))[0][0]
            info[i] = ' parallax {0:.2f}, para_err {1:.2f}, g band {2:.2f}'.format(sample['parallax'][k], sample['parallax_over_error'][k], sample['phot_g_mean_mag'][k])
    print('GAIA info: PS0:'+ info[0] +'\tPS1:' + info[1])
    
    plt.show() 
    sz = len(rgb_default)
    plt.text(sz/20,sz/20*18,ID,color='white',fontsize=15)
    plt.text(sz/20,sz/20*16.5,"I-band offset: {0:.2f} arcsec".format(offset_list[2]),fontsize=15,color='white')
    plt.text(sz/20,sz/20*15,"I-band mag: {0:.2f} {1:.2f}".format(AGN_mag_list[2][0], AGN_mag_list[2][1]),fontsize=15,color='white')
    plt.text(sz/20,sz/20*1, Band[0]+Band[1]+Band[2]+'-band Used',color='white',fontsize=15)
    plt.imshow(rgb_default, origin='lower')
    # plt.savefig('savefig/{0}_colorimag.png'.format(ID))
    plt.show()
    # plt.close()
