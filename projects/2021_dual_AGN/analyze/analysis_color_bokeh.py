#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:05:16 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob 
from ast import literal_eval
# Run_ID_list =  ['00:00:50.56 -01:30:55.2', '00:14:59.72 +00:23:19.2', '00:29:38.47 +00:20:39.7', 
#             '00:36:59.26 -00:19:22.7', '00:36:59.43 -00:18:50.2', '01:12:27.87 -00:31:51.6', 
#             '01:21:10.93 +01:07:03.3', '01:37:36.57 +00:57:42.3', '01:38:34.18 -00:05:09.3', 
#             '01:39:30.80 -00:21:31.6', '02:24:04.85 +01:49:41.9', '02:29:06.04 -05:14:28.9', 
#             '02:36:00.28 -01:04:32.3', '02:38:29.90 -01:12:24.2', '09:06:54.53 +02:13:15.2', 
#             '09:25:32.13 -02:08:06.1', '09:52:18.04 -00:04:59.1', '10:46:44.31 +00:03:29.7', 
#             '10:54:58.01 +04:33:10.6', '12:46:18.51 -00:17:50.2', '13:15:12.46 +01:50:21.6', 
#             '13:24:41.58 -01:54:01.8', '13:42:57.16 -01:39:12.9', '15:02:16.66 +02:57:19.8', 
#             '15:30:08.91 +42:56:34.8', '16:25:01.98 +43:09:31.6', '22:06:42.82 +00:30:16.2', 
#             '22:10:11.62 -00:16:54.9', '22:11:01.45 +00:14:49.0', '23:04:02.77 -00:38:55.4', 
#             '23:37:18.07 +00:25:50.6', '12:04:17.10 +00:36:53.7', '12:46:04.03 -01:09:54.6', 
#             '12:52:16.06 +00:31:41.1', '14:43:08.16 -00:49:13.4', '14:53:47.46 +00:39:27.0', 
#             '22:04:22.46 +07:31:38.3', '22:09:10.38 -00:16:01.5', '15:21:12.96 +44:14:52.5']

from ID_list import ID_list as Run_ID_list

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

ID_list, g_r_list, r_i_list, AGN_mag_list, AGN_offset_list, AGN_RA_DEC_list = [], [], [], [], [], []
for k in range(len(Run_ID_list)):
    string = Run_ID_list[k]
    string = string.replace(':', '')
    string = string.replace(' ', '')
    ID = [string][0]
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
    band_seq = ['G', 'R', 'I', 'Z', 'Y']
    AGN_mags = [[-99,-99]] * 5
    AGN_RA_DEC = [ [[-99,-99], [-99,-99]] ]* 5
    AGN_offset = [-99] * 5
    for i in range(len(band_seq)):
        file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format(band_seq[i]))
        if file != []:
            if file != []:
                f = open(file[0],"r")    
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            trust = 2
            if ID == '120417.10+003653.7':
                trust = 1            
            l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
            AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])
            AGN_mag = [AGN_dic[i]['magnitude'] for i in range(2)]
            AGN_mags[i] = AGN_mag
            AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
            AGN_offset[i] =  np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5 
            l0 = [i for i in range(len(lines)) if 'AGN0 position:' in lines[i]]            
            AGN_RA_DEC[i] = [[ float(lines[l0[trust-1]].split('RA: ')[1][:11]), float(lines[l0[trust-1]].split('DEC: ')[1][:11])  ],
                          [ float(lines[l0[trust-1]].split('RA: ')[2][:11]), float(lines[l0[trust-1]].split('DEC: ')[2][:11])  ]]
    g_r_0 = AGN_mags[0][0] - AGN_mags[1][0]
    r_i_0 = AGN_mags[1][0] - AGN_mags[2][0]
    g_r_1 = AGN_mags[0][1] - AGN_mags[1][1]
    r_i_1 = AGN_mags[1][1] - AGN_mags[2][1]
    if np.min(AGN_mags[:3]) != -99:
        g_r_list.append([g_r_0,g_r_1])
        r_i_list.append([r_i_0, r_i_1])
        ID_list.append(ID)
        AGN_mag_list.append(AGN_mags)
        AGN_offset_list.append(AGN_offset)
        AGN_RA_DEC_list.append(AGN_RA_DEC)

import pandas as pd
sample = pd.read_csv('material/Gaia_catalog.csv')

#%%
from bokeh.models import ColumnDataSource, Select, RangeSlider, Slider
import random
from bokeh.plotting import figure, output_file, show
output_file("color_fig.html")
blue_line = np.loadtxt('Shenli_materials/stellar_blueline.txt')
p = figure(plot_width=800, plot_height=800, x_axis_label = 'HSC g-r', y_axis_label = 'HSC r-i')
p.line(blue_line[:,0], blue_line[:,1], line_width=4)
random.seed(4)
for i in range(len(g_r_list)):
    color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    info = ['None'] * 2
    for j in range(2):
        _ra, _dec  = AGN_RA_DEC_list[i][2][j]
        _dis = [ np.sqrt((_ra - sample['ra'][k])**2 + (_dec - sample['dec'][k])**2) for k in range(len(sample))]
        if np.min(_dis)< 1.e-04:
            k = np.where(_dis == np.min(_dis))[0][0]
            info[j] = ' parallax {0:.2f}, para_err {1:.2f}'.format(sample['parallax'][k], sample['parallax_over_error'][k])
    qso_com = ColumnDataSource(data={
        'g_r_0' : [g_r_list[i][0], -99],
        'r_i_0' : [r_i_list[i][0], -99],
        'g_r_1' : [g_r_list[i][1], -99],
        'r_i_1' : [r_i_list[i][1], -99],
        'ID' : [ID_list[i]]*2,
        'imgs0' : ['savefig/{0}_colorimag.png'.format(ID_list[i])]*2,
        'imgs1' : ['savefig/{0}_fiveband.png'.format(ID_list[i])]*2,
        'offset': ["offset: {0:.2f}".format(AGN_offset_list[i][2])]*2,
        'RADEC': ['PS0 PS1 RA Dec {0} {1}'.format(AGN_RA_DEC_list[i][2][0], AGN_RA_DEC_list[i][2][1])]*2,
        'Gmag': ["G band mag {0:.2f} {1:.2f}".format(AGN_mag_list[i][0][0], AGN_mag_list[i][0][1])]*2,
        'Rmag': ["R band mag {0:.2f} {1:.2f}".format(AGN_mag_list[i][1][0], AGN_mag_list[i][1][1])]*2,
        'Imag': ["I band mag {0:.2f} {1:.2f}".format(AGN_mag_list[i][2][0], AGN_mag_list[i][2][1])]*2,
        'gaia': ['GAIA info:'+ info[0] + info[1]]*2
    })    
    p.circle('g_r_0', 'r_i_0', source=qso_com, size=17, color=color, alpha=0.5, name="foo")
    p.circle('g_r_1', 'r_i_1', source=qso_com, size=17, color=None, line_color=color, alpha=0.5, name="foo")
    p.line(g_r_list[i],  r_i_list[i], line_width=1, color = color)
from bokeh.models import HoverTool,Range1d
p.x_range = Range1d(-0.4, 1.6)
p.y_range = Range1d(-0.5, 3)
hover = HoverTool(names=["foo"],
        tooltips="""
        <div>
            <div>
                <img
                    src="@imgs1" height="450"
                    style="float: left; margin: 0px 15px 15px 0px;"
                    border="-10"
                ></img>
            </div>
            <div>
                <img
                    src="@imgs0" width="300"
                    style="float: left; margin: 0px 15px 15px 0px;"
                    border="0"
                ></img>
            </div>
            <div>
                <span style="font-size: 17px; font-weight: bold;">@ID</span>
            </div>           
            <div>
                <span style="font-size: 17px; font-weight: bold;">@Gmag</span>
            </div>    
            <div>
                <span style="font-size: 17px; font-weight: bold;">@Rmag</span>
            </div>    
            <div>
                <span style="font-size: 17px; font-weight: bold;">@Imag</span>
            </div>                
            <div>
                <span style="font-size: 17px; font-weight: bold;">@RADEC</span>
            </div>                  
            <div>
                <span style="font-size: 17px; font-weight: bold;">@offset</span>
            </div>       
            <div>
                <span style="font-size: 17px; font-weight: bold;">@gaia</span>
            </div>             
        </div>
        """
    )
p.add_tools(hover)
p.xaxis.axis_label_text_font_size = "25pt"
p.yaxis.axis_label_text_font_size = "25pt"
p.yaxis.major_label_text_font_size = "15pt"
p.xaxis.major_label_text_font_size = "15pt"
show(p)
    