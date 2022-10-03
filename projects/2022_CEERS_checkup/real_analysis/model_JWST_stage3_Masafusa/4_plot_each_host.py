#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:34:41 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from functions_for_result import esti_smass, load_prop, load_info
from scipy.ndimage import zoom
import copy, matplotlib
from matplotlib.colors import LogNorm

# ID, mags, z = 'idx0', 
# 1,2,0,51,35
idx = 51
# root_folder = '../*/*'  #Include HST
root_folder = './*'  #Not include HST
fit_run_dict = load_prop(idx, root_folder = root_folder, prop_name='fit_run')
# filt_list = list(fit_run_dict.keys())
use_filt = ['F200W', 'F444W']
image_list, label = [],[]
deltaPix_list = []
for filt in use_filt:
    fit_run = fit_run_dict[filt]
    data = fit_run.flux_2d_out['data']
    host = fit_run.flux_2d_out['data-Point Source']
    image_list.append(data)
    image_list.append(host)
    label.append(filt+'\ndata')
    label.append(filt+'\nhost (data-quasar)')
    deltaPix_list.append(fit_run.fitting_specify_class.deltaPix)
    deltaPix_list.append(fit_run.fitting_specify_class.deltaPix)

my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=15):
    if flipped:
        p0 = d - d / 15. - dist
        p1 = d / 15.
        ax.plot([p0, p0 + dist], [p1, p1], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    else:
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=2, color=color)
        ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')


def total_compare(flux_list_2d, label_list_2d,
                  deltaPix_list , target_ID = 'target_ID',
                  if_annuli=False, z=None ,
                  show_plot = True, cmap=None):
    """
    Make quick plots to compare the flux profiles in a list and show the normalized residual.
    
    Parameter
    --------
        flux_list_2d: 
            A list of 2D flux array, that will use plt.imshow() to plot and show.
            e.g., [data, pointsource_list, galaxy_model_list, normalized residual]
            
        label_list_2d: 
            A list of lables for flux_list_2d.
            e.g., ['data', 'model', 'point source(s)', 'galaxy(s)']
            
        flux_list_1d:  
            A list of 2D flux array, that will be plot as 1D profile in the very right panel.
        
        label_list_1d: 
            The labels for flux_list_1d.
        
        mask_image: 
            A 2D mask for the flux_list_2d image.
        
        arrows: bool. 
            If show the arrows for pointing the North and East.
        
        if_annuli: bool.
            If True, the 1D profile will show the surface brightness in the annuli apertures. 
    """
    # norm = LogNorm() #ImageNormalize(stretch=SqrtStretch())
    cl_num = len(flux_list_2d)
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    
    f = plt.figure(0, figsize=(6.5+ (cl_num-1)*3.5,4))    
    ax_l = [plt.subplot2grid((6,cl_num), (0,i), rowspan=6) for i in range(len(flux_list_2d))] #The image plot
    if cmap == None:
        cmap = my_cmap
    for i in range(len(flux_list_2d)):
        frame_size = len(flux_list_2d[i])
        if i >1:
            flux_list_2d[i] = flux_list_2d[i]
        if i == 0:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap, norm=LogNorm(vmax = flux_list_2d[0].max(), vmin = 1.e-4))
            clim=im_i.properties()['clim'] #To uniform the color bar scale.
            # ax_l[i].set_ylabel(target_ID, fontsize=15, weight='bold')
        else:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap, norm=LogNorm(), clim=clim)
            ax_l[i].get_yaxis().set_visible(False)
        ax_l[i].get_xaxis().set_visible(False)
        deltaPix = deltaPix_list[i]
        scale_bar(ax_l[i], frame_size, dist=0.5/deltaPix, text='0.5"', color = 'white')
        ticks= np.array([1.e-4, 1.e-3, 1.e-2,1.e-1,0, 10])
        cb_i = f.colorbar(im_i, ax=ax_l[i], shrink=0.48, pad=0.01,  orientation="horizontal", 
                          aspect=15, ticks=ticks)
        fontsize = 17
        ax_l[i].text(frame_size*0.05, frame_size*0.8, label_list_2d[i],fontsize=fontsize, 
                     weight='bold', color='white')
        
        if z != None:
            ax_l[i].text(frame_size*0.95, frame_size*0.05, 'z={0}'.format(z),fontsize=fontsize, 
                         ha='right', va='bottom',
                         weight='bold', color='white', bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
        
        if i == 0:
            ax_l[i].text(frame_size*0.95, frame_size*0.95, target_ID,fontsize=fontsize, 
                         ha='right', va='top',
                         weight='bold', color='white',  bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})

    #Plot normalized residual map:
    norm_residual = flux_list_2d[-1]
    plt.subplots_adjust(wspace=-0.5, hspace=0)
    plt.show()       
    return f

name_list = {1:'SDSS1420A', 2: 'SDSS1420B', 0: 'SDSS1419', 
             51: 'AEGIS 477', 35: 'AEGIS 482'}
ID, z = load_info(idx)
fig = total_compare(image_list, label, deltaPix_list, target_ID=name_list[idx], z=z)
fig.savefig('outcomes/datahost_{0}.pdf'.format(idx))
