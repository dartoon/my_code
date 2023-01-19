#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:58:31 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
import copy, matplotlib
from matplotlib.colors import LogNorm



def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=20):
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
        color = 'white'
        if i >1:
            flux_list_2d[i] = flux_list_2d[i]
        if i == 0:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap, norm=LogNorm(vmax = flux_list_2d[0].max(), vmin = 1.e-4))
            clim=im_i.properties()['clim'] #To uniform the color bar scale.
            ticks= np.array([1.e-4, 1.e-3, 1.e-2,1.e-1,0, 10])
        elif i == 3:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap='bwr', vmin=-6, vmax=6)
            ax_l[i].get_yaxis().set_visible(False)
            color = 'black'
            cb_i = f.colorbar(im_i, ax=ax_l[i], shrink=0.48, pad=0.01,  orientation="horizontal", 
                              aspect=15)
        else:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap, norm=LogNorm(), clim=clim)
            ax_l[i].get_yaxis().set_visible(False)
        ax_l[i].get_xaxis().set_visible(False)
        deltaPix = deltaPix_list[i]
        scale_bar(ax_l[i], frame_size, dist=0.5/deltaPix, text='0.5"', color = color)
        # if i == 0:
        #     ax_l[0].text(frame_size*0.25, frame_size*0.09, "~2.8kpc",fontsize=20, 
        #                  color='white')
        
        # if not (idx == 0 and filt == 'F150W'):
        #     ax_l[2].text(frame_size*0.40, frame_size*0.05, "Host", weight='bold', fontsize=20, 
        #                  color='white')
        fontsize = 20
        if i <3:
            cb_i = f.colorbar(im_i, ax=ax_l[i], shrink=0.48, pad=0.01,  orientation="horizontal", 
                              aspect=15, ticks=ticks)
        if i <2:
            ax_l[i].text(frame_size*0.05, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                         weight='bold', color='white')
        if i ==2:
            ax_l[i].text(frame_size*0.02, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                         weight='bold', color='white')
        if i == 3:
            ax_l[i].text(frame_size*0.02, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                         weight='bold', color='black')
        
        if z != None:
            ax_l[i].text(frame_size*0.95, frame_size*0.05, 'z={0}'.format(z),fontsize=fontsize, 
                         ha='right', va='bottom',
                         weight='bold', color='white', bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
        
        if i == 0:
            ax_l[i].text(frame_size*0.95, frame_size*0.95, target_ID,fontsize=fontsize, 
                         ha='right', va='top',
                         weight='bold', color='white',  bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
   
    ax_l[2].scatter(ps_x, ps_y, marker= 'x', c='lightblue', s = 80)
    # if not (filt =='F150W' and idx ==0):
    #     ax_l[2].scatter(host_x, host_y, marker= '+', c='lightblue', s= 100)
    # from matplotlib.patches import Rectangle
    # ax_l[2].add_patch( Rectangle((0, 0), frame_size, frame_size, fc='none', ec='r', lw=10) )
    
    # for axis in ['top', 'bottom', 'left', 'right']:
    #     ax_l[2].spines[axis].set_linewidth(5.5)  # change width
    #     ax_l[2].spines[axis].set_color('red')    # change color

    ax_l[0].set_ylabel(filt, fontsize=20, weight='bold')
    plt.subplots_adjust(wspace=-0.5, hspace=0)
    plt.show()       
    return f

idx = 0
# filters = ['F150W', 'F356W']
filters = ['F150W']
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
z_str = str(z)

import copy, matplotlib
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]
        
        if filt == 'F150W' :
            cmap = 'inferno'
        else:
            cmap = 'gist_heat'
        my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)) # copy the default cmap
        my_cmap.set_bad('black')

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        if idx !=1:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        elif idx ==1:
            fit_files = glob.glob(run_folder+'*fit_material*/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        # fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        # fit_run.plot_final_qso_fit(target_ID = filt, save_plot = True, cmap = my_cmap)
fit_run.cal_astrometry()
ps_x, ps_y = np.array(fit_run.final_result_ps[0]['position_xy']) + len(fit_run.image_host_list[0])/2
host_x, host_y = np.array(fit_run.final_result_galaxy[0]['position_xy']) + len(fit_run.image_host_list[0])/2

label = list(fit_run.flux_2d_out.keys())
label.sort()
label[1],label[2] = label[2], label[1]
image_list = [fit_run.flux_2d_out[label[i]] for i in range(len(label)) ]
label[2] = "data$-$point source"
fig = total_compare(image_list, label, [fit_run.fitting_specify_class.deltaPix]*4, 
                    target_ID=None, z=None,)
fig.savefig('figures/{1}_{0}_qso_final_plot.pdf'.format(filt,target_id))
