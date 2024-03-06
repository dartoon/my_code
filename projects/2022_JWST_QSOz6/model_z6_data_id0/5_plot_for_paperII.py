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
# import matplotlib as mat
# mat.rcParams['font.family'] = "sans-serif"

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
    plt.rcParams["font.family"] = "sans-serif"
    import matplotlib
    matplotlib.rcParams['pdf.fonttype'] = 42
    
    ax_l = [plt.subplot2grid((6,cl_num), (0,i), rowspan=6) for i in range(len(flux_list_2d))] #The image plot
    # if cmap == None:
    #     cmap = my_cmap
    for i in range(len(flux_list_2d)):
        frame_size = len(flux_list_2d[i])
        color = 'white'
        if i >1:
            flux_list_2d[i] = flux_list_2d[i]
        if i == 0:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap_list[i], norm=LogNorm(vmax = flux_list_2d[0].max(), vmin = 1.e-4))
            clim=im_i.properties()['clim'] #To uniform the color bar scale.
            ticks= np.array([1.e-4, 1.e-3, 1.e-2,1.e-1,0])
        # elif i == 3:
        #     im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap='bwr', vmin=-6, vmax=6)
        #     ax_l[i].get_yaxis().set_visible(False)
        #     color = 'black'
        #     cb_i = f.colorbar(im_i, ax=ax_l[i], shrink=0.48, pad=0.01,  orientation="horizontal", 
        #                       aspect=15)
        else:
            im_i = ax_l[i].imshow(flux_list_2d[i],origin='lower',cmap=cmap_list[i], norm=LogNorm(), clim=clim)
        ax_l[i].get_yaxis().set_visible(False)
        ax_l[i].get_xaxis().set_visible(False)
        deltaPix = deltaPix_list[i]
        if i == 0 or i == 2:
            scale_bar(ax_l[i], frame_size, dist=0.5/deltaPix, text='0.5"', color = color)
        # if i == 0:
        #     ax_l[0].text(frame_size*0.25, frame_size*0.09, "~2.8kpc",fontsize=20, 
        #                  color='white')
        
        # if not (idx == 0 and filt == 'F150W'):
        #     ax_l[2].text(frame_size*0.40, frame_size*0.05, "Host",  fontsize=20, 
        #                  color='white')
        
        # if filt != 'F150W':
            # theta_dict = {0: 137.475, 1: 139.818, }
            from astropy.coordinates import Angle
            # theta = Angle(theta_dict[idx], 'deg')
            if target_info[str(idx)]['theta'] != None:
                theta = Angle(target_info[str(idx)]['theta'], 'deg')
                f_center = len(flux_list_2d[0])/2
                w = 0.2 / deltaPix_list[0]
                # h = 0.6 / deltaPix_list[0]
                h = 3.2 / deltaPix_list[0]
                from photutils.aperture import RectangularAperture
                aper = RectangularAperture((ps_x, ps_y), w, h, theta=theta)
                aper.plot(color='white',
                          lw=0.4,axes=ax_l[0])
                h = 0.6 / deltaPix_list[0]
                from photutils.aperture import RectangularAperture
                aper = RectangularAperture((ps_x, ps_y), w, h, theta=theta)
                aper.plot(color='white',
                          lw=1,axes=ax_l[0])
                # axins.add_patch(aper)
            
        fontsize = 18
        cb_i = f.colorbar(im_i, ax=ax_l[i], shrink=0.48, pad=0.01,  orientation="horizontal", 
                          aspect=15, ticks=ticks)
        cb_i.ax.tick_params(labelsize=15) 
        if i <2:
            ax_l[i].text(frame_size*0.05, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                          color='white')
        if i ==2:
            ax_l[i].text(frame_size*0.02, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                          color='white')
        if i == 3:
            ax_l[i].text(frame_size*0.02, frame_size*0.87, label_list_2d[i],fontsize=fontsize, 
                          color='white')
        
        if z != None:
            ax_l[i].text(frame_size*0.95, frame_size*0.05, 'z={0}'.format(z),fontsize=fontsize, 
                         ha='right', va='bottom',
                          color='white', bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
        
        if i == 0:
            ax_l[i].text(frame_size*0.95, frame_size*0.95, target_ID,fontsize=fontsize, 
                         ha='right', va='top',
                          color='white',  bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
        
        if i == 0 or i ==2:
            ax_l[i].text(frame_size*0.95, frame_size*0.95, filt_list[i],fontsize=fontsize, 
                         ha='right', va='top',
                          color='white',  bbox={'facecolor': 'gold', 'alpha': 0.6, 'pad': 3})
   
    # ax_l[2].scatter(ps_x, ps_y, marker= 'x', c='lightblue', s = 80)
    
    ax_l[0].set_xlim([0, len(flux_list_2d[0])])
    ax_l[0].set_ylim([0, len(flux_list_2d[0])])
    
    ax_l[0].yaxis.set_tick_params(labelsize=15)
    # if not (filt =='F150W' and idx ==0):
    #     ax_l[2].scatter(host_x, host_y, marker= '+', c='lightblue', s= 100)
    # from matplotlib.patches import Rectangle
    # ax_l[2].add_patch( Rectangle((0, 0), frame_size, frame_size, fc='none', ec='r', lw=10) )
    
    # if filt != 'F150W':
    #     plt.title(target_id, fontsize=20,loc='left', x=-1.5,y=1.05)
    
    plt.subplots_adjust(wspace=-0.5, hspace=0)
    plt.show()       
    return f

idx = 1
# filters = ['F150W', 'F356W']
filters = ['F356W', 'F150W']
from target_info import target_info
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx) #!!!
z_str = str(z)

image_list = []
deltaPix_list = []
cmap_list = []
filt_list = []
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
        cmap_list.append(my_cmap)
        cmap_list.append(my_cmap)

        PSF_lib_files = glob.glob(run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        if idx !=1:
            if filt == 'F356W':
                fit_files = glob.glob(run_folder+'*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
            if filt == 'F150W':
                # fit_files = glob.glob(run_folder+'*fit_material_super2/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
                fit_files = glob.glob(run_folder+'*fit_material/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        elif idx ==1:
            if filt == 'F356W':
                fit_files = glob.glob(run_folder+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
            if filt == 'F150W':
                fit_files = glob.glob(run_folder+'*fit_material_super2/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
                # fit_files = glob.glob(run_folder+'*fit_material/fit_run*_fixn1_*idx{0}_{1}_*.pkl'.format(idx, filt))#+\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i],'rb')))
        chisqs = np.array([fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()  
        # print('idx', idx, filt, "Total PSF NO.", 'chisq',chisqs[sort_Chisq[top_psf_id]], len(sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        # fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        # fit_run.plot_final_qso_fit(target_ID = filt, save_plot = True, cmap = my_cmap)
        # folder = '../NIRCam_data/*/bkg_removed/SHELLQs_J0844m0132_F150W_Dec7_i2d_rmbkg.fits' 
        data_files = glob.glob('../NIRCam_data/*/bkg_removed/SHELLQs_{0}*_{1}_*_i2d_rmbkg.fits'.format(target_id[:5], filt))
        if len(data_files) > 1:
            if idx == 1:
                data_files = [data_files[1]]
            else:
                raise ValueError("error!")
        data_file = data_files[0]
        _fitsFile = pyfits.open(data_file)
        fov_image_org = _fitsFile[1].data # check the back grounp
        # fit_run.fitting_specify_class.data_process_class.fov_image = fov_image
        fov_image = copy.deepcopy(fov_image_org)
        fit_run.targets_subtraction(sub_qso_list=[0], org_fov_data=fov_image)
        
        radius = 63
        if filt == 'F150W':
            radius = radius * 2
            
        from galight.tools.cutout_tools import cutout
        image_list.append(cutout(image = fov_image_org, center = fit_run.fitting_specify_class.data_process_class.target_pos, 
                                 radius=radius))
        image_list.append(cutout(image = fit_run.fov_image_targets_sub, center = fit_run.fitting_specify_class.data_process_class.target_pos, 
                                 radius=radius))
        
        # image_list.append(fit_run.flux_2d_out['data'])
        # image_list.append(fit_run.flux_2d_out['data-point source'])
        
        deltaPix_list.append(fit_run.fitting_specify_class.deltaPix)
        deltaPix_list.append(fit_run.fitting_specify_class.deltaPix)
        #%%
        if count == 0:
            fit_run.cal_astrometry()
            ps_x, ps_y = np.array(fit_run.final_result_ps[0]['position_xy']) + len(fit_run.image_host_list[0])/2
            host_x, host_y = np.array(fit_run.final_result_galaxy[0]['position_xy']) + len(fit_run.image_host_list[0])/2
        
            host_x = host_x + (radius - (len(fit_run.flux_2d_out['data'])/2))
            host_y = host_y + (radius - (len(fit_run.flux_2d_out['data'])/2))
            ps_x = ps_x + (radius - (len(fit_run.flux_2d_out['data'])/2))
            ps_y = ps_y + (radius - (len(fit_run.flux_2d_out['data'])/2))
        
        filt_list.append(filt)
        filt_list.append(filt)
        # label = list(fit_run.flux_2d_out.keys())
        # label.sort()
        # label[1],label[2] = label[2], label[1]
        # image_list = [fit_run.flux_2d_out[label[i]] for i in range(len(label)) ]
        # label[2] = "data$-$point source"
label = ['data', "data$-$point source", 'data', "data$-$point source"]
fig = total_compare(image_list, label, deltaPix_list, 
                    target_ID=None, z=None,cmap =cmap_list)
fig.savefig('/Users/Dartoon/Downloads/{0}_result.pdf'.format(target_id))
print(target_id) 
        
    
