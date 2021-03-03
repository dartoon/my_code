#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:58:39 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
#from gen_fit_id import gen_fit_id
from photutils import make_source_mask
import os
#plt.ion()

import sys
sys.path.insert(0,'../../../py_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result
from mask_objects import detect_obj
from flux_profile import profiles_compare, flux_profile
from matplotlib.colors import LogNorm
import pickle
#import copy
#import time

#singel_exp = 625.0
deep_seed = False
pltshow = 1 #Note that setting plt.ion() in line27, the plot won't show anymore if running in terminal.
fixcenter = True
run_MCMC = False
use_true_PSF = 'same_drizzle' #'same_drizzle' #False
save_SNR = True
zp_dic = {'F160W':25.9463, 'F105W':26.2687} 

# seed_l, seed_h = [[201, 208], [208,215], [215, 221]][0]
seed_l, seed_h = [[0,1],[0, 25], [25, 50], [50, 75], [75, 99]][0]

for seed in range(seed_l, seed_h):
    for ID in [0, 1]:#range(1, 7):
        for filt in ['F160W', 'F105W']:
            folder = filt + '/sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed)
            if use_true_PSF == True:
                f = open(folder+"/sim_info.txt","r")
                string = f.read()
                lines = string.split('\n')   # Split in to \n            
                psf_id = int([lines[i] for i in range(len(lines)) if 'PSF' in lines[i]][0].split('\t')[1])
                save_name = 'fit_result_truePSF'
                psf = pyfits.getdata('webPSF/drizzle_PSF_{0}/Drz_PSF_id{1}.fits'.format(filt, psf_id))
            elif use_true_PSF == False:
                psf_id = 0
                save_name = 'fit_result_diffPSF'
                psf = pyfits.getdata('webPSF/drizzle_PSF_{0}/Drz_PSF_id{1}.fits'.format(filt, psf_id))
            elif use_true_PSF == 'same_drizzle':
                psf = pyfits.getdata(folder+'/Drz_POINTclean_image.fits')
                psf = psf/psf.sum()
                save_name = 'fit_result_samedrz'
            print("Fitting "+folder)                
            if glob.glob(folder+'/{0}.pkl'.format(save_name)) !=[]:
                continue
            else:
                #Setting the fitting condition:
                pix_org_scale = 0.13 #After dirzzled
                pix_scale = 0.08 #After dirzzled
                zp= zp_dic[filt]
                framesize = 61
                QSO_img_org = pyfits.getdata(folder+'/Drz_QSO_image.fits')
                rms_org = pyfits.getdata(folder+'/noise_map.fits')                    
                half_r = int(framesize/2)
                peak = np.where(QSO_img_org==QSO_img_org.max())
                peak = [peak[0][0], peak[1][0]]                
                QSO_img = QSO_img_org[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]                
                QSO_std = rms_org[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
                plt.imshow(QSO_img, origin='lower', cmap='gist_heat', norm=LogNorm(), vmin = 0.001, vmax = QSO_img.max())
                plt.colorbar()
                plt.show()   
                if framesize < len(psf)-10:
                    psf_framesize = len(QSO_img)
                    psf_half_r = int(psf_framesize/2)
                    psf_peak = np.where(psf==psf.max())
                    psf_peak = [psf_peak[0][0], psf_peak[1][0]]
                    psf = psf[psf_peak[0]-psf_half_r:psf_peak[0]+psf_half_r+1,psf_peak[1]-psf_half_r:psf_peak[1]+psf_half_r+1]
                fixed_source = []
                kwargs_source_init = []
                kwargs_source_sigma = []
                kwargs_lower_source = []
                kwargs_upper_source = []
                fixed_source.append({})  
                kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
                kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
                kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
                kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
                source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
                #
                #==============================================================================
                # to fit and save the inference
                #==============================================================================
                tag = folder+'/'+ save_name
                source_result, ps_result, image_ps, image_host, error_map=fit_qso(QSO_img, psf_ave=psf, psf_std = None,
                                                                                  source_params=source_params, fixcenter=fixcenter,
                                                                                  pix_sz = pix_scale, no_MCMC = (run_MCMC==False),
                                                                                  QSO_std =QSO_std, tag=tag, deep_seed= deep_seed, pltshow=pltshow,
                                                                                  corner_plot=False, flux_ratio_plot=True, dump_result=run_MCMC, pso_diag =True)
                if pltshow == 0:
                    plot_compare=False
                    fits_plot =False
                else:
                    plot_compare=True
                    fits_plot =True
                result = transfer_to_result(data=QSO_img, pix_sz = pix_scale,
                                            source_result=source_result, ps_result=ps_result, image_ps=image_ps[0], image_host=image_host, error_map=error_map,
                                            zp=zp, fixcenter=fixcenter,ID='ID'+repr(ID), tag=tag, plot_compare = plot_compare)
            
                print("Saving results:")
                pickle.dump([result, framesize], open(folder+'/{0}.pkl'.format(save_name), 'wb'))
                filtname = folder+'/{0}.txt'.format(save_name)
                result_file =  open(filtname,'w') 
                result_file.write(repr(result))
                result_file.close()         
            #    #%%
            #    print("The truth:")
            #    with open(folder+"/sim_info.txt") as f: # The with keyword automatically closes the file when you are done
            #        print(f.read())
            #                         
            #    print("The inferred results:")                     
            #    print("host_flux:",  result['host_amp'])
            #    print("host Reff:",  result['R_sersic'])
            #    print("host n:",  result['n_sersic'])
            #    print("host q:",  result['q'])
            #    print("AGN flux :",  result['QSO_amp'])
    #            #%%Print SNR map:
    #        #    print("Host galaxy SNR map:")
                if save_SNR == True:
                    host_clean = pyfits.getdata(folder+'/Drz_HOSTclean_image.fits')
                    host_clean = host_clean[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
                    host_SNR = host_clean/QSO_std
                    plt.imshow(host_SNR, origin='lower')#,cmap='gist_heat', norm=LogNorm())
                    plt.colorbar()
                    plt.savefig(folder+'/SNR_host.pdf')
                    plt.close()
                #    print("Point source SNR map:")
                    point_clean = pyfits.getdata(folder+'/Drz_POINTclean_image.fits')
                    point_clean = point_clean[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
                    point_SNR = point_clean/QSO_std
                    plt.imshow(point_SNR, origin='lower')#,cmap='gist_heat', norm=LogNorm())
                    plt.colorbar()
                    plt.savefig(folder+'/SNR_point.pdf')
                    plt.close()
                #    print("AGN (total) SNR map:")
                    AGN_clean = pyfits.getdata(folder+'/Drz_AGNclean_image.fits')
                    AGN_clean = AGN_clean[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
                    Total_SNR = AGN_clean/QSO_std
                    plt.imshow(Total_SNR, origin='lower')#,cmap='gist_heat', norm=LogNorm())
                    plt.colorbar()
                    plt.savefig(folder+'/SNR_Total.pdf')
                    plt.close()
