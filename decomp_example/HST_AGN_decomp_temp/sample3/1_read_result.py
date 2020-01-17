#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 17:55:30 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import astropy.io.fits as pyfits

import sys
sys.path.insert(0,'../../../py_tools/')

picklename = 'example.pkl'
result = pickle.load(open(picklename,'rb'))

#[best_fit,chain_r, trans_paras] = result      
best_fit, chain_result, trans_paras   = result         
source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, _ = chain_result
sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]

#chain_list, param_list, _ = pso_fit
#samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
_, _, mcmc_new_list, labels_new, _ = trans_paras            
mcmc_new_list = np.asarray(mcmc_new_list)    
#    [source_result, image_host, ps_result, image_ps, samples_mcmc, param_mcmc, paras, chain_list, param_list] = result
#    [source_params_2, ps_param_2, mcmc_new_list, labels_new] = paras
    #print "The fixed parameters in galaxy:", source_params_2
    #%%
#    print "plot the overall parameter contour:"
#    plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
#    plt.show()
    
print "plot the flux contour for all the hosts"
plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
plt.show()

#readout the parameter
print "=================================\noverall fitting parameters:", param_mcmc
idx = 0     #The inferred Reff of the host
v_l=np.percentile(samples_mcmc[:,idx],16,axis=0)
v_m=np.percentile(samples_mcmc[:,idx],50,axis=0)
v_h=np.percentile(samples_mcmc[:,idx],84,axis=0)
print "The inferred", param_mcmc[idx], ":"
print "lower limit:", v_l
print "The mid fit:", v_m
print "upper limit", v_h

##Readout the translated flux.
print "=================================\noverall flux components:", labels_new
mcmc_new_list = np.asarray(mcmc_new_list)
idx = 1     #The translated flux for the host
v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
zp=26.4524
#print labels_new[idx], ":", v_l, v_m, v_h
print "The inferred", labels_new[idx], "mag:"
print "lower limit:", -2.5 * np.log10(v_h) + zp
print "The mid fit:", -2.5 * np.log10(v_m) + zp
print "upper limit", -2.5 * np.log10(v_l) + zp

#%%
print "Check the convergency of the PSO chains:"
import lenstronomy.Plots.chain_plot as out_plot
for i in range(len(chain_list)):
    f, axes = out_plot.plot_chain_list(chain_list,0)
plt.show()

#%%test the MCMC chain convergency
#        
#import lenstronomy.Plots.output_plots as plot_mcmc_behaviour
fig = plt.figure(figsize=(20, 15))
ax = fig.add_subplot(111)
out_plot.plot_mcmc_behaviour(ax, samples_mcmc, param_mcmc, dist_mcmc)       
plt.show()

#%% Plot the image again:
pix_scale = 0.0642 
from flux_profile import total_compare
psf, QSO_img, QSO_std = pyfits.getdata('psf_16.fits'),  pyfits.getdata('CID454_cutout.fits'),  pyfits.getdata('wht_err.fits')

frame_size = 61
#frame = '{0}'.format(frame_size)
QSO_fm = len(QSO_img)
ct = (QSO_fm-frame_size)/2     # If want to cut to 61, i.e. (121-61)/2=30

QSO_img = QSO_img[ct:-ct,ct:-ct]
QSO_std = QSO_std[ct:-ct,ct:-ct]
psf = psf[ct:-ct,ct:-ct]   
if len(image_host) == 1:
    host = image_host[0]
    label = ['data', 'QSO', 'host', 'model', 'normalized residual']
elif len(image_host) >1:
    host = np.zeros_like(image_host[0])
    for i in range(len(image_host)):
        host += image_host[i]
    label = ['data', 'QSO', 'host as {0} components'.format(i+1), 'model', 'normalized residual']  #Print the numbers of objects
#error_map = err_map[ct:-ct,ct:-ct]
flux_list = [QSO_img, image_ps, host, QSO_std]

total_compare(label_list = label, flux_list = flux_list, target_ID = 'example', pix_sz=pix_scale, zp = zp,
                    plot_compare = False, msk_image = np.ones_like(QSO_img))

#fig.savefig("{0}_SB_profile.pdf".format(name_save), bbox_inches = 'tight')
from transfer_to_result import transfer_to_result
result = transfer_to_result(data=QSO_img, pix_sz = pix_scale,
                            source_result=source_result, ps_result=ps_result, image_ps=image_ps, image_host=image_host, error_map=QSO_std,
                            zp=zp, fixcenter=False,ID='Example', plot_compare = False)
