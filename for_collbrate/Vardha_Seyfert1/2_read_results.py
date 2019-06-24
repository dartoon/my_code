#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:48:10 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy

ID = 'l106'
picklename = ID+'_fitting.pkl'
result = pickle.load(open(picklename,'rb'))
[best_fit,pso_fit,mcmc_fit, trans_paras] = result

source_result, image_host, ps_result, image_ps, _ =best_fit
chain_list, param_list, _ = pso_fit
samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
_, _, mcmc_new_list, labels_new, _ = trans_paras

print "best-fit source_result:", source_result
print "best-fit ps_result:", ps_result
pix_sz = 0.04

#%%Replot the profiles
import sys
sys.path.insert(0,'./fitting_tools/')
from flux_profile import total_compare
zp = 25.0985

flux_hdulist = pyfits.open(ID+'_flux_list.fits')  # includes [agn_image,image_ps, extended_source, error_map, QSO_msk]
flux_list = [flux_hdulist[i].data for i in range(4)]
QSO_msk = flux_hdulist[4].data
label = ['data', 'QSO', 'extended sources', 'model', 'normalized residual']
fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz, zp = zp,
                    plot_compare = False, msk_image = QSO_msk)
plt.show()

print "plot individual frame:"
plt.imshow(flux_list[0] - flux_list[1], norm=LogNorm(),origin='low')
plt.colorbar()
frame_size = len(flux_list[0])
plt.text(frame_size*0.05, frame_size*0.9, 'data - Point Source', weight='bold',
         fontsize=17, color='white')
plt.show()

#%%diagnose the PSO chain convergency
import lenstronomy.Plots.output_plots as out_plot
for i in range(len(chain_list)):
    if len(param_list[i]) > 0:
        f, axes = out_plot.plot_chain(chain_list[i], param_list[i])

#%%test the MCMC chain convergency
#        
#import lenstronomy.Plots.output_plots as plot_mcmc_behaviour
fig = plt.figure(figsize=(20, 15))
ax = fig.add_subplot(111)
out_plot.plot_mcmc_behaviour(ax, samples_mcmc, param_mcmc, dist_mcmc)       

#%%Plot the transferred corner plots
# here the (non-converged) MCMC chain of the non-linear parameters
if not samples_mcmc == []:
#    plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
#    plt.show()
    plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
    plt.show()

#%% To readout the parameter

#Read the fitting parameter.
idx = 0
v_l=np.percentile(samples_mcmc[:,idx],16,axis=0)
v_m=np.percentile(samples_mcmc[:,idx],50,axis=0)
v_h=np.percentile(samples_mcmc[:,idx],84,axis=0)
print param_mcmc[idx], ":", v_l, v_m, v_h

#For the translated totol flux.
mcmc_new_list = np.asarray(mcmc_new_list)
idx = 2
v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
print labels_new[idx], ":", v_l, v_m, v_h

##%%To save inividual fits file (Take host 0 as example):
#host_id = 0
#center_QSO = np.array([2053, 2475])   #!!! The is the array of the QSO position that introduced in 0_cutout.py
#fitsFile = pyfits.open(ID+"_sci.fits")
#file_header = copy.deepcopy(fitsFile[0].header)
#file_header['CRPIX1'] = file_header['CRPIX1']-center_QSO[0]+len(flux_list[0])/2
#file_header['CRPIX2'] = file_header['CRPIX2']-center_QSO[1]+len(flux_list[0])/2
#pyfits.PrimaryHDU(image_host[host_id],header=file_header).writeto(ID+'fitted_host{0}.fits'.format(host_id),overwrite=True)
##thdu_fluxlist.writeto(ID+'_flux_list.fits', overwrite=True)

