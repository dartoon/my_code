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
import glob
import astropy.io.fits as pyfits
from gen_fit_id import gen_fit_id

#image_folder = '../images_directory/'
image_folder = './image/'

folder = 'RESULTS/'

galaxy_ID = '18184' 
galaxy_RA = 150.121811
galaxy_DEC = 2.394572

band = 'I'
zp = 27.0

picklename = folder+'fit_image_'+galaxy_ID+ '_HSC-{0}.pkl'.format(band)
if_file = glob.glob(picklename)

if if_file == []:
    print "MCMC for {0}-{1} not performed!".format(galaxy_ID,band)
else:
    result = pickle.load(open(picklename,'rb'))

    best_fit, chain_result, trans_paras   = result         
    source_result, image_host, _ =best_fit
    chain_list, _ = chain_result
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]    
    
    _, mcmc_new_list, labels_new, _ = trans_paras            
    mcmc_new_list = np.asarray(mcmc_new_list)    
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
    #print labels_new[idx], ":", v_l, v_m, v_h
    print "The inferred", labels_new[idx], "mag:"
    print "lower limit:", -2.5 * np.log10(v_h) + zp
    print "The mid fit:", -2.5 * np.log10(v_m) + zp
    print "upper limit", -2.5 * np.log10(v_l) + zp

#%%
print "Check the convergency of the PSO chains:"
import lenstronomy.Plots.output_plots as out_plot
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
band_seq = ['G', 'R', 'I', 'Z', 'Y']
filename = galaxy_ID+'_HSC-{0}.fits'.format(band)
sub_bkg = True
from flux_profile import galaxy_total_compare
galaxy_im, err_map, galaxy_bkg, PSF, pix_scale, zp, qso_fr_center, fr_c_RA_DEC = gen_fit_id(image_folder, galaxy_RA, galaxy_DEC, filename, cut_frame=120, subbkl=sub_bkg)
ct = (len(galaxy_im)-len(image_host[0]))/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

if len(image_host) == 1:
    host = image_host[0]
    label = ['data', 'QSO', 'host', 'model', 'normalized residual']
elif len(image_host) >1:
    host = np.zeros_like(image_host[0])
    for i in range(len(image_host)):
        host += image_host[i]
    label = ['data', 'QSO', 'host as {0} components'.format(i+1), 'model', 'normalized residual']  #Print the numbers of objects
agn_image = galaxy_im[ct:-ct,ct:-ct]
error_map = err_map[ct:-ct,ct:-ct]

objs_img = np.zeros_like(image_host[0])
if len(source_result)>1:
    for i in range(1,len(image_host)):
        objs_img += image_host[i]
flux_list = [agn_image, image_host[0]*0, image_host[0]+objs_img, error_map]

fig = galaxy_total_compare(label_list = ['data', 'model', 'normalized residual'], flux_list = flux_list, target_ID = galaxy_ID, pix_sz=pix_scale,
                  zp = zp, plot_compare=True, msk_image = np.ones_like(image_host[0]))
   

