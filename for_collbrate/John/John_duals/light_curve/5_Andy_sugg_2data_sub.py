#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:53:18 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from matplotlib.colors import LogNorm
import copy, matplotlib
from galight.tools.measure_tools import measure_FWHM
# name = 'PSPS_fixpos_result'
name = 'PSPS+sersic_fixpos_nearbyPSF_result'
import imageio

my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')

# from galight.tools.astro_tools import plt_fits
bands = 'GRIZY'
# band = 'G'
for band in bands:
    files = glob.glob(name+"*band{0}*pkl".format(band))
    files.sort()
    filename_list = []
    ct = 0
    for file in files:
        date = file.split('result-')[1][:10]
        filename = file
        fit_run = pickle.load(open(filename,'rb'))
        # fit_run.plot_final_qso_fit()
        ps_result = fit_run.final_result_ps
        if ct == 0:
            sub_data = fit_run.fitting_specify_class.kwargs_data['image_data']
        else:
            mag0, mag1 = ps_result[0]['magnitude'], ps_result[1]['magnitude']
            
            if ps_result[0]['ra_image']<0:
                print(file)
                print('Warning, miss PS0 PS1', ps_result[0]['ra_image'], ps_result[1]['ra_image'])
                fit_run.plot_final_qso_fit()
                fit_run.fitting_specify_class.plot_fitting_sets()
            # date = date[2:4]+date[5:7]+date[8:10]
            
            fig, ax = plt.subplots(figsize=None)
            # norm = LogNorm()#np.max(img[~np.isnan(img)]))
            plt.imshow(fit_run.fitting_specify_class.kwargs_data['image_data'] - sub_data, cmap=my_cmap, 
                       norm=LogNorm(vmax=100, vmin = 0.001), origin='lower') 
            # plt.colorbar()
            # plt.xticks([])
            # plt.yticks([])
            plt.text(1,1,date + '-'+band,
                      {'color': 'white', 'fontsize':25})
            # x0 = len(fit_run.fitting_specify_class.kwargs_data['image_data'])/2 -ps_result[0]['ra_image']/fit_run.fitting_specify_class.deltaPix
            # y0 = len(fit_run.fitting_specify_class.kwargs_data['image_data'])/2 +ps_result[0]['dec_image']/fit_run.fitting_specify_class.deltaPix
            # plt.text(x0-3,y0+2,round(mag0,2),
            #          {'color': 'black', 'fontsize':12})
            
            # x1 = len(fit_run.fitting_specify_class.kwargs_data['image_data'])/2 -ps_result[1]['ra_image']/fit_run.fitting_specify_class.deltaPix
            # y1 = len(fit_run.fitting_specify_class.kwargs_data['image_data'])/2 +ps_result[1]['dec_image']/fit_run.fitting_specify_class.deltaPix
            # plt.text(x1-3,y1+2,round(mag1,2),
            #          {'color': 'black', 'fontsize':12})
            
            # FWHM = np.mean(measure_FWHM(fit_run.fitting_specify_class.psf_class.kernel_pixel) )
            # plt.text(1,55,'PSF FWHM =' + str(round(FWHM,1))+' pix',
            #          {'color': 'white', 'fontsize':14})
            filename = 'figs/'+date+band+'_sub.png'
            plt.savefig(filename)
            filename_list.append(filename)
            plt.close()     
        ct = ct + 1

    with imageio.get_writer('figs/'+band+'_days_sub.gif', mode='I', fps=1) as writer:
        for filename in filename_list:
            image = imageio.imread(filename)
            writer.append_data(image)