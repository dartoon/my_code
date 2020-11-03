#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import copy

from lenstronomy.Util import constants as const
import lenstronomy.Util.param_util as param_util

#file name:
filt='f160w'

#==============================================================================
# # import main simulation class of lenstronomy
#==============================================================================
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import ImageData as Data
from lenstronomy.Data.psf import PSF

# data specifics
numPix = 241  #  pixel size  #!!!
deltaPix = 0.13/4 #  pixel size in arcsec (area per pixel = deltaPix**2)

kwargs_data_high_res = sim_util.data_configure_simple(numPix, deltaPix) #,inverse=True)
data_class = Data(**kwargs_data_high_res)
kwargs_numerics = {'supersampling_factor': 2, 'supersampling_convolution': False}

import sys
sys.path.insert(0,'../../share_tools/')
from gene_para import gene_para
for seed in range(700, 701):
    print(seed)
    para=gene_para(seed=seed,fixh0=102)
    #==============================================================================
    # ##### lens mass model 
    #==============================================================================
    from lenstronomy.LensModel.lens_model import LensModel
    kwargs_spemd = para.spemd()
    
    kwargs_spemd['q'] = np.random.uniform(0.6,0.9)
    kwargs_spemd['e1'], kwargs_spemd['e2'] = param_util.phi_q2_ellipticity(phi=kwargs_spemd['phi_G'], q=kwargs_spemd['q'])
    lens_model_list = ['PEMD','SHEAR']
    #kwargs_spemd['gamma'] = 2.
    kwargs_mass_copy = copy.deepcopy([kwargs_spemd])
    del kwargs_spemd['phi_G']
    del kwargs_spemd['q']    
    ex_shear = {'gamma1': para.shear()[0]['e1'], 'gamma2': para.shear()[0]['e2']}
    kwargs_lens_list = [kwargs_spemd, ex_shear]
    lens_model_class = LensModel(lens_model_list)
    
    # #==============================================================================
    # # #########source light
    # #==============================================================================
    source_pos=[-0.092, -0.0212]  #Seed 220  cored-Powerlaw
    #==============================================================================
    # #The information of the QSO as PSF
    #==============================================================================
    from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
    lensEquationSolver = LensEquationSolver(lens_model_class)
    x_image, y_image = lensEquationSolver.findBrightImage(source_pos[0], source_pos[1], kwargs_lens_list, numImages=4,
                                                          min_distance=deltaPix, search_window=numPix * deltaPix)
    # #==============================================================================
    # # Creat a folder save the fits file
    # #==============================================================================
    # import os
    # if os.path.exists(sim_folder_name)==True:
    #     import shutil 
    #     shutil.rmtree(sim_folder_name)
    # os.mkdir(sim_folder_name)
    #==============================================================================
    # Time delay
    #==============================================================================
    TD_distance=(1+para.z_lens)*para.D_l*para.D_s/para.D_ls
    fermat_po=lens_model_class.fermat_potential(x_image, y_image, x_source=source_pos[0], y_source=source_pos[1], kwargs_lens=kwargs_lens_list)
    pix_ra, pix_dec= x_image/deltaPix+numPix/2., y_image/deltaPix+numPix/2.
    abc_list = ['A', 'B', 'C', 'D']
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    for i in range(len(pix_ra)):
    	x_, y_ =pix_ra[i], pix_dec[i]
    	ax.plot(x_, y_, 'or')
    	ax.text(x_, y_, abc_list[i], fontsize=20, color='k') 
    ax.plot(numPix/2., numPix/2., 'ob')
    # ax.text(numPix/2., numPix/2., 'lens', fontsize=20, color='b') 
    plt.xlim(0, numPix)
    plt.ylim(0, numPix)
    # fig.savefig(sim_folder_name+'/ABCD.png')
    plt.show()  
    TD = TD_distance/const.c * fermat_po / const.day_s * const.arcsec**2 * const.Mpc
    print("fermat_po", fermat_po)
    print("Time delay of ABCD - A :\n\t", TD-TD[0])
    
    ###################Active this part if need to see the caustic and critical line
    from  lenstronomy.Plots import lens_plot
    f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=False, sharey=False)
    lens_plot.lens_model_plot(ax, lensModel=lens_model_class, kwargs_lens=kwargs_lens_list,
                              sourcePos_x=source_pos[0], sourcePos_y=source_pos[1],
                              point_source=True, with_caustics=True)
    f.show()

##%%
###################Plot the kappa profile
#import lenstronomy.Util.util as util
#import sys 
#sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
#from flux_profile import SB_profile
#sys.path.insert(0,'../../../submission_analysis/comparing_model/')
#from share_tools import cal_gamma
#
#x_grid, y_grid = data_class.pixel_coordinates
#x_grid1d = util.image2array(x_grid)
#y_grid1d = util.image2array(y_grid)
#kappa = lens_model_class.kappa(x_grid1d, y_grid1d, kwargs_lens_list)
#kappa = util.array2image(kappa)

##%%
#gridspace, y_log, if_annuli='log', True, True
#kappa_value, x = SB_profile(kappa, center = [len(kappa)/2]*2, radius=len(kappa)/2, start_p=3.5, grids=70,
#                             gridspace=gridspace,if_annuli=if_annuli, fits_plot=False)
#fig, ax = plt.subplots(figsize=(10,7))
#plt.plot(x*deltaPix, kappa_value, label='Kappa profile', c='blue')
#s, Rein_r = cal_gamma(kappa)
#
#print("effective R_ein", round(Rein_r*deltaPix,3))
#if len(kwargs_lens_list) ==2:
#    kappa_decomp_list = []
#    for i in range(len(kwargs_lens_list)):
#        lens_model_class_sub = LensModel([lens_model_list[i]])
#        kappa_i = lens_model_class_sub.kappa(x_grid1d, y_grid1d, [kwargs_lens_list[i]])    
#        kappa_i = util.array2image(kappa_i)
#        kappa_value_i, x = SB_profile(kappa_i, center = [len(kappa_i)/2]*2, radius=len(kappa_i)/2, start_p=3.5, grids=70,
#                                     gridspace=gridspace,if_annuli=if_annuli, fits_plot=False)        
#        plt.plot(x*deltaPix, kappa_value_i, label="{0} profile".format(lens_model_list[i]))
#x_p = np.linspace(-1,2)*0 + Rein_r *deltaPix
#y_reg = np.linspace(10**(-0.6), (kappa_value.max() )*1.2)
#plt.plot(x_p, y_reg, c='red',label='Einstein Radius')
#plt.legend(fontsize=15, loc=3)
##legend1 = plt.legend([l1[0]], ["Einstein Radius"], loc=3, fontsize=20)
#plt.tick_params(labelsize=25)
#ax.set_ylabel("profiles of Convergence map (log10 space)", fontsize=20)
#ax.set_xlabel("Radius (arcsec)",fontsize=20)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
