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

# psf = pyfits.open('psf_F160W_sub4.fits'.format(filt))
# psf_data = psf[0].data
# cut =25 
# psf_data = psf_data[1+cut:-cut,1+cut:-cut]+ 0  #shave PSF to singular size, as max in the middle
# plt.imshow(psf_data, origin='lower',cmap='gist_heat', norm=LogNorm())
# plt.colorbar()
# plt.close()

# kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': deltaPix}
# psf_class = PSF(**kwargs_psf_high_res)

zp= 25.9463
kwargs_data_high_res = sim_util.data_configure_simple(numPix, deltaPix) #,inverse=True)
data_class = Data(**kwargs_data_high_res)
kwargs_numerics = {'supersampling_factor': 3, 'supersampling_convolution': False}

import sys
sys.path.insert(0,'../../projects/Sim_HST_JWST/share_tools/')
from gene_para import gene_para
# for seed in range(702, 703):
seed = 502
    # print(seed)
para=gene_para(seed=seed,fixh0=102, fixz = [1, 2])
#==============================================================================
# #######lens light 
#==============================================================================
from lenstronomy.LightModel.light_model import LightModel
lens_light_para=para.lens_light()
np.random.seed(seed)

# import magVSamp as mva
# lens_amp=mva.getAmp(SERSIC_in_mag=lens_light_para,zp=zp,deltaPix=deltaPix)
# lens_light_para['amp_sersic']=lens_amp * 963.490605#!!!
# lens_light_para['q'] = 0.6 + np.random.normal(0,0.01)

# lens_light_tran_Reff=lens_light_para['R_sersic']/np.sqrt(lens_light_para['q'])  #!!!
# #lens_light_tran_Reff = 1.1
# lens_light_para['e1'], lens_light_para['e2'] = param_util.phi_q2_ellipticity(phi=lens_light_para['phi_G'], q=lens_light_para['q'])
# lens_light_model_list = ['SERSIC_ELLIPSE']
# kwargs_lens_light = {'amp':lens_light_para['amp_sersic'], 'R_sersic': lens_light_tran_Reff, 'n_sersic': lens_light_para['n_sersic'],
# 					  'center_x': 0.0, 'center_y': 0.0, 'e1':lens_light_para['e1'], 'e2':lens_light_para['e2']}
# kwargs_lens_light_copy = copy.deepcopy(kwargs_lens_light)
# kwargs_lens_light_copy['phi_G'] =  lens_light_para['phi_G']
# kwargs_lens_light_copy['q'] =  lens_light_para['q']

# kwargs_lens_light_list = [kwargs_lens_light]
# lens_light_model_class = LightModel(light_model_list=lens_light_model_list)

#==============================================================================
# ##### lens mass model 
#==============================================================================
from lenstronomy.LensModel.lens_model import LensModel
kwargs_spemd = para.spemd()

kwargs_spemd['q'] = 0.6 + np.random.normal(0,0.01)
kwargs_spemd['e1'], kwargs_spemd['e2'] = param_util.phi_q2_ellipticity(phi=kwargs_spemd['phi_G'], q=kwargs_spemd['q'])
lens_model_list = ['PEMD','SHEAR']
#kwargs_spemd['gamma'] = 2.
kwargs_mass_copy = copy.deepcopy([kwargs_spemd])
del kwargs_spemd['phi_G']
del kwargs_spemd['q']    
ex_shear = {'gamma1': para.shear()[0]['e1'], 'gamma2': para.shear()[0]['e2']}
kwargs_lens_list = [kwargs_spemd, ex_shear]
lens_model_class = LensModel(lens_model_list)
#==============================================================================
# #########source light
#==============================================================================
count = 0
# mag_diff_list = []
# TD_diff_list = []
# source_pos_list = []
for i in range(1):
    for j in range(1):
        #%%
        # source_pos=[-0.02, 0.2]  #Seed 220  cored-Powerlaw
        # source_pos=[-0.2+0.002*i, -0.2+0.002*j]  #Seed 220  cored-Powerlaw
        source_pos = source_pos_list
        source_para=para.source_light()
        #==============================================================================
        # #The information of the QSO as PSF
        #==============================================================================
        from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
        np.random.seed(seed)
        amp_qRh_s_plane=np.random.uniform(0.8,1./0.8)
        qso_amp= 10.**(-0.4*(source_para['mag_sersic']-zp))*amp_qRh_s_plane
        
        lensEquationSolver = LensEquationSolver(lens_model_class)
        x_image, y_image = lensEquationSolver.findBrightImage(source_pos[0], source_pos[1], kwargs_lens_list, numImages=4,
                                                              min_distance=deltaPix, search_window=numPix * deltaPix)
        mag = lens_model_class.magnification(x_image, y_image, kwargs=kwargs_lens_list)
        
        #==============================================================================
        # Time delay
        #==============================================================================
        TD_distance=(1+para.z_lens)*para.D_l*para.D_s/para.D_ls
        fermat_po=lens_model_class.fermat_potential(x_image, y_image, x_source=source_pos[0], y_source=source_pos[1], kwargs_lens=kwargs_lens_list)
        pix_ra, pix_dec= x_image/deltaPix+numPix/2., y_image/deltaPix+numPix/2.
        abc_list = ['A', 'B', 'C', 'D']
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(len(pix_ra)):
        	x_, y_ =pix_ra[i], pix_dec[i]
        	ax.plot(x_, y_, 'or')
        	ax.text(x_, y_, abc_list[i], fontsize=20, color='k') 
        # 	ax.matshow(np.log10(image_highres),origin='lower')
        # fig.savefig(sim_folder_name+'/ABCD.png')
        plt.close()  
        TD = TD_distance/const.c * fermat_po / const.day_s * const.arcsec**2 * const.Mpc
        TD_ = TD[TD.argsort()]
        mag_ = mag[TD.argsort()]
        TD_diff = TD_[1:] - TD_[:-1]
        # if np.min(TD_diff) < 5/24.:
        #     idx = np.where(TD_diff == np.min(TD_diff))[0][0]
        #     mag_diff = abs(mag_[idx]/mag_[idx+1])
        #     TD_diff_list.append(TD_diff)
        #     mag_diff_list.append(mag_diff)
        #     source_pos_list.append(source_pos)
        #     if mag_diff> 100:
        #         print(source_pos)
                
        # if count%1000 == 0:
        #     print(count)
        # count = count + 1
        
        # print("Time delay of ABCD - A :\n\t", TD-TD[0])
        # print("mag:", mag)
        # if 


#%%
##################Active this part if need to see the caustic and critical line
from  lenstronomy.Plots import lens_plot
f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=False, sharey=False)
lens_plot.lens_model_plot(ax, lensModel=lens_model_class, kwargs_lens=kwargs_lens_list,
                          sourcePos_x=source_pos[0], sourcePos_y=source_pos[1],
                          point_source=True, with_caustics=True)
f.show()
