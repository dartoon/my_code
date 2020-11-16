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

psf = pyfits.open('psf_F160W_sub4.fits'.format(filt))
psf_data = psf[0].data
cut =25 
psf_data = psf_data[1+cut:-cut,1+cut:-cut]+ 0  #shave PSF to singular size, as max in the middle
plt.imshow(psf_data, origin='lower',cmap='gist_heat', norm=LogNorm())
plt.colorbar()
plt.close()

kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': deltaPix}
psf_class = PSF(**kwargs_psf_high_res)

zp= 25.9463
kwargs_data_high_res = sim_util.data_configure_simple(numPix, deltaPix) #,inverse=True)
data_class = Data(**kwargs_data_high_res)
kwargs_numerics = {'supersampling_factor': 3, 'supersampling_convolution': False}

import sys
sys.path.insert(0,'../../share_tools/')

# [[718, [-0.07, 0.06]]]
source_pos_id = [[702, [0.07, -0.06]],  #c
                 [703, [0.07, -0.06]], #f
                 [705, [0.07, -0.06]], #c
                 [706, [0.06, -0.05]], #f
                 [707, [0.02, -0.04]], #c
                 [708, [0.02, -0.04]], #+
                 [709, [0.02, -0.04]], #+
                 [710, [0.02, -0.02]], #+
                 [711, [0.02, -0.02]], #F
                 [712, [0.07, -0.02]], #c
                 [713, [0.07, -0.02]], #F
                 [713, [0.07, -0.02]], #c
                 [714, [0.07, -0.02]], #F
                 [715, [0.07, -0.02]], #F
                 [717, [0.02, 0.02]], #+
                 [718, [0.02, 0.02]], #f
                 [719, [0.06, 0.05]], #c
                 [720, [0.06, 0.05]], #f
                 [721, [0.03, -0.0]], #c
                 [724, [0.06, 0.05]], #f
                 [727, [0.06, 0.05]], #f
                 [728, [0.06, 0.05]], #c
                 [729, [0.06, 0.05]], #c
                 [730, [0.06, 0.05]], #c
                 [731, [0.06, 0.05]], #+
                 [733, [0.045, 0.045]], #f
                 [734, [0.06, 0.05]], #c
                 [735, [0.06, 0.05]], #f
                 [736, [0.06, -0.02]], #c
                 [737, [0.06, -0.05]], #c
                 [738, [0.06, 0.05]], #c
                 [739, [0.06, -0.05]], #c
                 [740, [0.06, 0.05]], #f
                 [741, [0.01, 0.05]], #+
                 [742, [0.06, 0.05]], #f
                 [743, [0.06, -0.05]], #c
                 [744, [0.06, 0.05]], #f
                 [745, [0.02, 0.05]], #f
                 [747, [0.06, 0.05]], #f
                 [748, [0.06, 0.05]], #f
                 [749, [0.06, 0.05]], #+
                 [751, [0.06, 0.05]], #f c
                 [752, [0.06, -0.02]], #c
                 [753, [0.07, 0.06]], #c
                 [754, [0.06, 0.05]], #f
                 [755, [0.06, 0.05]], #f
                 [757, [0.06, 0.05]], #f
                 [759, [0.06, 0.05]], #f
                 [760, [0.06, 0.05]],  #f
                 ]

from gene_para import gene_para
for seed in range(source_pos_id[-1][0], source_pos_id[-1][0]+1):
    print(seed)
    para=gene_para(seed=seed,fixh0=102)
    #==============================================================================
    # #######lens light 
    #==============================================================================
    from lenstronomy.LightModel.light_model import LightModel
    lens_light_para=para.lens_light()
    np.random.seed(seed)
    
    import magVSamp as mva
    lens_amp=mva.getAmp(SERSIC_in_mag=lens_light_para,zp=zp,deltaPix=deltaPix)
    lens_light_para['amp_sersic']=lens_amp * 963.490605#!!!
    # lens_light_para['q'] = 0.9 + np.random.normal(0,0.01)
    
    lens_light_tran_Reff=lens_light_para['R_sersic']/np.sqrt(lens_light_para['q'])  #!!!
    #lens_light_tran_Reff = 1.1
    lens_light_para['e1'], lens_light_para['e2'] = param_util.phi_q2_ellipticity(phi=lens_light_para['phi_G'], q=lens_light_para['q'])
    lens_light_model_list = ['SERSIC_ELLIPSE']
    kwargs_lens_light = {'amp':lens_light_para['amp_sersic'], 'R_sersic': lens_light_tran_Reff, 'n_sersic': lens_light_para['n_sersic'],
    					  'center_x': 0.0, 'center_y': 0.0, 'e1':lens_light_para['e1'], 'e2':lens_light_para['e2']}
    kwargs_lens_light_copy = copy.deepcopy(kwargs_lens_light)
    kwargs_lens_light_copy['phi_G'] =  lens_light_para['phi_G']
    kwargs_lens_light_copy['q'] =  lens_light_para['q']
    
    kwargs_lens_light_list = [kwargs_lens_light]
    lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
    
    #==============================================================================
    # ##### lens mass model 
    #==============================================================================
    from lenstronomy.LensModel.lens_model import LensModel
    kwargs_spemd = para.spemd()
    
    # kwargs_spemd['q'] = 0.9 + np.random.normal(0,0.01)
    kwargs_spemd['e1'], kwargs_spemd['e2'] = param_util.phi_q2_ellipticity(phi=kwargs_spemd['phi_G'], q=kwargs_spemd['q'])
    lens_model_list = ['PEMD','SHEAR']
    #kwargs_spemd['gamma'] = 2.
    kwargs_mass_copy = copy.deepcopy([kwargs_spemd])
    print(kwargs_spemd['q'] )
    del kwargs_spemd['phi_G']
    del kwargs_spemd['q']    
    ex_shear = {'gamma1': para.shear()[0]['e1'], 'gamma2': para.shear()[0]['e2']}
    kwargs_lens_list = [kwargs_spemd, ex_shear]
    lens_model_class = LensModel(lens_model_list)
    #==============================================================================
    # #########source light
    #==============================================================================
    source_pos=source_pos_id[-1][1]
  #Seed 220  cored-Powerlaw
    source_model_list = ['SERSIC_ELLIPSE']
    source_para=para.source_light()
    source_amp= mva.getAmp(SERSIC_in_mag=source_para,zp=zp,deltaPix=deltaPix)
    source_para['amp_sersic']=source_amp * 963.490605   #!!!
    source_light_tran_Reff=source_para['R_sersic']/np.sqrt(source_para['q'])  #To unfiy the defintion R_eff in the paper
    kwargs_source_light = {'amp': source_para['amp_sersic'], 'R_sersic': source_light_tran_Reff, 'n_sersic': source_para['n_sersic'],
                             'center_x': source_pos[0], 'center_y': source_pos[1], 'phi_G': source_para['phi_G'], 'q': source_para['q']} 
    kwargs_source_light['e1'], kwargs_source_light['e2'] = param_util.phi_q2_ellipticity(phi=kwargs_source_light['phi_G'], q=kwargs_source_light['q'])
    kwargs_source_light_copy = copy.deepcopy(kwargs_source_light)
    #kwargs_source_light_copy['mag'] = source_para['mag_sersic']
    del kwargs_source_light['phi_G']
    del kwargs_source_light['q']
    kwargs_source_list = [kwargs_source_light]
    source_model_class = LightModel(light_model_list=source_model_list)
    
    #==============================================================================
    # Setup the simulating
    #==============================================================================
    from lenstronomy.ImSim.image_model import ImageModel
    imageModel_arc = ImageModel(data_class, psf_class, lens_model_class=lens_model_class,
                            source_model_class=source_model_class, kwargs_numerics=kwargs_numerics)
    
    # generate the arc image
    image_arc = imageModel_arc.image(kwargs_lens_list, kwargs_source_list, unconvolved=True)
    
    src_mag=-2.5*np.log10(np.sum(image_arc))+zp
    pc = int(numPix/2)
    if image_arc[pc,pc]/image_arc[pc,pc+1]>20:
        print("Average center arc light for seed:", seed)
        image_arc[pc,pc]= (image_arc[pc,pc-2] + image_arc[pc,pc+2] + image_arc[pc-2,pc] + image_arc[pc+2,pc])/4  #This makes no sense
        
    import scipy.signal as signal
    image_arc_conv= signal.fftconvolve(image_arc, psf_class.kernel_point_source, mode='same')  #convolve the image
        
    #==============================================================================
    # #The information of the QSO as PSF
    #==============================================================================
    from lenstronomy.PointSource.point_source import PointSource
    from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
    np.random.seed(seed)
    amp_qRh_s_plane=np.random.uniform(0.8,1./0.8)
    qso_amp= 10.**(-0.4*(source_para['mag_sersic']-zp))*amp_qRh_s_plane
    
#    add_qso = int(input("add QSO?:\n input 0 no qso, others add qso:\t"))
    add_qso= 1
#    add_qso= 0
    
    if add_qso == 0:
    	qso_amp = 0
                  
    lensEquationSolver = LensEquationSolver(lens_model_class)
    x_image, y_image = lensEquationSolver.findBrightImage(source_pos[0], source_pos[1], kwargs_lens_list, numImages=4,
                                                          min_distance=deltaPix, search_window=numPix * deltaPix)
    mag = lens_model_class.magnification(x_image, y_image, kwargs=kwargs_lens_list)
    kwargs_ps = {'ra_image': x_image, 'dec_image': y_image, 'point_amp': np.abs(mag) * qso_amp}  # quasar point source position in the source plane and intrinsic brightness
    point_source_list = ['LENSED_POSITION']
    point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])
    
    imageModel_without_arc = ImageModel(data_class, psf_class,
                                   lens_model_class=lens_model_class, source_model_class=None,    # No arc, i.e. source_model_class=None
                                   lens_light_model_class= lens_light_model_class, point_source_class= point_source_class,
                                   kwargs_numerics=kwargs_numerics)
    image_without_arc = imageModel_without_arc.image(kwargs_lens = kwargs_lens_list, kwargs_source = [{}],
                                             kwargs_lens_light = kwargs_lens_light_list, kwargs_ps = [kwargs_ps])
    
    image_deflector = imageModel_without_arc.image(kwargs_lens = kwargs_lens_list, kwargs_source = [{}],
                                             kwargs_lens_light = kwargs_lens_light_list, kwargs_ps = [kwargs_ps], point_source_add=False, unconvolved=True)
    imageModel_source = ImageModel(data_class, psf_class, lens_model_class=None,
                            source_model_class=source_model_class, kwargs_numerics=kwargs_numerics)
    image_source = imageModel_source.image([{}], kwargs_source_list, unconvolved=True)
    kwargs_lens_light_copy['mag'] = -2.5*np.log10(np.sum(image_deflector))+zp
    kwargs_source_light_copy['mag'] = -2.5*np.log10(np.sum(image_source))+zp
    
    image_highres = image_without_arc + image_arc_conv
    print('total flux:', np.sum(image_highres))
    #print('total_mag:', -2.5*np.log10(np.sum(image_highres))+zp)
    plt.imshow(np.log10(image_highres),origin='lower',vmin=-4.5)
    plt.colorbar()
    plt.show()
    
    #%%
###################Active this part if need to see the caustic and critical line
    from  lenstronomy.Plots import lens_plot
    f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=False, sharey=False)
    lens_plot.lens_model_plot(ax, lensModel=lens_model_class, kwargs_lens=kwargs_lens_list,
                              sourcePos_x=source_pos[0], sourcePos_y=source_pos[1],
                              point_source=True, with_caustics=True)
    f.show()

    