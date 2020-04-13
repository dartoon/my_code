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
sys.path.insert(0,'../share_tools/')
from gene_para import gene_para
for seed in range(622, 623):
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
    lens_light_para['q'] = 0.9 + np.random.normal(0,0.01)
    
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
    
    kwargs_spemd['q'] = 0.9 + np.random.normal(0,0.01)
    kwargs_spemd['e1'], kwargs_spemd['e2'] = param_util.phi_q2_ellipticity(phi=kwargs_spemd['phi_G'], q=kwargs_spemd['q'])
    lens_model_list = ['SPEMD','SHEAR']
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
    source_pos=[-0.02, 0.02]  #Seed 220  cored-Powerlaw
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
    
    #print(image_arc.sum())
    #print(10.**(-0.4*(source_para['mag_sersic']-zp)))
    plt.imshow(np.log10(image_arc),origin='lower')
    plt.colorbar()
    plt.show()
    src_mag=-2.5*np.log10(np.sum(image_arc))+zp
    pc = int(numPix/2)
    if image_arc[pc,pc]/image_arc[pc,pc+1]>20:
        print("Average center arc light for seed:", seed)
        image_arc[pc,pc]= (image_arc[pc,pc-2] + image_arc[pc,pc+2] + image_arc[pc-2,pc] + image_arc[pc+2,pc])/4  #This makes no sense
        
    import scipy.signal as signal
    image_arc_conv= signal.fftconvolve(image_arc, psf_class.kernel_point_source, mode='same')  #convolve the image
    #image_arc_conv = imageModel_arc.image(kwargs_lens_list, kwargs_source_list, unconvolved=False)
    ###The place to derive the truncated lensed arc
    plt.imshow(np.log10(image_arc_conv),origin='lower',vmin=-4.5)
    plt.colorbar()
    plt.show()
    
    #==============================================================================
    # #The information of the QSO as PSF
    #==============================================================================
    from lenstronomy.PointSource.point_source import PointSource
    from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
    np.random.seed(seed)
    amp_qRh_s_plane=np.random.uniform(0.8,1./0.8)
    qso_amp= 10.**(-0.4*(source_para['mag_sersic']-zp))*amp_qRh_s_plane
    
#    add_qso = int(input("add QSO?:\n input 0 no qso, others add qso:\t"))
#    add_qso= 1
    add_qso= 0
    
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
    
    plt.imshow(np.log10(image_without_arc),origin='lower',vmin=-4.5)
    plt.colorbar()
    plt.show()
    image_highres = image_without_arc + image_arc_conv
    print('total flux:', np.sum(image_highres))
    #print('total_mag:', -2.5*np.log10(np.sum(image_highres))+zp)
    plt.imshow(np.log10(image_highres),origin='lower',vmin=-4.5)
    plt.colorbar()
    plt.show()
    
    if add_qso == 0:
    	sim_folder_name = 'sim_lens_noqso_ID_'+repr(seed)
    else:
    	sim_folder_name = 'sim_lens_ID_'+repr(seed)
    #==============================================================================
    # Creat a folder save the fits file
    #==============================================================================
    import os
    if os.path.exists(sim_folder_name)==True:
        import shutil 
        shutil.rmtree(sim_folder_name)
    os.mkdir(sim_folder_name)
    
    ##==============================================================================
    ## #Bin the image res. from high to low. 
    ##==============================================================================
    sys.path.insert(0, '../share_tools')
    import rebin
    factor=4
    pattern_x=[0,2,0,2,1,3,1,3]
    pattern_y=[0,0,2,2,3,3,1,1]      #from the info. given by observation
    ################Bin the lensed image################
    exp_grid=rebin.expend_grid(image_highres)
    cut_out=np.zeros([len(pattern_x),image_highres.shape[0]-5,image_highres.shape[1]-5])
    image_bin =np.zeros([len(pattern_x),int(image_highres.shape[0]/factor)-1,int(image_highres.shape[1]/factor)-1])
    for i in range(len(pattern_x)):
        cut_out[i]=exp_grid[pattern_x[i]:(numPix-5)+pattern_x[i],pattern_y[i]:(numPix-5)+pattern_y[i]]   #the size before bin
        image_bin[i]=rebin.block(cut_out[i],(int(numPix/factor)-1,int(numPix/factor)-1),factor=factor)
    plt.imshow(image_bin[0], origin='lower',cmap='gist_heat', norm=LogNorm())
    plt.colorbar()
    plt.show()
    ################Bin the PSF and save it################
    #exp_psf=rebin.expend_grid(psf_pixel_high_res)
    cut_fd=int((len(psf_data)-((int(len(psf_data)/8*2)-1)*4+3))/2)
    exp_psf_o=psf_data[1+cut_fd:-cut_fd,1+cut_fd:-cut_fd]+ 0  # To change it from 251 to 247.
    exp_psf=rebin.expend_grid(exp_psf_o)
    cut_len=int(round(len(exp_psf_o)/factor)*factor)
    cut_out_psf=np.zeros([len(pattern_x),cut_len,cut_len])
    image_bin_psf=np.zeros([len(pattern_x),int(cut_len/factor),int(cut_len/factor)])
    for i in range(len(pattern_x)):
        cut_out_psf[i]=exp_psf[pattern_x[i]:cut_len+pattern_x[i],pattern_y[i]:cut_len+pattern_y[i]]   #the size before bin
        image_bin_psf[i]=rebin.block(cut_out_psf[i],(int(cut_len/factor),int(cut_len/factor)),factor=factor)
        image_bin_psf[i] /= np.sum(image_bin_psf[i])  #unify the psf value
        pyfits.PrimaryHDU(image_bin_psf[i]).writeto(sim_folder_name+'/non_drizzled_psf-{0}.fits'.format(i+1),overwrite=False)
    #==============================================================================
    # Add the noise same as Ding et al. 2017a 
    ######Since two long pics only ###########
    #==============================================================================
    bf_noz = image_bin#input simulate data to bf_noz
    rms = np.zeros_like(image_bin) #input rms
    noiz = np.zeros_like(image_bin) #input noiz
    image_data_noz=np.zeros_like(image_bin) #image after noiz
    stddlong=0.016
    #stddshort=0.265
    explong=599.
    #expshort=43.98
    #exp_tot=2*explong+expshort
    for i in range(len(pattern_x)):
        rms[i]=(bf_noz[i]/(2*explong)+1/2.*stddlong**2)**0.5
        bkg_noise=(1/2.*stddlong**2)**0.5
        noiz[i]=np.random.normal(0, bkg_noise, size=rms[i].shape)
        image_data_noz[i]=noiz[i]+np.random.poisson(lam=bf_noz[i]*2*explong)/(2*explong)
        pyfits.PrimaryHDU(image_data_noz[i]).writeto(sim_folder_name+'/non_drizzled-image-{0}.fits'.format(i+1),overwrite=False)
    plt.matshow(np.log10(image_data_noz[0]),origin='lower')
    plt.colorbar()
    plt.show()
    
    #%%
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
    	ax.matshow(np.log10(image_highres),origin='lower')
    fig.savefig(sim_folder_name+'/ABCD.png')
    plt.show()  
    TD = TD_distance/const.c * fermat_po / const.day_s * const.arcsec**2 * const.Mpc
    print("Time delay of ABCD - A :\n\t", TD-TD[0])
    
    #%%Save infor
    qso_mag=-2.5*np.log10(np.sum(kwargs_ps['point_amp']))+zp
    from roundme import roundme
    #Write information to record true information
    lens_info = open(sim_folder_name+'/lens_all_info.txt','w') 
    lens_info.write("Unit:\n\tLength in arcsecond scale, 'phi_G' is Angle in radian system start from x axis anticlockwise.")
    lens_info.write("\nCosmological para\n\tFlatLambdaCDM, with Om=0.27 and H0: "+repr(roundme(para.H0))+ "km/s/Mpc")
    lens_info.write("\nPixel size is 0.13'' and 0.08'' before and after drizzle")
    lens_info.write("\nTime delay distance: TD_distance=(1+z_l)*D_l*D_s/D_ls:"+ repr(roundme(TD_distance)) + "Mpc")
    if len(TD)==4:
    	lens_info.write("\nTime delay of BCD - A :\n\t" + repr(roundme(TD-TD[0])[1:])  + "days")
    if len(TD)==2:
    	lens_info.write("\nTime delay of B - A :\n\t" + repr(roundme(TD-TD[0])[1:])  + "days")
    lens_info.write("\nZeropoint of filter (AB system): \t"+ repr(zp))
    lens_info.write("\nLens/Source redshift:\t"+ repr([roundme(para.z_lens), roundme(para.z_source)]))
    re_shear=para.shear()
    re_shear[1]['pa'] *= (np.pi)/180
    re_shear[1]['phi_G']=re_shear[1].pop('pa')
    lens_info.write("\nLens mass model:\n"+"\tSPEMD:"+repr(roundme(para.spemd()))+"\n\tShear: \t"+repr(roundme(re_shear)))
    lens_info.write("'\t Note that e1=-b*cos(2*phi_G); e2=-b*sin(2*phi_G)'")
    lens_info.write("\nLens light: \n\t"+repr(roundme(kwargs_lens_light_copy)))
    lens_info.write("\nSource light in source plane:\n\t"+repr(roundme(kwargs_source_light_copy)))
    lens_info.write("\nAGN light:" +\
                    "\n\tAGN position in source plane:\t" + repr(roundme(source_pos[0]))+","+repr(roundme(source_pos[1]))+\
                    "\n\tAGN amplitude in source plane:\t"+repr(roundme(qso_amp)) +\
                    "\n\tAGN position in image plane:\t" + "x: " +repr(roundme(kwargs_ps['ra_image'])) + "\ty: " + repr(roundme(kwargs_ps['dec_image']))+\
                    "\n\tAGN amplitude in image plane:\t"+repr(roundme(kwargs_ps['point_amp'])))
    lens_info.write("\nHost galaxy mag in the image plane: "+ repr(roundme(src_mag))+" mag, AGN total mag in the image plane: "+ repr(roundme(qso_mag))+" mag")
    np.random.seed(seed)
    k_ext= 0 #np.random.normal(scale=0.025)
    TD_obs=TD*(1-k_ext)
    lens_info.write("\nkappa_ext: \n\t"+repr(roundme(k_ext, prec=4)))
    lens_info.write("\nTime delay with external kappa i.e. TD_obs=TD*(1-k_ext):\n\t" + repr(roundme(TD_obs-TD_obs[0])[1:])  + "days")
    lens_info.close()
    
    #Write information for good team
    lens_info_4Goodteam = open(sim_folder_name+'/lens_info_for_Good_team.txt','w') 
    lens_info_4Goodteam.write("Pixel size is 0.13'' and 0.08'' before and after drizzle")
    lens_info_4Goodteam.write("\nZeropoint of filter (AB system): \t"+ repr(zp))
    lens_info_4Goodteam.write("\nLens/Source redshift:\t"+ repr([roundme(para.z_lens), roundme(para.z_source)]))
    lens_info_4Goodteam.write("\nExternal Convergence: Kext= 0 +/- 0.025")
    TD_obs_t=TD_obs-TD_obs[0]
    TD_err_l=abs(TD_obs_t*0.01)
    TD_err_l[1:][TD_err_l[1:]<0.25]=0.25
    TD_err=np.random.normal(0, TD_err_l, size=TD_err_l.shape)
    TD_obs_err=TD_obs_t+TD_err
    if qso_amp != 0:
        if len(TD)==4:
            lens_info_4Goodteam.write("\nTime delay of BCD - A :\n\t" + repr(roundme(TD_obs_err)[1:])  + "days, error level: " + repr(roundme(TD_err_l, prec=2)[1:]) + "days")
        if len(TD)==2:
            lens_info_4Goodteam.write("\nTime delay of B - A :\n\t" + repr(roundme(TD_obs_err)[1:])  + "days, error level: " + repr(roundme(TD_err_l, prec=2)[1:]) + "days")
    lens_info_4Goodteam.close()
    
    import pickle
    picklename = sim_folder_name + '/sim_kwargs.pkl'
    model_lists = [lens_model_list, lens_light_model_list, source_model_list, point_source_list]
    para_s = [kwargs_lens_list, kwargs_lens_light_list, kwargs_source_list, kwargs_ps]
    lens_info = [para.z_lens, para.z_source, TD_distance, TD_obs_t, TD_obs_err, TD_err_l]
    pickle.dump([model_lists, para_s, lens_info], open(picklename, 'wb'))

#%%
###################Active this part if need to see the caustic and critical line
#from  lenstronomy.Plots import lens_plot
#f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=False, sharey=False)
#lens_plot.lens_model_plot(ax, lensModel=lens_model_class, kwargs_lens=kwargs_lens_list,
#                          sourcePos_x=source_pos[0], sourcePos_y=source_pos[1],
#                          point_source=True, with_caustics=True)
#f.show()

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
