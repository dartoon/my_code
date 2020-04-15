#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:11:45 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from astropy.cosmology import FlatLambdaCDM
import glob

from lenstronomy.ImSim.image_model import ImageModel
import lenstronomy.Util.param_util as param_util
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.psf import PSF
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LightModel.light_model import LightModel

sys.path.insert(0, '../share_tools')
#%%Set up data basic information
from quasar_info import qso_info
import rebin #From my share_tools    
oversample = 4
seed = 100
ID = 1 
for filt_i in range(4): #int(input("which filter 0: 'F444W', 1: 'F356W', 2: 'F200W', 3: 'F150W':\n"))
    filt  = ['F444W', 'F356W', 'F200W', 'F150W'][filt_i]
    numPix = [341, 341, 645, 645][filt_i]  # total frame pixel size #!!!Need to be changed for different filter
#    zp = [28., 27.9, 26.7, 27.75][filt_i]   #Using ETC, (616-556) total flux for 23.5 ab mag objects.  #Need to check
    zp = [27.3012, 27.1841, 27.0383, 26.8627][filt_i]   #Using mirage
    pix_scale = [0.063, 0.063, 0.031, 0.031][filt_i] #After dirzzled
    #properties:
    z_s = qso_info['ID'+repr(ID)]['z']       #AGN redshift
    point_mag = [28, 28.7, 29.0, 28.9][filt_i]
    host_mag = qso_info['ID'+repr(ID)]['galaxy_{0}_mag'.format(filt)]
    host_n = np.random.uniform(2,4)   #Host effective radius, unit: Kpc
    host_Reff_kpc = np.random.uniform(2,3)   #Host effective radius, unit: Kpc
    np.random.seed(seed = seed)
    #host_ratio = np.random.uniform(0.4, 0.7) #Set the random host flux ratio [40% - 70%].
    host_flux = 1
    point_flux = 10**(0.4*(zp - point_mag))
    total_flux =  point_flux + host_flux
    host_ratio = host_flux/total_flux
    if host_ratio< 0.1:
        host_ratio =np.random.uniform(0.1,0.2)
        total_flux = point_flux/(1-host_ratio)
        host_flux = total_flux - point_flux
    elif host_ratio > 0.95:
        host_ratio =np.random.uniform(0.8,0.95)
        total_flux = point_flux/(1-host_ratio)
        host_flux = total_flux - point_flux
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    #Take the PSF: Use PSF 1 - 8 to simulate, but use PSF0 to model the data
    psf_take_id = np.random.randint(1,8)
    psf_take_name = 'webPSF/highres_PSF_'+filt+'/PSF_id{0}.fits'.format(psf_take_id)
    psf = pyfits.open(psf_take_name)
    #Build up the simulation:
    pix_s = psf[0].header['PIXELSCL'] #* oversample
    host_Reff = host_Reff_kpc/scale_relation   #In arcsec
    psf_data = psf[0].data
    psf_data = psf_data[1:,1:]
    psf_data /= psf_data.sum()
    kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': pix_s}
    kwargs_data_high_res = sim_util.data_configure_simple(numPix, pix_s)
    data_class = ImageData(**kwargs_data_high_res)
    psf_class = PSF(**kwargs_psf_high_res)
    center_x, center_y = np.random.uniform(-5,5) * pix_s*oversample, np.random.uniform(-5,5)* pix_s*oversample
    point_amp = point_flux
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
    light_model_list = ['SERSIC_ELLIPSE']
    lightModel = LightModel(light_model_list=light_model_list)
    q = np.random.uniform(0.5,0.9)
    phi = np.random.uniform(0.,2*np.pi)
    e1, e2 = param_util.phi_q2_ellipticity(phi=phi, q=q)
    kwargs_numerics = {'supersampling_factor': 2, 'supersampling_convolution': False}
    kwargs_sersic_medi = {'amp': 1. , 'n_sersic': host_n, 'R_sersic': host_Reff/np.sqrt(q), 'e1': e1, 'e2': e2,
                     'center_x': center_x, 'center_y': center_y}
    kwargs_host_medi = [kwargs_sersic_medi]
    imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
    medi_host_flux = np.sum(imageModel.image(kwargs_lens_light=kwargs_host_medi, unconvolved=True))
    amp = 1. / medi_host_flux * host_flux
    kwargs_sersic = {'amp': amp, 'n_sersic': host_n, 'R_sersic': host_Reff/np.sqrt(q), 'e1': e1, 'e2': e2,
                     'center_x': center_x, 'center_y': center_y}
    kwargs_host = [kwargs_sersic]
    
    ## simulate image with the parameters we have defined above #
    total_highres = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=kwargs_ps, unconvolved=False)
    host_highres = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=[{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [0]}], unconvolved=False)
    point_highres = total_highres - host_highres
#        print("AGN image:")
#        plt.imshow(total_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
#        plt.colorbar()
#        plt.show()
#        print("HOST image:")    
#        plt.imshow(host_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
#        plt.colorbar()
#        plt.show()
#        print("Point image:")    
#        plt.imshow(point_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
#        plt.colorbar()
#        plt.show()
    #%%
    ##==============================================================================
    ## #Bin the image res. from high to low. 
    ##==============================================================================
    factor=oversample
    pattern_x=[0]
    pattern_y=[0]      #from the info. given by observation
    ################Bin the lensed image################
    point_exp_grid=rebin.expend_grid(point_highres)
    point_cut_out=np.zeros([len(pattern_x),point_highres.shape[0]-5,point_highres.shape[1]-5])
    point_image_bin =np.zeros([len(pattern_x),int(point_highres.shape[0]/factor)-1,int(point_highres.shape[1]/factor)-1])
            
            
    for i in range(len(pattern_x)):
        point_cut_out[i]=point_exp_grid[pattern_x[i]:(numPix-5)+pattern_x[i],pattern_y[i]:(numPix-5)+pattern_y[i]]   #the size before bin
        point_image_bin[i]=rebin.block(point_cut_out[i],(int(numPix/factor)-1,int(numPix/factor)-1),factor=factor)
#        print("Rebin image:")    
#        plt.imshow(total_image_bin[0], origin='lower',cmap='gist_heat', norm=LogNorm())
#        plt.colorbar()
#        plt.show()
    #==============================================================================
    # Add the noise same as Ding et al. 2017a 
    ######Since two long pics only ###########
    #==============================================================================
    rms = np.zeros_like(point_image_bin) #input rms
    noiz = np.zeros_like(point_image_bin) #input noiz
    image_data_noise=np.zeros_like(point_image_bin) #image after noiz
    exptim= 10000.   #units of seconds. 
    stdd = [1.6, 1.06, 0.77, 0.75][filt_i] /np.sqrt(exptim)   #An empirical formula from ETC
    noise = (abs(point_image_bin[0]/exptim)+stdd**2)**0.5
    Total_SNR = point_image_bin[0] / noise

    r = [0.16, 0.16, 0.08, 0.08][filt_i]
    framesize = r / pix_scale * 2

    half_r = int(framesize/2)
    peak = np.where(Total_SNR==Total_SNR.max())
    peak = [peak[0][0], peak[1][0]]
    Total_SNR = Total_SNR[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
    point_image = point_image_bin[0][peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
    noise = noise[peak[0]-half_r:peak[0]+half_r+1,peak[1]-half_r:peak[1]+half_r+1]
    print("Pring {0} SNR map:".format(filt))#, "Total SNR:", Total_SNR.sum())
    print("Total SNR:", np.sum(point_image) / np.mean(noise))
    plt.imshow(Total_SNR, origin='lower')#,cmap='gist_heat', norm=LogNorm())
    plt.colorbar()
    plt.show()
    