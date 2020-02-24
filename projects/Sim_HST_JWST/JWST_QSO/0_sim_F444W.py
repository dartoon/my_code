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
#import sys
from astropy.cosmology import FlatLambdaCDM

#%%

#%%Generate PSF:
import webbpsf
nc = webbpsf.NIRCam()
nc.filter =  'F444W'
print("Generate PSF higher resolution by factor of 4:")
oversample = 4
psf = nc.calc_psf(oversample=oversample)     # returns an astropy.io.fits.HDUlist containing PSF and header
plt.imshow(psf[0].data, origin='lower',cmap='gist_heat', norm=LogNorm())
plt.colorbar()
plt.show()
print("Done.")

#%%Set up the simulation parameters:
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
pix_s = psf[0].header['PIXELSCL'] #* oversample
#host_total_flux = 27.0    #AB magnitude.
#host_ratio = 0.5     #Host to total flux ratio. 

z_s =6.0        #AGN redshift
scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc

host_flux = 100
point_flux = 50

host_Reff_kpc = 2.5   #Host effective radius, unit: Kpc
host_Reff = host_Reff_kpc/scale_relation   #In arcsec
host_n = 2.5   #Host effective radius, unit: Kpc

import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.psf import PSF
from lenstronomy.Data.imaging_data import ImageData

numPix = 321  #  pixel size

psf_data = psf[0].data
psf_data = psf_data[1:,1:]
psf_data /= psf_data.sum()
kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': pix_s}
kwargs_data_high_res = sim_util.data_configure_simple(numPix, pix_s)
data_class = ImageData(**kwargs_data_high_res)

psf_class = PSF(**kwargs_psf_high_res)
center_x = 0.02
center_y = 0.01

point_amp = point_flux
from lenstronomy.PointSource.point_source import PointSource
point_source_list = ['UNLENSED']
pointSource = PointSource(point_source_type_list=point_source_list)
kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]

from lenstronomy.LightModel.light_model import LightModel
light_model_list = ['SERSIC_ELLIPSE']
lightModel = LightModel(light_model_list=light_model_list)
import lenstronomy.Util.param_util as param_util
e1, e2 = param_util.phi_q2_ellipticity(phi=0.3, q=0.6)

kwargs_sersic_init = {'amp': 1, 'n_sersic': host_n, 'R_sersic': host_Reff, 'e1': e1, 'e2': e2,
                 'center_x': center_x, 'center_y': center_y}
kwargs_host_ini = [kwargs_sersic_init]
from lenstronomy.ImSim.image_model import ImageModel
kwargs_numerics = {'supersampling_factor': 3, 'supersampling_convolution': False}
imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
# simulate image with the parameters we have defined above #
image_host = imageModel.image(kwargs_lens_light=kwargs_host_ini, unconvolved=True)
amp = host_flux/image_host.sum()

kwargs_sersic = {'amp': amp, 'n_sersic': host_n, 'R_sersic': host_Reff, 'e1': e1, 'e2': e2,
                 'center_x': center_x, 'center_y': center_y}
kwargs_host = [kwargs_sersic]
image_host = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=kwargs_ps, unconvolved=False)

plt.imshow(image_host, origin='lower',cmap='gist_heat', norm=LogNorm())
plt.colorbar()
plt.show()

#%%
##==============================================================================
## #Bin the image res. from high to low. 
##==============================================================================
import rebin
factor=oversample
pattern_x=[0,2,0,2,1,3,1,3]
pattern_y=[0,0,2,2,3,3,1,1]      #from the info. given by observation
################Bin the lensed image################
exp_grid=rebin.expend_grid(image_host)
cut_out=np.zeros([len(pattern_x),image_host.shape[0]-5,image_host.shape[1]-5])
image_bin =np.zeros([len(pattern_x),int(image_host.shape[0]/factor)-1,int(image_host.shape[1]/factor)-1])
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
#    pyfits.PrimaryHDU(image_bin_psf[i]).writeto('../../../../TDLMC_material/mock_data/{2}/{3}/{4}-seed{0}/non_drizzled_psf-{1}.fits'.format(seed,i+1,rung,code,filt),overwrite=False)
plt.imshow(image_bin_psf[0], origin='lower',cmap='gist_heat', norm=LogNorm())
plt.colorbar()
plt.show()

#==============================================================================
# Add the noise same as Ding et al. 2017a 
######Since two long pics only ###########
#==============================================================================
bf_noz = image_bin#input simulate data to bf_noz
rms = np.zeros_like(image_bin) #input rms
noiz = np.zeros_like(image_bin) #input noiz
image_data_noz=np.zeros_like(image_bin) #image after noiz
#stddlong=0.016  #!!! Need to be confirmed For 10000s, 0.016. 
stddlong=0.042  #!!! Need to be confirmed For 1250s, 0.042. 
explong=1250.  #units of seconds.
for i in range(len(pattern_x)):
    rms[i]=(bf_noz[i]/(explong)+stddlong**2)**0.5
    bkg_noise=(1/2.*stddlong**2)**0.5
    noiz[i]=np.random.normal(0, bkg_noise, size=rms[i].shape)
    image_data_noz[i]=noiz[i]+np.random.poisson(lam=bf_noz[i]*2*explong)/(2*explong)
#    pyfits.PrimaryHDU(image_data_noz[i]).writeto('../../../../TDLMC_material/mock_data/{2}/{3}/{4}-seed{0}/non_drizzled-lens-image-{1}.fits'.format(seed,i+1,rung,code,filt),overwrite=False)
plt.imshow(image_data_noz[0], origin='lower',cmap='gist_heat', norm=LogNorm())
plt.colorbar()
plt.show()

