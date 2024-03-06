#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 14:58:51 2023

@author: Dartoon
"""
#Simulation based on the realistic JWST imaging data
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.astro_tools import read_pixel_scale

filt = 'F356W'
#%%
folder = '../../NIRCam_data/Nov14/bkg_removed'
files = glob.glob(folder+'/SHELLQs_J2236p0032_{0}_Nov14_i2d_rmbkg.fits'.format(filt))  #!!! Load in your own JWST data
file = files[0]
im = pyfits.open(file)
data = im[1].data
header = im[1].header

zp = -2.5*np.log10(2.350443 * 10**(-5) *read_pixel_scale(header)**2/3631) #- 2.5*np.log10(flux_mjsr)  #zp for flux

#Load in the PSF information:
PSF_lib_files = glob.glob('../../model_z6_data_id0/stage3_all/material/*'+filt[:-1]+'*_PSF_Library_idx0.pkl')[0]
_, PSF_list_clean, _, _ = pickle.load(open(PSF_lib_files,'rb'))  #A list of PSF



#%% Use grab a empty sky
rad = 50                #The radius to cutout, final simulation in 2*rad+1
for i in range(1000):
    pos = [int(np.random.uniform(400, len(data)-400)), int(np.random.uniform(400, len(data)-400))]  #!!!
    cut1 = data[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1]
    if np.sum(abs(cut1)>0.02) <10 and filt == 'F356W' or np.sum(abs(cut1)>0.1) <10 and filt == 'F150W':
        break
from galight.tools.astro_tools import plt_fits

print("The region will be added with our our simulation using 'pos' values")
plt_fits(cut1)  


#%% Start to make the mocks:
import warnings
warnings.filterwarnings("ignore")
from lenstronomy.ImSim.image_model import ImageModel
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LightModel.light_model import LightModel
import lenstronomy.Util.param_util as param_util
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
    
PSF1, PSF2 = PSF_list_clean[0],  PSF_list_clean[2]  #Use PSF1 to simulation and PSF2 to model the mocks

seed = 0
psf = PSF1 #!!! Use PSF 
psf[psf<0] = 0.  #!!! Only allow positive values in our PSF
pix_s = read_pixel_scale(header)

#Use Lenstronomy to make the simulation

kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': psf, 'pixel_size': pix_s}
psf_class = PSF(**kwargs_psf)
kwargs_data = sim_util.data_configure_simple(2*rad+1, pix_s, inverse=True)
data_class = ImageData(**kwargs_data)
point_source_list = ['UNLENSED']


kwargs_numerics = {'supersampling_factor': 5, 'supersampling_convolution': False} # 'point_source_supersampling_factor': 2}

kwargs_sersic = {'R_sersic': 0.3,   #The Reff size of the Sersic to use
                 'n_sersic': 3}

light_model_list = ['SERSIC_ELLIPSE']
lightModel = LightModel(light_model_list=light_model_list,sersic_major_axis=True) #sersic_major_axis = True in Galight and Galfit

host_mag = 23  #The host magnitude #!!!
_host_flux = 10**(0.4*(zp-host_mag))  
q = np.random.uniform(0.5,0.9)
phi = np.random.uniform(0.,2*np.pi)   
e1, e2 = param_util.phi_q2_ellipticity(phi=phi, q=q)
kwargs_sersic['e1'] = e1
kwargs_sersic['e2'] = e2
kwargs_sersic['amp'] = 1
center_x, center_y = np.random.uniform(-1.5, 1.5) * pix_s, np.random.uniform(-1.5,1.5)* pix_s
kwargs_sersic['center_x'] = center_x
kwargs_sersic['center_y'] = center_y
kwargs_host = [kwargs_sersic]
imageModel_host = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                point_source_class=None, kwargs_numerics=kwargs_numerics)
medi_host_flux = np.sum(imageModel_host.image(kwargs_lens_light=kwargs_host, unconvolved=True))
amp = 1. / medi_host_flux * _host_flux        
kwargs_sersic['amp'] = amp
kwargs_host = [kwargs_sersic]
host_highres = imageModel_host.image(kwargs_lens_light=kwargs_host, unconvolved=False)
if im[0].header['CHANNEL'] == 'LONG':
    gain_value = 2
else:
    gain_value = 1.8

exp = header['XPOSURE']
flux_mjsr = header['PHOTMJSR']
#Adding Possion noise to the host mock
host_highres_noise = np.random.poisson(lam=host_highres/header['PHOTMJSR']*exp*gain_value) / (1/header['PHOTMJSR']*exp*gain_value)  #Non-drizzled imaged

pointSource = PointSource(point_source_type_list=point_source_list)
imageModel_ps = ImageModel(data_class, psf_class, lens_light_model_class=None,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
ps_mag = 21.5   #The PS magnitude #!!!
_ps_amp = 10**(0.4*(zp-ps_mag))  
PS_pos = np.random.uniform(-1.5, 1.5) * pix_s, np.random.uniform(-1.5,1.5)* pix_s
kwargs_pointsource= {'ra_image': np.array([PS_pos[0]]),
                     'dec_image': np.array([PS_pos[1]]),
                     'point_amp': np.array([_ps_amp])}
kwargs_ps = [kwargs_pointsource]
ps_highres = imageModel_ps.image(kwargs_ps=kwargs_ps, unconvolved=False)

print("The simulation for Sersic host, without adding noise ")
plt_fits(host_highres_noise)  
print("The simulation for PS, without adding noise ")
plt_fits(ps_highres)  

print("The simulation for Sersic+PS, without adding noise ")
plt_fits(host_highres_noise+ps_highres)  

#%% Adding simulated Host and PS to the mock data
import copy 
data_mock = copy.deepcopy(data)
data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += host_highres_noise
data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += ps_highres

#%%Fitting the mock data
wht = im[4].data # The WHT map
exp = im[0].header['EFFEXPTM']
if im[0].header['CHANNEL'] == 'LONG':
    gain_value = 2
    expsize = 1
else:
    gain_value = 1.8
    expsize = 1 
flux_mjsr = header['PHOTMJSR']
# exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
exp_map = exp * wht/wht.max() / flux_mjsr * gain_value
data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header=header,
                            rm_bkglight = False, if_plot=False, zp = zp, exptime=exp_map)  #Gain value assuming as 1
data_process.generate_target_materials(radius=rad, 
                                       create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=False)
data_process.apertures = [data_process.apertures[0]]
use_psf = PSF2
use_psf[use_psf<0] = 0.
data_process.PSF_list = [use_psf]
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
fit_sepc.build_fitting_seq()
# fit_sepc.kwargs_params['lens_light_model'][0] = [fit_run_.kwargs_result['kwargs_lens_light'][0]]
plot_fit_name = 'sim_result_{1}_seed{0}'.format(seed,filt)
fit_run = FittingProcess(fit_sepc, savename = plot_fit_name, fitting_level='norm')
fit_run.run(algorithm_list = ['PSO','PSO'], fitting_level=['norm','deep'])
fit_run.plot_final_qso_fit()
# fit_run.dump_result()

