#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:59:42 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#%%Load image and define a empty place

filt = 'f356w'
folder = 'JWST_CEERS/'
file = 'ceers5_{filt}_i2d.fits'.format(filt=filt)

im = pyfits.open(folder+file)
data_sb = im[1].data
header = im[1].header

# print('For flux value in unit of MJy/sr.') #https://en.wikipedia.org/wiki/AB_magnitude
value_unit = header['BUNIT']
# flux(Mjy/sr) * 2.350443 * 10**(-5) *0.03**2   #Flux to Jy  https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/18/
print("Data unit:", value_unit)
# data_sb = data_sb[2000:4000,2000:4000]
data_sb = data_sb[100:2100,100:2100]
zp = -2.5*np.log10(2.350443 * 10**(-5) *0.03**2/3631)

header0 = im[0].header
img_filter = header0['FILTER']
img_cam = header0['APERNAME'] #In JDAT'simulation it is 'DETECTOR'
exptime = header0['TEXPTIME'] #The assumed exp time.
from galight.tools.astro_tools import plt_fits
plt_fits(data_sb)
#%%Obtain PSF stars:
from galight.data_process import DataProcess
# zp = 31.4 - 2.5*np.log10(header['PHOTMJSR'])  #Calculate the correspondingly zp as DN/S #This is wrong! See 3_...
data_process_ = DataProcess(fov_image = data_sb, target_pos = [1170., 940.], pos_type = 'pixel', header = header,
                           rm_bkglight = True, exptime = np.ones_like(data_sb)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
# data_process_.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
#                                       exp_sz= 1.2, npixels = 15, if_plot=True)
data_process_.find_PSF(radius = 60, user_option = True, psf_edge=10)
data_process_.plot_overview(label = 'Example', target_label = None)
#%%

psf_ids = [3, 8]  #use 5 to do mock and use 10 to do fitting.
# psf_ids = [5,9]  #use 5 to do mock and use 10 to do fitting.
psf_list = [data_process_.PSF_list[i]/np.sum(data_process_.PSF_list[i]) for i in psf_ids]

# from galight.tools.cutout_tools import psf_clean
# psf_list[0] = psf_clean(psf_list[0], if_plot=True)
# plt_fits(psf)
#%%
# data_process_.find_PSF(radius = 60, PSF_pos_list = [[341., 444.], [1925., 1127.]])
# data_process_.plot_overview(label = 'Example', target_label = None)
# psf_list = [data_process_.PSF_list[i]/np.sum(data_process_.PSF_list[i]) for i in [0,1]]


#%%For mock galaxy
#Generate the QSO galaxy info
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/my_code/projects/Sim_HST_JWST/JWST_QSO')
from quasar_info import qso_info
from astropy.cosmology import FlatLambdaCDM
z_s = qso_info['ID3']['z']
qso_mag = qso_info['ID3']['AGN_{}_mag'.format(img_filter)]
qso_flux = 10**(-0.4*(qso_mag-zp))
galaxy_mag = qso_info['ID3']['galaxy_{}_mag'.format(img_filter)]
galaxy_flux = 10**(-0.4*(galaxy_mag-zp))
host_Reff_kpc = np.random.uniform(1,3)   #Host effective radius, unit: Kpc
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
host_Reff = host_Reff_kpc/scale_relation   #In arcsec
host_Reff_pix = host_Reff/0.03

#Load mock galaxy
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source')
from source_info import source_list #[0]: file name [1]: total size [2]: galfit R_e [3]:R_e/totalzise
source_id = 4
source_name=source_list[0][source_id]
print(source_name)
Re = source_list[2][source_id]
hdu = pyfits.open('/Users/Dartoon/Astro/Projects/Lens_Model_challenge/TDSLMC/simulating/material/real_source/fix/{0}_fix.fits'.format(source_name))
hd_gal_img = hdu[0].data
# hd_gal_img = hd_gal_img

hdu.close()
print("source name:", source_name)
from galight.tools.astro_tools import plt_fits
plt_fits(hd_gal_img)
#project:
from scipy.ndimage import zoom
project_gal_img = zoom(hd_gal_img, host_Reff_pix/Re)
project_gal_img = project_gal_img/np.sum(project_gal_img)*galaxy_flux
plt_fits(project_gal_img)
#convolve:
from scipy import signal
conv_project_gal_img = signal.fftconvolve(project_gal_img, psf_list[0], mode='full')
plt_fits(conv_project_gal_img)

import copy
conv_project_qso_img = copy.deepcopy(conv_project_gal_img)
leng = len(conv_project_gal_img)

cut = (len(conv_project_gal_img) - len(psf_list[0]))/2
if cut == int(cut):
    cut = int(cut)
    conv_project_qso_img[cut:-cut,cut:-cut] = conv_project_qso_img[cut:-cut,cut:-cut] + qso_flux*psf_list[0]
else:
    cut = int(cut)
    conv_project_qso_img[cut:-cut-1,cut:-cut-1] = conv_project_qso_img[cut:-cut-1,cut:-cut-1] + qso_flux*psf_list[0]
plt_fits(conv_project_qso_img)


#Add Noise:
conv_project_qso_img=np.random.poisson(lam=abs(conv_project_qso_img)*exptime)/(exptime)
plt_fits(conv_project_qso_img)

#Add to fov image
# pos = [1270, 940]
pos = [750, 500]
rad = len(conv_project_qso_img)/2
data_mock = copy.deepcopy(data_sb)
if rad != int(rad):
    rad = int(rad)
    data_mock[ pos[1]-rad:pos[1]+rad+1, pos[0]-rad:pos[0]+rad+1] += conv_project_qso_img
elif rad == int(rad):
    rad = int(rad)
    data_mock[ pos[1]-rad:pos[1]+rad, pos[0]-rad:pos[0]+rad] += conv_project_qso_img

#%%Obtain PSF stars:
data_process = DataProcess(fov_image = data_mock, target_pos = pos, pos_type = 'pixel', header = header,
                            rm_bkglight = True, exptime = np.ones_like(data_sb)*exptime, if_plot=False, zp = zp)  #Gain value assuming as 1
data_process.generate_target_materials(radius=65, create_mask = False, nsigma=2.8, if_select_obj=False,
                                      exp_sz= 1.2, npixels = 15, if_plot=True)
# data_process.find_PSF(radius = 50, user_option = True, psf_edge=10)
data_process.plot_overview(label = 'Example', target_label = None)
data_process.PSF_list = [psf_list[0]]

#Start to produce the class and params for lens fitting.
from galight.fitting_specify import FittingSpecify

# data_process.apertures = []
fit_sepc = FittingSpecify(data_process)
fit_sepc.prepare_fitting_seq(point_source_num = 1) #, fix_n_list= [[0,4],[1,1]])
fit_sepc.build_fitting_seq()

#Plot the initial settings for fittings. 
fit_sepc.plot_fitting_sets()

#Setting the fitting method and run.
from galight.fitting_process import FittingProcess
fit_run = FittingProcess(fit_sepc, savename = 'savename', fitting_level='norm')
fit_run.run(algorithm_list = ['PSO', 'PSO'])
fit_run.plot_final_qso_fit()
print('inferred galaxy flux, mag, Re:\n\t', round(fit_run.final_result_galaxy[0]['flux_within_frame'],2), 
      round(fit_run.final_result_galaxy[0]['magnitude'],2), 
      round(fit_run.final_result_galaxy[0]['R_sersic'],2))
      # round(fit_run.final_result_galaxy[0]['n_sersic'],2),)
print('True galaxy flux, mag, Re:\n\t', round(galaxy_flux,2), round(galaxy_mag,2),
      round(host_Reff,2))

# print(fit_run.final_result_ps[0])