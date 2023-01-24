#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:26:38 2023

@author: Dartoon
"""

import matplotlib
import copy
from scipy.optimize import curve_fit
from photutils.aperture import aperture_photometry
from matplotlib.colors import LogNorm
from photutils.aperture import RectangularAperture
from astropy.coordinates import Angle
from target_info import target_info
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import pickle

idx = 1
info = target_info[str(idx)]
target_id, RA, Dec, z = info['target_id'], info['RA'], info['Dec'], info['z']

# files = glob.glob('./*fit_material*sm*/data_process_idx{0}_*_psf*.pkl'.format(idx))
# files.sort()
run_folder = '../model_z6_data_id{0}/stage3_all/'.format(idx)  # !!!
# run_folder = 'stage3_second_half/' #!!!
z_str = str(z)

# filters = ['F150W', 'F356W']
filters = ['F356W']
for top_psf_id in [0]:
    for count in range(len(filters)):
        fit_run_list = []
        # idx = idx_info
        filt = filters[count]
        if filt == 'F150W':
            cmap = 'inferno'
        else:
            cmap = 'gist_heat'
        my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap)
                            )  # copy the default cmap
        my_cmap.set_bad('black')

        PSF_lib_files = glob.glob(
            run_folder+'material/*'+filt[:-1]+'*_PSF_Library_idx{0}.pkl'.format(idx))[0]
        # idx, filt= item
        # fit_files = glob.glob(run_folder+'*fit_material*/fit_run_withcentralMask_idx{0}_{1}_FOV*.pkl'.format(idx, filt))#+\
        fit_files = glob.glob(
            run_folder+'*fit_material*/fit_run_idx{0}_{1}_*.pkl'.format(idx, filt))  # +\
        if idx == 1:
            fit_files = glob.glob(
                run_folder+'*fit_material*/fit_run_fixn1__idx{0}_{1}_*.pkl'.format(idx, filt))  # +\
        fit_files.sort()
        for i in range(len(fit_files)):
            fit_run_list.append(pickle.load(open(fit_files[i], 'rb')))
        chisqs = np.array(
            [fit_run_list[i].reduced_Chisq for i in range(len(fit_run_list))])
        sort_Chisq = chisqs.argsort()
        print('idx', idx, filt, "Total PSF NO.", 'chisq', chisqs[sort_Chisq[top_psf_id]], len(
            sort_Chisq), fit_files[sort_Chisq[top_psf_id]])
        fit_run = fit_run_list[sort_Chisq[top_psf_id]]
        fit_run.savename = 'figures/' + fit_run.savename+'_'+filt
        fit_run.plot_final_qso_fit(
            target_ID=target_id+'$-$'+filt, save_plot=True, cmap=my_cmap)

# %%
# Creat the slit
theta_dict = {0: 137.475, 1: 139.818}

theta = Angle(theta_dict[idx], 'deg')

image_host = fit_run.flux_2d_out['data-point source']
# image_host = fit_run.image_host_list[0]
image_ps = fit_run.image_ps_list[0]
image_total = fit_run.flux_2d_out['data']
fit_run.cal_astrometry()
f_center = len(image_host)/2
deltaPix = fit_run.fitting_specify_class.deltaPix
w = 0.2 / deltaPix
h = 0.6 / deltaPix

aper = RectangularAperture((f_center+fit_run.final_result_ps[0]['position_xy'][0],
                            f_center++fit_run.final_result_ps[0]['position_xy'][1]), w, h, theta=theta)

my_cmap = copy.copy(matplotlib.cm.get_cmap(
    'gist_heat'))  # copy the default cmap
my_cmap.set_bad('black')
vmin = 1.e-3
image = image_host
vmax = image.max()
plt.imshow(image, origin='lower', cmap=my_cmap, norm=LogNorm(
    vmin=vmin, vmax=vmax))  # , vmin=vmin, vmax=vmax)
aper.plot(color='blue',
          lw=3.5)
plt.show()

phot_table_host = aperture_photometry(image_host, aper)
phot_table_ps = aperture_photometry(image_ps, aper)
phot_table_total = aperture_photometry(image_total, aper)
print('slit flux of host:', round(phot_table_host['aperture_sum'].value[0], 3),
      '\nslit flux of qso:', round(phot_table_ps['aperture_sum'].value[0], 3),
      '\nslit flux of total,', round(phot_table_total['aperture_sum'].value[0], 3))
# %%
splt = 12
plt.imshow(image, origin='lower', cmap=my_cmap, norm=LogNorm(
    vmin=vmin, vmax=vmax))  # , vmin=vmin, vmax=vmax)
image_host = fit_run.flux_2d_out['data-point source']
image_host_sersic = fit_run.image_host_list[0]
image_ps = fit_run.image_ps_list[0]
fluxes = []
fluxes_sersic = []
fluxes_ps = []
for i in range(splt):
    l = h/2
    xdis = np.abs(l * np.cos((theta.value)/180*np.pi))
    ydis = np.abs(l * np.sin((theta.value)/180*np.pi))
    x_c = f_center+fit_run.final_result_ps[0]['position_xy'][0]
    y_c = f_center++fit_run.final_result_ps[0]['position_xy'][1]
    x_ci = x_c - xdis + xdis/(splt/2)*i
    y_ci = y_c - ydis + ydis/(splt/2)*i
    aper = RectangularAperture((x_ci, y_ci), w, h/splt, theta=theta)
    if i >0: #%2 == 1:
        aper.plot(color='blue',
                  lw=1.5, label='comp {0}'.format(i))
    fluxes.append(aperture_photometry(image_host, aper)
                  ['aperture_sum'].value[0])
    fluxes_sersic.append(aperture_photometry(
        image_host_sersic, aper)['aperture_sum'].value[0])
    fluxes_ps.append(aperture_photometry(
        image_ps, aper)['aperture_sum'].value[0])
plt.show()
fluxes = np.array(fluxes)
fluxes_sersic = np.array(fluxes_sersic)
fluxes_ps = np.array(fluxes_ps)
x_data = np.linspace(0, len(fluxes)-1, len(fluxes)) * \
    (h/splt)*deltaPix  # Pixel grid in x, arcsec
plt.plot(x_data, fluxes/fluxes.max(),label='data-ps (host) profile')
plt.plot(x_data, fluxes_sersic/fluxes_sersic.max(), label='host sersic profile')
plt.plot(x_data, fluxes_ps/fluxes_ps.max(), label='Point source profile')
plt.xlabel('arcsec',  fontsize=16)
plt.ylabel('flux',  fontsize=16)
plt.legend(fontsize=10)
plt.show()


def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

parameters, covariance = curve_fit(func, x_data, fluxes)
sigma = parameters[-1]
print(sigma)
parameters, covariance = curve_fit(func, x_data, fluxes_sersic)
sigma = parameters[-1]
print(sigma)
parameters, covariance = curve_fit(func, x_data, fluxes_ps)
sigma = parameters[-1]
print(sigma)

read = fluxes