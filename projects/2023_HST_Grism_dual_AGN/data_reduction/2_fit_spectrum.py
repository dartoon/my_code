#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:32:11 2023

@author: Dartoon
"""
import astropy.io.fits as pyfits
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
import glob
from astropy.wcs import WCS
# %matplotlib inline/''

# from 1_analysis_spectrum import z_dic
z_dic = {"SDSSJ1246-0017": 2.559, "SDSSJ1502+0257": 1.484, "SDSSJ1625+4309":1.647, 
            "SDSSJ2304-0038": 2.775, "SDSSJ2206+0030": 1.444, "SDSSJ0224+0149":2.721}


# from hstaxe import axetasks

cwd = os.getcwd()
# print("We are in %s" % (cwd))

#%%
name = "SDSSJ1246-0017"  #G102   #!!!
# name = "SDSSJ1625+4309"    #G141   #!!!
# name = "SDSSJ1502+0257"    #G141
# name = "SDSSJ2304-0038"    #G102
# name = "SDSSJ2206+0030"    #G141
# name = "SDSSJ0224+0149"    #G102

ID = 2   #Dual ID
slit_width = 3

z = z_dic[name]

CIV, MgII, Hb, OIII, Halpha = 1549, 2798, 4861, 5007, 6563

os.chdir(cwd+'/'+name)
try:
    filt_spec = glob.glob('G102')[0]
    filt_image = glob.glob('F105W')[0]
except:
    filt_spec, filt_image = 'G141', 'F140W'

from matplotlib.colors import LogNorm
import copy, matplotlib
if filt_spec == 'G141':
    x1, x2 = 11500, 16500
else:
    x1, x2 = 8000, 11500
    
#%%
# =============================================================================
# For calibrate the spectrum. The black line is the data and the orange line is the aXe inferred flux
# The ratio will be used to recalibrate the spectrum of each dual
# =============================================================================
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')

fits = pyfits.open("./DRIZZLE/aXeWFC3_{0}_mef_ID{1}.fits".format(filt_spec, ID))  #This file provide the header info
header = fits['SCI'].header

#Show spectrum
d = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.STP.fits".format(filt_spec))["BEAM_%dA" % (ID)].data #This file is use to infer the flux in SPC (checked manually)
wcs = WCS(header, naxis=1, relax=False, fix=False)

lam = wcs.wcs_pix2world(np.arange(len(d.T)), 0)[0]
plt.plot((lam*10**10)[(lam*10**10>x1) & (lam*10**10<x2)], 
         (np.sum(d,axis=0)/(lam*10**24))[(lam*10**10>x1) & (lam*10**10<x2)], c= 'k', zorder =100, linewidth = 4)

fin = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.SPC.fits".format(filt_spec))
tdata = fin["BEAM_%dA" % (ID)].data
x = tdata["LAMBDA"]
f = tdata["FLUX"]   #!!! This is the total flux, for dual together.
e = tdata["FERROR"]
c = tdata["CONTAM"]
vg = (x>x1) & (x<x2)
plt.plot(x[vg],f[vg])
plt.errorbar(x[vg],f[vg],e[vg])   #The total spectrum
# plt.plot(x[vg],c[vg])
plt.xlim([x1-100,x2+100])
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux ($erg/s/cm^2/\AA/s$)');
plt.close()  #Compare 1D spec from 2D grism data VS. aXe SPC to get a unit_ratio

unit_ratio = f[vg]/(np.sum(d,axis=0)/(lam*10**24))[(lam*10**10>x1) & (lam*10**10<x2)]

#%% Show each line's spectrunm

print("The selected region to generate the SPC are not the same for the DRZZLE and OUTPUT")

for s in glob.glob("OUTPUT/*2.SPC.fits"):
    # print (s)
    d1 = pyfits.open(s)["BEAM_%dA" % (ID)].data
    w = d1["LAMBDA"]
    f = d1["FLUX"]
    e = d1["FERROR"]
    c = d1["CONTAM"]
    wg = (w>x1) & (w<x2)
    plt.errorbar(w[wg],f[wg],e[wg])
    # plt.plot(w[wg],c[wg])
fin = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.SPC.fits".format(filt_spec))
tdata = fin["BEAM_%dA" % (ID)].data
x = tdata["LAMBDA"]
f = tdata["FLUX"]
e = tdata["FERROR"]
c = tdata["CONTAM"]
vg = (x>x1) & (x<x2)
#plt.errorbar(x[vg],y[vg],e[vg])
plt.plot(x[vg],f[vg],color='k',lw=5)
plt.errorbar(x[vg],f[vg],e[vg],color='k',lw=2)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux ($erg/s/cm^2/\AA/s$)');
plt.xlim([x1-100,x2+100])
plt.close() #Compare 1D spec from aXe SPC for DRZZILE and flt


#%% Show 2D image
from scipy.signal import argrelextrema
print(name, 'redshift:', z, 'filters:', filt_image, filt_spec)
image_file = glob.glob(filt_image+'/*_drz.fits')[0]
image = pyfits.open(image_file)['SCI'].data
pos = np.loadtxt(filt_image+'/cookbook.cat')[ID-1][:2]
plt.imshow(image[int(pos[1])-20:int(pos[1])+20,int(pos[0])-20:int(pos[0])+20], 
          norm=LogNorm(), origin='lower',cmap = my_cmap)

l_stk = np.sum(image[int(pos[1])-8:int(pos[1])+8,int(pos[0])-3:int(pos[0])+3], axis = 1)
l_stk[l_stk<5] = 0
id_y_pos = argrelextrema( l_stk, np.greater)[0]
if len(id_y_pos) != 2:
    raise ValueError("The NO. of PS is not 2.")
elif abs(id_y_pos - 8)[0]<abs(id_y_pos - 8)[1]:
    ID_pos = 0
else:
    ID_pos = 1
    
plt.axhline(id_y_pos[ID_pos]+20-8-int((slit_width)/2),c='white')
plt.axhline(id_y_pos[ID_pos]+20-8+int((slit_width)/2),c='white')
plt.show()

#%% Show 1D spec from Grism image, for DRZZILE and flt
fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
fits = pyfits.open("./DRIZZLE/aXeWFC3_{0}_mef_ID{1}.fits".format(filt_spec, ID))  #This file provide the header info
header = fits['SCI'].header
#Show spectrum
d = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.STP.fits".format(filt_spec))["BEAM_%dA" % (ID)].data #This file is use to infer the flux in SPC (checked manually)
l_stk = np.sum(d[:,75:150], axis=1)
l_stk = np.nan_to_num(l_stk)
l_stk[l_stk<np.max(l_stk)/10] = 0
id_y_pos = argrelextrema(l_stk, np.greater)[0]
if len(id_y_pos) == 2:
    id_y_pos_ = id_y_pos[ID_pos]
elif len(id_y_pos) == 1:
    id_y_pos_ = id_y_pos[0]
else:
    plt.imshow(d[:,75:150], norm=LogNorm(), origin='lower',cmap = my_cmap)
    plt.show()
    raise ValueError("The NO. of PS is not 2.")
y = id_y_pos_ - int(slit_width/2)
d_sub = d[y:y+slit_width,:]

wcs = WCS(header, naxis=1, relax=False, fix=False)
lam = wcs.wcs_pix2world(np.arange(len(d.T)), 0)[0]

driz_spec = unit_ratio*(np.sum(d_sub,axis=0)/(lam*10**24))[(lam*10**10>x1) & (lam*10**10<x2)]
plt.plot((lam*10**10)[(lam*10**10>x1) & (lam*10**10<x2)], driz_spec, c= 'k', zorder =100, linewidth = 4)

#For the fitting later
lam = (lam*10**10)[(lam*10**10>x1) & (lam*10**10<x2)]
flux = driz_spec


s_list = glob.glob("OUTPUT/*2.STP.fits")
d_flt_list = []
for i, s in enumerate(s_list):
    # print (s)
    d_flt = pyfits.open(s)["BEAM_%dA" % (ID)].data
    header_flt = pyfits.open(s)["BEAM_%dA" % (ID)].header
    # plt.imshow(d1, norm=LogNorm(), origin='lower',cmap = my_cmap)
    # y = int(len(d_flt)/2)-int(slit_width/2)

    l_stk = np.sum(d_flt[:,50:150], axis=1)
    l_stk = np.nan_to_num(l_stk)
    l_stk[l_stk<np.max(l_stk)/10] = 0
    id_y_pos = argrelextrema(l_stk, np.greater)[0]
    if len(id_y_pos) == 2:
        id_y_pos_ = id_y_pos[ID_pos]
    elif len(id_y_pos) == 1:
        id_y_pos_ = id_y_pos[0]
    else:
        id_y_pos_ = id_y_pos[:2][ID_pos]
        # plt.imshow(d[:,75:150], norm=LogNorm(), origin='lower',cmap = my_cmap)
        # plt.show()
        # raise ValueError("The NO. of PS is not 2.")
                
    y = id_y_pos_ - int(slit_width/2)
    d_flt = d_flt[y:y+slit_width,:]
    d_flt_list.append(d_flt)
    wcs = WCS(header_flt, naxis=1, relax=False, fix=False)
    lam_flt = wcs.wcs_pix2world(np.arange(len(d_flt.T)), 0)[0]
    # if i == 1:
    plt.plot((lam_flt*10**10)[(lam_flt*10**10>x1) & (lam_flt*10**10<x2)],unit_ratio*(np.sum(d_flt,axis=0)/(lam_flt*10**24))[(lam_flt*10**10>x1) & (lam_flt*10**10<x2)])
lines = CIV, MgII, Hb, OIII, Halpha
plt.xlim([x1-100,x2+100])
plt.ylim([0.8*np.min(driz_spec), 1.2*np.max(driz_spec[10:])])
for line in lines:
    if line * (1+z) > x1 and line * (1+z) < x2:
        plt.axvline(x = line * (1+z), color = 'b', label = str(line))
        y_max = ax.get_ylim()[1]
        plt.text(line * (1+z), y_max*0.9 , str(line))
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux ($erg/s/cm^2/\AA/s$)')
plt.show()


#%% Show entrie 2D spectrum
xticks = np.int0(ax.get_xticks())
wcs = WCS(header, naxis=1, relax=False, fix=False)
print('Final drizzled Spectrum:')
fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
ax.imshow(d, norm=LogNorm(), origin='lower',cmap = my_cmap)
l_stk = np.sum(d[:,50:150], axis=1)
l_stk[l_stk<3] = 0
id_y_pos = argrelextrema(l_stk, np.greater)[0]
if len(id_y_pos) == 2:
    id_y_pos_ = id_y_pos[ID_pos]
elif len(id_y_pos) == 1:
    id_y_pos_ = id_y_pos[0]
else:
    plt.imshow(d[:,50:150], norm=LogNorm(), origin='lower',cmap = my_cmap)
    plt.show()
    raise ValueError("The NO. of PS is not 2.")
y = id_y_pos_ - int(slit_width/2)
plt.axhline(y,c='white')
plt.axhline(y+slit_width,c='white')
pos = [np.sum((lam*10**10 - xticks[i])<0) for i in range(len(xticks))]
ax.set_xticks(pos)
ax.set_xticklabels(xticks)
plt.show()

#%%
# =============================================================================
# #Each dither's 2D's case
# =============================================================================
for i,s in enumerate(s_list):
    print(s)
    d_flt = pyfits.open(s)["BEAM_%dA" % (ID)].data
    header_flt = pyfits.open(s)["BEAM_%dA" % (ID)].header
    l_stk = np.sum(d_flt[:,50:150], axis=1)
    l_stk = np.nan_to_num(l_stk)
    l_stk[l_stk<np.max(l_stk)/10] = 0
    id_y_pos = argrelextrema(l_stk, np.greater)[0]
    if len(id_y_pos) == 2:
        id_y_pos_ = id_y_pos[ID_pos]
    elif len(id_y_pos) == 1:
        id_y_pos_ = id_y_pos[0]
    else:
        # plt.imshow(d[:,75:150], norm=LogNorm(), origin='lower',cmap = my_cmap)
        # plt.show()
        # raise ValueError("The NO. of PS is not 2.")
        id_y_pos_ = id_y_pos[:2][ID_pos]
        print(id_y_pos_)
    y = id_y_pos_ - int(slit_width/2)
    fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
    ax.imshow(d_flt, norm=LogNorm(), origin='lower',cmap = my_cmap)    
    plt.axhline(y,c='white')
    plt.axhline(y+slit_width,c='white')
    plt.show()
    lam_flt = wcs.wcs_pix2world(np.arange(len(d_flt.T)), 0)[0]
    y = id_y_pos_ - int(slit_width/2)
    d_flt = d_flt[y:y+slit_width,:]
    flux_flt = unit_ratio*(np.sum(d_flt,axis=0)/(lam_flt*10**24))[(lam_flt*10**10>x1) & (lam_flt*10**10<x2)]
    # if i == 0:
    #     wavelength = (lam_flt*10**10)[(lam_flt*10**10>x1) & (lam_flt*10**10<x2)]/(1+z)
    #     F_lambda = flux_flt

# fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
# plt.plot(wavelength, F_lambda)
# # plt.xlim([x1-100,x2+100])
# plt.ylim([-4*np.min(driz_spec), 0.5*np.max(driz_spec[10:])])
# plt.show()

#%% The 1D spectrum
# fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
# plt.plot(lam/(1+z), flux)
# # plt.xlim([x1-100,x2+100])
# plt.ylim([-4*np.min(driz_spec), 0.5*np.max(driz_spec[10:])])
# plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import emcee

# # Assume data is loaded into arrays called wavelength and F_lambda
wavelength = lam/(1+z)
F_lambda = flux

F_lambda_  = F_lambda[wavelength>2500] * 1e18
wavelength_  = wavelength[wavelength>2500]

# Define the Gaussian function
def gaussian(wavelength, amp, mu, sigma):
    return amp * np.exp(-(wavelength - mu)**2 / (2. * sigma**2))

# Define the model: a power-law plus a sum of Gaussians
def model(wavelength, amp_power, index, gaussian_params):
    power_law = amp_power * (wavelength / 3000.) ** index
    total_gaussian = sum(gaussian(wavelength, *params) for params in gaussian_params)
    return power_law + total_gaussian

# Define the objective function for minimization (chi-squared)
def log_likelihood(params):
    amp_power, index, *gaussian_params_flat = params
    gaussian_params = [(gaussian_params_flat[i], gaussian_params_flat[i+1], gaussian_params_flat[i+2])
                       for i in range(0, len(gaussian_params_flat), 3)]
    model_eval = model(wavelength_, amp_power, index, gaussian_params)
    # chi_squared = np.sum((F_lambda - model_eval)**2 / model_eval)
    chi_squared = - 0.5* np.sum( (F_lambda_ - model_eval) ** 2 / np.std(F_lambda_)**2 )
    # print(chi_squared, amp_power, index, gaussian_params)
    return chi_squared



# # Define the objective function for minimization (chi-squared)
# def log_likelihood(params):
#     model_eval = model(wavelength_, params)
#     # chi_squared = np.sum((F_lambda_ - model_eval)**2 / model_eval)
#     chi_squared =  - 0.5* np.sum( (F_lambda_ - model_eval) ** 2 / np.std(F_lambda_)**2 )
#     # print(chi_squared, params)
#     return chi_squared

# Initial guess
params_initial = [np.mean(F_lambda_), -1.5, 1., 2798, 20.]
# params_initial = [np.mean(F_lambda_), -1.5, 1., 4861, 20., 1., 5007, 10.]
# plt.plot(wavelength, F_lambda)
# plt.plot(wavelength, model(wavelength, params_initial))
# plt.show()

nll = lambda *args: -log_likelihood(*args)

# Minimization
result_min = minimize(nll, params_initial,args=())
params_best_fit = result_min.x

amp_power_best, index_best, *gaussian_params_best = params_best_fit

gaussian_params_best_split = [(gaussian_params_best[i], gaussian_params_best[i+1], gaussian_params_best[i+2])
                   for i in range(0, len(gaussian_params_best), 3)]

plt.plot(wavelength, F_lambda * 1e18)
plt.plot(wavelength_, model(wavelength_, amp_power_best, index_best, gaussian_params_best_split))
# plt.ylim(0,2*np.max(model(wavelength_, params_best_fit)))
plt.show()
chi_squared = np.sum( ((F_lambda_ - model(wavelength_, amp_power_best, index_best, gaussian_params_best_split)) / np.std(F_lambda_) )**2 )
print(params_best_fit, chi_squared)

# #%% MCMC is not needed
# def log_prior(params):
#     amp_power, index, amp_gauss, sigma = params
#     if 0.1 < amp_power < 50 and 0.5 < index < 3 and 0.1 < amp_gauss < 4.0 and 0<sigma<200:
#         return 0.0
#     return -np.inf
# def log_probability(params):
#     lp = log_prior(params)
#     if not np.isfinite(lp):
#         return -np.inf
#     return lp + log_likelihood(params)

# # Now set up the MCMC
# ndim, nwalkers = len(params_initial), 100
# sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)

# # Initial positions of walkers
# pos = [params_best_fit + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# # Run MCMC
# sampler.run_mcmc(pos, 10000, progress = True)

# # Optionally, plot the results
# samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

# # for i in range(ndim):
# #     plt.figure()
# #     plt.hist(samples[:, i], 100, color="k", histtype="step")
# #     plt.title("Dimension {0:d}".format(i))
# # plt.show()

# import corner
# flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
# fig = corner.corner(
#     flat_samples)

# for i in range(ndim):
#     mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
#     print(mcmc[1])
#     # q = np.diff(mcmc)
#     # txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
#     # txt = txt.format(mcmc[1], q[0], q[1])
#     # print(txt)

