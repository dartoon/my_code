#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:32:44 2020

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xy 
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

# Create x and y indices
x = np.linspace(0, 200, 201)
y = np.linspace(0, 200, 201)
x, y = np.meshgrid(x, y)

#create data
xy = (x, y)
data = twoD_Gaussian(xy, 3, 100, 100, 20, 40, 0, 10)

# plot twoD_Gaussian data generated above
plt.figure()
plt.imshow(data.reshape(201, 201))
plt.colorbar()

# add some noise to the data and try to fit the data generated beforehand
initial_guess = (3,100,100,20,40,0,10)

data_noisy = data + 0.2*np.random.normal(size=data.shape)
import scipy.optimize as opt
popt, pcov = opt.curve_fit(twoD_Gaussian, xy, data_noisy, p0=initial_guess)
data_fitted = twoD_Gaussian((x, y), *popt)
fig, ax = plt.subplots(1, 1)
ax.imshow(data_noisy.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    extent=(x.min(), x.max(), y.min(), y.max()))
ax.contour(x, y, data_fitted.reshape(201, 201), 8, colors='w')
plt.show()

#%%Test fit data:
data = pyfits.getdata('PSF_data.fits')
data = pyfits.getdata('QSO_img.fits')

# plt.plot(range(len(data)**2), data.ravel())
x = np.linspace(0, len(data)-1, len(data))
y = np.linspace(0, len(data)-1, len(data))
x, y = np.meshgrid(x, y)
xy = (x, y)
popt, pcov = opt.curve_fit(twoD_Gaussian, xy, data.ravel(), p0=(3,len(data)/2,len(data)/2,2,2,0,10))
popt[3] = abs(popt[3])
popt[4] = abs(popt[4])
data_fitted = twoD_Gaussian(xy, *popt)
fig, ax = plt.subplots(1, 1)
ax.imshow(data, cmap=plt.cm.jet, origin='bottom',
    extent=(x.min(), x.max(), y.min(), y.max()))
ax.contour(x, y, data_fitted.reshape(len(data), len(data)), 0.2, colors='w')
plt.show()
print(popt[4]/popt[3], popt[4])