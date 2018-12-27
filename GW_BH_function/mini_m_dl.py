#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Generate the simulating data
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from BH_mass_function import gene_BHBH, dl, solve_z


test = gene_BHBH(h0=70)
dis_dotN, year_rate, Ctheta2 = test.num_year_rate(ave_chirp_mass=6.7)
print "The yearly rate using numrical calculation:", year_rate

event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=2.35, mbh_max=80., mbh_min=5.)
print "The yearly rate using MC simulating:", event_rate

#Test the data
#event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=2.35, mbh_max=7.72, mbh_min=7.71)
#print "The yearly rate using MC simulating:", event_rate

plt.hist(masses[:,:2])
plt.show()

plt.hist(zs_detected)
plt.show()

#averaged Chirp mass
np.average(masses[:,0])
lum_dis = dl(zs_detected)

#The simulated data:
zs, chirp_mass, m1, m2, lumi_dis = zs_detected, masses[:,0], masses[:,1], masses[:,2], lum_dis
  
plt.plot(m1,m2,'.')
plt.show()

plt.hist(m1)
plt.yscale('log', nonposy='clip')
plt.show()

plt.hist(rhos_detected)
plt.yscale('log', nonposy='clip')
plt.show()


#import corner
#fig = corner.corner(np.column_stack((m1,m2)), labels=["$m1$", "$m2$"],title_kwargs={"fontsize": 12})
##                    quantiles=[0.16, 0.84],show_titles=True,
##                    title_kwargs={"fontsize": 12})
##                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
##                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2))
##                    range=[(0.15,0.45),(60,90)] )
#plt.show()