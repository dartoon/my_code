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
import pickle
import glob


test = gene_BHBH(h0=70)
dis_dotN, year_rate, Ctheta2 = test.num_year_rate(ave_chirp_mass=6.7)
print "The yearly rate using numrical calculation:", year_rate

a, mbh_max, mbh_min = 2.35, 80., 5.
filename = 'sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
if_file = glob.glob(filename)  
if if_file==[]:
    event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min)
    sim_data = [event_rate, zs_detected, masses, rhos_detected]
    pickle.dump(sim_data, open(filename, 'wb'))
else:
    event_rate, zs_detected, masses, rhos_detected=pickle.load(open(filename,'rb'))
print "The yearly rate using MC simulating:", event_rate

#averaged Chirp mass
print 'averaged Chirp mass, m1, m2: ', np.average(masses[:,0]), np.average(masses[:,1]), np.average(masses[:,2])
lum_dis = dl(zs_detected)
zs, chirp_mass, m1, m2, lumi_dis = zs_detected, masses[:,0], masses[:,1], masses[:,2], lum_dis                                                      
                                                      
#plt.hist(masses[:,:2])
#plt.hist(zs_detected)
#plt.hist(m1)
#plt.hist(rhos_detected)
#plt.plot(m1,m2,'.')
##plt.yscale('log', nonposy='clip')
#plt.show()
#import corner
#fig = corner.corner(np.column_stack((m1,m2)), labels=["$m1$", "$m2$"],title_kwargs={"fontsize": 12})
##                    quantiles=[0.16, 0.84],show_titles=True,
##                    title_kwargs={"fontsize": 12})
##                    plot_datapoints=True,smooth=1.0,smooth1d=1.0,
##                    levels=1.0 - np.exp(-0.5 * np.arange(1, 2.1, 1) ** 2))
##                    range=[(0.15,0.45),(60,90)] )
#plt.show()

#==============================================================================
#Setting the noise level 
#==============================================================================
#dl_nlevel = 0.3
m_noise_level = 0.20
m1_nlevel= m1 * m_noise_level
m1_noise = m1 + np.random.normal(0, m1_nlevel, size=m1_nlevel.shape)
plt.hist([m1[1000:],m1_noise[1000:]])
plt.yscale('log', nonposy='clip')
plt.show()
