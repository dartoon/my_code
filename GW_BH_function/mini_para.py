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
import random
import scipy.optimize as op
import time

a, mbh_max, mbh_min = 2.35, 80., 5.
filename = 'sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
if_file = glob.glob(filename)  
if if_file==[]:
    event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min)
    dl_zs = dl(zs_detected)
    sim_data = [event_rate, zs_detected, masses, rhos_detected, dl_zs]
    pickle.dump(sim_data, open(filename, 'wb'))
else:
    event_rate, zs_detected, masses, rhos_detected, dl_zs=pickle.load(open(filename,'rb'))
#print "The yearly rate using MC simulating:", event_rate
#
##averaged Chirp mass
#print 'averaged Chirp mass, m1, m2: ', np.average(masses[:,0]), np.average(masses[:,1]), np.average(masses[:,2])
zs, chirp_mass, m1, m2, lumi_dis = zs_detected, masses[:,0], masses[:,1], masses[:,2], dl_zs                                                      
                                                      
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
m_noise_level = 0.20
m1_nlevel= m1 * m_noise_level
m1_obs = m1 + np.random.normal(0, m1_nlevel, size=m1_nlevel.shape)

dl_l = 0.3
dl_sig = lumi_dis * dl_l
dl_obs =lumi_dis + np.random.normal(0, dl_sig, size=dl_sig.shape)

chirp_m_l = 0.08
chirp_m_sig = chirp_mass*chirp_m_l
chirp_m_obs = chirp_mass + np.random.normal(0, chirp_m_sig, size=chirp_m_sig.shape)

from cal_likelihood import posterior
import time
#def posterior(m1_obs, dm1_obs, dl, dl_sig,mass_Chirp, mass_Chirp_sig, a=2.35, mbh_max=80., mbh_min=5.,thetas=None, r0 = 1527):

bool_pos =  [dl_obs>0]
#Remove all the negative luminosities
m1_obs, m1_nlevel, dl_obs,dl_sig, chirp_m_obs, chirp_m_sig = m1_obs[bool_pos], m1_nlevel[bool_pos], dl_obs[bool_pos],dl_sig[bool_pos], chirp_m_obs[bool_pos], chirp_m_sig[bool_pos]

def lnprob(para, sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig):
#    a=2.35, mbh_max=80., mbh_min=5.
    ######################MCMC#########################
    a, mbh_max, mbh_min = para
    if 1.1 < a < 3 and 79 < mbh_max < 81 and 4.5 < mbh_min < 5.5:
        vec_posterior = np.vectorize(posterior)
        post = vec_posterior(sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
        return -np.sum(np.log(post))
    else:
        return np.inf

nll = lambda *args: lnprob(*args)
bnds = ((1.2, 3), (79.9, 80.1),(4.9, 5.1))

from scipy.optimize import fmin
para = (2.35, 80, 5)

shuffe_list = range(len(m1_obs))
vol = 50
t1 = time.time()
for loop in range(1):
    random.shuffle(shuffe_list)
    sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig =\
m1_obs[shuffe_list][:vol], m1_nlevel[shuffe_list][:vol], dl_obs[shuffe_list][:vol],dl_sig[shuffe_list][:vol], chirp_m_obs[shuffe_list][:vol], chirp_m_sig[shuffe_list][:vol]
#    print 'lnprob', lnprob((2.25, 80., 5.), sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig)
#    print 'lnprob', lnprob((2.35, 80., 5.), sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig)
#    print 'lnprob', lnprob((2.45, 80., 5.), sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig)
    result = op.minimize(nll, para, method='SLSQP', bounds=bnds,args=(sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig))
    print result["x"]
#    mini=fmin(lnprob,para,maxiter=10000, args=(sam_m1_obs, sam_dm1_obs, sam_dl, sam_dl_sig, sam_mass_Chirp, sam_mass_Chirp_sig))
t2 = time.time()
print t2-t1, 's'