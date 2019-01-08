#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:57:36 2018

@author: Dartoon

Calculate the likelihood for by given the m1
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from BH_mass_function import Theta, solve_z

def mass_fun_i(m1, a=2.35, mbh_max=80., mbh_min=5.):
    if m1>=mbh_min and m1<=mbh_max:
        N = (mbh_max)**(-a+1)/(1-a) - (mbh_min)**(-a+1)/(1-a)             # The intergral function of m**a for normalize the power law function
        return m1**(-a)/N
    else:
        return 0.
poss_mass_fun=np.vectorize(mass_fun_i)

def poss_gaussian(m1, mu, sigma):
    poss =  1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(m1 - mu)**2/(2*sigma**2))
    return poss

def likelihood(m1_obs, dm1_obs, a=2.35, mbh_max=80., mbh_min=5., bins = None):
#    m1_samp_grid = np.logspace(np.log(mbh_min-0.1),np.log(mbh_max+1),800)
    if bins==None:
        if abs(m1_obs-mbh_min)>5:
            bins = 51
        else:
            max_bins = 5001
            min_bins = 51
            gap = 5   # from 5-5 to 10-5
            bins = int(max_bins-(m1_obs-mbh_min)*(max_bins-min_bins)/gap)     #Set the max bins as 1001
    m1_samp_grid = np.linspace(m1_obs-5*dm1_obs,m1_obs+5*dm1_obs,bins)
    likeli = (poss_gaussian(m1_samp_grid, m1_obs, dm1_obs) * poss_mass_fun(m1_samp_grid, a=a, mbh_max=mbh_max, mbh_min=mbh_min))[1:] * (m1_samp_grid[1:]-m1_samp_grid[:-1])
    return likeli.sum()

def random_Theta(seed_vol=100000):
    '''
    Purpose:
    Use MC to realize the selection effect given the: Chripmass, dl
    Parameter
    --------
        a: The blash of blash
        b: The blash of blash
        c: The blash of blash
    Return
    --------
        A sth sth
    '''
    theta_class = Theta(vol = seed_vol)
    thetas = theta_class.gene_theta()
    print "This means the random_Theta calculation is on, which shouldn't appear more than once"
    return thetas

def fac_s_eff_s(dl,mass_Chirp,thetas=None,r0 = 1527):
    '''
    Purpose:
    Use the MC result of Theta's to evalute the selection effect factor given the specific value
    of dl and mass_Chrip.
    Parameter
    --------
        dl: The luminosity distance
        mass_Chirp: The Chirp mass
        theta: The result of the Thetas distribution. If None, will MC by it self.
    Return
    --------
        The possibility
    
    Note
    --------
        The zs is calculated correspondingly from the dl, based on Om=0.3 and H0=70.
    '''    
    if thetas is None:
        thetas = random_Theta()
    zs = solve_z(dl)
    rhos = 8.*thetas * r0/dl * ((1+zs)*mass_Chirp/1.2)**(5/6.)
    ratio = len(rhos[rhos>8])/float(len(rhos))
    return ratio

def count_poss(need_bin, a):
    '''
    Purpose
    Calculate the 'possibility' in b which exceed need_bin
    '''
    return 1-a[:need_bin].sum()/a.sum()

def twod_fac_s_eff_s(v_dl,v_mass_Chirp,thetas=None,r0 = 1527):
    '''
    same as fac_s_eff_s but in 2-dimensional. dl as first D and mass as second D.
    '''
    if thetas is None:
        thetas = random_Theta()
    vec_solve_z = np.vectorize(solve_z)
    v_zs = vec_solve_z(v_dl)
    need_theta = 8./(8. * r0/v_dl[:,None] * ((1+v_zs[:,None])*v_mass_Chirp/1.2)**(5/6.))
#    hist_a,hist_b,_ = plt.hist(thetas, bins=201)
#    plt.close()
#    need_bins = np.sum((need_theta>hist_b[:,None,None]),axis=0)
#    v_ratio = np.zeros((11,10))
#    for i in range(len(v_ratio)):
#        for j in range(len(v_ratio.T)):
#            v_ratio[i][j] = count_poss(need_bins[i][j], hist_a)
##            print i, j, count_poss(need_bins[i][j], hist_a)
    v_ratio = np.sum(need_theta<thetas[:,None,None],axis=0)/float(len(thetas))
    return v_ratio
#Test twod_fac_s_eff_s
#thetas = random_Theta()
#v_ratio = twod_fac_s_eff_s(v_dl = np.linspace(7000,8000,11), v_mass_Chirp = np.linspace(7,8,10), thetas=thetas)
#print v_ratio

#def twod_fac_s_eff_s(v_dl,v_mass_Chirp,thetas=None,r0 = 1527):
#    '''
#    same as fac_s_eff_s but in 2-dimensional. dl as first D and mass as second D.
#    '''
#    if thetas is None:
#        thetas = random_Theta()
#    vec_solve_z = np.vectorize(solve_z)
#    v_zs = vec_solve_z(v_dl)
#    v_rhos = 8.*thetas * r0/v_dl[:,None,None] * ((1+v_zs[:,None,None])*v_mass_Chirp[:,None]/1.2)**(5/6.)  # This is supposed in 3D...
#    v_ratio = np.sum(v_rhos>8,axis=2)/float(len(thetas))
#    return v_ratio
##Test twod_fac_s_eff_s
#thetas = random_Theta()
#v_ratio = twod_fac_s_eff_s(v_dl = np.linspace(7000,8000,11), v_mass_Chirp = np.linspace(7,8,10), thetas=thetas)
#print v_ratio

def fac_s_eff_un(dl, dl_sig,mass_Chirp, mass_Chirp_sig,thetas=None,r0 = 1527):
    '''
    Purpose:
    Similar to fac_s_eff_s, but the data are provided with uncertainty. The possiblity
    is calculated by marginazed over the space.
    
    Parameter
    --------
        dl, dl_sig: The luminosity distance and the uncertainty level
        mass_Chirp, mass_Chirp_sig: The Chirp mass and the uncertainty level
        theta: The result of the Thetas distribution. If None, will MC by it self.
    Return
    --------
        The possibility
    
    Note
    --------
        The zs is calculated correspondingly from the dl, based on Om=0.3 and H0=70.
    '''
    dl_low = dl-3*dl_sig
    if dl_low<0:
        dl_low = 1
    dl_grid = np.linspace(dl_low,dl+3*dl_sig,num=11)
    chirpm_grid = np.linspace(mass_Chirp-5*mass_Chirp_sig,mass_Chirp+5*mass_Chirp_sig,num=10)
    twod_fac = twod_fac_s_eff_s(v_dl=dl_grid,v_mass_Chirp=chirpm_grid,thetas=thetas,r0 = r0)
    marginalize_aix0 = np.sum( (twod_fac * poss_gaussian(dl_grid,dl,dl_sig)[:,None])[1:] *(dl_grid[1:]-dl_grid[:-1]), axis=0)
    marginalize_aix1 = np.sum( (marginalize_aix0 * poss_gaussian(chirpm_grid,mass_Chirp,mass_Chirp_sig))[1:] * (chirpm_grid[1:]- chirpm_grid[:-1]))
    return marginalize_aix1

thetas = random_Theta()
def posterior(m1_obs, dm1_obs, dl, dl_sig,mass_Chirp, mass_Chirp_sig, a=2.35, mbh_max=80., mbh_min=5., r0 = 1527):
    '''
    Purpose:
    Calculate the posterior given the data include the m1, dl, chirpmass and their uncertainties.

    Parameter
    --------
    m1_obs, dm1_obs: m1 and m1's uncertainty

    Return
    --------
        The possibility
    '''
    prior = fac_s_eff_un(dl=dl, dl_sig=dl_sig,mass_Chirp=mass_Chirp, mass_Chirp_sig=mass_Chirp_sig,thetas=thetas,r0 = r0)
    likeli = likelihood(m1_obs, dm1_obs, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
    poster = prior * likeli
    return poster


#==============================================================================
#For the test: 
#==============================================================================
### Test if poss_mass_fun is well normizled:
#a, mbh_max, mbh_min = 2.35, 80., 5.
#m = np.logspace(np.log10(mbh_min-0.1),np.log10(mbh_max+1),801)        #It truns out that log space is more effective to simpling
##m = np.linspace(5,100,800)
#poss = poss_mass_fun(m,a=a, mbh_max=mbh_max, mbh_min=mbh_min)
#plt.plot(m, poss)
#plt.show()
#print 'if == 1:', np.sum(poss[1:] * (m[1:]-m[:-1]))

## Test if poss_gaussian is well normizled:
#mu, sigma= 40, 2
##m = np.logspace(np.log10(mbh_min),np.log10(mbh_max),800)        #It truns out that log space is more effective to simpling
#m = np.linspace(mu-5*sigma,mu+5*sigma,41)
#num = poss_gaussian(m, mu=mu, sigma=sigma)
#plt.plot(m, num)
#plt.show()
#print 'if == 1:', np.sum(num[1:] * (m[1:]-m[:-1]))

##Test if likelihood is well defined and the difference by changing the bins.
#m, dm=5, 1
#print mass_fun_i(m), likelihood(m,dm)
#print (likelihood(m,dm)-likelihood(m,dm, bins=100002))/likelihood(m,dm, bins=100002) * 100, '%'
##print (likelihood(6,1, bins=400).sum()-likelihood(6,1).sum())/likelihood(6,1).sum()

##Test the func: fac_s_eff_s
#thetas = random_Theta()
#print 'the possibility of detecting:', fac_s_eff_s(dl=7849.2, mass_Chirp=5.14, thetas=thetas)

##Test twod_fac_s_eff_s
#thetas = random_Theta()
#v_ratio = twod_fac_s_eff_s(v_dl = np.linspace(7000,8000,11), v_mass_Chirp = np.linspace(7,8,10), thetas=thetas)
#print v_ratio.shape

##Test func fac_s_eff_un
#thetas = random_Theta()
#fac_s_eff_un = fac_s_eff_un(7000, 1000,9, 1,thetas=thetas)
#print fac_s_eff_s(7000,9,thetas=thetas), fac_s_eff_un

##Test func posterior
#print posterior(6,1,7000,1000,5.6,1)
