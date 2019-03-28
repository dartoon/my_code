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
from BH_mass_function import Theta, solve_z, random_z, cal_chirpmass
from BH_mass_function import dl as cal_dl
from scipy.integrate import quad
from scipy import interpolate
import pickle

def mass_fun_i(m1, a, mbh_max, mbh_min):
    if m1>=mbh_min and m1<=mbh_max:
        N = (mbh_max)**(-a+1)/(1-a) - (mbh_min)**(-a+1)/(1-a)             # The intergral function of m**a for normalize the power law function
        return m1**(-a)/N
    else:
        return 0.
poss_mass_fun=np.vectorize(mass_fun_i)

def poss_gaussian(m1, mu, sigma):
    poss =  1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(m1 - mu)**2/(2*sigma**2))
    return poss

def poss_ln_gaussian(m1, mu, sigma):
    poss =  1/(m1 * sigma * np.sqrt(2 * np.pi)) * np.exp(-(np.log(m1) - mu)**2/(2*sigma**2))
    return poss
"""
However, it turns out the likelihood is not "necessary":
#def likelihood(m1_obs, dm1_obs, a=2.35, mbh_max=80., mbh_min=5., bins = None):
##    m1_samp_grid = np.logspace(np.log(mbh_min-0.1),np.log(mbh_max+1),800)
#    if bins==None:
#        if abs(m1_obs-mbh_min)>5:
#            bins = 51
#        else:
#            max_bins = 161
#            min_bins = 51
#            gap = 5   # from 5-5 to 10-5
#            bins = int(max_bins-(m1_obs-mbh_min)*(max_bins-min_bins)/gap)     #Set the max bins as 1001
#        m1_samp_grid = np.linspace(m1_obs-5*dm1_obs,m1_obs+5*dm1_obs,bins)   #!!! This is wrong!!!
#        m1_samp_grid = m1_samp_grid[m1_samp_grid>0]
#        if (m1_obs-5*dm1_obs)<mbh_min<(m1_obs+5*dm1_obs): 
#            m1_samp_grid = np.append(m1_samp_grid,mbh_min)
#            m1_samp_grid = m1_samp_grid[m1_samp_grid.argsort()]
#        if (m1_obs-5*dm1_obs)<mbh_max<(m1_obs+5*dm1_obs):
#            m1_samp_grid = np.append(m1_samp_grid,mbh_max)
#            m1_samp_grid = m1_samp_grid[m1_samp_grid.argsort()]
#    likeli = (poss_gaussian(m1_samp_grid, mu=m1_obs, sigma=dm1_obs) * poss_mass_fun(m1_samp_grid, a=a, mbh_max=mbh_max, mbh_min=mbh_min))[1:] * (m1_samp_grid[1:]-m1_samp_grid[:-1])
#    return likeli.sum()
#def likelihood_lognorm(m1_obs, dm1_obs, a=2.35, mbh_max=80., mbh_min=5., use_method='med', f = None):
#    '''
#    The likelihood for log-norm
#    The way to summarize the likelihood should be re-write.
#    '''
##    if abs(m1_obs-mbh_min)>5:
##        bins = 51
##    else:
##        max_bins = 51
##        min_bins = 51
##        gap = 5   # from 5-5 to 10-5
##        bins = int(max_bins-(m1_obs-mbh_min)*(max_bins-min_bins)/gap)     #Set the max bins as 1001
##    bins = 41
##    sigma = dm1_obs/m1_obs  # Recover the noise level
##    m1_samp_grid = np.linspace(m1_obs/np.exp(sigma)**5,m1_obs*np.exp(sigma)**5,bins)
##    if use_method == 'med':
##        mu = np.log(m1_obs)             #0_med, i.e. directly use the med as generated
##    elif use_method == 'exp':
##        mu = np.log(m1_obs) - sigma**2/2.        #1_trans_exp np.log(mu) = np.log(exp-sig**2/2)
##    if (m1_obs/np.exp(sigma)**5)<mbh_min<(m1_obs*np.exp(sigma)**5):
##        m1_samp_grid = np.linspace(m1_obs/np.exp(sigma)**5,m1_obs*np.exp(sigma)**5,bins)
##        m1_samp_grid = np.append(mbh_min,m1_samp_grid)
##        m1_samp_grid = m1_samp_grid[m1_samp_grid.argsort()]
##    if (m1_obs/np.exp(sigma)**5)<mbh_max<(m1_obs*np.exp(sigma)**5):
##        m1_samp_grid = np.linspace(m1_obs/np.exp(sigma)**5,m1_obs*np.exp(sigma)**5,bins)
##        m1_samp_grid = np.append(m1_samp_grid,mbh_max)
##        m1_samp_grid = m1_samp_grid[m1_samp_grid.argsort()]
##    print m1_samp_grid.min(), m1_samp_grid.max()
##    if f == None:
##        x = np.logspace(np.log10(2/1.1), np.log10(400*1.1),300)
##        y = cov_twofun(x, a=a, mbh_max=mbh_max, mbh_min=mbh_min,sigma=sigma)
##        f = interpolate.interp1d(x, y)  # Use interpolate to make cal faster
#    joint = poss_ln_gaussian(m1_samp_grid, mu=mu, sigma=sigma) * f(m1_samp_grid)
#    likeli = (joint[1:]+joint[:-1])/2 * (m1_samp_grid[1:]-m1_samp_grid[:-1])
#    return likeli.sum()
##Test if likelihood is well defined and the difference by changing the bins.
##print mass_fun_i(79,a = 2.35, mbh_max=80, mbh_min=5),
##print likelihood_lognorm( 5.09 , 5*0.2 ,a = 2.35, mbh_max=80, mbh_min=5)
##print (likelihood(m,dm)-likelihood(m,dm, bins=100002))/likelihood(m,dm, bins=100002) * 100, '%'
##print (likelihood(6,1, bins=400).sum()-likelihood(6,1).sum())/likelihood(6,1).sum()
"""

def random_Theta(seed_vol=100000):
    '''
    Purpose:
    Use MC to realize the selection effect given the: Chripmass, dl
    Return
    --------
        A sth sth
    '''
    theta_class = Theta(vol = seed_vol)
    thetas = theta_class.gene_theta()
    print "This means the random_Theta calculation is on, which shouldn't appear more than once"
    return thetas

def m1_fac_select_effect(m1, mbh_min=5, thetas=None, itera=5000, r0 = 1527.):
    '''
    Purpose:
    Use MC to evaluate the possibility of detecting m1
    Return
    --------
        The possibility of detecting m1
    '''    
    if thetas is None:
        thetas = random_Theta()
    m2 = np.random.uniform(mbh_min, m1, size = itera)
    mass_Chirp = cal_chirpmass(m1,m2)
    zs = random_z(itera=itera)
    dl = cal_dl(zs)
    over_rho0 = []
    for i in range(itera):
        rhos = 8.*thetas * r0/dl[i] * ((1+zs[i])*mass_Chirp[i]/1.2)**(5/6.)
        over_rho0.append(np.sum([rhos>8])/float(len(rhos)))
    over_rho0 = np.asarray(over_rho0)
    return over_rho0
#thetas = random_Theta()
#over_rho0 = m1_fac_select_effect(m1=6, thetas=thetas)
#print 'm1_fac_select_effect', np.average(over_rho0)
def m1_array_select_effect(mbh_min=5., thetas=None, itera=10000, r0 = 1527., bins=250):
    '''
    Purpose:
    based on def m1_fac_select_effect to derive a m1_array prior
    '''    
    m1s = np.logspace(np.log10(mbh_min+0.01), np.log10(200), bins)
    priors = np.zeros(bins)
    for i in range(bins):
        m1 = m1s[i]
        priors[i] = np.average(m1_fac_select_effect(m1=m1,mbh_min=mbh_min, thetas=thetas, itera=itera, r0 = r0))
        if i/20 > (i-1)/20:
            print "Total bins:", bins, "; Finish i:", i
    return m1s, priors
#thetas = random_Theta()
#import time
#t1= time.time()
#m1s, priors = m1_array_select_effect(mbh_min=5, thetas=thetas)
#plt.plot(m1s,priors)
#plt.xlim(0,150)
#plt.show()
#t2 = time.time()
#print (t2-t1)/60. , "mins"    
#pickle.dump([m1s, priors], open('select_effect_MBHmin5', 'wb'))
def select_effect(m1, fname='select_effect_MBHmin5'):
    x, y = pickle.load(open(fname,'rb'))
    f = interpolate.interp1d(x, y)
    if x[0]<m1<x[-1]:
        return f(m1)
    elif m1:
        return 0
def joint_prior_lognorm(tau, t, sigma =0.2):
    return select_effect(tau) * poss_ln_gaussian(t, mu=np.log(tau), sigma = sigma)
def cov_prior_lognorm(t, sigma =0.2):
    inter = quad(joint_prior_lognorm, 0, 200, args=(t,sigma))[0]
    return inter
#cov_prior_lognorm=np.vectorize(cov_prior_lognorm)   
#m1_obs = np.logspace(np.log10(0.01), np.log10(150), 200)
#cov_priors = cov_prior_lognorm(m1_obs)
#plt.plot(m1_obs,cov_priors)
#plt.xlim(0,150)
#pickle.dump([m1_obs, cov_priors], open('select_effect_MBHmin5_cov_lognorm0.2', 'wb'))
#plt.show()


def joint_twofun(tau, t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    return mass_fun_i(tau, a=a, mbh_max=mbh_max, mbh_min=mbh_min) * poss_ln_gaussian(t, mu=np.log(tau), sigma = sigma)
def cov_twofun(t, a=2.35, mbh_max=80, mbh_min=5, sigma =0.2):
    inter = quad(joint_twofun, 0, 120, args=(t, a, mbh_max, mbh_min,sigma))[0]
    return inter
cov_twofun=np.vectorize(cov_twofun)   

#x = np.logspace(np.log10(3/1.1), np.log10(80*1.1),50)  #50 could be smaller?
#y = cov_twofun(x, a=2.35, mbh_max=80,mbh_min=5,sigma=0.2) #
#plt.plot(x,y)
#plt.show()


def fac_s_eff_s(dl,mass_Chirp,thetas=None,r0 = 1527.):
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
##Test the func: fac_s_eff_s
#thetas = random_Theta()
#print 'the possibility of detecting:', fac_s_eff_s(dl=5000, mass_Chirp=7, thetas=thetas)

def fac_s_eff_v(dl,mass_Chirp,thetas=None,r0 = 1527.):
    '''
    Purpose:
        Same as fac_s_eff_s but vectorize
    '''    
    if thetas is None:
        thetas = random_Theta()
    vec_solve_z = np.vectorize(solve_z)
    zs = vec_solve_z(dl)
    need_theta = 8./(8. * r0/dl * ((1+zs)*mass_Chirp/1.2)**(5/6.))
#    print need_theta
    ratios = np.sum(need_theta<thetas[:,None],axis=0)/float(len(thetas))
    return ratios
#thetas = random_Theta()
#print 'the possibility of detecting:', fac_s_eff_v(dl = np.linspace(5000,8000,3), mass_Chirp = np.linspace(7,8,3), thetas=thetas)


def twod_fac_s_eff_s(dl,mass_Chirp,thetas=None,r0 = 1527.):
    '''
    same as fac_s_eff_s but in 2-dimensional. dl as first D and mass as second D.
    '''
    if thetas is None:
        thetas = random_Theta()
    vec_solve_z = np.vectorize(solve_z)
    v_zs = vec_solve_z(dl)
    need_theta = 8./(8. * r0/dl[:,None] * ((1+v_zs[:,None])*mass_Chirp/1.2)**(5/6.))
    v_ratio = np.sum(need_theta<thetas[:,None,None],axis=0)/float(len(thetas))
    return v_ratio
##Test twod_fac_s_eff_s
#thetas = random_Theta()
#v_ratio = twod_fac_s_eff_s(dl = np.linspace(7000,8000,11), mass_Chirp = np.linspace(7,8,10), thetas=thetas)
#print v_ratio.shape


def fac_s_eff_un(dl, dl_sig, mass_Chirp, mass_Chirp_sig,thetas=None,r0 = 1527):
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
    twod_fac = twod_fac_s_eff_s(dl=dl_grid,mass_Chirp=chirpm_grid,thetas=thetas,r0 = r0)
    marginalize_aix0 = np.sum( (twod_fac * poss_gaussian(dl_grid,dl,dl_sig)[:,None])[1:] *(dl_grid[1:]-dl_grid[:-1]), axis=0)
    marginalize_aix1 = np.sum( (marginalize_aix0 * poss_gaussian(chirpm_grid,mass_Chirp,mass_Chirp_sig))[1:] * (chirpm_grid[1:]- chirpm_grid[:-1]))
    return marginalize_aix1
##Test func fac_s_eff_un
#thetas = random_Theta()
#fac_s_eff_un = fac_s_eff_un(7000, 1000,9, 1,thetas=thetas)
#print fac_s_eff_s(7000,9,thetas=thetas), fac_s_eff_un

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

## Test if poss_ln_gaussian is well normizled:
#mu, sigma= 5., 1.
#sigma_star = np.exp(sigma/mu)
#m = np.linspace(mu-5*sigma,mu+5*sigma,41)
#num = poss_ln_gaussian(m, mu=np.log(mu), sigma=sigma/mu)
#plt.plot(m, num)
#plt.show()
#print mu/sigma_star,mu ,mu*sigma_star
#print 'if == 1:', np.sum(num[1:] * (m[1:]-m[:-1]))

##Test if likelihood is well defined and the difference by changing the bins.
#m, dm=5, 1
#print mass_fun_i(m), likelihood(m,dm)
#print (likelihood(m,dm)-likelihood(m,dm, bins=100002))/likelihood(m,dm, bins=100002) * 100, '%'
##print (likelihood(6,1, bins=400).sum()-likelihood(6,1).sum())/likelihood(6,1).sum()

##Test if likelihood is well defined and the difference by changing the bins.
##print mass_fun_i(79,a = 2.35, mbh_max=80, mbh_min=5),
#print likelihood_lognorm(2.50183031, 0.50036606,a = 2.35, mbh_max=80, mbh_min=5)
##print (likelihood(m,dm)-likelihood(m,dm, bins=100002))/likelihood(m,dm, bins=100002) * 100, '%'
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
