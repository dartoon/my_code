#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 16:05:10 2020

@author: dartoon
"""

import numpy as np
import matplotlib.pyplot as plt

infers  = np.loadtxt('./local_Vardha_data.txt', dtype=str)
z = infers[:, 3].astype(np.float)
MBH = infers[:, 6].astype(np.float)
Mstar = infers[:, 11]

z = z[Mstar!= '...']
MBH = MBH[Mstar!= '...']
Mstar = Mstar[Mstar!= '...']

Mstar = [Mstar[i].split('$\\pm$') for i in range(len(Mstar))]
Mstar = np.array(Mstar)
Mstar = Mstar.astype(np.float)

x=Mstar[:,0]
y = MBH
yerr = (Mstar[:,1] **2  + 0.4 **2)**0.5

# MBH = [[MBH[i], 0.4] for i in range(len(MBH))]
# y=np.array(MBH)

plt.errorbar(Mstar[:,0],MBH, xerr=Mstar[:,1] ,yerr=MBH*0 + 0.4,
             fmt='.',color='black',markersize=15)

def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
#    yerr=(m*(np.append(bloc[:,2], hloc[:,2]))**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # right error level
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.58, 1.55, 0.35], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]

xp = np.array([5, 13])
#plt.plot(xp, m_ml*xp+b_ml, 'r-')
def lnprior(theta):
    m, b, sint	 = theta
    if -5.0 < m < 5 and -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 1000)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
xl = np.linspace(5, 13, 100)
m, b, sint =np.percentile(samples, 50,axis=0)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)

def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
m=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    m=np.percentile(samples,posi,axis=0)[0]
    b=samples[:,1][samples[:,0]==find_n(samples[:,0],m)][0]   #may find out many numbers
    plt.plot(xl, m*xl+b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1000)
#plt.text(9.3, 6.24, r"log(M$_{\rm BH}/$10$^{7}$M$_{\odot}$)=%s+%slog(M$_*/$10$^{10}$M$_{\odot}$)"%(round(b_ml+m_ml*10-7,2),round(m_ml,2)),color='blue',fontsize=25)
