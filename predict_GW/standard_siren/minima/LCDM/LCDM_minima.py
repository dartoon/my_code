#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 19:25:55 2018

@author: dartoon
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../')
from pz import pz
from gene_data import gene_data
##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).
from likelihoods import lnprob

# =============================================================================
# ###############minmize lnprob#######################
# =============================================================================
############to generate the likehood with zs############
error_l=input("which error level?:\n")
writefile=open('test_{0}_boost_16'.format(error_l),'w')
writefile.write("# Omage H0"+"\n")
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
bnds = ((0, None), (0, None))
import time
points=input("how many minima points?:\n")  #default as 10,000
for loop in range(points):
    ticks1=time.time()
    dl = gene_data(10000,error_l)
    z=dl[:,0]
    y=dl[:,2]  
    err=dl[:,3] 
    result = op.minimize(nll, (0.30, 70), method='SLSQP', bounds=bnds,args=(y, err))
    ticks2=time.time()
    print "the precentage:", float(loop)/points*100, "%;", "time remain:", round((ticks2-ticks1)*(points-loop)/60,2), "mins;", "para:", result["x"], "\n\r",
    writefile.write(repr(result["x"][0]) + " " + repr(result["x"][1])+ "\n")
writefile.close()

#################perfrom MCMC##################
#import emcee
#z=z
#y=y           #set to the 1 truth, 2 biased
#err=err
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(y, err))  ###the input args() should change
#sampler.run_mcmc(pos, 500)
#
################save the data###############
#
#pickle.dump(sampler.chain, open("nozs_lcdm__%s"%(value),'wb'))
#
################load the data##############
#samplerchain=pickle.load(open("nozs_lcdm__%s"%(value),'rb'))
#burn=samplerchain[:,:,:].T
#plt.plot(burn[0,20:,:], '-', color='k', alpha=0.3)  #show the chain after 50 steps 
#samples = samplerchain[:, 40:, :].reshape((-1, ndim))
#import corner
#fig = corner.corner(samples, labels=["$om$", "$h0$"],
#                       quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
#                     #,  plot_datapoints=False,smooth=2.0,smooth1d=2.0,plot_density=False,levels=(0.6826, 0.9544),\
#                     #  color='#0000ff',show_titles=True, title_fmt='.3f',title_kwargs={"fontsize": 16} )
##fig.savefig("triangle.png")
######
#plt.show()   

