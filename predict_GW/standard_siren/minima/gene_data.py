#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 20:33:16 2018

@author: dartoon

a function to generate the data
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from pz import pz

def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z,om):
      return quad(EZ, 0, z , args=(om))[0]

def gene_data(vol,error_ls=(5,10),H0=70., om=0.3, reverse=False):
    """
    A function to generate data understand LCDM model with Om=0.3, H0=70:
        input:
            vol=the volume of the data size
            error_l= error level, time 100, i.e. the number before '%'
        output:
            z, dl, dl_noised, errorbar
    """
    R=pz(om,H0)
    zs=np.zeros(vol)
    error_l = np.random.uniform(error_ls[0],error_ls[1],vol)
    for i in range(vol):
        idx = int(np.sum(np.random.random()>R[:,2]))-1
        zs[i] = R[idx, 0] #np.random.uniform(R[idx, 0],R[idx+1, 0])
    vec_EE=np.vectorize(EE)
    c=299790.        #speed of light [km/s]
    dl=(1+zs)*c*vec_EE(zs,om)/H0            # unit in Mpc 
    DL=np.zeros([len(zs),5])
    if reverse==False:
        DL[:,0]=zs            #zs
        DL[:,1]=dl            #true dl distance 
        erbar=dl*error_l/100.   #1% level as the based
        DL[:,3]=erbar        #error
        DL[:,2]=np.random.normal(dl,erbar)       #after noise
        DL[:,4]=DL[:,2]*error_l/100.  #Realisic noise level?
    elif reverse==True:
        DL[:,0]=zs            #zs
        DL[:,1]=dl            #true dl distance 
        erbar=1/dl*error_l/100.   #1% level as the based
        DL[:,3]=dl*error_l/100.        #error
        bias=np.random.normal(0,erbar)   ##include system error to the data, still with seed(1), normal and random is different type, it's linear to dl*error/100
        dl_nois=(1/dl+bias)**-1
        DL[:,2]=dl_nois       #after noise
        DL[:,4]=dl_nois*error_l/100.  #Realisic noise level?
    return DL
