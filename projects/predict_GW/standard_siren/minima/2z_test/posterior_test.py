#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:09:04 2018

@author: Dartoon

Use 2 redshift distribution to test whether the posterior could present the data.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

H0=70
Om=0.3
c=299790.

def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z,om):
      return quad(EZ, 0, z , args=(om))[0]
 
z1=1.
z2=3.

pz = 0.7
def gene_data(number,error_l): 
    seed = np.random.random(number)    
    zs = np.zeros(number) + z1
    zs[seed>pz]=z2
    dl = np.zeros(number)
    dl[seed<pz]= (1+z1)*c*EE(z1,Om)/H0 #assume zs =3
    dl[seed>pz]= (1+z2)*c*EE(z2,Om)/H0 #assume zs =1
#    error_l = 10 #Set error bar as 10%
    erbar = dl * error_l/100.
    bias=np.random.normal(0,erbar)
    dl_real =dl+bias
    return dl_real, erbar, zs
    #bias=np.random.normal(0,erbar)

def lnlike(theta, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta
    zs = np.array([z1,z2])
    vec_EE=np.vectorize(EE)
    rzs=vec_EE(zs,om)
    model = (1+zs)*c*rzs/h0             #133 numbers of the models
#    print model.shape
    likh = np.exp(-0.5*(pow((y[:,None]-model),2.)/pow(err[:,None],2.)))*np.array([pz,1-pz])
#    print likh.shape
    int_likh=np.sum(likh,axis=1)
    return np.sum(np.log(int_likh))

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
bnds = ((0, 1), (0, None))
error_level =  input("how do you want to set you error level? 10? 15? 20?:\n")  #default as 10,000
points=10 # input("how many minima points?:\n")  #default as 10,000
for loop in range(points):
    y,err,zs = gene_data(10000,error_level) 
    result = op.minimize(nll, (0.35, 65), method='SLSQP', bounds=bnds,args=(y, err))
    print "para:", result["x"], "\n\r",
    like_r= lnlike((0.3,70),y, err)
    like_w= lnlike((0.25,72),y, err)
    print like_r, like_w, like_r-like_w

#d1 = (1+z1)*c*EE(z1,Om)/H0
#d2 = (1+z2)*c*EE(z2,Om)/H0
#k = d2 - d1
#sig = error_level/100.
#print "exp wrong:,", 0.7*10000*np.log(0.7*np.exp(-0.5)+0.3*np.exp(-0.5*(1+k**2/(d1*sig)**2)))+\
#0.3*10000*np.log(0.3*np.exp(-0.5)+0.7*np.exp(-0.5*(1+k**2/(d2*sig)**2)))

#==============================================================================
# Calculate the expected value
#==============================================================================
def post(d, om, h0, err_l, dtype=1):
    z1, z2 = 1., 3.
    zs = np.array([z1,z2])
    vec_EE=np.vectorize(EE)
    rzs=vec_EE(zs,om)
    D1,D2 = (1+zs)*c*rzs/h0
#    print D1, D2, err_l
#    f = (d-mu)**2  #Test the Sigma
    om_fid, h0_fid = 0.3, 70
    rzs_fid=vec_EE(zs,om_fid)
    D1_fid,D2_fid = (1+zs)*c*rzs_fid/h0_fid
    if dtype==1:
        mu, sig = D1_fid, D1_fid* err_l
    elif dtype==2:
        mu, sig = D2_fid, D2_fid* err_l
    f = np.log(0.7*np.exp(-0.5*((d-D1)/sig)**2)+0.3*np.exp(-0.5*((d-D2)/sig)**2))
    f_q = f * 1. /(np.sqrt(2*np.pi)*sig)*np.exp(-0.5*((d-mu)/sig)**2)    
    return f_q

def exp_post(om, h0, err_l, dtype):
    exp = quad(post,0, 25422*1.5, args=(om, h0, err_l, dtype))
    return exp

err_l = error_level/100.

om, h0=0.3,70
exp_right = 7000*exp_post(om=om, h0=h0, err_l=err_l,dtype=1)[0]+3000*exp_post(om=om, h0=h0, err_l=err_l,dtype=2)[0]
print "exp",exp_right
om_w, h0_w=0.25,72
exp_wrong = 7000*exp_post(om=om_w, h0=h0_w, err_l=err_l,dtype=1)[0]+3000*exp_post(om=om_w, h0=h0_w, err_l=err_l,dtype=2)[0]
print "exp",exp_right,exp_wrong, exp_right-exp_wrong
 

    